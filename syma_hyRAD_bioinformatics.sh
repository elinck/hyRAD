################################################################################
# hyRAD data bioinformatics pipeline (used for S. torotoro; Aves: Halcyonidae) #
# based around Berkeley QB3's denovoTargetedCapturePopGen wrapper scripts      #
# developed w/ Zach Hanna (https://github.com/calacademy-research              #
################################################################################

#!/bin/bash
# requires: pyrad, perl, samtools, bedtools blastn, angsd, QB3 denovoTargetCapturePopGen scripts (https://github.com/CGRL-QB3-UCBerkeley/denovoTargetCapturePopGen.git)

########################################
# copy scripts over                    #
########################################

mkdir /data/elinck/syma_assembly
cd /data/elinck/syma_assembly/
git clone https://github.com/CGRL-QB3-UCBerkeley/denovoTargetCapturePopGen.git 

########################################
# count the number of reads per sample #
########################################

mkdir /data/elinck/syma_assembly/probes
mkdir /data/elinck/syma_assembly/captures

for i in *; do echo $i; gunzip -c $i | echo $((`wc -l`/4)); done > hiseq_readcounts.txt

##########################################################
# run pyRAD on probe library to assembly probe sequences #
##########################################################

cd /data/elinck/syma_assembly/probes
python /home/elinck/bin/pyrad/pyrad/pyRAD.py -p /home/elinck/syma_assembly/probes/syma_params.txt -s 12345

########################################
# Cleaning raw reads                   #
########################################

mkdir /data/elinck/syma_assembly/pipeline_07212016/
mkdir /data/elinck/syma_assembly/pipeline_07212016/raw/ #obtained raw reads from QB3 via FTP, put here

# obtain ecoli reference genome over to test or contamination
cp -r /data/zhanna/bin/denovoTargetCapturePhylogenomics/ecoli/ /data/elinck/syma_assembly/pipeline_07212016/

# run 1-ScrubReads_a
perl /data/elinck/syma_assembly/1-ScrubReads_a cleanPE -t /home/elinck/bin/Trimmomatic0.36/trimmomatic-0.36.jar -f /data/elinck/syma_assembly/pipeline_07212016/raw/ -k 12 -z -c /data/elinck/syma_assembly/pipeline_07212016/ecoli/e_coli_K12.fasta -o /data/elinck/syma_assembly/pipeline_07212016/cleaned/ 1>/data/elinck/syma_assembly/pipeline_07212016/2016Jul21_scrubreads_run1.log 2>/data/elinck/2016Jul21_scrubreads_run1.err

#how many duplicates?
cat *.fq | echo $((`wc -l`/4)) > reads_original.txt

cat *.fq | echo $((`wc -l`/4)) > reads_precelean.txt

# 1 - (reads_preclean / sample count in reads_original) = % duplicates

########################## 
# Generating Assemblies  #
##########################

mkdir /data/elinck/syma_assembly/pipeline_07212016/assemblies

# run 2-GenerateAssembliesPhylo
perl /data/elinck/syma_assembly/2-GenerateAssembliesPhylo spades -reads /data/elinck/syma_assembly/pipeline_07212016/cleaned -out /data/elinck/syma_assembly/pipeline_07212016/assemblies -np 16 1 > assembly_log1.txt 2 > assembly_error1.txt

##########################
# Finding Targets        #
##########################

#make directories to find targets separately by sample type (modern or ancient)
mkdir /data/elinck/syma_assembly/pipeline_07212016/assemblies/raw/reference

mkdir /data/elinck/syma_assembly/pipeline_07212016/assemblies/raw/ancient

# make a blast database from rad library
makeblastdb -in /data/elinck/syma_assembly/probes/clust.97/t_ochr_29.consens -parse_seqids -dbtype nucl

# blast rad library against itself to determine repeats
blastn -db /data/elinck/syma_assembly/probes/clust.97/t_ochr_29.consens  -query /data/elinck/syma_assembly/probes/clust.97/t_ochr_29.consens -evalue 0.00001 -outfmt 6 -out radprobe_selfblast_25june.txt

cut -f 1,1 radprobe_selfblast_25june.txt | sort | uniq -u | wc -l # 168814 unique results

cut -f 1,1 radprobe_selfblast_25june.txt | sort | uniq -u > rad_probe_unique.txt

blastdbcmd -db /data/elinck/syma_assembly/probes/clust.97/t_ochr_29.consens -entry_batch /data/elinck/syma_assembly/probes/clust.97/rad_probe_unique.txt > rad_probe_unique.fasta

# run kebi's script to find targets

perl /data/elinck/syma_assembly/3-FindingTargetsV8 combineExon -t /data/elinck/syma_assembly/probes/clust.97/rad_probe_unique.fasta -a /data/elinck/syma_assembly/pipeline_07212016/assemblies/raw/reference -p 0.95 -b 1 -e 4 l

perl 3-FindingTargetsV8 combineExon -t /data/elinck/syma_assembly/probes/clust.97/rad_probe_unique.fasta -a /data/elinck/syma_assembly/pipeline_07212016/assemblies/raw/ancient -p 0.95 -b 1 -e 4

# how many captured sequences?
grep ">" -c *.fa

##########################
# Evaluate Assemblies    #
##########################

#abyss  assembly stats evaluator (not used)
for i in /data/elinck/syma_assembly/pipeline_07212016/assemblies/raw/*fasta; do abyss-fac $i; echo $i; done > readcounts.txt;

##historic samples, in target regions

#rename files in /data/elinck/syma_assembly/pipeline_07212016/assemblies/raw/ancient/In_target/
find -type f -name "EL_hyRAD*" -exec rename 's/.intargetPremasked.fa/.fasta/' \{\} \;

#assembly eval 
perl /data/elinck/syma_assembly/AssemblyEvaluation BASIC -a /data/elinck/syma_assembly/pipeline_07212016/assemblies/raw/ancient/In_target/premasked

## historic samples, all assembled regions 

#rename files /data/elinck/syma_assembly/pipeline_07212016/assemblies/raw/ancient/original/
find -type f -name "EL_hyRAD*" -exec rename 's/.fasta.original/.fasta/' \{\} \;

#assemby eval 
perl /data/elinck/syma_assembly/AssemblyEvaluation BASIC -a /data/elinck/syma_assembly/pipeline_07212016/assemblies/raw/ancient/original

## modern samples, in target regions
find -type f -name "EL_hyRAD*" -exec rename 's/.intargetPremasked.fa/.fasta/' \{\} \;

#assembly eval 
perl /data/elinck/syma_assembly/AssemblyEvaluation BASIC -a /data/elinck/syma_assembly/pipeline_07212016/assemblies/raw/reference/In_target/premasked

## modern samples, all assembled regions
find -type f -name "EL_hyRAD*" -exec rename 's/.fasta.original/.fasta/' \{\} \;

#assembly eval 
perl /data/elinck/syma_assembly/AssemblyEvaluation BASIC -a /data/elinck/syma_assembly/pipeline_07212016/assemblies/raw/reference/original

#########################################
# GC content filtering (not run for ms) #
#########################################

# assembly eval reveals high GC content! use biopieces to filter raw reads for those 3 standard deviations (3 x 0.27 = 0.83) above mean for modern samples (44.57 + 0.83 = 45.4)
# run prinseq with 3 stdev GC cutoff for one sample (test)
perl /home/elinck/bin/prinseq-lite-0.20.4/prinseq-lite.pl -fastq EL_hyRAD_001T_S29_1_final.fq -max_gc 45.4

# run prinseq for all files, GC cutoff of 50
for i in *.fq; do perl /home/elinck/bin/prinseq-lite-0.20.4/prinseq-lite.pl -fastq $i -max_gc 50; done

# check stats -- how does GC content compare?
perl /home/elinck/bin/prinseq-lite-0.20.4/prinseq-lite.pl -fastq EL_hyRAD_001A_S29_1_final.fq -stats_all

##################################################
# Blast reference genome for mtDNA contamination #
##################################################

#upload t_sanctus_ref.fasta from desktop to /data/elinck/syma_assembly/pipeline_07212016/assemblies/raw/reference/In_target/Final/

# make a blast database from published kingfisher mtDNA (Todiramphus sanctus)
makeblastdb -in /data/elinck/syma_assembly/pipeline_07212016/assemblies/raw/reference/In_target/Final/t_sanctus_ref.fasta -parse_seqids -dbtype nucl

# blast reference genome against kingfisher mtDNA blast database
blastn -db t_sanctus_ref.fasta -query combined_targetedRegionAndFlanking.fasta -evalue 0.00001 -outfmt 6 -out mtDNA_in_reference_blast.txt

####################################################
# Blast reference genome for nonvert contamination #
####################################################

#blast reference genome against local blast nt database
blastn -db ncbi_nt -query combined_targetedRegionAndFlanking.fasta -evalue 0.00001  -outfmt 6 -out nonvertcontamination_in_reference_blast.txt

#directions for further processing (run by Z. Hanna; may be updated w/ actual commands)

#output top five hits for each contig

#take subset that align to nonvert ref genomes

#blast again

#evaluate results, enter both mtDNA hit contig IDs and nonvert contig IDs in one row per ID text file, feed to cutcontigs.R (see README.md)

##########################
# Alignment              #
##########################

mkdir /data/elinck/syma_assembly/pipeline_07212016/alignment

# align to extended reference genome
perl 5-Alignment -f /data/elinck/syma_assembly/pipeline_07212016/assemblies/raw/reference/In_target/Final/combined_targetedRegionAndFlanking.fasta -r /data/elinck/syma_assembly/pipeline_07212016/cleaned -o /data/elinck/syma_assembly/pipeline_07212016/alignment -i 235 -v 24 -l -t 90 -P /home/elinck/bin/picard-tools-2.4.1/picard.jar -G /home/elinck/bin/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -p 16

#evaluate exon capture performance (not working)
perl /data/elinck/syma_assembly/6-ExonCaptureEvaluation Evaluation -genome /data/elinck/syma_assembly/pipeline_07212016/assemblies/raw/reference/In_target/Final/combined_targetedRegionAndFlanking.fasta -cleanDir /data/elinck/syma_assembly/pipeline_07212016/cleaned -rawDir /data/elinck/syma_assembly/pipeline_07212016/raw/pre-clean -bamDir /data/elinck/syma_assembly/pipeline_07212016/alignment -InstrID HS -readLen 100 -resDir /data/elinck/syma_assembly/pipeline_07212016/evaluation_20160808 -bedFile /data/elinck/syma_assembly/pipeline_07212016/assemblies/raw/ancient/In_target/Final/combined_targetedRegionforExonCapEval.bed 1>/data/elinck/syma_assembly/20160808_eval1.log 2>/data/elinck/syma_assembly/20160808_eval1.err

# subtract preclean from original, tk

# determine specificity per sample 
for bam_file in *.bam
do
total=$(samtools view -c $bam_file)
mapped=$(samtools view -c -F 4 $bam_file)
unmapped=$(samtools view -c -f 4 $bam_file)
echo "$bam_file $total $mapped $unmapped"
done

# determine average depth per sample
for bam_file in *.bam
do
average_depth=$(samtools depth $bam_file |  awk '{sum+=$3} END { print "Average = ",sum/NR}')
echo "$bam_file $average_depth"
done

##########################
# SNP calling            #
##########################

mkdir /data/elinck/syma_assembly/pipeline_07212016/alignment/angsd/
mkdir /data/elinck/syma_assembly/pipeline_07212016/alignment/angsd/historic
mkdir /data/elinck/syma_assembly/pipeline_07212016/alignment/angsd/modern

#move bams to subfolders
mv EL_hyRAD_001A* EL_hyRAD_001B* EL_hyRAD_001C* EL_hyRAD_001D* EL_hyRAD_001E EL_hyRAD_001F* ./modern
mv EL_hyRAD* ./historic

# prepare merged .bam file for variant discovery, modern samples
samtools merge merge.bam /data/elinck/syma_assembly/pipeline_07212016/alignment/angsd/historic/*.bam

# prepare merged .bam file for variant discovery, historic samples
samtools merge merge.bam /data/elinck/syma_assembly/pipeline_07212016/alignment/angsd/modern/*.bam

# sort merged .bam for modern
mkdir /data/elinck/syma_assembly/pipeline_07212016/alignment/angsd/modern/merged/

mv merge.bam merged
 
samtools sort /data/elinck/syma_assembly/pipeline_07212016/alignment/angsd/modern/merged/merge.bam merge_sorted

# sort merged .bam for historic
mkdir /data/elinck/syma_assembly/pipeline_07212016/alignment/angsd/modern/historic/merged

mv merge.bam merged
 
samtools sort /data/elinck/syma_assembly/pipeline_07212016/alignment/angsd/historic/merged/merge.bam merge_sorted

# prefiltering steps before calling variants
perl /data/elinck/syma_assembly/9-preFiltering percentile -b /data/elinck/syma_assembly/pipeline_07212016/alignment/angsd/historic/merged/merge_sorted.bam -o s_toro_historic

perl /data/elinck/syma_assembly/9-preFiltering percentile -b /data/elinck/syma_assembly/pipeline_07212016/alignment/angsd/modern/merged/merge_sorted.bam -o s_toro_modern

#make historic vcf
samtools mpileup -B -D -I -S -uf /data/elinck/syma_assembly/pipeline_07212016/assemblies/raw/reference/In_target/Final/combined_targetedRegionAndFlanking.fasta /data/elinck/syma_assembly/pipeline_07212016/alignment/angsd/historic/*sorted.bam | bcftools view -cg - > raw_historic_extended.vcf

#make modern vcf
samtools mpileup -B -D -I -S -uf /data/elinck/syma_assembly/pipeline_07212016/assemblies/raw/reference/In_target/Final/combined_targetedRegionAndFlanking.fasta /data/elinck/syma_assembly/pipeline_07212016/alignment/angsd/modern/*sorted.bam | bcftools view -cg - > raw_modern_extended.vcf

#SNP cleaning historic
perl /data/elinck/syma_assembly/10-SNPcleaner.pl -A /data/elinck/syma_assembly/pipeline_07212016/assemblies/raw/reference/In_target/Final/combined_targetedRegionAndFlanking.fasta -M CT_GA -d 5 -k 7 -u 3 -a 0 -r /data/elinck/syma_assembly/pipeline_07212016/alignment/angsd/historic/s_toro_gene_depth_percentile.txt -B syma_torotoro_historic.bed  -p syma_torotoro_dropped_historic -v /data/elinck/syma_assembly/pipeline_07212016/SNPs/raw_historic_extended.vcf > out_historic_cleaned.vcf

#SNP cleaning modern 
perl /data/elinck/syma_assembly/10-SNPcleaner.pl -A /data/elinck/syma_assembly/pipeline_07212016/assemblies/raw/reference/In_target/Final/combined_targetedRegionAndFlanking.fasta -M CT_GA -d 5 -k 3 -u 3 -a 0 -r /data/elinck/syma_assembly/pipeline_07212016/alignment/angsd/modern/s_toro_gene_depth_percentile.txt -B syma_torotoro_modern.bed  -p syma_torotoro_dropped_modern -v /data/elinck/syma_assembly/pipeline_07212016/SNPs/raw_modern_extended.vcf > out_modern_cleaned.vcf

#report count of alleles for all bed files
wc -l /data/elinck/syma_assembly/pipeline_07212016/alignment/angsd/*/*.bed

#find shared sites among extended ref genome
bedtools intersect -a /data/elinck/syma_assembly/pipeline_07212016/angsd/historic/syma_torotoro_historic.bed -b /data/elinck/syma_assembly/pipeline_07212016/alignment/angsd/modern/syma_torotoro_modern.bed > /data/elinck/syma_assembly/pipeline_07212016/alignment/angsd/shared_sites.bed

#make keep file for angsd
cut -f 1,2 /data/elinck/syma_assembly/pipeline_07212016/alignment/angsd/shared_sites.bed  > shared_sites.keep

#index sites for angsd
angsd sites index shared_sites.keep

#ANGSD, full data
angsd -bam /data/elinck/syma_assembly/pipeline_07212016/angsd/bam.filelist.txt -sites /data/elinck/syma_assembly/pipeline_07212016/angsd/shared_sites.keep -anc /data/elinck/syma_assembly/pipeline_07212016/assemblies/raw/reference/In_target/Final/combined_targetedRegionAndFlanking.fasta -minQ 20 -fold 1 -out snps_extended  -GL 1  -doGeno 2  -doPost 1  -postCutoff 0.95 -doCounts 1 -doGlf 2 -geno_minDepth 6 -SNP_pval 0.05 -minMaf 0.0125 -doMaf 2 -doMajorMinor 1 -doSaf 1 -doPlink 2

#ANGSD, no MD
angsd -bam /data/elinck/syma_assembly/pipeline_07212016/angsd/bam.filelist.txt -sites /data/elinck/syma_assembly/pipeline_07212016/angsd/shared_sites.keep -anc /data/elinck/syma_assembly/pipeline_07212016/assemblies/raw/reference/In_target/Final/combined_targetedRegionAndFlanking.fasta -minQ 20 -fold 1 -out snps_extended_all  -GL 1  -doGeno 2  -doPost 1  -postCutoff 0.95 -doCounts 1 -doGlf 2 -geno_minDepth 6 -SNP_pval 0.05 -minMaf 0.0125 -doMaf 2 -doMajorMinor 1 -doSaf 1 -doPlink 2 -minInd 19

#report number of loci
wc -l *geno 

#output 100% complete matrix genind for adegenet
perl /data/elinck/syma_assembly/11-PopGenTools Adegenet -g snps_extended.geno -n 19 -s 2750 -o snps_exended_adg

#output full data matrix genind for adegenet
perl /data/elinck/syma_assembly/11-PopGenTools Adegenet -g snps_extended_all.geno -n 19 -s 68545 -o snps_extended_all_adg

#############################################
# Downstream analyses (not reported in ms)  #
#############################################

# run NGSadmix on 75% matrix
for i in 1 2 3 4 5 6 7 8 9 10; do NGSadmix -likes snps_15.beagle -K $i -P 4 -o ngsadmix_${i}_15 -minMaf 0.0125

# run admixture on 75% matrix
for K in 1 2 3 4 5 6 7 8 9 10; do admixture --cv hyrad_admixture.bed $K | tee log${K}.out; done

# choose K with cross validation error
grep -h CV log*.out



