################################################################################
# hyRAD data bioinformatics pipeline (used for S. torotoro; Aves: Halcyonidae) #
# based around Berkeley QB3's denovoTargetedCapturePopGen wrapper scripts      #
# developed w/ Zach Hanna (https://github.com/calacademy-research              #
################################################################################

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

# how many duplicates reads?
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

# make directories to find targets separately by sample type (modern or ancient)
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

perl /data/elinck/syma_assembly/3-FindingTargetsV8 combineExon -t /data/elinck/syma_assembly/probes/clust.97/rad_probe_unique.fasta -a /data/elinck/syma_assembly/pipeline_07212016/assemblies/raw/ancient -p 0.95 -b 1 -e 4

# how many captured sequences?
grep ">" -c *.fa

##################################
# Evaluate and clean Assemblies  #
##################################

# find mtDNA seqs in assemblies
for i in /data/elinck/syma_assembly/pipeline_07212016/assemblies/raw/ancient/fasta/*.fasta; do blastn -db /data/elinck/syma_assembly/pipeline_07212016/assemblies/raw/reference/In_target/Final/t_sanctus_ref.fasta -query $i -evalue 0.00001 -outfmt 6 -out ${i}_mtDNAhits.txt; done

for i in /data/elinck/syma_assembly/pipeline_07212016/assemblies/raw/reference/In_target/premasked/test/*.fasta; do blastn -db /data/elinck/syma_assembly/pipeline_07212016/assemblies/raw/reference/In_target/Final/t_sanctus_ref.fasta -query $i -evalue 0.00001 -outfmt 6 -out ${i}_mtDNAhits.txt; done

# eliminate line breaks in fasta files for R treatment
for i in /data/elinck/syma_assembly/pipeline_07212016/assemblies/raw/ancient/intarget/.*fasta; do awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' $i > ${i}_oneline; done

# run functions in extractcontigIDs.R and cutcontigsbatch.R to remove mtDNA contaminant assemblies

# abyss  assembly stats evaluator (not used)
for i in /data/elinck/syma_assembly/pipeline_07212016/assemblies/raw/*fasta; do abyss-fac $i; echo $i; done > readcounts.txt;

## historic samples, in target regions

# rename files in /data/elinck/syma_assembly/pipeline_07212016/assemblies/raw/ancient/In_target/
find -type f -name "EL_hyRAD*" -exec rename 's/.intargetPremasked.fa/.fasta/' \{\} \;

# assembly eval
perl /data/elinck/syma_assembly/AssemblyEvaluation BASIC -a /data/elinck/syma_assembly/pipeline_07212016/assemblies/raw/ancient/In_target/premasked/test/

## historic samples, all assembled regions 

# rename files /data/elinck/syma_assembly/pipeline_07212016/assemblies/raw/ancient/original/
find -type f -name "EL_hyRAD*" -exec rename 's/.fasta.original/.fasta/' \{\} \;

# assemby eval 
perl /data/elinck/syma_assembly/AssemblyEvaluation BASIC -a /data/elinck/syma_assembly/pipeline_07212016/assemblies/raw/ancient/original

## modern samples, in target regions
find -type f -name "EL_hyRAD*" -exec rename 's/.intargetPremasked.fa/.fasta/' \{\} \;

# assembly eval 
perl /data/elinck/syma_assembly/AssemblyEvaluation BASIC -a /data/elinck/syma_assembly/pipeline_07212016/assemblies/raw/reference/In_target/premasked/test/eval/

## modern samples, all assembled regions
find -type f -name "EL_hyRAD*" -exec rename 's/.fasta.original/.fasta/' \{\} \;

# assembly eval 
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

# upload t_sanctus_ref.fasta from desktop to /data/elinck/syma_assembly/pipeline_07212016/assemblies/raw/reference/In_target/Final/

# make a blast database from published kingfisher mtDNA (Todiramphus sanctus)
makeblastdb -in /data/elinck/syma_assembly/pipeline_07212016/assemblies/raw/reference/In_target/Final/t_sanctus_ref.fasta -parse_seqids -dbtype nucl

# blast reference genome against kingfisher mtDNA blast database
blastn -db t_sanctus_ref.fasta -query combined_targetedRegionAndFlanking.fasta -evalue 0.00001 -outfmt 6 -out mtDNA_in_reference_blast.txt

# contigs blasting to  mtDNA: combined_Contig2441, combined_Contig70086

####################################################
# Blast reference genome for nonvert contamination #
####################################################

# blast reference genome against NCBI nt dastabase
blastn -db nt -query combined_targetedRegionAndFlanking.fasta -outfmt 10 -num_alignments 5 -max_hsps 1 -out blast-output.m6

# get non-vert hits
GItaxidIsVert.py blast-output.m6 -n -a #Henderson, James B., Hanna, Zachary R. 2016. GItaxidIsVert. Version 1.0.0. DOI: 10.5281/zenodo.163737

# pulled potentially contaminated scaffolds with bioawk version 1.0 (Li, 2013)

# blasted pulled scaffolds
blastn -db nt -query potential-nonVert-contigs.fasta -html -out blast-output.html

# contigs w/ first or only blast hit as nonvert: combined_Contig118715, combined_Contig131378, combined_Contig136283, combined_Contig150008, combined_Contig30796, combined_Contig52365, combined_Contig58599, combined_Contig73564

##########################
# Alignment              #
##########################

mkdir /data/elinck/syma_assembly/pipeline_07212016/alignment

# align to extended reference genome
perl 5-Alignment -f /data/elinck/syma_assembly/pipeline_07212016/assemblies/raw/reference/In_target/Final/combined_targetedRegionAndFlanking.fasta -r /data/elinck/syma_assembly/pipeline_07212016/cleaned -o /data/elinck/syma_assembly/pipeline_07212016/alignment -i 235 -v 24 -l -t 90 -P /home/elinck/bin/picard-tools-2.4.1/picard.jar -G /home/elinck/bin/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -p 16

# align to extended reference genome
perl /data/elinck/syma_assembly/5-Alignment -f /data/elinck/syma_assembly/pipeline_07212016/alignment/cleaned/pseudo_ref_genome_cleaned.fa -r /data/elinck/syma_assembly/pipeline_07212016/cleaned -o /data/elinck/syma_assembly/pipeline_07212016/alignment/cleaned -i 235 -v 24 -l -t 90 -P /home/elinck/bin/picard-tools-2.4.1/picard.jar -G /home/elinck/bin/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -p 24

# evaluate exon capture performance (not working)
perl /data/elinck/syma_assembly/6-ExonCaptureEvaluation Evaluation -genome /data/elinck/syma_assembly/pipeline_07212016/assemblies/raw/reference/In_target/Final/combined_targetedRegionAndFlanking.fasta -cleanDir /data/elinck/syma_assembly/pipeline_07212016/cleaned -rawDir /data/elinck/syma_assembly/pipeline_07212016/raw/pre-clean -bamDir /data/elinck/syma_assembly/pipeline_07212016/alignment -InstrID HS -readLen 100 -resDir /data/elinck/syma_assembly/pipeline_07212016/evaluation_20160808 -bedFile /data/elinck/syma_assembly/pipeline_07212016/assemblies/raw/ancient/In_target/Final/combined_targetedRegionforExonCapEval.bed 1>/data/elinck/syma_assembly/20160808_eval1.log 2>/data/elinck/syma_assembly/20160808_eval1.err

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

# move bams to subfolders
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

# make historic vcf
samtools mpileup -B -D -I -S -uf /data/elinck/syma_assembly/pipeline_07212016/assemblies/raw/reference/In_target/Final/combined_targetedRegionAndFlanking.fasta /data/elinck/syma_assembly/pipeline_07212016/alignment/angsd/historic/*sorted.bam | bcftools view -cg - > raw_historic.vcf

# make modern vcf
samtools mpileup -B -D -I -S -uf /data/elinck/syma_assembly/pipeline_07212016/assemblies/raw/reference/In_target/Final/combined_targetedRegionAndFlanking.fasta /data/elinck/syma_assembly/pipeline_07212016/alignment/angsd/modern/*sorted.bam | bcftools view -cg - > raw_modern.vcf

# remove contaminant contigs from historic vcf
vcftools --vcf raw_historic_extended.vcf --not-chr combined_Contig2441 --not-chr combined_Contig70086 --not-chr combined_Contig118715 --not-chr combined_Contig131378 --not-chr combined_Contig136283 --not-chr combined_Contig150008 --not-chr combined_Contig30796 --not-chr combined_Contig52365 --not-chr combined_Contig58599 --not-chr combined_Contig73564 --remove-indels --recode --recode-INFO-all --out historic_nocontam_nomito

# remove contaminant contigs from modern vcf
vcftools --vcf raw_modern_extended.vcf --not-chr combined_Contig2441 --not-chr combined_Contig70086 --not-chr combined_Contig118715 --not-chr combined_Contig131378 --not-chr combined_Contig136283 --not-chr combined_Contig150008 --not-chr combined_Contig30796 --not-chr combined_Contig52365 --not-chr combined_Contig58599 --not-chr combined_Contig73564 --remove-indels --recode --recode-INFO-all --out modern_nocontam_nomito

# SNP cleaning historic
perl /data/elinck/syma_assembly/10-SNPcleaner.pl -A /data/elinck/syma_assembly/pipeline_07212016/assemblies/raw/reference/In_target/Final/combined_targetedRegionAndFlanking.fasta -M CT_GA -d 5 -k 7 -u 3 -a 0 -r /data/elinck/syma_assembly/pipeline_07212016/alignment/angsd/historic/s_toro_gene_depth_percentile.txt -B syma_torotoro_historic_cleaned.bed  -p syma_torotoro_dropped_historic -v /data/elinck/syma_assembly/pipeline_07212016/SNPs/historic_nocontam_nomito.recode.vcf > out_historic_nocontam.vcf

# SNP cleaning modern 
perl /data/elinck/syma_assembly/10-SNPcleaner.pl -A /data/elinck/syma_assembly/pipeline_07212016/assemblies/raw/reference/In_target/Final/combined_targetedRegionAndFlanking.fasta -M CT_GA -d 5 -k 3 -u 3 -a 0 -r /data/elinck/syma_assembly/pipeline_07212016/alignment/angsd/modern/s_toro_gene_depth_percentile.txt -B syma_torotoro_modern_cleaned.bed  -p syma_torotoro_dropped_modern -v /data/elinck/syma_assembly/pipeline_07212016/SNPs/modern_nocontam_nomito.recode.vcf > out_modern_nocontam.vcf

# report count of alleles for all bed files
wc -l /data/elinck/syma_assembly/pipeline_07212016/SNPs/*bed

# find shared sites among extended ref genome
bedtools intersect -a /data/elinck/syma_assembly/pipeline_07212016/SNPs/syma_torotoro_historic_cleaned.bed -b /data/elinck/syma_assembly/pipeline_07212016/SNPs/syma_torotoro_modern_cleaned.bed > /data/elinck/syma_assembly/pipeline_07212016/alignment/angsd/shared_sites_cleaned.bed

# make keep file for angsd
cut -f 1,2 /data/elinck/syma_assembly/pipeline_07212016/SNPs/shared_sites_cleaned.bed  > shared_sites_cleaned.keep

# index sites for angsd
angsd sites index shared_sites_cleaned.keep

# run ANGSD, full data
angsd -bam /data/elinck/syma_assembly/pipeline_07212016/SNPs/extended/bam.filelist.txt -sites /data/elinck/syma_assembly/pipeline_07212016/SNPs/shared_sites_cleaned.keep -anc /data/elinck/syma_assembly/pipeline_07212016/assemblies/raw/reference/In_target/Final/combined_targetedRegionAndFlanking.fasta -minQ 20 -fold 1 -out snps_cleaned_05 -GL 1  -doGeno 2  -doPost 1  -postCutoff 0.95 -doCounts 1 -doGlf 2 -geno_minDepth 6 -SNP_pval 0.05 -minMaf 0.05 -doMaf 2 -doMajorMinor 1 -doSaf 1 -doPlink 2

# run ANGSD, no MD
angsd -bam /data/elinck/syma_assembly/pipeline_07212016/SNPs/extended/bam.filelist.txt -sites /data/elinck/syma_assembly/pipeline_07212016/SNPs/shared_sites_cleaned.keep -anc /data/elinck/syma_assembly/pipeline_07212016/assemblies/raw/reference/In_target/Final/combined_targetedRegionAndFlanking.fasta -minQ 20 -fold 1 -out snps_cleaned_all_025  -GL 1  -doGeno 2  -doPost 1  -postCutoff 0.95 -doCounts 1 -doGlf 2 -geno_minDepth 6 -SNP_pval 0.05 -minMaf 0.025 -doMaf 2 -doMajorMinor 1 -doSaf 1 -doPlink 2 -minInd 19

mkdir SNPcounts

cd SNPcounts

# obtain threader script
git clone https://github.com/slager/threader.git

# open threadme file
nano threadme.txt

# to count SNPs per missing individual individual req, copy above command 19 times, tweaking -minInd flag from 1-19, paste in threadme.txt

# select 19 threads

python threader.py

gunzip *gz

# report number of loci
wc -l *geno 

cd ../

# output 100% complete matrix genind for adegenet
perl /data/elinck/syma_assembly/11-PopGenTools Adegenet -g snps_cleaned_all_05.geno -n 19 -s 1690 -o snps_cleaned_complete_maf05_adg

# output full data matrix genind for adegenet
perl /data/elinck/syma_assembly/11-PopGenTools Adegenet -g snps_cleaned.geno -n 19 -s 39105 -o snps_cleaned_all_adg

#############################################
# Downstream analyses (not reported in ms)  #
#############################################

# run NGSadmix on 75% matrix
for i in 1 2 3 4 5 6 7 8 9 10; do NGSadmix -likes snps_15.beagle -K $i -P 4 -o ngsadmix_${i}_15 -minMaf 0.0125; done

# run admixture on 75% matrix
for K in 1 2 3 4 5 6 7 8 9 10; do admixture --cv hyrad_admixture.bed $K | tee log${K}.out; done

# choose K with cross validation error
grep -h CV log*.out



