# hyRAD Bioinformatics Tutorial  

This tutorial describes a bioinformatics pipeline for cleaning, assembling reads, and performaing variant discovery with hyRAD data, starting with raw Illumina fastq reads. Our pipeline is primarily based on Derek Eaton's pyRAD pipeline and Berkeley QB3's denovoTargetCapturePopGen wrapper scripts for handling target capture data.     

## Getting Started

First, you'll need [Python](https://www.python.org/) and [Perl](https://www.perl.org/) installed; these languages form the basis of the Pyrad and QB3 denovoTargetCapturePopGen wrapper scripts we'll rely on heavily throughout the pipeline. Additionally, you'll need to separately install a variety of genomic software tools for individual steps, indicated below.  

Then, establish your working directory:

```{sh}
mkdir /home/user/project
mkdir /home/user/project/raw_reads #FTP or SCP sequence data into this directory
```
...clone QB3's denovoTargetCapturePopGen perl wrapper scripts...

```{sh}
cd /home/user/project
git clone https://github.com/CGRL-QB3-UCBerkeley/denovoTargetCapturePopGen.git 
```
...clone QB3's denovoTargetCapturePhylogenomics perl wrapper scripts (you'll need these in a few cases, too)...

```{sh}
git clone https://github.com/CGRL-QB3-UCBerkeley/denovoTargetCapturePhylogenomics.git
```

...count the number of raw reads and save it as file "hiseq_readcounts.txt"...

```{sh}
cd /home/user/project/raw_reads
for i in *; do echo $i; gunzip -c $i | echo $((`wc -l`/4)); done > hiseq_readcounts.txt
```

...and install pyRAD.

```{sh}
mkdir /home/user/project/probes
cd /home/user/project/probes
git clone https://github.com/dereneaton/pyrad.git #clones latest development version
```

## Assemble Probe Library / Psuedo-Reference Genome

To create a psuedo-reference genome, you'll want to run pyrad steps 1 through 5 on raw reads from your sequenced RAD (probe) library. This requires a "params.txt" file set to user specifications, which is explained in detail over at [Deren Eaton's website](http://dereneaton.com/software/pyrad/). We suggest setting a clustering threshold of 0.97 (## 10.), and keeping trimmed reads with a minimum length of 70 (## 32.opt). 

```{sh}
cd /home/user/project/probes
nano syma_params.txt #edit with preferred specs
python /home/user/project/probes/pyRAD.py -p /home/user/project/probes/syma_params.txt -s 12345
```

After creating a raw this, you'll need to remove repetitive genome regions by creating a BLAST database from your consensus fasta and then using the same original fasta as the query sequence to search for alignments. We can then use unix tools to identify duplicates, remove them, and save only unique contigs. Note the e-value and BLAST out format tags; you may want to change these based on your system / preferences.

```{sh}
makeblastdb -in /home/user/project/probes/outfiles/samplename.consens -parse_seqids -dbtype nucl
blastn -db /home/user/project/probes/outfiles/samplename.consens -query /home/user/project/probes/clust.97/samplename.consens -evalue 0.00001 -outfmt 6 -out samplename_selfblast.txt
cut -f 1,1 samplename_selfblast.txt | sort | uniq -u | wc -l #counts unique BLAST hits
cut -f 1,1 samplename_selfblast.txt | sort | uniq -u > rad_probe_unique.txt #saves unique hits
blastdbcmd -db /home/user/project/probes/clust.97/samplename.consens -entry_batch /home/user/project/probes/clust.97/rad_probe_unique.txt > rad_probe_unique.fasta
```

## Cleaning hyRAD library Reads

Now, we'll turn to the captured (or hyRAD) libraries. First, you'll need to tun QB3's 1-ScrubReads wrapper script. This step cleans up raw data, which includes trimming for quality, removing adapters, merging overlapping reads, removing duplicates and reads sourced from contamination. The flag -c should point to a fasta file with reference genome for any putative contaminants, e.g. *E. coli*. It also requires [Trimmomatic](www.usadellab.org/cms/?page=trimmomatic), which you should now install. (The flag "-t" indicates a path to the trimmomatic .jar executable). 

```{sh}
perl /home/user/project/1-ScrubReads_a cleanPE -t /home/user/bin/Trimmomatic0.36/trimmomatic-0.36.jar -f /home/user/project/raw_reads/ -k 12 -z -c /home/user/project/e_coli_K12.fasta -o /home/user/project/cleaned/ 1>/home/user/project/cleaned/scrubreads_run1.err 2>/home/user/project/cleaned/scrubreads_run1.log
```

To determine the number of duplicate reads in your libraries, subtract the number of scrubbed reads in the "preclean" folder from the original. Expect >50%.

```{sh}
cd home/user/project/raw_reads
cat *.fq | echo $((`wc -l`/4)) > reads_original.txt
cd home/user/project/cleaned/preclean
cat *.fq | echo $((`wc -l`/4)) > reads_precelean.txt # 1 - (reads_preclean / sample count in reads_original) = % duplicates
```

## Generating hyRAD library assemblies

Once you've done a preliminary clean of your reads, it's time to generate assemblies (clusters). The QB3 donovoTargetCapture pipeline provides support for several different assemblers, discussed in detail over at their [wiki](https://cgrlucb.wikispaces.com/Bioinformatics+Pipeline+for+de+novo+Targeted+Capture). We chose to use [Spades](http://bioinf.spbau.ru/spades), which automatically selects an optimal K-mer (substring) value with minimal futzing. (This is only an option with the wrapper 2-GenerateAssembliesPhylo from the "Phylogenomics" package.) Note the flag "-np", which determines how many threads to use. We chose 16, as we were running analyses on a cluster via a websever, but you'll obviously want to work within the bounds of the computer. (This step can be time consuming.)    

```{sh}
mkdir assemblies
perl home/user/project/2-GenerateAssembliesPhylo spades -reads home/user/project/cleaned -out /home/user/project/assemblies -np 16 1 > assembly_log1.txt 2 > assembly_error1.txt
```

## Finding targets

After generating your assemblies, it's time to determine which represent "on target" (or captured) loci in each sample. If you have both modern and historic samples (e.g., from museum specimens or other "ancient" DNA source), you'll want to run this step separately for each, as flanking regions to probe sequences will be merged to extend the pseudo-reference genome and only modern samples (presumably contaminant free) should be included. 

* The "-p" flag selects a clustering threshhold -- we reccomend 0.95
* The -b" flag should be set to 1 to merge modern DNA assemblies into probe sequences to extend the pseudo-reference genome. (We'll use this downstream for SNP calling.) 
* The "-e" option should be set to 4 to note were using "random" probe sequences, not exons or other known genomic regions. Note we're again using a script (3-FindingTargetsV8) from the "Phylogenomics" wrapper package, with the "combineExon" sub.

```{sh}
cd assemblies
mkdir reference
mkdir ancient
perl home/user/project/3-FindingTargetsV8 combineExon -t /home/user/project/probes/clust.97/rad_probe_unique.fasta -a /home/user/project/assemblies/reference -p 0.95 -b 1 -e 4 
perl home/user/project/3-FindingTargetsV8 combineExon -t /home/user/project/probes/clust.97rad_probe_unique.fasta -a /home/user/project/assemblies/ancient -p 0.95 -b 1 -e 4
```

If you're curious about how many captured sequences you have in each sample -- which you should be -- you can print them with the following command:

```{sh}
grep ">" -c *.fa
```

## Evaluating and cleaning assemblies

After determining which assemblies are on target, you'll want to filter for accidental mitochondrial DNA assemblies. You can do this with BLAST, a fasta file of a full mitochondrial genome of your organism (or a close relative; "organism.fast" here), and the R script [cutcontigsbatch.R](https://github.com/elinck/hyRAD/blob/master/cutcontigsbatch.R). Additionally, it's a good time to install [Samtools](http://www.htslib.org/).

First, make your organim.fasta mtDNA file a blast database.

```{sh}
makeblastdb -in /home/user/project/assemblies/organism.fasta -parse_seqids -dbtype nucl
```

Then, run a loop in each of your sets of samples (e.g., modern and ancient DNA) to query your assemblies for alignments:

```{sh}
for i in /home/user/project/assemblies/reference/fasta/*.fasta; do blastn -db /home/user/project/assemblies/organism.fasta -query $i -evalue 0.00001 -outfmt 6 -out ${i}_mtDNAhits.txt; done
for i in /home/user/project/assemblies/ancient/fasta/*.fasta; do blastn -db /home/user/project/organism.fasta -query $i -evalue 0.00001 -outfmt 6 -out ${i}_mtDNAhits.txt; done
```

You may need to eliminate line breaks in your fasta file prior to using the cutcontigsbatch.R scripts. To do so:

```{sh}
for i in /home/user/project/assemblies/reference/fasta/intarget/.*fasta; do awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' $i > ${i}_oneline; done
```

You can then use [extractcontigID.R](https://github.com/elinck/hyRAD/blob/master/extractcontigIDs.R) to collect the contig ("chromosome") ID of each potential mtDNA contaminant sequence across all samples, and then use and [cutcontigsbatch.R](https://github.com/elinck/hyRAD/blob/master/extractcontigIDs.R) to export fasta files filtered for these sequences. (See documentation in R scripts for details.)  

Once you have filtered, final versions of your assemblies, you probably will want to generate some basic summary statistics, using the QB3 wrapper AssemblyEvaluation with option BASIC. Note the paths here are not necessarily what you'll use: the assembly wrapper will generate a confusing number of different versions of each assemblies, and it's important to investigate each one and determine which is the most "final." 

```{sh}
perl /data/elinck/syma_assembly/AssemblyEvaluation BASIC -a /home/user/project/assemblies/reference/In_target/premasked/test/
perl /data/elinck/syma_assembly/AssemblyEvaluation BASIC -a /home/user/project/assemblies/ancient/In_target/premasked/test/
```

It's possible the AssemblyEvaluation script will throw an error if fasta file names are incorrect. If so, this line of code will ensure each file (sampleID*) ends with a clean ".fasta":

```{sh}
find -type f -name "sampleID*" -exec rename 's/.fasta.original/.fasta/' \{\} \;
```

## BLAST pseudo-reference genome for contamination

At this point, you have cleaned and filtered assemblies for all your samples, and a pseudo-reference genome composed of unique original probe sequences and their flanking regions. If you're working at deeper timescales, you can take these consensus assemblies and align them (possibly with [MAFFT](http://mafft.cbrc.jp/alignment/software/) for export for downstream phylogenetic analyses. But if you're working with population level samples, as we are here, you'll need to return to align raw reads to the psuedo-reference genome for phasing of alleles and variant (SNP) discovery. And while you've presumably filtered contaminant reads from all samples, and mtDNA content from consensus *assemblies*, you still need to make sure the pseudo-reference genome contains no contaminant sequence. We'll use the BLAST database of your organism's full mtDNA genome and your full pseudo-reference genome as a query: 

```{sh}
blastn -db t_sanctus_ref.fasta -query combined_targetedRegionAndFlanking.fasta -evalue 0.00001 -outfmt 6 -out mtDNA_in_reference_blast.txt
```

Record any matching contig IDs. Next -- if you are working with vertebrates and if you are working with ancient DNA or have another reason to suspect exogenous microbial contamination -- you can additionally filter the pseudo-reference genome for any contigs that primarily align to non-vertebrate taxa. You'll need a local copy of the [NCBI Genbank nucleotide database](https://blast.ncbi.nlm.nih.gov/Blast.cgi), and a copy of Zach Hanna and Jim Henderson's [GItaxidIsvert program](https://github.com/calacademy-research/GItaxidIsVert). 

```{sh}
git clone https://github.com/calacademy-research/GItaxidIsVert.git
```

Use your pseudo-reference genome (should be named "combined_targetedRegionAndFlanking.fasta" now) as a query in a search for alignments against the full NCBI database"

```{sh}
blastn -db nt -query combined_targetedRegionAndFlanking.fasta -outfmt 10 -num_alignments 5 -max_hsps 1 -out blast-output.m6
```

Then, use the GItacidIsVert.py script to identify non-vertebrate sequences. (See [that program's documentation](https://github.com/calacademy-research/GItaxidIsVert) for details.)

```{sh}
python GItaxidIsVert.py blast-output.m6 -n -a #Henderson, James B., Hanna, Zachary R. 2016. GItaxidIsVert. Version 1.0.0. DOI: 10.5281/zenodo.163737
```

Pull putative contaminant contigs with whatever method you see fit -- we used bioawk version 1.0 (Li, 2013). Finally, take these contigs and perform a second BLAST search against the full nucleotide database, exporting output as an .html.

```{sh}
blastn -db nt -query potential-nonVert-contigs.fasta -html -out blast-output.html
```

Within this html, identify contigs with their first / only alignment to non-vertebrate taxa, and add their IDs to the list that already includes putative mtDNA contamination. We'll return to this during SNP calling. 

## Alignment

Now that you are reasonably sure your pseudo-reference genome is free of contamination and repetitive regions, it's time to align your cleaned reads against it for allele phasing and SNP calling with QB3's 5-Alignment wrapper script. You'll need to do this for each folder of samples, e.g., both modern and ancient DNA. For this step, you'll need both [Picard Tools](https://broadinstitute.github.io/picard/) and the [GATK](https://software.broadinstitute.org/gatk/) installed in your bin, as well as [NovoAlign](http://www.novocraft.com/products/novoalign/). There are also several flags that must be set to your unique specifications: 
* "-i" indicates average insert size 
* "-v" indicates the standard deviation of insert size 
* "-t" indicates maximum alignment score (should be 90 for pop gen studies) 
* "-p" indicates number of threads to call
* "-P" and "-G" are paths to the Picard and GATK executables, respectively. 

```{sh}
mkdir /home/user/bin/alignment/
perl 5-Alignment -f /home/user/project/assemblies/reference/In_target/Final/combined_targetedRegionAndFlanking.fasta -r /home/user/project/raw_reads/cleaned -o /home/user/project/alignment -i 235 -v 24 -t 90 -P /home/user/bin/picard-tools-2.4.1/picard.jar -G /home/user/bin/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -p 16
```

You can then determine the specificity, or percentage or reads that mapped to the reference for each sample...

```{sh}
for bam_file in *.bam
do
total=$(samtools view -c $bam_file)
mapped=$(samtools view -c -F 4 $bam_file)
unmapped=$(samtools view -c -f 4 $bam_file)
echo "$bam_file $total $mapped $unmapped"
done
```
...as well as the average depth per sample:

```{sh}
for bam_file in *.bam
do
average_depth=$(samtools depth $bam_file |  awk '{sum+=$3} END { print "Average = ",sum/NR}')
echo "$bam_file $average_depth"
done
```

## Data filtering

Before moving on to SNP calling, we need to perform some preliminary filtering steps on our alignment, which we'll do separately for each sample type. Begin by setting up subdirectories:

```{sh}
mkdir /home/user/project/SNPs
mkdir /home/user/project/SNPs/reference
mkdir /home/user/project/SNPs/ancient
```

Move all .bam files from the alignment step into these folders by sample type: 

```{sh}
mv sampleID1* sampleID2* ../SNPs/reference
mv sampleID3* sampleID4* ../SNPs/ancient
```

Next, you need to create merged .bam files for each folder of samples:

```{sh}
samtools merge merge.bam /home/user/project/SNPs/reference/*.bam
samtools merge merge.bam /home/user/project/SNPs/ancient/*.bam
```
Create a merged file for each, and move this to its own subdirectory. (Getting it out of the way of files with the same extension will make generating a .vcf file easier downstream.)

```{sh}
mkdir /home/user/project/SNPs/reference/merged/
mv merge.bam merged
samtools sort /home/user/project/SNPs/reference/merged/merge.bam merge_sorted
mkdir /home/user/project/SNPs/ancient/merged/
mv merge.bam merged
samtools sort /home/user/project/SNPs/ancient/merged/merge.bam merge_sorted
```
We'll then input these merged .bam files to QB3's 9-Prefiltering wrapper with both its "percentile" option. This will generate empirical site and gene coverage distributions, which will then be used for variant site filtering. 

```{sh}
perl 9-preFiltering percentile -b /home/user/project/SNPs/reference/merged/merge_sorted.bam -o samples_reference
perl 9-preFiltering percentile -b /home/user/project/SNPs/ancient/merged/merge_sorted.bam -o samples_modern
```

We also need to generate a .bed file, indicating the alignment positions of our pseudo-reference genome relative to the original target. The relevant output will be labeled "All_contigs.bed," and will be useful to us later on.

```{sh}
perl 9-preFiltering bed /home/user/project/probes/clust.97/rad_probe_unique.txt /home/user/project/assemblies/reference/In_target/Final/combined_targetedRegionAndFlanking.fasta
```

Next, we'll create a "raw" .vcf file for each folder, using samtools mpileup function. ("Raw" here meaning retaining all variant sites, regardless of quality.) Note we are using .bam files for individual samples, not the merged file, which we indicate with the regex "*sorted.bam." (This is why we placed the merged file in a seperate subdirectory above.)

```{sh}
samtools mpileup -B -D -I -S -uf /home/user/project/alignment/reference/In_target/Final/combined_targetedRegionAndFlanking.fasta /home/user/project/SNPs/reference/*sorted.bam | bcftools view -cg - > raw_reference.vcf
samtools mpileup -B -D -I -S -uf /home/user/project/alignment/reference/In_target/Final/combined_targetedRegionAndFlanking.fasta /home/user/project/SNPs/ancient/*sorted.bam | bcftools view -cg - > raw_ancient.vcf
```

Now, it's time to remove the contigs you identified as mtDNA or microbial contamination in the pseudo-reference genome, using [vcftools](http://vcftools.sourceforge.net/) --not-chr flag. (You'll need to install this program if you don't have it already.) By waiting until this stage to remove contaminant sites, we're reducing the probability of mismatches during alignment. 

```{sh}
vcftools --vcf raw_reference.vcf --not-chr combined_Contig2441 --not-chr combined_Contig70086 --remove-indels --recode --recode-INFO-all --out reference_nocontam
vcftools --vcf raw_ancient.vcf --not-chr combined_Contig2441 --not-chr combined_Contig70086 --remove-indels --recode --recode-INFO-all --out ancient_nocontam
```

We can then use QB3's SNPcleaner script to perform several additional filtering steps. These include:

1. Removing contigs that show extremely low or high coverage based on the empirical coverage distribution across all contigs
2. Removing contigs with at least one SNP having allele frequencies highly deviating from Hardy–Weinberg equilibrium expectations 
3. Removing sites with excessively low or high coverage based on the empirical coverage distribution
4. Removing sites having allele frequencies highly deviating from Hardy–Weinberg equilibrium expectations
5. Removing sites with biases associated with reference and alternative allele Phred quality, mapping quality and distance of alleles from the ends of reads, and removing sites that show a bias towards SNPs coming from the forward or reverse strand
6. Removing sites for which there are not at least M of the individuals sequenced at N coverage each
7. Removing sites with a root mean square (RMS) mapping quality for SNPs across all samples below a certain threshold
8. Removing C to T and G to A SNPs from the dataset if working with ancient DNA with a high level of deamination / base misincorporation.

As always, you'll want to tweak these commands for your specific needs. 
* The "-M" flag here indicates which mutations to remove (C to T and G to A for our purposes, which is indicated CT_GA). 
* The "-d" flag indicates minimum site depth
* The "-k" flag indicates the minimum number of individuals to include with less that "-u" converage
* The "-a" flag  indicates the minimum number of alternate alleles.
* The "-v" glag indicates your contaminant free .vcf; 
* The "-p" flag indicates an (optional) file to log dropped sites with "-p";
* The "-X" flag inputs your "All_contigs.bed" file 
* The "-B" flag indicates a name for the cleaned .bed file for each sample type.

```{sh}
perl 10-SNPcleaner.pl -A /home/user/project/assemblies/reference/In_target/Final/combined_targetedRegionAndFlanking.fasta -M CT_GA -d 5 -k 7 -u 3 -a 0 -r /home/user/project/SNPs/reference/reference_gene_depth_percentile.txt -B samples_reference.bed  -p dropped_reference -X All_contigs.bed -v /home/user/project/SNPs/reference/reference_nocontam.vcf > reference_cleaned.vcf
perl 10-SNPcleaner.pl -A /home/user/project/assemblies/reference/In_target/Final/combined_targetedRegionAndFlanking.fasta -M CT_GA -d 5 -k 3 -u 3 -a 0 -r /home/user/project/SNPs/ancient/ancient_gene_depth_percentile.txt -B samples_ancient.bed  -p dropped_ancient -X All_contigs.bed -v /home/user/project/SNPs/ancient/ancient_nocontam.vcf > ancient_cleaned.vcf
```

You can report the number of sites present in each .bed file with this command:

```{sh}
wc -l /home/user/project/SNPs/*/*bed
```

Next, we'll restrict our final matrix to only those sites passing filters for BOTH sample types (again, modern and ancient DNA in our example). We'll do this using [bedtools](http://bedtools.readthedocs.io/en/latest/) "intersect" function. (You'll need to install this program).

```{sh}
bedtools intersect -a /home/user/project/SNPs/reference/samples_reference.bed -b /home/user/project/SNPs/reference/samples_ancient.bed > /home/user/project/SNPs/shared_sites.bed
```

Lastly, we need to transform this file (shared_sites.bed) into a "keep" file for [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD), the program we'll use for SNP calling, and then use ANGSD itself to "index" this file. (Install the program now.)

```{sh}
cut -f 1,2 /home/user/project/SNPs/shared_sites.bed > /home/user/project/SNPs/shared_sites.keep
angsd sites index shared_sites.keep
```

## SNP calling

Your data should now be ready for SNP calling. As mentioned above, we'll be using ANGSD, which is well-suited to low-coverage sequence data as it calls SNPs and estimates allele frequencies using an empirical Bayesian framework. We recommend reading the program's extensive documentation [here](http://www.popgen.dk/angsd/index.php/ANGSD). A quick-start guide follows. 
* You'll need a list of the full paths of the .bam files of all samples you plan on including in your SNP calling in a text file (one path per line), named "bam.filelist.txt", indicated with the flag "-bam." 
* You'll also need your .keep file ("-sites"), and an "ancestral" sequence (use your pseudo-reference genome if you don't have one). 
* The "-minQ" flag sets a minimum acceptable Phred-scaled base quality 
* "-fold 1" makes the program estimate the folded site frequency spectrum and then uses your supplied reference as ancestral the ancestral sequence
* "-out" indicates outfile prefixes
* "-doGeno" selects a format for writing genotypes from ANGSD's list of options (see documentation for details; "2" will work for downstream adegent analyses)
* "-doPost 1" indicates the posterior genotype probability will be estimated with allele frequency as a prior
* "-postCutoff" uses only genotypes with a posterior probability above the given value
* "-SNP_pval" sets the P-value for a site being variable required to call a SNP
* "-doCounts 1" calculates the frequency of different bases
* "-genoMinDepth" indicates the minimum per site depth
* "-doMajorMinor" infers the major and minor alleles "-doMaf 1" calculates minor allele frequency, while "-minMaf" selects a minumum minor allele frequency
* "-doSaf" will calculate the likelihood of the sample allele frequency for each site and dump these into a ".saf" file
* "-minInd" sets the minimum number of individuals a site must be present in to be called as a SNP.  

```{sh}
angsd -bam /home/user/project/SNPs/bam.filelist.txt -sites /home/user/project/SNPs/shared_sites.keep -anc /home/user/project/assemblies/reference/In_target/Final/combined_targetedRegionAndFlanking.fasta -minQ 20 -fold 1 -out snps_cleaned -GL 1 -doGeno 2  -doPost 1  -postCutoff 0.95 -doCounts 1 -geno_minDepth 6 -SNP_pval 0.05 -minMaf 0.05 -doMaf 2 -doMajorMinor 1 -doSaf 1 -minInd 10
```
Based on your choice for "-doGeno," genotypes of SNPs passing all filters will be printed in any number of ways that will be difficult to use outside of ANGSD's sister program, [NGSadmix](http://www.popgen.dk/software/index.php/NgsAdmix). As you're likely to want to perform standard population genetic or phylogenetic analyses, we can again use a QB3 wrapper script to transform this output into the input for a number of popular programs. Described at length in the useage description for [11-PopGenTools](https://github.com/CGRL-QB3-UCBerkeley/denovoTargetCapturePopGen/blob/master/11-PopGenTools), the following command generates a .genind file for the R package adegenet, using a .geno file from angsd with 20 individuals, 2000 sites, and a custom output prefix.

```{sh}
perl 11-PopGenTools Adegenet -g snps_cleaned.geno -n 20 -s 2000 -o snps_cleaned_adg
```
Other conversions are similarly simple.

## Authorship
[Ethan Linck](https://github.com/elinck/) and [Zach Hanna](https://github.com/calacademy-research) cobbled together this pipeline from a number of different tools. The QB3 denovoTargetCapture wrapper scripts were written by Ke Bi, Sonal Singhal, and Tyler Linderoth. Pyrad was written by Deren Eaton. Ethan Linck wrote this tutorial and converted it to markdown. 





















