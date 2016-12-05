# Analysis of *Syma torotoro* (Aves: Halcyonidae) hyRAD data

Working repository for analysis of *S. torotoro* hyRAD data associated with the following manuscript in preparation:

Linck, E., Hanna, Z., Sellas, A., Dumbacher, J. In prep. Evaluating hybridization capture with RAD probes as a tool for museum genomics with historical bird specimens. 

### Bioinformatics

**Script:** syma_hyRAD_bioinformatics.sh  
**Description:** Shell commands for sequence data manipulation running QB3's [denovoTargetCapturePopGen wrapper scripts](https://github.com/CGRL-QB3-UCBerkeley/denovoTargetCapturePopGen), and running external bioinformatics tools (e.g., samtools)   
**Notes:** Don't try and execute this, obviously. Data stored at NCBI's SRA (### pending); Dryad (### pending). 

**Script:** cutcontigs.R  
**Description:** Simple R function for removing specific contigs from a fasta file. Used to remove contigs matching mitochondrial DNA sequences and non-vertebrate BLAST search matches (e.g., contaminant DNA) from pseudo-reference genome.  
**Notes:** Credit to [R. Harris] (https://github.com/rebzzy/) for the lapply() solution.  

### Population genetic analyses

**Script:** syma_pca.R  
**Description:** PCA / DAPC analysis of 100% complete SNP matrix implemented in adegenet  
**Notes:** Requires dapcplot.R. Data: elinck/hyRAD/data/snps_extended_adg  

### Statistical analysis

**Script:** hyrad_performance.R  
**Description:** Basic statistical analyses of hyRAD's sequencing and assembly performance as relates to various input DNA variables.   
**Notes:** Data: elinck/hyRAD/data/syma_torotoro_hyRAD_performance.csv

### Plotting tools

**Script:** dapcplot.R  
**Description:** An alternative to adegenet's "compoplot()" function, for producing STRUCTURE-like plots with DAPC population membership probabilities. Extends compoplots functionality by allowing consistent(ish) color assignment across individuals when increasing values of K, and allows sorting of individuals (bars) by population assignment.   
**Notes:** Adapted from C.J. Battey's [similar function for STRUCTURE plots] (https://github.com/cjbattey/RADplots/blob/master/structurePlot.R)

