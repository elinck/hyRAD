### quick plot of site frequency spectrum from angsd

#library(ggplot2);library(plyr)
#setwd("~/Dropbox/Syma/hyRAD/")

derived <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)
numloci <- c(1.005061,0.000000,1488.243402,3917.796615,1205.680385,774.548338,973.329280,595.955838,427.955476,502.444814,349.693775,489.239975,410.210897,218.298125,422.634176,350.344884,199.571866,437.537579,255.943029,580.566484)
proportion <- numloci/sum(numloci)
sfs <- data.frame(derived,numloci,proportion)
ggplot(sfs, aes(x=sfs$derived,y=sfs$proportion)) +
  geom_bar(stat="identity") +
  xlab("Derived Allele Frequency") +
  ylab("Proportion of SNPs")

