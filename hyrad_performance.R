setwd("~/Dropbox/Syma/hyRAD/")
#install.packages("gridExtra"); install.packages("grid");install.packages("ggpmisc");install.packages("lme4");install.packages("AICcmodavg")
library(ggplot2);library(ggpmisc);library(plyr);library(dplyr);library(reshape2);library(gridExtra);library(grid);library(lme4);library(AICcmodavg)

### analyses exploring hyRAD's sequencing / assembly performance by various input DNA variables

#read in assembly performance data, subset by extract type
data <- read.csv("./data/syma_torotoro_hyRAD_performance.csv")
modern <- subset(data, extract_type=="modern")
historic <- subset(data, extract_type=="historic")

### explore differences between extract types ###

# calculate means and standard deviations by extract type
mpd <- ddply(data, c("extract_type"), summarise,
             mean = mean(percent_duplicates),
             sd = sd(percent_duplicates)
)
mnc <- ddply(data, c("extract_type"), summarise,
             mean = mean(num_loci_on_target),
             sd = sd(num_loci_on_target)
)
mcov <- ddply(data, c("extract_type"), summarise,
               mean = mean(eval_coverage),
               sd = sd(eval_coverage)
)
mmean <- ddply(data, c("extract_type"), summarise,
               mean = mean(mean),
               sd = sd(mean)
)
mgc <- ddply(data, c("extract_type"), summarise,
             mean = mean(percent_gc_extended_ref),
             sd = sd(percent_gc_extended_ref)
)
msp <- ddply(data, c("extract_type"), summarise,
             mean = mean(specificity),
             sd = sd(specificity)
)
msen <- ddply(data, c("extract_type"), summarise,
           mean = mean(sensitivity),
           sd = sd(sensitivity)
)
menr <- ddply(data, c("extract_type"), summarise,
          mean = mean(enrichment),
          sd = sd(enrichment)
)

#generate plots for means / sds
p1 <- ggplot(mpd, aes(x = factor(extract_type), y = mean)) + 
  geom_bar(position=position_dodge(), width=0.5, stat="identity", aes(fill = factor(extract_type))) +
  coord_cartesian(ylim=c(0.50,0.95)) +
  guides(fill=guide_legend(title=NULL)) +
  theme(legend.text=element_text(size=12)) +
  scale_fill_discrete(breaks=c("historic","modern"),labels=c("Historic","Modern")) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1) +
  xlab("") +
  ylab("Mean percent duplicate reads") +
  theme(axis.title.y = element_text(size=12)) +
  theme(axis.text.x = element_blank())

p2 <- ggplot(mnc, aes(x = factor(extract_type), y = mean)) + 
  geom_bar(position=position_dodge(), width=0.5, stat="identity", aes(fill = factor(extract_type))) +  
  guides(fill=guide_legend(title=NULL)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1) +
  xlab("") +
  ylab("Mean number of on-target contigs") +
  theme(axis.title.y = element_text(size=12)) +
  theme(axis.text.x = element_blank())

p3 <- ggplot(mcov, aes(x = factor(extract_type), y = mean)) + 
  geom_bar(position=position_dodge(), width=0.5, stat="identity", aes(fill = factor(extract_type))) +
  coord_cartesian(ylim=c(0,16)) +
  guides(fill=guide_legend(title=NULL)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1) +
  xlab("") +
  ylab("Mean coverage") +
  theme(axis.title.y = element_text(size=12)) +
  theme(axis.text.x = element_blank())

p4 <- ggplot(mmean, aes(x = factor(extract_type), y = mean)) + 
  geom_bar(position=position_dodge(), width=0.5, stat="identity", aes(fill = factor(extract_type))) +
  guides(fill=guide_legend(title=NULL)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1) + 
  xlab("") +
  ylab("Mean contig size") +
  theme(axis.title.y = element_text(size=12)) +
  theme(axis.text.x = element_blank())

p5 <- ggplot(mgc, aes(x = factor(extract_type), y = mean)) + 
  geom_bar(position=position_dodge(), width=0.5, stat="identity", aes(fill = factor(extract_type))) +
  coord_cartesian(ylim=c(40,55)) +
  guides(fill=guide_legend(title=NULL)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1) + 
  xlab("") +
  ylab("Mean percent GC content") +
  theme(axis.title.y = element_text(size=12)) +
  theme(axis.text.x = element_blank())

p6 <- ggplot(msp, aes(x = factor(extract_type), y = mean)) + 
  geom_bar(position=position_dodge(), width=0.5, stat="identity", aes(fill = factor(extract_type))) +
  coord_cartesian(ylim=c(0,25)) +
  guides(fill=guide_legend(title=NULL)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1) + 
  xlab("") +
  ylab("Mean specificity") +
  theme(axis.title.y = element_text(size=12)) +
  theme(axis.text.x = element_blank())

p7 <- ggplot(msen, aes(x = factor(extract_type), y = mean)) + 
  geom_bar(position=position_dodge(), width=0.5, stat="identity", aes(fill = factor(extract_type))) +
  coord_cartesian(ylim=c(0,100)) +
  guides(fill=guide_legend(title=NULL)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1) + 
  xlab("") +
  ylab("Mean sensitivity") +
  theme(axis.title.y = element_text(size=12)) +
  theme(axis.text.x = element_blank())

p8 <- ggplot(menr, aes(x = factor(extract_type), y = mean)) + 
  geom_bar(position=position_dodge(), width=0.5, stat="identity", aes(fill = factor(extract_type))) +
  coord_cartesian(ylim=c(0,20)) +
  guides(fill=guide_legend(title=NULL)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1) + 
  xlab("") +
  ylab("Mean fold enrichment") +
  theme(axis.title.y = element_text(size=12)) +
  theme(axis.text.x = element_blank())

#load grid / shared legend function
grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  
  plots <- list(p1,p4,p5,p2,p3,p6,p7,p8)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  grid.newpage()
  grid.draw(combined)
  
}

# plot all together w/ shared legend
grid_arrange_shared_legend(ncol=4,nrow=2,p1,p4,p5,p2,p3,p6,p7,p8)

### test for significant differences between historic / modern samples ###

t.test(historic$read_count_original, modern$read_count_original)
t.test(historic$read_count_preclean, modern$read_count_preclean)
t.test(historic$coverage, modern$coverage)
t.test(historic$on_target_loci_1kb, modern$on_target_loci_1kb)
t.test(historic$percent_duplicates, modern$percent_duplicates)
t.test(historic$num_loci_on_target, modern$num_loci_on_target)
t.test(historic$specificity_extended_ref, modern$specificity_extended_ref)
t.test(historic$mean, modern$mean)
t.test(historic$percent_gc_extended_ref, modern$percent_gc_extended_ref)
t.test(historic$missing_snps, modern$missing_snps)
t.test(historic$specificity, modern$specificity)
t.test(historic$sensitivity, modern$sensitivity)
t.test(historic$enrichment, modern$enrichment)
t.test(historic$eval_coverage, modern$eval_coverage)

### explore relationship between num_loci and multiple variables ###

# ended up dropping this analysis -- too few samples
# modern samples, drivers of number of loci on target
# only include interactions w/ plausible impact
# mnl1 <- lm(num_loci_on_target ~ read_count_preclean + ng_uL_qubit + read_count_preclean*ng_uL_qubit, modern)
# summary(mnl1) 

# mnl2 <- lm(num_loci_on_target ~ read_count_preclean + ng_uL_qubit, modern)
# summary(mnl2) 

# mnl3 <- lm(num_loci_on_target ~ read_count_preclean, modern)
# summary(mnl3) 

# AICc(mnl1)
# AICc(mnl2)
# AICc(mnl3) # selects mln3, significant

# historic samples, drivers of number of loci on target
hnl1 <- lm(num_loci_on_target ~ age + input_quantity + read_count_preclean, historic)
summary(hnl1) 

hnl2 <- lm(num_loci_on_target ~ input_quantity + read_count_preclean, historic)
summary(hnl2)

hnl3 <- lm(num_loci_on_target ~ read_count_preclean, historic)
summary(hnl3)

# corrected AIC (small sample size) to choose model
AICc(hnl1)
AICc(hnl2)
AICc(hnl3)

# plot sig. relationship between concentation, loci on target
p6 <- ggplot(historic, aes(x=read_count_preclean, color=extract_type, y=num_loci_on_target)) + 
  guides(color=FALSE) +
  geom_point() + 
  xlab("Read count") +
  ylab("Number of captured loci") +
  theme(axis.title.y = element_text(size=12)) +
  theme(axis.title.x = element_text(size=12)) +
  stat_smooth(method="lm", fullrange=TRUE, alpha=0.1) 

# historic samples, mean locus length
hmll1 <- lm(mean ~ age + ng_uL_qubit + read_count_preclean, historic)
summary(hmll1) 

hmll2 <- lm(mean ~ age + read_count_preclean, historic)
summary(hmll2) 

hmll3 <- lm(mean ~ read_count_preclean, historic)
summary(hmll3) 

#historical samples, gc content
hgc1 <- lm(percent_gc_extended_ref ~ age + input_quantity + read_count_preclean, historic)
summary(hgc1)

hgc2 <- lm(percent_gc_extended_ref ~ age + input_quantity, historic)
summary(hgc2)

hgc3 <- lm(percent_gc_extended_ref ~ input_quantity, historic)
summary(hgc3)

#historical samples, specificity
hsp1 <- lm(specificity ~ age + input_quantity + read_count_preclean, historic)
summary(hgc1)

hsp2 <- lm(specificity ~ age + input_quantity, historic)
summary(hgc2)

hsp3 <- lm(specificity ~ input_quantity, historic)
summary(hgc3) # significant, p=0.0161

#historical samples, sensitivity

hse1 <- lm(sensitivity ~ age + input_quantity + read_count_preclean, historic)
summary(hse1)

hse2 <- lm(specificity ~ age + read_count_preclean, historic)
summary(hse2)

hse3 <- lm(specificity ~ age, historic)
summary(hse3) #not significant

#historical samples, enrichment

hen1 <- lm(enrichment ~ age + input_quantity + read_count_preclean, historic)
summary(hen1)

hen2 <- lm(enrichment ~ age + read_count_preclean, historic)
summary(hen2)

hen3 <- lm(enrichment ~ age, historic)
summary(hen3)

# plot sig. relationship between concentration, gc content in pseudo ref genome
p7 <- ggplot(historic, aes(x=input_quantity, color=extract_type, y=percent_gc_extended_ref)) + 
  guides(color=FALSE) +
  geom_point() + 
  xlab("Input DNA quantity (ng)") +
  ylab("Percent GC content") +
  theme(axis.title.y = element_text(size=12)) +
  theme(axis.title.x = element_text(size=12)) +
  stat_smooth(method="lm", fullrange=TRUE, alpha=0.1) 

# explore num of loci on target, historic samples
mh1 <- lm(num_loci_on_target ~ age + input_quantity + read_count_preclean, historic)
summary(mh1)

mh2 <- lm(num_loci_on_target ~ age + read_count_preclean, historic)
summary(mh2)

mh3 <- lm(num_loci_on_target ~ read_count_preclean, historic)
summary(mh3)

# explore drivers of percent duplicate reads, historic samples
hpd1 <- lm(percent_duplicates ~ age + input_quantity + read_count_preclean, historic)
summary(hpd1)

hpd2 <- lm(percent_duplicates ~ age + input_quantity, historic)
summary(hpd2)

hpd3 <- lm(percent_duplicates ~ input_quantity, historic)
summary(hpd3)

# explore drivers of specificity, historic samples
hsp1 <- lm(specificity_extended_ref ~ age + input_quantity + read_count_preclean, historic)
summary(hsp1)

hsp2 <- lm(specificity_extended_ref ~ age + read_count_preclean, historic)
summary(hsp2)

hsp3 <- lm(specificity_extended_ref ~ read_count_preclean, historic)
summary(hsp3)

# plot significant relationships, same grid
grid.arrange(p6,p7,ncol=2)

# make dataframe number of SNPs by matrix completeness
# entered manually because I'm lazy
count <- c(39105,39073,38929,38642,38077,37167,35697,33159,29818,26392,23175,20278,17570,14935,12172,9401,6825,4196,1690)
minInd <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)
snps <- as.data.frame(cbind(minInd,count))

# plot 
p8 <- ggplot(snps, aes(x=minInd, y=count)) + 
  guides(color=FALSE) +
  geom_point() + 
  xlab("Minimum individuals required") +
  ylab("Total number of SNPs") +
  theme(axis.title.y = element_text(size=12)) +
  theme(axis.title.x = element_text(size=12))



