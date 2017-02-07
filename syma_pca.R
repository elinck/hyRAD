setwd("/Users/ethanlinck/Dropbox/Syma/hyRAD/data")
#install.packages("data.table");install.packages("ggplot2");install.packages("adegenet");install.packages("StAMPP");install.packages("data.table");install.packages("ggmap");install.packages("mapdata")
library(ggplot2);library(adegenet);library(reshape2);library(StAMPP);library(data.table);library(pegas);library(plyr);library(RColorBrewer);library(ggmap);library(mapdata);library(maps)
source("/Users/ethanlinck/Dropbox/Syma/hyRAD/dapcplot.R")

### PCA / DAPC analysis of population genetic structure in S. torotoro

### upload replicates of the same dataset for testing K1 - K8
a <- read.delim("snps_cleaned_complete_maf05_adg", sep="\t")
b <- read.delim("snps_cleaned_complete_maf05_adg", sep="\t")
c <- read.delim("snps_cleaned_complete_maf05_adg", sep="\t")
d <- read.delim("snps_cleaned_complete_maf05_adg", sep="\t")
e <- read.delim("snps_cleaned_complete_maf05_adg", sep="\t")
f <- read.delim("snps_cleaned_complete_maf05_adg", sep="\t")
g <- read.delim("snps_cleaned_complete_maf05_adg", sep="\t")
h <- read.delim("snps_cleaned_complete_maf05_adg", sep="\t")

### turn replicate datasets into genind files
torotoro_a <- df2genind(a, NA.char="NA", ploidy=2, ncode=2)
torotoro_b <- df2genind(b, NA.char="NA", ploidy=2, ncode=2)
torotoro_c <- df2genind(c, NA.char="NA", ploidy=2, ncode=2)
torotoro_d <- df2genind(d, NA.char="NA", ploidy=2, ncode=2)
torotoro_e <- df2genind(e, NA.char="NA", ploidy=2, ncode=2)
torotoro_f <- df2genind(f, NA.char="NA", ploidy=2, ncode=2)
torotoro_g <- df2genind(g, NA.char="NA", ploidy=2, ncode=2)
torotoro_h <- df2genind(h, NA.char="NA", ploidy=2, ncode=2)

### identify sample names, putative pops (=subspecies)
torotoro_a@pop <- factor(c("pseuestes","pseuestes","meeki","pseuestes","meeki","torotoro","tentelare","pseuestes","pseuestes","torotoro","ochracea","torotoro","torotoro","torotoro","torotoro","torotoro","torotoro","torotoro","torotoro"))
indNames(torotoro_a) <- c("t_pseu_31", "t_pseu_32", "t_meek_33", "t_pseu_34", "t_meek_26", "t_toro_36","t_tent_04", "t_pseu_12", "t_pseu_13", "t_ochr_03","t_toro_14","t_toro_02","t_toro_06", "t_toro_07", "t_toro_08", "t_toro_11", "t_toro_15", "t_toro_18", "t_toro_19")
torotoro_b@pop <- factor(c("pseuestes","pseuestes","meeki","pseuestes","meeki","torotoro","tentelare","pseuestes","pseuestes","torotoro","ochracea","torotoro","torotoro","torotoro","torotoro","torotoro","torotoro","torotoro","torotoro"))
indNames(torotoro_b) <- c("t_pseu_31", "t_pseu_32", "t_meek_33", "t_pseu_34", "t_meek_26", "t_toro_36","t_tent_04", "t_pseu_12", "t_pseu_13", "t_ochr_03","t_toro_14","t_toro_02","t_toro_06", "t_toro_07", "t_toro_08", "t_toro_11", "t_toro_15", "t_toro_18", "t_toro_19")
torotoro_c@pop <- factor(c("pseuestes","pseuestes","meeki","pseuestes","meeki","torotoro","tentelare","pseuestes","pseuestes","torotoro","ochracea","torotoro","torotoro","torotoro","torotoro","torotoro","torotoro","torotoro","torotoro"))
indNames(torotoro_c) <- c("t_pseu_31", "t_pseu_32", "t_meek_33", "t_pseu_34", "t_meek_26", "t_toro_36","t_tent_04", "t_pseu_12", "t_pseu_13", "t_ochr_03","t_toro_14","t_toro_14","t_toro_06", "t_toro_07", "t_toro_08", "t_toro_11", "t_toro_15", "t_toro_18", "t_toro_19")
torotoro_d@pop <- factor(c("pseuestes","pseuestes","meeki","pseuestes","meeki","torotoro","tentelare","pseuestes","pseuestes","torotoro","ochracea","torotoro","torotoro","torotoro","torotoro","torotoro","torotoro","torotoro","torotoro"))
indNames(torotoro_d) <- c("t_pseu_31", "t_pseu_32", "t_meek_33", "t_pseu_34", "t_meek_26", "t_toro_36","t_tent_04", "t_pseu_12", "t_pseu_13", "t_ochr_03","t_toro_14","t_toro_02","t_toro_06", "t_toro_07", "t_toro_08", "t_toro_11", "t_toro_15", "t_toro_18", "t_toro_19")
torotoro_e@pop <- factor(c("pseuestes","pseuestes","meeki","pseuestes","meeki","torotoro","tentelare","pseuestes","pseuestes","torotoro","ochracea","torotoro","torotoro","torotoro","torotoro","torotoro","torotoro","torotoro","torotoro"))
indNames(torotoro_e) <- c("t_pseu_31", "t_pseu_32", "t_meek_33", "t_pseu_34", "t_meek_26", "t_toro_36","t_tent_04", "t_pseu_12", "t_pseu_13", "t_ochr_03","t_toro_14","t_toro_02","t_toro_06", "t_toro_07", "t_toro_08", "t_toro_11", "t_toro_15", "t_toro_18", "t_toro_19")
torotoro_f@pop <- factor(c("pseuestes","pseuestes","meeki","pseuestes","meeki","torotoro","tentelare","pseuestes","pseuestes","torotoro","ochracea","torotoro","torotoro","torotoro","torotoro","torotoro","torotoro","torotoro","torotoro"))
indNames(torotoro_f) <- c("t_pseu_31", "t_pseu_32", "t_meek_33", "t_pseu_34", "t_meek_26", "t_toro_36","t_tent_04", "t_pseu_12", "t_pseu_13", "t_ochr_03","t_toro_14","t_toro_02","t_toro_06", "t_toro_07", "t_toro_08", "t_toro_11", "t_toro_15", "t_toro_18", "t_toro_19")
torotoro_g@pop <- factor(c("pseuestes","pseuestes","meeki","pseuestes","meeki","torotoro","tentelare","pseuestes","pseuestes","torotoro","ochracea","torotoro","torotoro","torotoro","torotoro","torotoro","torotoro","torotoro","torotoro"))
indNames(torotoro_g) <- c("t_pseu_31", "t_pseu_32", "t_meek_33", "t_pseu_34", "t_meek_26", "t_toro_36","t_tent_04", "t_pseu_12", "t_pseu_13", "t_ochr_03","t_toro_14","t_toro_02","t_toro_06", "t_toro_07", "t_toro_08", "t_toro_11", "t_toro_15", "t_toro_18", "t_toro_19")
torotoro_h@pop <- factor(c("pseuestes","pseuestes","meeki","pseuestes","meeki","torotoro","tentelare","pseuestes","pseuestes","torotoro","ochracea","torotoro","torotoro","torotoro","torotoro","torotoro","torotoro","torotoro","torotoro"))
indNames(torotoro_h) <- c("t_pseu_31", "t_pseu_32", "t_meek_33", "t_pseu_34", "t_meek_26", "t_toro_36","t_tent_04", "t_pseu_12", "t_pseu_13", "t_ochr_03","t_toro_14","t_toro_02","t_toro_06", "t_toro_07", "t_toro_08", "t_toro_11", "t_toro_15", "t_toro_18", "t_toro_19")

#pca plot
scaled.data_a <- scaleGen(torotoro_a, NA.method = c("zero"))
torotoro.pca_a <- dudi.pca(scaled.data_a,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)
c <- torotoro.pca_a$li
c$pop <- torotoro_a@pop
ggplot()+
  geom_point(data=c, aes(x=Axis1,y=Axis2,col=pop), size=3)

### investigate correlation between PCs and gc content, age, and initial conc. ###
pcas <- torotoro.pca_a$li
pc1 <- pcas[,1]
pc2 <- pcas[,2]
pc3 <- pcas[,3]

#variables input by hand because I'm lazy -- will clean up
data <- read.csv("syma_torotoro_hyRAD_performance.csv")
sensitivity <- data$sensitivity
specificity <- data$specificity
enrichment <- data$enrichment

#linear models to test significance of correlations
pca_lm <- lm(pc1 ~ sensitivity)
summary(pca_lm) #significant, p=7.11e-09
pca_lm2 <- lm(pc2 ~ sensitivity)
summary(pca_lm2) 
pca_lm3 <- lm(pc3 ~ sensitivity)
summary(pca_lm3) 
pca_lm4 <- lm(pc1 ~ specificity)
summary(pca_lm4) #significant, p=7.11e-09
pca_lm5 <- lm(pc2 ~ specificity)
summary(pca_lm5) 
pca_lm6 <- lm(pc3 ~ specificity)
summary(pca_lm6) 
pca_lm7 <- lm(pc1 ~ enrichment)
summary(pca_lm7) #significant, 7.74e-08 ***
pca_lm8 <- lm(pc2 ~ enrichment)
summary(pca_lm8) 
pca_lm9 <- lm(pc3 ~ enrichment)
summary(pca_lm9) 
pca_lm10 <- lm(pc1 ~ age)
summary(pca_lm10) #significant, 5.55e-09 ***
pca_lm11 <- lm(pc2 ~ age)
summary(pca_lm11)
pca_lm12 <- lm(pc3 ~ age)
summary(pca_lm12)
pca_lm13 <- lm(pc1 ~ conc)
summary(pca_lm13) #significant, 2.43e-05 ***
pca_lm14 <- lm(pc2 ~ conc)
summary(pca_lm14)
pca_lm15 <- lm(pc3 ~ conc)
summary(pca_lm15)

# plot GC by pca 1
# not run for ms
#p1 <- ggplot(data, aes(x=gc, y=pc1)) + 
#  guides(color=FALSE) +
#  geom_point() + 
#  xlab("%gc") +
#  ylab("pca 1 score") +
#  theme(axis.title.y = element_text(size=16)) +
#  theme(axis.title.x = element_text(size=16)) +
#  stat_smooth(method="lm", fullrange=TRUE, alpha=0.1) 

# plot GC across individuals (for figure)
# not run for ms
# names <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)
# gc <- c(44.84,49.22,49.21,45.55,49.42,48.43,45.49,47.61,51.46,45.99,48.8,52.11,47.69,48.58,49.32,50.87,49.75,51.65,55.08)
# data2 <- as.data.frame(cbind(gc,names))

# plot 
# p2 <- ggplot(data2, aes(x=names, y=gc, group=1)) + 
#   guides(color=FALSE) +
#  geom_line() +
#  geom_point() + 
#  xlab("") +
#  ylab("%GC Content") +
#  theme(axis.title.y = element_text(size=14)) +
#  theme(axis.title.x = element_text(size=14))

###dapc w/max #PC's, then test for opt # to keep w/alpha scores
temp_a <- dapc(torotoro_a, n.pca=6, n.da=6)
opt.pc_a <- optim.a.score(temp_a) #selects 6 PCAs

###find clusters w/o using pop priors
clust_a <- find.clusters(torotoro_a, n.pca=6, n.clust=1, max.n.clust=18)
clust_b <- find.clusters(torotoro_b, n.pca=6, n.clust=2, max.n.clust=18)
clust_c <- find.clusters(torotoro_c, n.pca=6, n.clust=3, max.n.clust=18)
clust_d <- find.clusters(torotoro_d, n.pca=6, n.clust=4, max.n.clust=18)
clust_e <- find.clusters(torotoro_e, n.pca=6, n.clust=5, max.n.clust=18)
clust_f <- find.clusters(torotoro_f, n.pca=6, n.clust=6, max.n.clust=18)
clust_g <- find.clusters(torotoro_g, n.pca=6, n.clust=7, max.n.clust=18)
clust_h <- find.clusters(torotoro_h, n.pca=6, n.clust=8, max.n.clust=18)

###dapc and loadings
torotoro_a@pop <- clust_a$grp
torotoro_b@pop <- clust_b$grp
torotoro_c@pop <- clust_c$grp
torotoro_d@pop <- clust_d$grp
torotoro_e@pop <- clust_e$grp
torotoro_f@pop <- clust_f$grp
torotoro_g@pop <- clust_g$grp
torotoro_h@pop <- clust_h$grp

# run dapc for different cluster assignments K1-K8
torotoro.dapc_a <- dapc(torotoro_a,n.pca=6,n.da=2)
torotoro.dapc_b <- dapc(torotoro_b,n.pca=6,n.da=2)
torotoro.dapc_c <- dapc(torotoro_c,n.pca=6,n.da=2)
torotoro.dapc_d <- dapc(torotoro_d,n.pca=6,n.da=2)
torotoro.dapc_e <- dapc(torotoro_e,n.pca=6,n.da=2)
torotoro.dapc_f <- dapc(torotoro_f,n.pca=6,n.da=2)
torotoro.dapc_g <- dapc(torotoro_g,n.pca=6,n.da=2)
torotoro.dapc_h <- dapc(torotoro_h,n.pca=6,n.da=2)

#identify majority pop membership visually
compoplot(torotoro.dapc_a)
compoplot(torotoro.dapc_b)
compoplot(torotoro.dapc_c)
compoplot(torotoro.dapc_d)
compoplot(torotoro.dapc_e)
compoplot(torotoro.dapc_f)
compoplot(torotoro.dapc_g)
compoplot(torotoro.dapc_h)

###recode pops based on compoplot results
pb <- c(1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2)
pc <- c(1,1,1,1,1,1,2,2,2,3,3,3,3,3,3,3,3,3,3)
pd <- c(1,1,1,1,1,1,2,2,2,4,3,3,3,3,3,3,3,3,3)
pe <- c(1,1,1,1,1,1,2,2,2,4,5,3,3,3,3,3,3,3,3)
pf <- c(1,1,1,1,6,6,2,2,2,4,5,3,3,3,3,3,3,3,3)


#define colors
set <- c("gold","forestgreen","magenta3","orangered","cornflowerblue","orange","sienna","dodgerblue4")

#plot 
dapcplot(torotoro.dapc_b,pb,colors=set,sort=FALSE)
dapcplot(torotoro.dapc_c,pc,colors=set,sort=FALSE)
dapcplot(torotoro.dapc_d,pd,colors=set,sort=FALSE)

#map k2
localities <- read.csv("/Users/ethanlinck/Dropbox/Syma/figures/syma_localities.csv")
localities1 <- cbind(localities,pb)
map <- map_data("worldHires", xlim=c(130,154), ylim=c(-10,1))
ggplot() + coord_map()+
  geom_path(data=map, aes(x=long, y=lat, group=group)) +
  geom_point(data=localities1, size=5, aes(x=long, y=lat, col=as.factor(pb))) +
  scale_colour_manual(values = set,name="DAPC Cluster") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Longitude") +
  ylab("Latitude") +
  guides(fill=guide_legend(title=NULL))

#map k3
localities2 <- cbind(localities,pc)
map <- map_data("worldHires", xlim=c(130,154), ylim=c(-10,1))
ggplot() + coord_map()+
  geom_path(data=map, aes(x=long, y=lat, group=group)) +
  geom_point(data=localities2, size=5, aes(x=long, y=lat, col=as.factor(pc))) +
  scale_colour_manual(values = set,name="DAPC Cluster") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Longitude") +
  ylab("Latitude") +
  guides(fill=guide_legend(title=NULL))

#map k4
localities3 <- cbind(localities,pd)
map <- map_data("worldHires", xlim=c(130,154), ylim=c(-10,1))
ggplot() + coord_map()+
  geom_path(data=map, aes(x=long, y=lat, group=group)) +
  geom_point(data=localities2, size=5, aes(x=long, y=lat, col=as.factor(pd))) +
  scale_colour_manual(values = set,name="DAPC Cluster") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Longitude") +
  ylab("Latitude") +
  guides(fill=guide_legend(title=NULL))


pg <- c(1,1,2,1,3,3,1,3,3,5,6,6,3,6,7,6,6,3,6)

#map mtDNA
localities4 <- cbind(localities,pg)
map <- map_data("worldHires", xlim=c(130,154), ylim=c(-10,1))
ggplot() + coord_map()+
  geom_path(data=map, aes(x=long, y=lat, group=group)) +
  geom_point(data=localities2, size=5, aes(x=long, y=lat, col=as.factor(pg))) +
  scale_colour_manual(values = set,name="DAPC Cluster") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Longitude") +
  ylab("Latitude") +
  guides(fill=guide_legend(title=NULL))

pg <- c(1,1,2,1,3,3,1,3,3,5,6,6,3,6,7,6,6,3,6)
ph <- c(1,1,1,1,3,3,1,1,1,5,3,3,3,3,7,3,3,3,3)


#map nuDNA
localities5 <- cbind(localities,ph)
map <- map_data("worldHires", xlim=c(130,154), ylim=c(-10,1))
ggplot() + coord_map()+
  geom_path(data=map, aes(x=long, y=lat, group=group)) +
  geom_point(data=localities2, size=5, aes(x=long, y=lat, col=as.factor(ph))) +
  scale_colour_manual(values = set,name="DAPC Cluster") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Longitude") +
  ylab("Latitude") +
  guides(fill=guide_legend(title=NULL))

### mantel test for IBD
z <- read.delim("snps_cleaned_complete_maf05_adg", sep="\t")
torotoro_z <- df2genind(z, NA.char="NA", ploidy=2, ncode=2)
torotoro_z@pop <- factor(c("t_pseu_31", "t_pseu_32", "t_meek_33", "t_pseu_34", "t_meek_26", "t_toro_36","t_tent_04", "t_pseu_12", "t_pseu_13", "t_ochr_03","t_toro_14","t_toro_02","t_toro_06", "t_toro_07", "t_toro_08", "t_toro_11", "t_toro_15", "t_toro_18", "t_toro_19"))
toro <- genind2genpop(torotoro_z)
xy <- cbind(localities$lat, localities$long)
distxy <- dist(xy)
distgen <- dist.genpop(toro,method=2)
ibd <- mantel.randtest(distgen,distxy) #not significant
plot(ibd)

