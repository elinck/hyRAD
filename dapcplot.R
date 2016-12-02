#### DAPCplot: an R function for plotting DAPC membership probabilities from adegenet ####
# usage: dapcplot(x, pop, colors), where
#        x = .dapc object  
#        pop = vector of population assignments (e.g., pop_a <- c(1,1,1,2,2,2))
#        colors = vector of at least k colors (e.g., c("gold","forestgreen","magenta3","orangered"))
#        sort = TRUE will sort bars in plot by pop membership
# improves on adegenet's native compoplot() by holding colors constant in individuals and offering ability to sort by K (=pop membership)
# inspired by cjbattey's structurePlot.R function (https://github.com/cjbattey/RADplots/blob/master/structurePlot.R)
# to make plots that look like the standard STRUCTURE plots from the GUI, set spacing=0 and outline='black'
# good colors: c("gold","forestgreen","magenta3","orangered","cornflowerblue","orange","sienna","dodgerblue4"), or
# c(colors=brewer.pal(12,"Set3"))

dapcplot <- function(x, pop, colors, sort=FALSE){
  df1 <- as.data.frame(x$posterior)
  df1 <- df1[,1:ncol(df1)]
  k <- ncol(df1)
  df1 <- cbind(df1,pop)
  colnames(df1) <- c(1:k,"pop")
  if (sort=="TRUE"){
    df1 <- arrange(df1,pop)
  }
  n.pops <- max(as.numeric(as.character(df1$pop)))
  palette <- c(rep("NA",k))
  for(i in 1:n.pops){
    d <- subset(df1,pop==i)
    e <- colMeans(d[,1:k])
    f <- as.numeric(names(e[which(e==max(e))]))
    if (palette[f] == "NA"){
      palette[f] <- colors[i]
    }
  }
  palette[which(palette == "NA")] <- colors[which(colors%in%palette == FALSE)]
  colors <- palette
  barplot(t(df1[,1:k]),border = NA,col=colors,cex.names = 1,legend = FALSE)
}
