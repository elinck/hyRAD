#### cutcontigs: an R function for removing specific contigs from a fasta file  ####
# credit to R. Harris (https://github.com/rebzzy/) for the lapply() solution
# useage: cutcontigs(data,contigs,output) where
#         data = fasta file (e.g., data <- readLines("data.fasta"))    
#         contigs = text file with contig names of interest, 1 per row (contigs <-readLines("test_contigs.txt"))
#         output = outfile name

cutcontigs <- function(data,contigs,output){
  matchingLines <- unlist(lapply(contigs, function(x) {
    pattern <- paste0("^>.*", x)
    grep(pattern, data)
  }))
  toGrab <- sort(c(matchingLines, matchingLines + 1))
  write(data[-toGrab], file=paste(output,".fasta",sep=""))
}

setwd("/Users/ethanlinck/Dropbox/Syma/hyRAD/")
data <- readLines("combined_targetedRegionAndFlanking.fasta")
contigs <- readLines("contam_contigs.txt")
  
cutcontigsbatch <- function(contigs){
  files <- list.files(".fasta")
  for(i in 1:length(files)) {
    lst[[i]] <- read.csv(files[i])
  }
  lapply(lst, function(x) {
    matchingLines <- unlist(lapply(contigs, function(x) {
      pattern <- paste0("^>.*", x)
      grep(pattern, data)
    }))
    toGrab <- sort(c(matchingLines, matchingLines + 1))
    write(data[-toGrab], "output.fa")
  })
}








