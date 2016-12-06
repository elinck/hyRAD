#### extractcontigIDs: an R function for identifying matching contig IDs from multiple BLAST (m6) results ####
# usage: extractcontigIDs(directory, outfile_name), where
#        directory = path to files, e.g. "/User/elinck/data/"
#        outfile_name = name for resulting list of BLAST matches, e.g. "contigs.txt"

extractcontigIDs <- function(directory, outfile_name){
  files <- list.files(directory, pattern = ".txt")
  lst <- vector("list", length(files))
  for(i in 1:length(files)) {
    lst[[i]] <- as.data.frame(read.table(files[i],sep="\t"))
  }
  IDs <- unlist(lapply(lst, function(x){
    as.vector(x$V1)
  }))
  contigs <- sort(c(IDs))
  write(contigs, outfile_name)
}