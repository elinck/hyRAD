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
  write(data[-toGrab], "output.fa")
}






