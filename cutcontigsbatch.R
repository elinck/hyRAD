#### cutcontigsbatch: an R function for removing specific contigs from multiple fasta file  ####
# credit to R. Harris (https://github.com/rebzzy/) for the lapply() solution
# useage: cutcontigsbatch(directory, contigs) where
#         directory = path to directory that only contains fasta files of interest
#         contigs = text file with contig names of interest, 1 per row (contigs <-readLines("test_contigs.txt"))
# notes: use extractcontigIDs.R to generate contigs input text file from BLAST results (e.g., to remove contamination)

cutcontigsbatch <- function(directory, contigs){
  files <- list.files(directory, pattern = ".fasta")
  for (i in files){
    fasta <- readLines(i)
    matchingLines <- unlist(lapply(contigs, function(y) {
      pattern <- paste0("^>.*", y)
      grep(pattern, fasta)
    }))
    toGrab <- sort(c(matchingLines, matchingLines + 1))
    print(fasta[-toGrab])
    write(fasta[-toGrab],file=paste("cleaned_",i,sep=""))
    }
  }
  