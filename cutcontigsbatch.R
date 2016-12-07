

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
  