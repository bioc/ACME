"read.resultsGFF" <-
  function (fnames, path = ".", samples = NULL, notes = NULL, skip = 0, 
            sep = "\t", quote = "\"", ...) 
{
  if (is.null(path)) 
    fullfnames <- fnames
  else fullfnames <- file.path(path, fnames)
  fname <- fullfnames[1]
  Chromosome <- Source <- Type <- Location <- End <- Score <- Phast <- Strand <- Comments <- NULL
  first.pass <- TRUE
  for (f in fullfnames) {
    print(paste("Reading", f))
    tmp <- read.delim(f, header = F, skip = skip, sep = sep, 
                      quote = quote, ...)
    if (first.pass) {
      Annotation <- data.frame(tmp[, c(1:5, 7:9)])
      colnames(Annotation) <- c("Chromosome", "Source", 
                                "Type", "Location", "End", "Phase", "Strand", "Comment")
      first.pass <- FALSE
    }
    Score <- cbind(Score, as.numeric(tmp[, 6]))
  }
  fnames <- gsub(".gff", "", fnames)
  colnames(Score) <- fnames
  ord <- order(Annotation$Chromosome, Annotation$Location)
  Score <- as.matrix(Score[ord,])
  ret <- new("aGFF", annotation = Annotation[ord, ],
             data = Score, samples = data.frame(samples))
  return(ret)
}
