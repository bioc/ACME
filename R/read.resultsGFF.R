"read.resultsGFF" <-
  function (fnames, path = NULL, samples = NULL, notes = NULL, skip = 0, 
            sep = "\t", quote = "\"", comment.char="#", ...) 
{
  if (is.null(path)) 
    fullfnames <- fnames
  else fullfnames <- file.path(path, fnames)
  fname <- fullfnames[1]
  Chromosome <- Source <- Type <- Location <- End <- Score <- Phast <- Strand <- Comments <- NULL
  first.pass <- TRUE
  for (f in fullfnames) {
    print(paste("Reading", f))
    if(first.pass) {
      tmp <- read.delim(f, header = F, skip = skip, sep = sep, comment.char=comment.char, 
                        quote = quote, colClasses=c('character','character','character','integer','integer','numeric','character','character','character'),
                        ...)
      Annotation <- data.frame(tmp[, c(1:5, 7:9)])
      colnames(Annotation) <- c("Chromosome", "Source", 
                                "Type", "Position", "End", "Phase", "Strand", "Comment")
      first.pass <- FALSE
      Score <- cbind(Score, as.numeric(tmp[, 6]))
    } else {
      tmp <- read.delim(f, header = F, skip = skip, sep = sep, comment.char=comment.char, 
                        quote = quote, colClasses=c("NULL","NULL","NULL","NULL","NULL",'numeric',"NULL","NULL","NULL"),
                        ...)
      Score <- cbind(Score,as.numeric(tmp[,1]))
    }
  }
  fnames <- gsub(".gff", "", fnames)
  colnames(Score) <- fnames
  ord <- order(Annotation$Chromosome, Annotation$Position)
  rownames(Score) <- make.unique(as.character(1:nrow(Score)))
  rownames(Annotation) <- rownames(Score)
  Score <- as.matrix(Score[ord,])
  return(new("ACME",exprs=Score,featureData=as(Annotation,"AnnotatedDataFrame")))
#  ret <- list(annotation = Annotation[ord, ],
#             data = Score, samples = data.frame(samples))
  return(ret)
}
