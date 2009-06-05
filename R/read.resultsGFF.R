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
      colnames(Annotation) <- c("chromosome", "source", 
                                "type", "start", "end", "phase", "strand", "comment")
      first.pass <- FALSE
    }
    Score <- cbind(Score, as.matrix(as.numeric(tmp[, 6]),ncol=1))
  }
  fnames <- gsub(".gff", "", fnames)
  colnames(Score) <- fnames
  ret <- new("ACMESet",exprs=Score,featureData=
           as(Annotation,"AnnotatedDataFrame"),phenoData=
           as(data.frame(fullfnames,row.names=fnames),"AnnotatedDataFrame"))
  sampleNames(ret) <- fnames
  ret <- ret[order(chromosome(ret),start(ret),end(ret)),]
  return(ret)
}
