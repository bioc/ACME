write.sgr <- function(agff,raw=TRUE,vals=TRUE,directory=NULL) {
  if (!(class(agff) %in% c('ACME','ACMECalc'))) {
    stop('ACME or derived class object. Got an object of class: ',class(agff))
  }
  sampnames <- sampleNames(agff)
  sub1 <- (!is.na(Chromosome(agff))) & (!is.na(Position(agff)))
  agff <- agff[sub1,]
  for(i in 1:ncol(agff)) {
    if (class(agff)=='ACMECalc'){
      if(vals) {
        if(!is.null(directory)) {
          filename <- file.path(directory,sprintf("%s_thresh%3.2f.sgr",sampnames[i],threshold(agff)))
        } else {
          filename <- sprintf("%s_thresh%3.2f.sgr",sampnames[i],threshold(agff))
        }
        cat(filename,"\n")
        write.table(data.frame(Chromosome(agff),Position(agff),-log10(pvals(agff)[,i]+min(pvals(agff)[pvals(agff)[,i]>0,i]))),file=filename,
                    sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
      }
    }
    if(raw) {
      if(!is.null(directory)) {
        filename <- file.path(directory,sprintf("%s_raw.sgr",sampnames[i]))
      } else {
        filename <- sprintf("%s_raw.sgr",sampnames[i])
      }
      cat(filename,"\n")
      write.table(data.frame(Chromosome(agff),Position(agff),exprs(agff)[,i]),file=filename,
                  sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
    }
  }
}
