write.sgr <- function(x,raw=TRUE,vals=TRUE,directory='.') {
  for(i in 1:ncol(x)) {
    if (class(x)=='ACMECalcSet'){
      if(vals) {
        filename <- file.path(directory,sprintf("%s_thresh%3.2f.sgr",sampleNames(x)[i],threshold(x)))
        cat(filename,"\n")
        write.table(data.frame(chromosome=chromosome(x),start=start(x),-log10(vals(x)[,i]+min(vals(x)[vals(x)[,i]>0,i]))),
                    file=filename,sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
      }
    }
    if(raw) {
      filename <- file.path(directory,sprintf("%s_raw.sgr",sampleNames(x)[i]))
      cat(filename,"\n")
      write.table(data.frame(chromosome=chromosome(x),start=start(x),exprs(x)[,i]),
                  file=filename,sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
    }
  }
}
