write.sgr <- function(agff,raw=TRUE,vals=TRUE,directory='.') {
  if (!(class(agff) %in% c('aGFF','aGFFCalc'))) {
    stop('Need agff to be an aGFF or aGFFCalc object')
  }
  for(i in 1:ncol(agff@data)) {
    if (class(agff)=='aGFFCalc'){
      if(vals) {
        filename <- file.path(directory,sprintf("%s_thresh%3.2f.sgr",colnames(agff@vals)[i],agff@threshold))
        cat(filename,"\n")
        write.table(data.frame(agff@annotation[,c('Chromosome','Location')],-log10(agff@vals[,i]+min(agff@vals[agff@vals[,i]>0,i]))),file=filename,
                    sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
      }
    }
    if(raw) {
      filename <- file.path(directory,sprintf("%s_raw.sgr",colnames(agff@data)[i]))
      cat(filename,"\n")
      write.table(data.frame(agff@annotation[,c('Chromosome','Location')],agff@data[,i]),file=filename,
                  sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
    }
  }
}
