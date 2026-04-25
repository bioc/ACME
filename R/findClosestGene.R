"findClosestGene" <-
  function(chrom,pos,genome="hg17",position='txStart') {
    if (!exists(genome, envir=.acmeCache)) {
      assign(genome, getRefflat(genome), envir=.acmeCache)
    }
    rf <- get(genome, envir=.acmeCache)
    chromsub <- rf$chrom==chrom
    diffdist <- rf[chromsub,position]-pos
    sub <- which(abs(diffdist)==min(abs(diffdist)))
    rf <- rf[chromsub,1:9][sub,]
    return(data.frame(rf,Distance=diffdist[sub]))
  }
