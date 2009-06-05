"findClosestGene" <-
  function(chrom,pos,genome="hg17",position='txStart') {
    if (!exists('refflat')) {
      refflat <<- list()
      refflat[[genome]] <<- getRefflat(genome)
    } else if (!(genome %in% names(refflat))) {
      refflat[[genome]] <<- getRefflat(genome)
    }
    rf <- refflat[[genome]]
    chromsub <- rf$chrom==chrom
    diffdist <- rf[chromsub,position]-pos
    sub <- which(abs(diffdist)==min(abs(diffdist)))
    rf <- rf[chromsub,1:9][sub,]
    return(data.frame(rf,Distance=diffdist[sub]))
  }
