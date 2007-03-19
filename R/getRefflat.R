"getRefflat" <-
  function(genome="hg17") {
    tmpfile <- tempfile()
    download.file(paste('http://hgdownload.cse.ucsc.edu/goldenPath/',
                        genome,'/database/refFlat.txt.gz',sep=""),
                  tmpfile,mode='wb')
    rf <- read.delim(gzfile(tmpfile),header=FALSE,sep="\t")
    colnames(rf) <- c('geneName','name','chrom','strand','txStart','txEnd',
                      'cdsStart','cdsEnd','exonCount','exonStarts','exonEnds')
    txEndNeg <- rf$txStart
    txStartNeg <- rf$txEnd
    cdsStartNeg <- rf$cdsEnd
    cdsEndNeg <- rf$cdsStart
    NegStrand <- rf$strand=='-'
    ## Fix negative strand stuff
    rf[NegStrand,'cdsEnd'] <- cdsEndNeg[NegStrand]
    rf[NegStrand,'cdsStart'] <- cdsStartNeg[NegStrand]
    rf[NegStrand,'txEnd'] <- txEndNeg[NegStrand]
    rf[NegStrand,'txStart'] <- txStartNeg[NegStrand]
    return(rf)
  }

