"do.aGFF.calc" <-
function (x, window, thresh) 
{
    chroms <- unique(as.character(Chromosome(x)))
    chroms <- chroms[!is.na(chroms)]
    nsamps <- dim(x)[2]
    ngenes <- dim(x)[1]
    y <- matrix(NA,nr=nrow(exprs(x)),nc=ncol(exprs(x)))
    cutpoints <- vector()
    for (i in 1:nsamps) {
      writeLines(paste("Working on sample", i))
      cutpoints[i] <- quantile(exprs(x)[,i], probs = thresh)
      vals <- exprs(x)[,i]>cutpoints[i]
      vals[vals==TRUE] <- 1
      vals[vals==FALSE] <- 0
      positive.count <- sum(vals)
      for (j in chroms) {
        writeLines(paste("Working on chromosome", j))
        sub <- Chromosome(x) == j
        sub[is.na(sub)] <- FALSE
        z <- windowChisq(Position(x)[sub],
                         vals[sub],
                         window,
                         length(vals),
                         positive.count)
        y[sub, i] <- z$p.vals;
      }
    }
    colnames(y) <- sampleNames(x)
    names(cutpoints) <- sampleNames(x)
    ret <- new("ACMECalc", featureData=featureData(x), phenoData=phenoData(x),pvals=y,exprs=exprs(x),
               threshold = thresh, windowsize=window)
    return(ret)
  }

windowChisq <- function(locations,ratios,windowsize,totprobes,posprobes) {
  ret <- .Call('windowChisq',locations,ratios,windowsize,totprobes,posprobes,PACKAGE='ACME')
  return(list(posProbes=ret[[1]],
              nProbes=ret[[2]],
              chivals=ret[[3]],
              values =ret[[4]],
              p.vals=1-pchisq(ret[[3]],1)))
}
