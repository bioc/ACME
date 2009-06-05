"do.aGFF.calc" <-
function (x, window, thresh) 
{
  chroms <- unique(chromosome(x))
  nsamps <- ncol(x)
  ngenes <- nrow(x)
  y <- matrix(NA,ncol=nsamps,nrow=ngenes)
  cutpoints <- vector()
  for (i in 1:nsamps) {
    writeLines(paste("Working on sample", i))
    cutpoints[i] <- quantile(exprs(x)[, i], probs = thresh, na.rm=TRUE)
    vals <- exprs(x)[,i]>cutpoints[i]
    ## Convert TRUE/FALSE to numeric
    vals <- as.numeric(vals)
    vals[vals==TRUE] <- 1
    vals[vals==FALSE] <- 0
    positive.count <- sum(vals)
    writeLines("Working on chromosome: ")
    for (j in chroms) {
      cat(j," ")
      sub <- chromosome(x) == j
      z <- windowChisq(start(x)[sub],
                       vals[sub],
                       window,
                       length(vals),
                       positive.count)
      y[sub, i] <- z$p.vals;
    }
  }
  colnames(y) <- sampleNames(x)
  names(cutpoints) <- sampleNames(x)
  ret <- new("ACMECalcSet",phenoData=phenoData(x),featureData=featureData(x),
             experimentData=experimentData(x),annotation=annotation(x),
             threshold=thresh,cutpoints=cutpoints,exprs=exprs(x),vals=y)
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
