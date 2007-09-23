"do.aGFF.calc" <-
function (x, window, thresh) 
{
    chroms <- unique(as.character(x@annotation$Chromosome))
    nsamps <- dim(x@data)[2]
    ngenes <- dim(x@data)[1]
    y <- x@data
    cutpoints <- vector()
    for (i in 1:nsamps) {
      writeLines(paste("Working on sample", i))
      cutpoints[i] <- quantile(x@data[, i], probs = thresh)
      vals <- x@data[,i]>cutpoints[i]
      vals[vals==TRUE] <- 1
      vals[vals==FALSE] <- 0
      positive.count <- sum(vals)
      for (j in chroms) {
        writeLines(paste("Working on chromosome", j))
        sub <- x@annotation$Chromosome == j
        z <- windowChisq(x@annotation$Location[sub],
                         vals[sub],
                         window,
                         length(vals),
                         positive.count)
        y[sub, i] <- z$p.vals;
      }
    }
    colnames(y) <- colnames(x@data)
    names(cutpoints) <- colnames(x@data)
    ret <- new("aGFFCalc", vals = y, threshold = thresh, cutpoints = cutpoints, 
               data = x@data, annotation = x@annotation, samples = x@samples, 
               call = match.call())
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
