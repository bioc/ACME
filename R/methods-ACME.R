setMethod("initialize",signature(.Object="ACME"),
          function(.Object,
                   Chromosome,
                   Position,
                   featureData=as(data.frame(Chromosome=Chromosome,Position=Position),"AnnotatedDataFrame"),
                   ...) {
            print('ACMEInit')
            print(match.call())
            if (missing(featureData)) {
              if (missing(Chromosome) || missing(Position)) {
                stop("One of featureData or Chromosome and Position must be specified")
              }
              featureData=as(data.frame(Chromosome=Chromosome,Position=Position),"AnnotatedDataFrame")
            }
            .Object <- callNextMethod(.Object,
                                      featureData=featureData,
                           ...)
            .Object <- .Object[order(pData(featureData(.Object))$Chromosome,pData(featureData(.Object))$Position),]
            return(.Object)
          })

setMethod("initialize",signature(.Object="ACMECalc"),
          function(.Object,
                   Chromosome,
                   Position,
                   featureData=as(data.frame(Chromosome=Chromosome,Position=Position),"AnnotatedDataFrame"),
                   windowsize,
                   threshold,
                   ...) {
            print("ACMECalcInit")
            print(match.call())
            .Object <- callNextMethod(.Object,featureData=featureData,...)
            .Object@windowsize=as.integer(windowsize)
            .Object@threshold=threshold
            return(.Object)
          })

setValidity("ACME",function(object) {
  if(!assayDataValidMembers(assayData(object),c("exprs"))) {
    stop('assayData does not contain an "exprs" member')
  }
  if(any(is.na(match(c("Chromosome","Position"),varLabels(featureData(object)))))) {
    stop('featureData must contain both a "Chromosome" and a "Position" column')
  }
  return(TRUE)
})

setValidity("ACMECalc",function(object) {
  callNextMethod()
})

setMethod("windowsize","ACMECalc",function(object) {
  return(object@windowsize)
})

setMethod("threshold","ACMECalc",function(object) {
  return(object@threshold)
})

setMethod("show","ACMECalc",function(object) {
  callNextMethod(object)
  cat("Calculation Protocol Parameters:\n")
  cat(sprintf("  threshold:  %f\n  windowsize: %d base pairs\n",threshold(object),windowsize(object)))
})

setMethod("Position","ACME",function(object) {
  return(pData(featureData(object))$Position)
})

setMethod("Chromosome","ACME",function(object) {
  return(pData(featureData(object))$Chromosome)
})


setMethod("Pos","ACME",function(object) {
  return('Testing')
})


setMethod("Position","ACMECalc",function(object) {
  return(pData(featureData(object))$Position)
})

setMethod("Chromosome","ACMECalc",function(object) {
  return(pData(featureData(object))$Chromosome)
})

setMethod("exprs","ACME",function(object) {
  return(assayData(object)$exprs)
})

setMethod("pvals","ACMECalc",function(object) {
  return(assayData(object)$pvals)
})

setMethod("QCData","ACME",function(object) {
  return(summary(exprs(object)))
})

setMethod("QCData","ACMECalc",function(object) {
  return(callNextMethod(object))
})

setMethod("combine",c("ACMECalc","ACMECalc"),function(x,y,...) {
  if(class(x) != class(y))
    stop(paste("objects must be the same class, but are ",
               class(x), ", ", class(y), sep=""))
  if((window(x)!=window(y)) || (threshold(x)!=threshold(y))) {
    stop("Window and threshold parameters do not match")
  }
})

plotgffcalc <- function(x,y='missing',chrom,samples=NULL,...) {
  nsamps <- 0
  if (is.null(samples)) {
    samples <- 1:ncol(x)
    nsamps <- length(samples)
  }
  if (!is.numeric(samples))
    samples <- which(sampleNames(x) %in% samples)
  sub <- Chromosome(x)==chrom
  sub[is.na(sub)] <- FALSE
  if (sum(sub)==0)
    stop('No matching chromosome in the data')
  if (nsamps>1) par(ask=T)
  par(mar=c(4,5,4,5))
  for (i in samples) {
    plot(Position(x)[sub],exprs(x)[sub,i],
         axes=F,col='gray85',pch=20,xlab="",ylab="",main="",...)
    axis(side=4)
    abline(h=0,col='gray85')
    abline(h=quantile(exprs(x)[,i],threshold(x)),col='gray85',lty=2)
    par(new=T)
    plot(Position(x)[sub],-log10(pvals(x)[sub,i]),
         main=paste('Chromosome:',chrom,', Sample:',sampleNames(x)[i]),
         ylab='-log10(p-value)',xlab='Chromosome Position',type='b',pch=20,
         col='red',...)
  }
  par(ask=F)
}
