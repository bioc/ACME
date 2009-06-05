### The ACMESet class is a simple container
### that must have at least the columns
### chromosome, start, and end (with those names)
### in the fData.
setClass('ACMESet',
         contains="ExpressionSet",
         validity=function(object) {
           msg <- NULL
           cnames <- colnames(fData(object))
           chk <- match(c('chromosome','start','end'),cnames)
           if(any(is.na(chk)) & nrow(fData(object))>0) {
             msg <- "Column names of featuredata must include 'chromosome','start', and 'end'"
           }
           if(is.null(msg)) return(TRUE)
           return(msg)
         })
setClass('ACMECalcSet',
         representation(cutpoints="numeric",
                        threshold="numeric"),
         contains='ACMESet',
         validity=function(object) {
           msg <- NULL
           if(length(object@cutpoints)!=ncol(object)) {
             msg <- "cutpoints should be the same length as the number of samples"
           }
           if(!is.null(msg)) return(msg)
           return(TRUE)
         })
setMethod("initialize", "ACMECalcSet",
          function(.Object,
                   assayData,
                   phenoData = annotatedDataFrameFrom(assayData, byrow=FALSE),
                   featureData = annotatedDataFrameFrom(assayData, byrow=TRUE),
                   experimentData = new("MIAME"),
                   annotation = character(),
                   cutpoints=numeric(),
                   threshold=numeric(),
                   exprs = new("matrix"),
                   vals = new("matrix"),
                   ... ) {
            .Object@cutpoints <- cutpoints
            .Object@threshold <- threshold
            if (missing(assayData)) {
              if (missing(phenoData))
                phenoData <- annotatedDataFrameFrom(exprs, byrow=FALSE)
              if (missing(featureData))
                featureData <- annotatedDataFrameFrom(exprs, byrow=TRUE)
              .Object <- callNextMethod(.Object,
                                        phenoData=phenoData,
                                        featureData=featureData,
                                        experimentData=experimentData,
                                        annotation=annotation,
                                        exprs=exprs,
                                        vals=vals,
                                        ...)
            } else if(missing(exprs) & missing(vals)) {
              .Object <- callNextMethod(.Object,
                                        assayData = assayData,
                                        phenoData = phenoData,
                                        featureData = featureData,
                                        experimentData = experimentData,
                                        annotation = annotation,
                                        ...)
            } else stop("provide at most one of 'assayData' or c('exprs','vals') to initialize ACMECalcSet",
                        call.=FALSE)
            return(.Object)
          })
setValidity("ACMECalcSet", function(object) {
  msg <- NULL
  msg <- validMsg(msg, assayDataValidMembers(assayData(object), c("exprs","vals")))
  msg <- validMsg(msg, if(length(object@cutpoints)!=ncol(object))
                  return(c("cutpoints should be same length as number of samples")) else TRUE)
  msg <- validMsg(msg, if(length(object@threshold)!=1) return(c("threshold should be a single number"))
                  else return(TRUE)
                )
  if (is.null(msg)) TRUE else msg
})

### Chromosome
if(!isGeneric("chromosome")) {
  setGeneric("chromosome",function(object,...) standardGeneric("chromosome"))
}
setMethod("chromosome",c(object="ACMESet"),function(object,...) {
  return(fData(object)$chromosome)
})
### Start
if(!isGeneric("start")) {
  setGeneric("start",function(x,...) standardGeneric("start"))
}
setMethod("start",c(x="ACMESet"),function(x,...) {
  return(fData(x)$start)
})
### End
if(!isGeneric("end")) {
  setGeneric("end",function(x,...) standardGeneric("end"))
}
setMethod("end",c(x="ACMESet"),function(x,...) {
  return(fData(x)$end)
})
### cutpoints
if(!isGeneric("cutpoints")) {
  setGeneric("cutpoints",function(x,...) standardGeneric("cutpoints"))
}
setMethod("cutpoints",c(x="ACMECalcSet"),function(x,...) {
  return(x@cutpoints)
})
### threshold
if(!isGeneric("threshold")) {
  setGeneric("threshold",function(x,...) standardGeneric("threshold"))
}
setMethod("threshold",c(x="ACMECalcSet"),function(x,...) {
  return(x@threshold)
})
### vals
if(!isGeneric("vals")) {
  setGeneric("vals",function(x,...) standardGeneric("vals"))
}
setMethod("vals",c(x="ACMECalcSet"),function(x,...) {
  return(assayDataElement(x,'vals'))
})

### show for ACMECalcSet
setMethod("show","ACMECalcSet",
          function(object) {
            callNextMethod(object)
            cat(paste("Threshold:",threshold(object),"\n"))
          })

### Plotting
## ACMESet
if(!isGeneric("plot")) {
  setGeneric("plot",function(x,y,...) standardGeneric("plot"))
}
.plotACMESet <- function(x,y='missing',chrom,samples=NULL,...) {
  nsamps <- 0
  if (is.null(samples)) {
    samples <- 1:ncol(x)
    nsamps <- length(samples)
  }
  if (!is.numeric(samples))
    samples <- which(sampleNames(x) %in% samples)
  sub <- chromosome(x)==chrom
  if (sum(sub)==0)
    stop('No matching chromosome in the data')
  if (nsamps>1 & interactive()) par(ask=TRUE)
  for (i in samples)
    plot(start(x)[sub],exprs(x)[sub,i],...,
         main=paste('Chromosome:',chrom,', Sample:',sampleNames(x)[i]),
         ylab='Data',xlab='Chromosome Position')
  par(ask=FALSE)
}
setMethod("plot","ACMESet",.plotACMESet)
## ACMECalcSet
.plotACMECalcSet <- function(x,y='missing',chrom,samples=NULL,...) {
  nsamps <- 0
  if (is.null(samples)) {
    samples <- 1:ncol(x)
    nsamps <- length(samples)
  }
  if (!is.numeric(samples))
    samples <- which(sampleNames(x) %in% samples)
  sub <- chromosome(x)==chrom
  if (sum(sub)==0)
    stop('No matching chromosome in the data')
  if (nsamps>1 & interactive()) par(ask=TRUE)
  par(mar=c(4,5,4,5))
  for (i in samples) {
    plot(start(x)[sub],exprs(x)[sub,i],
         axes=F,col='gray85',pch=20,xlab="",ylab="",main="",...)
    axis(side=4)
    abline(h=0,col='gray85')
    abline(h=cutpoints(x)[i],col='gray85',lty=2)
    par(new=T)
    plot(start(x)[sub],-log10(vals(x)[sub,i]),
         main=paste('Chromosome:',chrom,', Sample:',sampleNames(x)[i]),
         ylab='-log10(p-value)',xlab='Chromosome Position',type='b',pch=20,
         col='red',...)
  }
  par(ask=FALSE)
}
setMethod("plot","ACMECalcSet",.plotACMECalcSet)


setClass('aGFF',
         representation(annotation='data.frame',
                        data='matrix',
                        samples='data.frame'),
         prototype=list(annotation=NULL,
           data=matrix(nrow=0,ncol=0),
           samples=NULL)
         )
setClass('aGFFCalc',
         representation(call='call',
                        threshold='numeric',
                        cutpoints='numeric',
                        vals='matrix'),
         contains='aGFF')
printgff <- function(x) {
  numgenes <- dim(x@data)[1]
  numsamples <- dim(x@data)[2]
  cat("\n")
  cat("Array GFF object")
  cat("\n")
  cat(paste("Number of Samples:",numsamples))
  cat("\n")
  cat(paste("Number of probes: ",numgenes))
  cat("\n")
  cat("\n")
  cat("===== Data =====")
  cat("\n")
  print(cbind(x@annotation[1:5,c('Chromosome','Location')],x@data[1:5,]))
  cat(paste("With",numgenes-5,"more rows..."))
  cat("\n")
  cat("\n")
  cat("=== Samples ===")
  cat("\n")
  if (is.null(x@samples)) print ('No sample information')
  else print(x@samples)
  cat("\n")
  cat("\n")
  cat("== Annotation ==")
  cat("\n")
  print(x@annotation[1:5,])
  cat(paste("With",numgenes-5,"more rows..."))
  cat("\n")
}
setMethod("print","aGFF",printgff)
setMethod("show","aGFF",function(object) print(object))
printgffcalc <- function(x) {
  numgenes <- dim(x@data)[1]
  numsamples <- dim(x@data)[2]
  cat("\n")
  cat("Array GFF Calculation object")
  cat("\n")
  cat(paste("Number of Samples:",numsamples))
  cat("\n")
  cat(paste("Threshold: ",x@threshold))
  cat("\n")
  cat("Cutpoints:\n")
  print(x@cutpoints)
  cat("\n")
  cat("Call:\n")
  print(x@call)
  cat("\n")
  cat("===== Data =====")
  cat("\n")
  print(cbind(x@annotation[1:5,c('Chromosome','Location')],x@data[1:5,]))
  cat(paste("With",numgenes-5,"more rows..."))
  cat("\n")
  cat("\n")
  cat("=== Samples ===")
  cat("\n")
  if (is.null(x@samples)) print ('No sample information')
  else print(x@samples)
  cat("\n")
  cat("\n")
  cat("== Annotation ==")
  cat("\n")
  print(x@annotation[1:5,])
  cat(paste("With",numgenes-5,"more rows..."))
  cat("\n")
  cat("\n")
  cat("== Values ==")
  cat("\n")
  print(x@vals[1:5,])
  cat(paste("With",numgenes-5,"more rows..."))
  cat("\n")
}
setMethod("print","aGFFCalc",printgffcalc)
setMethod("show","aGFFCalc",function(object) print(object))
plotgff <- function(x,y='missing',chrom,samples=NULL,ask=FALSE,...) {
  nsamps <- 0
  if (is.null(samples)) {
    samples <- 1:dim(x@data)[2]
    nsamps <- length(samples)
  }
  if (!is.numeric(samples))
    samples <- which(colnames(x@data) %in% samples)
  sub <- x@annotation$Chromosome==chrom
  if (sum(sub)==0)
    stop('No matching chromosome in the data')
  if (nsamps>1) par(ask=ask)
  for (i in samples)
    plot(x@annotation$Location[sub],x@data[sub,i],...,
         main=paste('Chromosome:',chrom,', Sample:',colnames(x@data)[i]),
         ylab='Data',xlab='Chromosome Position')
  par(ask=F)
}
setMethod("plot","aGFF",plotgff)
plotgffcalc <- function(x,y='missing',chrom,samples=NULL,...) {
  nsamps <- 0
  if (is.null(samples)) {
    samples <- 1:dim(x@data)[2]
    nsamps <- length(samples)
  }
  if (!is.numeric(samples))
    samples <- which(colnames(x@data) %in% samples)
  sub <- x@annotation$Chromosome==chrom
  if (sum(sub)==0)
    stop('No matching chromosome in the data')
  if (nsamps>1) par(ask=T)
  par(mar=c(4,5,4,5))
  for (i in samples) {
    plot(x@annotation$Location[sub],x@data[sub,i],
         axes=F,col='gray85',pch=20,xlab="",ylab="",main="",...)
    axis(side=4)
    abline(h=0,col='gray85')
    abline(h=x@cutpoints[i],col='gray85',lty=2)
    par(new=T)
    plot(x@annotation$Location[sub],-log10(x@vals[sub,i]),
         main=paste('Chromosome:',chrom,', Sample:',colnames(x@data)[i]),
         ylab='-log10(p-value)',xlab='Chromosome Position',type='b',pch=20,
         col='red',...)
  }
  par(ask=F)
}
setMethod("plot","aGFFCalc",plotgffcalc)
assign('[.aGFF',function(object,i,j) {
  if (nargs()!=3) stop("two subscripts required",call.=FALSE)
  if(missing(i)) {
    if(missing(j)) {
      return(object)
    } else {
      object@data <- object@data[,j,drop=FALSE]
      object@annotation <- object@annotation
      if (!is.null(object@samples)) object@samples <- object@samples[j,]
      return(object)
    }
  } else {
    if (missing(j)) {
      object@data <- object@data[i,,drop=FALSE]
      object@annotation <- object@annotation[i,,drop=FALSE]
      return(object)
    }
  }
  object@data <- object@data[i,j,drop=FALSE]
  object@annotation <- object@annotation[i,,drop=FALSE]
  if (!is.null(object@samples)) object@samples <- object@samples[j,]
  return(object)
}
       )
