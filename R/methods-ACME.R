setMethod("initialize","ACME",
          function(.Object,
                   assayData,
                   Chromosome,
                   Position,
                   phenoData        = annotatedDataFrameFrom(assayData,byrow=FALSE),
                   experimentData   = new("MIAME"),
                   annotation       = character(),
                   exprs            = new("matrix"),
                   featureData,
                   ...) {
            if (missing(featureData)) {
              if (missing(Chromosome) || missing(Position)) {
                stop("One of featureData or Chromosome and Position must be specified")
              }
              featureData=as(data.frame(Chromosome=Chromosome,Position=Position),"AnnotatedDataFrame")
            }
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
                             ...)
            } else if (missing(exprs)) {
              .Object <- callNextMethod(.Object,
                             assayData = assayData,
                             phenoData = phenoData,
                             featureData = featureData,
                             experimentData = experimentData,
                             annotation = annotation,
                             ...)
            } else stop("provide at most one of 'assayData' or 'exprs' to initialize ExpressionSet",
                        call.=FALSE)
            validObject(.Object)
            .Object <- .Object[order(pData(featureData(.Object))$Chromosome,pData(featureData(.Object))$Position),]
            return(.Object)
          })

#.Object <- callNextMethod(.Object,
#                                      assayData = assayData,
#                                      phenoData = phenoData,
#                                      experimentData = experimentData,
#                                      annotation = annotation,
#                                      featureData = featureData)
#            validObject(.Object)
#            return(.Object)
#          })

setMethod("initialize","ACMECalc",
          function(.Object,
                   assayData,
                   Chromosome,
                   Position,
                   phenoData        = annotatedDataFrameFrom(assayData,byrow=FALSE),
                   experimentData   = new("MIAME"),
                   annotation       = character(),
                   exprs            = new("matrix"),
                   featureData,
                   windowsize,
                   threshold,
                   ...) {
            .Object@windowsize <- as.integer(windowsize)
            .Object@threshold <- as.numeric(threshold)
            browser()
            if (missing(featureData)) {
              if (missing(Chromosome) || missing(Position)) {
                stop("One of featureData or Chromosome and Position must be specified")
              }
              featureData=as(data.frame(Chromosome=Chromosome,Position=Position),"AnnotatedDataFrame")
            }
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
                             ...)
            } else if (missing(exprs)) {
              .Object <- callNextMethod(.Object,
                             assayData = assayData,
                             phenoData = phenoData,
                             featureData = featureData,
                             experimentData = experimentData,
                             annotation = annotation,
                             ...)
            } else stop("provide at most one of 'assayData' or 'exprs' to initialize ExpressionSet",
                        call.=FALSE)
            validObject(.Object)
            .Object <- .Object[order(pData(featureData(.Object))$Chromosome,pData(featureData(.Object))$Position),]
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
  if(!assayDataValidMembers(assayData(object),c("exprs"))) {
    stop('assayData does not contain an "exprs" member')
  }
  if(is.na(windowsize(object))) {
    stop('windowsize not specified')
  }
  if(is.na(threshold(object))) {
    stop('threshold not specified')
  }
#  if(any(is.na(match(c("Chromosome","Position"),varLabels(featureData(object)))))) {
#    stop('featureData must contain both a "Chromosome" and a "Position" column')
#  }
  return(TRUE)
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

