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
