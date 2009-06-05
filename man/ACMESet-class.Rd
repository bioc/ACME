\name{ACMESet-class}
\Rdversion{1.1}
\docType{class}
\alias{ACMESet-class}
\alias{chromosome,ACMESet-method}
\alias{end,ACMESet-method}
\alias{plot,ACMESet-method}
\alias{start,ACMESet-method}
\title{Class "ACMESet" }
\description{An extension of ExpressionSet to deal with ACME data
  including chromosome locations}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ACMESet",
  assayData, phenoData, featureData, experimentData, annotation, exprs,
  ...)}.  The exprs assayDataElement stores the data.  The featureData
slot stores the chromosome location.  In practice, the data.frame
underlying the featureData MUST contain three columns named chromosome,
start, and end; this is enforced by the class validity method. 
}
\section{Slots}{
	 \describe{
    \item{\code{assayData}:}{Object of class \code{"AssayData"}.  See
      \linkS4class{ExpressionSet} for details. }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"} See
      \linkS4class{ExpressionSet} for details.}
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"} See
      \linkS4class{ExpressionSet} for details.}
    \item{\code{experimentData}:}{Object of class \code{"MIAME"} See
      \linkS4class{ExpressionSet} for details.}
    \item{\code{annotation}:}{Object of class \code{"character"} See
      \linkS4class{ExpressionSet} for details.}
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"} See
      \linkS4class{ExpressionSet} for details.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{ExpressionSet}"}, directly.
Class \code{"\linkS4class{eSet}"}, by class "ExpressionSet", distance 2.
Class \code{"\linkS4class{VersionedBiobase}"}, by class "ExpressionSet", distance 3.
Class \code{"\linkS4class{Versioned}"}, by class "ExpressionSet", distance 4.
}
\section{Methods}{
  \describe{
    \item{chromosome}{\code{signature(object = "ACMESet")}: Accessor for
    the chromosome.  Returns a vector of chromosomes.}
    \item{end}{\code{signature(x = "ACMESet")}: Accessor for the end
      location for a probe.  If that is not known, this could be set to
      the same value as the start location.  }
    \item{plot}{\code{signature(x = "ACMESet")}: A convenience plotting
      method that takes a sample name and chrom as well. }
    \item{start}{\code{signature(x = "ACMESet")}: Accessor for the start
    location for a probe.  }
	 }
}
\author{ Sean Davis <sdavis2@mail.nih.gov>}
\seealso{
	\code{\linkS4class{ExpressionSet}}, \code{\linkS4class{ACMECalcSet}} 

}
\examples{
showClass("ACMESet")
data(example.agff)
example.agff
head(chromosome(example.agff))
head(start(example.agff))
head(end(example.agff))
}
\keyword{classes}
