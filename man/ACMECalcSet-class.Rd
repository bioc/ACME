\name{ACMECalcSet-class}
\Rdversion{1.1}
\docType{class}
\alias{ACMECalcSet-class}
\alias{cutpoints,ACMECalcSet-method}
\alias{plot,ACMECalcSet-method}
\alias{show,ACMECalcSet-method}
\alias{threshold,ACMECalcSet-method}
\alias{vals,ACMECalcSet-method}

\title{Class "ACMECalcSet" }
\description{A subclass of \linkS4class{ACMESet} that can also store the
parameters and results of an ACME calculation}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ACMECalcSet",
  assayData, phenoData, featureData, experimentData, annotation,
  cutpoints, threshold, exprs, vals, ...)}.  In addition to the
constraints defined by \linkS4class{ACMESet}, this class can also hold
the results (in the assayDataElement \code{vals}) and the
\code{threshold} and \code{cutpoints} from an ACME do.aGFF.calc run}
\section{Slots}{
	 \describe{
    \item{\code{cutpoints}:}{Object of class \code{"numeric"} The values
    of the cutpoints used in an analysis by do.aGFF.calc, one per sample.}
    \item{\code{threshold}:}{Object of class \code{"numeric"} The
      threshold used in an analysis.}
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
Class \code{"\linkS4class{ACMESet}"}, directly.
Class \code{"\linkS4class{ExpressionSet}"}, by class "ACMESet", distance 2.
Class \code{"\linkS4class{eSet}"}, by class "ACMESet", distance 3.
Class \code{"\linkS4class{VersionedBiobase}"}, by class "ACMESet", distance 4.
Class \code{"\linkS4class{Versioned}"}, by class "ACMESet", distance 5.
}
\section{Methods}{
  \describe{
    \item{cutpoints}{\code{signature(x = "ACMECalcSet")}: A simple
      getter for the cutpoints. }
    \item{plot}{\code{signature(x = "ACMECalcSet")}: A convenience
      plotting method that also takes sample and chrom }
    \item{show}{\code{signature(object = "ACMECalcSet")}: A show method }
    \item{threshold}{\code{signature(x = "ACMECalcSet")}: A simple
      getter for the threshold}
    \item{vals}{\code{signature(x = "ACMECalcSet")}: an accessor for the
    p-values from a run of do.aGFF.calc.  Returns a matrix with samples
    in columns and probes in rows.}
	 }
}
\author{ Sean Davis <sdavis2@mail.nih.gov> }
\seealso{
  \code{\linkS4class{ACMESet}}
}
\examples{
showClass("ACMECalcSet")
data(example.agff)
b <- do.aGFF.calc(example.agff,thresh=0.95,window=1000)
b
head(vals(b))
threshold(b)
cutpoints(b)
}
\keyword{classes}
