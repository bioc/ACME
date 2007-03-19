\name{aGFFCalc-class}
\docType{class}
\alias{aGFFCalc-class}
\alias{plot,aGFFCalc-method}
\alias{print,aGFFCalc-method}
\alias{show,aGFFCalc-method}

\title{Class "aGFFCalc"}
\description{Store results of ACME calculations}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("aGFFCalc", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{call}:}{Object of class \code{"call"}, contains the
      exact call to do.aGFF.calc, for historical purposes}
    \item{\code{threshold}:}{Object of class \code{"numeric"}, the
      threshold used in the calculation}
    \item{\code{cutpoints}:}{Object of class \code{"numeric"}, the data
      value above which probes were considered positive}
    \item{\code{vals}:}{Object of class \code{"matrix"}, equivalent in
      size to the original data matrix, containing the calculated
      p-values from the ACME algorithm}
    \item{\code{annotation}:}{Object of class \code{"data.frame"},
      currently a copy of the original annotation, possibly reordered in
      chromosome order }
    \item{\code{data}:}{Object of class \code{"matrix"}, the original
      data, possibly reordered}
    \item{\code{samples}:}{Object of class \code{"data.frame"}, sample
      metadata}
  }
}
\section{Extends}{
Class \code{"aGFF"}, directly.
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "aGFFCalc", ask=FALSE)}: plot the results of an
    ACME calculation}
    \item{print}{\code{signature(x = "aGFFCalc")}: brief overview of the
    object}
    \item{show}{\code{signature(object = "aGFFCalc")}: brief overview of
    the object}
	 }
}
\author{Sean Davis <sdavis2@mail.nih.gov>}
\seealso{
  \code{\link{do.aGFF.calc}}, \code{\link{aGFF-class}}
}
\examples{
data(example.agff)
example.agffcalc <- do.aGFF.calc(example.agff,window=1000,thresh=0.9)
example.agffcalc
}
\keyword{classes}
