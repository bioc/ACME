\name{aGFF-class}
\docType{class}
\alias{aGFF-class}
\alias{plot,aGFF-method}
\alias{print,aGFF-method}
\alias{show,aGFF-method}

\title{Class for storing GFF-like data}
\description{The GFF format is quite versatile while remaining simple.
  This class simply stores the annotation associated with a set of GFF
  files from the same regions of the genome along with some information
  about the samples from which the data came and the data (from the
  "score" column of the GFF file) themselves.}

\section{Objects from the Class}{Objects can be created by calls of the
  form \code{new("aGFF", ...)}.  Also, the \code{read.resultsGFF()}
  function returns aGFF objects.}

\section{Slots}{
	 \describe{
    \item{\code{annotation}:}{Object of class \code{"data.frame"} with
      two columns absolutely necessary, "Chromosome" and "Location".
      Other columns can be included.}
    \item{\code{data}:}{Object of class \code{"matrix"} of the same
      number of rows as the annotation slot and the same number of
      columns as the number of rows in the samples slot, containing
      data for later analysis}
    \item{\code{samples}:}{Object of class \code{"data.frame"} for
      describing the samples, one row per sample}
  }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "aGFF")}: to plot a region along the
    genome.}
    \item{print}{\code{signature(x = "aGFF")}: simple method to display
      summary of aGFF object}
    \item{show}{\code{signature(object = "aGFF")}: simple method to display
      summary of aGFF object}
	 }
}
\author{Sean Davis}
\seealso{
  \code{\link{read.resultsGFF}} and\code{\link{aGFFCalc-class}}
}
\examples{
# Load an example
data(example.agff)
example.agff
}
\keyword{classes}
