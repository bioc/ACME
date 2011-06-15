\name{write.bedGraph}
\alias{write.bedGraph}
\title{Write bedGraph format tracks for UCSC genome browser}
\description{
  Generate bedGraph format files for the UCSC genome browser.  This function will
  write the bedGraph files associated with a aGFFcalc object.  There will be
  either one or two files (default two) representing the raw data and
  the calculated data (which is output as -log10(val) for visualization
  purposes for EACH sample).  
}
\usage{
write.bedGraph(x, raw = TRUE, vals = TRUE, directory = ".")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An ACMESet or ACMECalcSet object}
  \item{raw}{Boolean.  Create a file for the raw data?}
  \item{vals}{Boolean.  Create a file for the calculated p-values?}
  \item{directory}{Give a directory for storing the files}
}
\author{Sean Davis}
\examples{
data(example.agff)
write.bedGraph(example.agff)
}
\keyword{IO}
