\name{write.sgr}
\alias{write.sgr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Write Affy IGB .sgr format files}
\description{
  The affy Integrated Genome Browser (IGB) is a powerful, fast browser
  for genomic data.  The file format is simple (three columns:
  chromosome, location, and score) to generate.  This function will
  write the sgr files associated with a aGFFcalc object.  There will be
  either one or two files (default two) representing the raw data and
  the calculated data (which is output as -log10(val) for visualization
  purposes).  
}
\usage{
write.sgr(agff, raw = TRUE, vals = TRUE, directory = ".")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{agff}{An aGFFCalc object obtained after running do.aGFF.calc}
  \item{raw}{Create a file for the raw data?}
  \item{vals}{Create a file for the calculated p-values?}
  \item{directory}{Give a directory for storing the files}
}
\author{Sean Davis}
\examples{
data(example.agff)
write.sgr(example.agff)
}
\keyword{IO}
