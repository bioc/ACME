\name{read.resultsGFF}
\alias{read.resultsGFF}
\title{Read Nimblegen GFF files}
\description{
  A GFF format file is a quite flexible format for storing genomic
  data.  Nimblegen uses these format files as one format for making
  chip-chip data available.  This function reads these files, one per
  experiment and creates a resulting aGFF-class object.
}
\usage{
read.resultsGFF(fnames, path = NULL, samples = NULL, notes = NULL, skip = 0, sep = "\t", quote = "\"", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fnames}{A vector of filenames}
  \item{path}{The path to the filenames or \code{NULL} for relative to
    the local directory}
  \item{samples}{A data.frame containing sample information, one row per
    sample, in the same order as the files in fnames}
  \item{notes}{A character vector for notes--not currently stored}
  \item{skip}{Number of lines to skip if the file contains a header}
  \item{sep}{The field separator--should be a tab character for gff
    files, but can be set if necessary.}
  \item{quote}{The text quote character--again not used for gff file, typically}
  \item{\dots}{...}
}
\details{
  The output is an aGFF object.
}
\value{
  A single aGFF object.
}
\references{http://www.sanger.ac.uk/Software/formats/GFF/}
\author{Sean Davis <sdavis2@mail.nih.gov>}
\seealso{\code{\link{aGFF-class}}}
\examples{
datdir <- system.file('extdata',package='ACME')
fnames <- dir(datdir)
example.agff <- read.resultsGFF(fnames,path=datdir)
}
\keyword{IO}
\keyword{manip}
