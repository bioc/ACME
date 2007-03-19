\name{getRefflat}
\alias{getRefflat}
\title{ Get the refflat table from ucsc for the given genome}
\description{
  Fetches the refflat table from ucsc, stores in temp dir and then
  gunzips it and reads it in.
}
\usage{
getRefflat(genome = "hg17")
}
\arguments{
  \item{genome}{ The genome code from ucsc, like 'hg16', 'mm6', etc.}
}
\value{
  A data frame mirroring the UCSC table structure.
}
\references{ http://genome.ucsc.edu }
\author{Sean Davis <sdavis2@mail.nih.gov>}
\seealso{\code{\link{findClosestGene}}}
\examples{
rf <- getRefflat('hg17')
}
\keyword{IO}
