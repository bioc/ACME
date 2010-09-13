\name{findRegions}
\alias{findRegions}
\title{Find all regions in data above p-value threshold}
\description{
  After the ACME calculation, each probe is associated with a
  p-value of enrichment.  However, one often wants the contiguous
  regions associated with runs of p-values above a given p-value
  threshold.
}
\usage{
findRegions(x, thresh = 1e-04)
}
\arguments{
  \item{x}{An \code{ACMESetCalc} object}
  \item{thresh}{The p-value threshold}
}
\details{
  Runs of p-values above the p-value threshold will be reported as one
  "region".  These can be used for downstream analyses, export to
  browsers, submitted for transcription factor binding enrichment
  analyses, etc.
}
\value{
  A data frame with these columns:
  \item{Length}{The length of the region in probes}
  \item{TF}{Either TRUE or FALSE; TRUE regions represent regions of
    enrichment while FALSE regions are the regions between the TRUE regions}
  \item{StartInd}{The starting Index of the region}
  \item{EndInd}{The ending Index of the region}
  \item{Sample}{The sample containing the region}
  \item{Chromosome}{The Chromosome of the region}
  \item{Start}{The starting basepairof the region}
  \item{End}{The ending basepair of the region}
  \item{Median}{The median p-value in the region}
  \item{Mean}{The mean p-value in the region}
}
\author{Sean Davis <sdavis2@mail.nih.gov>}
\seealso{ \code{\link{do.aGFF.calc}}, \code{\link{findClosestGene}} }
\examples{
data(example.agff)
example.agffcalc <- do.aGFF.calc(example.agff,window=1000,thresh=0.9)
foundregions <- findRegions(example.agffcalc,thresh=0.001)
foundregions[1:6,]
}
\keyword{ manip }
