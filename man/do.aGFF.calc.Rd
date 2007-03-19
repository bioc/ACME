\name{do.aGFF.calc}
\alias{do.aGFF.calc}
\title{Perform ACME calculation}
\description{
  This function performs the moving window chi-square calculation.  It
  is written in C, so is quite fast.
}
\usage{
do.aGFF.calc(x, window, thresh)
}
\arguments{
  \item{x}{An \code{aGFF} class object}
  \item{window}{An integer value, representing the number of basepairs
    to include in the windowed chi-square calculation}
  \item{thresh}{The quantile of the data distribution for each sample
    that will be used to classify a probe as positive}
}
\details{
  A window size on the order of 2-3 times the average size of fragments
  from sonication, digestion, etc. and containing at least 8-10 probes
  is the recommended size.  Larger size windows are probably more
  sensitive, but obviously reduce the accuracy with which boundaries of
  signal can be called.

  A threshold of between 0.9 and 0.99 seems empirically to be adequate.
  If one plots the histogram of data values and there is an obvious
  better choice (such as a bimodal distribution, with one peak
  representing enrichment), a more data-driven approach may yield better
  results.  
}
\value{
  An object of class \code{aGFFCalc}
}
\author{Sean Davis <sdavis2@mail.nih.gov>}
\examples{
data(example.agff)
example.agffcalc <- do.aGFF.calc(example.agff,window=1000,thresh=0.9)
example.agffcalc
}
\keyword{htest}
