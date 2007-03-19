\name{findClosestGene}
\alias{findClosestGene}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Find closest refseq gene}
\description{
  This function is used to find the nearest refseq transcript(s) to a
  point in the genome specified.  Note that it is limited to the refseq
  transcripts listed at genome.ucsc.edu, where this function goes for
  information.  
}
\usage{
findClosestGene(chrom, pos, genome = "hg17", position = "txStart")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chrom}{ Usually specified like 'chr1', 'chr2', etc. }
  \item{pos}{ A position in base pairs in the genome }
  \item{genome}{ Something like 'hg16', 'hg17', 'mm6', etc. }
  \item{position}{ The location to measure distance from: one of
    'txStart', 'txEnd', 'cdsStart', 'cdsEnd'}
}
\details{
  The first time the function is run, it checks to see if the refflat
  table for the given genome is present in the package environment.  If
  not, it downloads it to the /tmp directory and gunzips it (using
  \link{getRefflat}.  It is then stored so that in future calls, there
  is no re-download required.
}
\value{
  A data frame with the gene name, refseq id(s), txStart, txEnd,
  cdsStart, cdsEnd, exon count, and distance.  Note that distance is
  measured as pos-position, so negative values mean that the point in
  the gene is to the left of the point specified in the function call
  (with the p-tel on the left).
}
\author{Sean Davis <sdavis2@mail.nih.gov>}
\note{ The function may return more than one transcript, as several
  transcripts may have the same start site }

\examples{
findClosestGene('chr1',100000000,'hg17')
}
\keyword{ IO }
