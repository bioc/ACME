\name{generics}
\Rdversion{1.1}
\alias{vals}
\alias{chromosome}
\alias{threshold}
\alias{end}
\alias{plot}
\alias{start}
\alias{cutpoints}
\title{Generics defined within ACME
}
\description{
 See methods descriptions for details. 
}
\usage{
vals(x, ...)
chromosome(object, ...)
end(x, ...)
start(x, ...)
plot(x, y, ...)
cutpoints(x, ...)
threshold(x, ...)
}
\arguments{
  \item{x}{An \code{ACMESet} or \code{ACMECalcSet} object (for cutpoints
    and threshold)}
  \item{object}{An \code{ACMESet} or \code{ACMECalcSet} object (for cutpoints
    and threshold)}
  \item{y}{Treated as missing for plotting these types of objects}
  \item{\dots}{Passed into method
}
}
\details{
  These are all getters for \code{ACMESet} and \code{ACMECalcSet} objects.
}
\value{
  See methods descriptions for details
}
\author{
  Sean Davis <sdavis2@mail.nih.gov>
}
\seealso{
  \linkS4class{ACMESet},   \linkS4class{ACMECalcSet}
}
\examples{
data(example.agff)
head(chromosome(example.agff))
head(end(example.agff))
head(start(example.agff))
}
\keyword{ manip }
