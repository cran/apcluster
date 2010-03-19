\name{show-methods}
\docType{methods}
\alias{show}
\alias{show-methods}
\alias{show,APResult-method}
\title{Display APResult objects}
\description{
 Display method for S4 class \code{\linkS4class{APResult}}
}
\usage{
\S4method{show}{APResult}(object)
}
\arguments{
  \item{object}{An \code{\linkS4class{APResult}} object created by
                affinity propagation clustering}
}
\details{
  \code{show} displays the most important information stored in
  \code{object},
  namely number of data samples, number of iterations, input
  preference, final objective
  function values, exemplars, and clusters. For accessing more
  detailed information, it is necessary to directly the corresponding
  slots (see  \linkS4class{APResult}).
}
\value{
  \code{show} returns an invisible \code{NULL}
}
\author{Ulrich Bodenhofer & Andreas Kothmeier
\email{apcluster@bioinf.jku.at}}
\references{\url{http://www.bioinf.jku.at/software/apcluster}

Frey, B. J. and Dueck, D. (2007) Clustering by passing messages
between data points. \emph{Science} \bold{315}, 972-976.
}
\seealso{\code{\link{APResult-class}}, \code{\link{apcluster}}}
\examples{
## create two Gaussian clouds
cl1 <- cbind(rnorm(100,0.2,0.05),rnorm(100,0.8,0.06))
cl2 <- cbind(rnorm(50,0.7,0.08),rnorm(50,0.3,0.05))
x <- rbind(cl1,cl2)

## compute similarity matrix (negative squared Euclidean)
sim <- negDistMat(x, r=2)

## run affinity propagation
apres <- apcluster(sim)

## show details of clustering results
show(apres)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}
\keyword{methods}
