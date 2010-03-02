\name{apcluster-package}
\alias{apcluster-package}
\docType{package}
\title{Affinity Propagation Package}
\description{
The apcluster package implements affinity propagation according to
Frey and Dueck. It further offers visualization of the algorithm's
progress and clustering results (limited to 2D data sets).
}
\details{
\tabular{ll}{
Package: \tab apcluster\cr
Type: \tab Package\cr
Version: \tab 1.0.1\cr
Date: \tab 2010-03-02\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}
The central function is \code{\link{apcluster}}. It runs affinity
propagation on a given similarity matrix. The function returns
an \code{\linkS4class{APResult}} object from which the clustering
itself and information about the affinity propagation run can be
obtained. Further functions are implemented to visualize the
results and to create distance matrices.
}
\author{Ulrich Bodenhofer & Andreas Kothmeier
\email{apcluster@bioinf.jku.at}}
\references{\url{http://www.bioinf.jku.at/software/apcluster}

Frey, B. J. and Dueck, D. (2007) Clustering by passing messages
between data points. \emph{Science} \bold{315}, 972-976.
}
\keyword{package}
\examples{
## create two Gaussian clouds
cl1 <- cbind(rnorm(100,0.2,0.05),rnorm(100,0.8,0.06))
cl2 <- cbind(rnorm(50,0.7,0.08),rnorm(50,0.3,0.05))
x <- rbind(cl1,cl2)

## compute similarity matrix (negative squared Euclidean)
sim <- negDistMat(x, r=2)

## run affinity propagation
apres <- apcluster(sim, details=TRUE)

## show details of clustering results
show(apres)

## plot information about clustering run
plot(apres)

## plot clustering result
plot(apres, x)
}
