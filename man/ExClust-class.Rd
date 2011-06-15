\name{ExClust-class}
\docType{class}
\alias{ExClust-class}
\alias{ExClust}
\alias{exclust}

\title{Class "ExClust"}
\description{S4 class for storing exemplar-based clusterings}
\section{Objects}{
  Objects of this class can be created by calling \code{\link{cutree}}
  to cut out a clustering level from a cluster hierarchy
  of class \code{\linkS4class{AggExResult}}. Moreover,
  \code{\link{cutree}} can also be used to convert an object of
  class \code{\linkS4class{APResult}} to class \code{ExClust}.
}
\section{Slots}{
The following slots are defined for \link{ExClust} objects:
    \describe{
    \item{\code{l}:}{number of samples in the data set}
    \item{\code{exemplars}:}{vector containing indices of exemplars}
    \item{\code{clusters}:}{list containing the clusters; the i-th
                            component is a vector of indices of
                            data points belonging to the i-th
                            exemplar (including the exemplar itself)}
    \item{\code{idx}:}{vector of length \code{l} realizing a
                       sample-to-exemplar mapping; the i-th entry
                       contains the index of the exemplar the i-th
                       sample belongs to}
  }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "ExClust", y = "matrix")}: see
       \code{\link{plot-methods}}}
    \item{show}{\code{signature(object = "ExClust")}: see
       \code{\link{show-methods}}}
    \item{labels}{\code{signature(object = "ExClust")}: see
       \code{\link{labels-methods}}}
    \item{cutree}{\code{signature(object = "ExClust", k="ANY", h="ANY")}: see
       \code{\link{cutree-methods}}}}
}
\author{Ulrich Bodenhofer & Andreas Kothmeier
\email{apcluster@bioinf.jku.at}}
\references{\url{http://www.bioinf.jku.at/software/apcluster}}
\seealso{\code{\link{aggExCluster}}, \code{\link{show-methods}},
  \code{\link{plot-methods}}, \code{\link{labels-methods}},
  \code{\link{cutree-methods}}, \code{\linkS4class{AggExResult}},
  \code{\linkS4class{APResult}}}
\examples{
## create two Gaussian clouds
cl1 <- cbind(rnorm(20,0.2,0.05),rnorm(20,0.8,0.06))
cl2 <- cbind(rnorm(25,0.7,0.08),rnorm(25,0.3,0.05))
x <- rbind(cl1,cl2)

## compute similarity matrix (negative squared Euclidean)
sim <- negDistMat(x, r=2)

## run affinity propagation
aggres <- aggExCluster(sim)

## extract level with two clusters
excl <- cutree(aggres, k=2)

## show details of clustering results
show(excl)

## plot information about clustering run
plot(excl, x)
}
\keyword{classes}
