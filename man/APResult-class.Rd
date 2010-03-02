\name{APResult-class}
\docType{class}
\alias{APResult-class}
\alias{APResult}
\alias{apresult}

\title{Class "APResult"}
\description{S4 class for storing results of affinity propagation
clustering}
\section{Objects from the Class}{
Objects of this class can be created by calling \code{\link{apcluster}}
for a given similarity matrix.
}
\section{Slots}{
The following slots are defined for \link{APResult} objects. Most names
are taken from Frey's and Dueck's original Matlab package:
    \describe{
    \item{\code{l}:}{number of samples in the data set}
    \item{\code{it}:}{number of iterations the algorithm ran}
    \item{\code{p}:}{input preference (either set by user or
                     computed by \code{\link{apcluster}})}
    \item{\code{netsim}:}{final total net similarity, defined as the
                          sum of \code{expref} and \code{dpsim}
                          (see below)}
    \item{\code{dpsim}:}{final sum of similarities of data points to
                         exemplars}
    \item{\code{expref}:}{final sum of preferences of the identified
                          exemplars}
    \item{\code{exemplars}:}{array containing indices of exemplars}
    \item{\code{clusters}:}{list containing the clusters; the i-th
                            component is a array of indices of
                            data points belonging to the i-th
                            exemplar (including the exemplar itself)}
    \item{\code{idx}:}{array of length \code{n} realizing a
                       sample-to-exemplar mapping; the i-th entry
                       contains the index of the exemplar the i-th
                       sample belongs to}
    \item{\code{netsimAll}:}{array containing the total net similarity
                             for each iteration; only available if
                             \code{\link{apcluster}} was called with
                             \code{details=TRUE}}
    \item{\code{dpsimAll}:}{array containing the sum of similarities
                            of data points to exemplars
                            for each iteration; only available if
                            \code{\link{apcluster}} was called with
                            \code{details=TRUE}}
    \item{\code{exprefAll}:}{array containing the sum of preferences
                             of the identified exemplars
                             for each iteration; only available if
                             \code{\link{apcluster}} was called with
                             \code{details=TRUE}}
    \item{\code{idxAll}:}{matrix with sample-to-exemplar indices
                          for each iteration; only available if
                          \code{\link{apcluster}} was called with
                          \code{details=TRUE}}
  }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "APResult")}: see
       \code{\link{APResult-methods}}}
    \item{plot}{\code{signature(x = "APResult", y = "matrix")}: see
       \code{\link{APResult-methods}}}
    \item{show}{\code{signature(object = "APResult")}: see
       \code{\link{APResult-methods}}}}
}
\author{Ulrich Bodenhofer & Andreas Kothmeier
\email{apcluster@bioinf.jku.at}}
\references{\url{http://www.bioinf.jku.at/software/apcluster}

Frey, B. J. and Dueck, D. (2007) Clustering by passing messages
between data points. \emph{Science} \bold{315}, 972-976.
}
\seealso{\link{APResult-methods}, \link{apcluster}}
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
\keyword{classes}
