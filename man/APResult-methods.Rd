\name{APResult-methods}
\alias{APResult-methods}
\alias{plot,APResult,missing-method}
\alias{plot,APResult,matrix-method}
\alias{show,APResult-method}
\title{Functions for APResult class}
\description{
  S4 methods for \code{\linkS4class{APResult}} class
}
\usage{
\S4method{show}{APResult}(object)
\S4method{plot}{APResult,missing}(x,
             type=c("netsim", "dpsim", "expref"),
             xlab="# Iterations",
             ylab="Similarity", ...)
\S4method{plot}{APResult,matrix}(x, y,
              connect=TRUE, xlab="", ylab="", ...)
}
\arguments{
  \item{object}{An \code{\linkS4class{APResult}} object created by
                affinity propagation clustering}
  \item{x}{an \code{\linkS4class{APResult}} object created by
           affinity propagation clustering}
  \item{y}{a matrix with the original data (must be 2D to work)}
  \item{type}{a string or array of strings indicating which
              performance measures should be plotted; valid values are
              \code{"netsim"}, \code{"dpsim"}, and \code{"expref"}
              which can be used in any combination or order;
              all other strings
              are ignored (for the meaning see
              \code{\linkS4class{APResult}})}
  \item{xlab}{label for x axis of plot}
  \item{ylab}{label for y axis of plot}
  \item{connect}{if \code{TRUE}, lines are drawn that connect every
                 exemplar with all its cluster members}
  \item{...}{all other arguments are passed to the standard
             \code{\link[graphics:plot]{plot}} command that is called
             internally to display the graphics window}
}
\details{
  \code{show} displays the most important information stored in
  \code{object},
  namely number of data samples, number of iterations, input
  preference, final objective
  function values, exemplars, and clusters. For accessing more
  detailed information, it is necessary to directly the corresponding
  slots (see  \linkS4class{APResult}).

  Provided that \code{\link{apcluster}} was called with
  \code{details=TRUE}, \code{plot} without a \code{y} argument
  displays graphs of performance
  measures over execution time.

  If \code{plot} is called with the second argument, \code{y}, being
  the original data matrix, a plot of clustering results superimposed
  on the original data set is created. This only works for 2D data,
  i.e. if \code{y} has 2 columns.
}
\value{
  \code{show} and both variants of \code{plot} return an
  invisible \code{NULL}
}
\author{Ulrich Bodenhofer & Andreas Kothmeier
\email{apcluster@bioinf.jku.at}}
\references{\url{http://www.bioinf.jku.at/software/apcluster}

Frey, B. J. and Dueck, D. (2007) Clustering by passing messages
between data points. \emph{Science} \bold{315}, 972-976.
}
\seealso{\link{APResult-methods}}
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}
\keyword{methods}
