\name{apclusterK}
\alias{apclusterK}
\title{Affinity Propagation for Pre-defined Number of Clusters}
\description{
  Runs affinity propagation clustering for a given similarity matrix
  adjusting input preferences iteratively in order to adjust number of
  clusters to given value
}
\usage{
apclusterK(s, K, prc=10, bimaxit=20, exact=FALSE,
           nonoise=FALSE, seed=NA, ...)
}
\arguments{
  \item{s}{an \eqn{l\times l}{lxl} similarity matrix}
  \item{K}{desired number of clusters; must be between 2 and l-1}
  \item{prc}{the algorithm stops if the number of clusters does not
    deviate more than prc percent from desired value K; set to 0 if you
    want to have exactly K clusters}
  \item{bimaxit}{maximum number of bisection steps to perform; note that
    no warning is issued if if the number of clusters is still not in
    the desired range}
  \item{exact}{flag indicating whether or not to compute the initial
    preference range exactly (see \code{\link{preferenceRange}})}
  \item{nonoise}{\code{apcluster} adds a small amount of noise to
                 \code{s} to prevent degenerate cases; if \code{TRUE},
                 this is disabled}
  \item{seed}{for reproducibility, the seed of the random number
              generator can be set to a fixed value before
              adding noise (see above), if \code{NA}, the seed remains
              unchanged}
  \item{...}{all other arguments are passed to
             \code{\link{apcluster}} as they are}
}
\details{Affinity Propagation clusters data, using a set of
real-valued pairwise data point similarities as input. Clusters 
are each represented by a cluster center data point (the exemplar). 
The method is iterative and searches for clusters so as to maximize 
an objective function, called net similarity.

Apart from minor adaptations and optimizations, the implementation is
largely analogous to Frey's and Dueck's Matlab code
(see \url{http://www.psi.toronto.edu/affinitypropagation/}).
}
\value{
  Upon successful completion, the function returns a
  \code{\linkS4class{APResult}} object.
}
\author{Ulrich Bodenhofer & Andreas Kothmeier
\email{apcluster@bioinf.jku.at}}
\references{\url{http://www.bioinf.jku.at/software/apcluster}

Frey, B. J. and Dueck, D. (2007) Clustering by passing messages
between data points. \emph{Science} \bold{315}, 972-976.
}
\seealso{\code{\link{apcluster}}, \code{\link{preferenceRange}},
         \code{\link{APResult-class}}, \code{\link{APResult-methods}}}
\examples{
## create three Gaussian clouds
cl1 <- cbind(rnorm(70, 0.2, 0.05), rnorm(70, 0.8, 0.06))
cl2 <- cbind(rnorm(50, 0.7, 0.08), rnorm(50, 0.3, 0.05))
cl3 <- cbind(rnorm(60, 0.8, 0.04), rnorm(60, 0.8, 0.05))
x <- rbind(cl1, cl2, cl3)

## create similarity matrix
sim <- negDistMat(x, r=2)

## run affinity propagation such that 3 clusters are obtained
apres <- apclusterK(sim, 3)

## show details of clustering results
show(apres)

## plot clustering result
plot(apres, x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}

