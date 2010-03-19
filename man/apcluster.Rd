\name{apcluster}
\alias{apcluster}
\alias{apclusterLM}
\title{Affinity Propagation}
\description{
   Runs affinity propagation clustering for a given similarity matrix
}
\usage{
apcluster(s, p=NA, q=NA, maxits=1000, convits=100, lam=0.9,
          details=FALSE, nonoise=FALSE, seed=NA)
apclusterLM(s, p=NA, q=NA, maxits=1000, convits=100, lam=0.9,
          details=FALSE, nonoise=FALSE, seed=NA)
}
\arguments{
  \item{s}{an \eqn{l\times l}{lxl} similarity matrix}
  \item{p}{input preference; can be a vector that specifies
           individual preferences for each data point. If scalar,
           the same value is used for all data points. If \code{NA},
           exemplar preferences are initialized according to the
           distribution of non-Inf values in \code{s}. How this
           is done is controlled by the parameter \code{q}.}
  \item{q}{If \code{p=NA}, exemplar preferences are initialized
           according to the distribution of non-Inf values in \code{s}.
           If \code{q=NA}, exemplar preferences are set to the median
           of non-Inf values in \code{s}. If \code{q} is a value
           between 0 and 1, the sample quantile with threshold
           \code{q} is used, whereas \code{q=0.5} again results in
           the median.}
  \item{maxits}{maximal number of iterations that should be executed}
  \item{convits}{the algorithm terminates if the examplars have not
                 changed for \code{convits} iterations}
  \item{lam}{damping factor; should be a value in the range [0.5, 1);
             higher values correspond to heavy damping which may be
             needed if oscillations occur}
  \item{details}{if \code{TRUE}, more detailed information about the
                 algorithm's progess is stored in the output object
                 (see \code{\linkS4class{APResult}})}
  \item{nonoise}{\code{apcluster} adds a small amount of noise to
                 \code{s} to prevent degenerate cases; if \code{TRUE},
                 this is disabled}
  \item{seed}{for reproducibility, the seed of the random number
              generator can be set to a fixed value before
              adding noise (see above), if \code{NA}, the seed remains
              unchanged}
}
\details{Affinity Propagation clusters data, using a set of
real-valued pairwise data point similarities as input. Clusters 
are each represented by a cluster center data point (the exemplar). 
The method is iterative and searches for clusters so as to maximize 
an objective function, called net similarity.

Apart from minor adaptations and optimizations, the
implementation of the function \code{apclusterLM} is 
largely analogous to Frey's and Dueck's Matlab code
(see \url{http://www.psi.toronto.edu/affinitypropagation/}). The
function \code{apcluster} uses the same ideas, but replaces the loops
in the computations of responsibilities and availabilities by pure
matrix operations. For moderate data sets, the variant \code{apcluster}
is approximately 60\% faster than \code{apclusterLM}. For large
data sets (several thousands of data samples), the use of
\code{apclusterLM} (LM = Less Memory) may be advantageous, since this
function requires less temporal storage (LM = Less Memory). For at most
5000 samples, we recommend to use \code{apcluster} (on up-to-date
systems that are not too tight with memory).

The new argument \code{q} allows for better controlling the number of
clusters without knowing the distribution of similarity
values. A meaningful range for the parameter \code{p} can be determined
using the function \code{\link{preferenceRange}}. Alternatively, a
certain fixed number of clusters may be desirable. For this purpose,
the function \code{\link{apclusterK}} is available.
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
\seealso{\code{\link{APResult-class}}, \code{\link{show-methods}},
  \code{\link{plot-methods}}, \code{\link{labels-methods}},
  \code{\link{preferenceRange}}, \code{\link{apclusterK}}}
\examples{
## create two Gaussian clouds
cl1 <- cbind(rnorm(100,0.2,0.05),rnorm(100,0.8,0.06))
cl2 <- cbind(rnorm(50,0.7,0.08),rnorm(50,0.3,0.05))
x <- rbind(cl1,cl2)

## create similarity matrix
sim <- negDistMat(x, r=2)

## run affinity propagation (p defaults to median of similarity)
apres <- apcluster(sim)

## show details of clustering results
show(apres)

## plot clustering result
plot(apres, x)

## run affinity propagation with default preference of 10% quantile
## of similarities; this should lead to a smaller number of clusters
apres <- apcluster(sim, q=0.1)
show(apres)
plot(apres, x)

## now try the same with RBF kernel
sim <- expSimMat(x, r=2)
apres <- apcluster(sim, q=0.2)
show(apres)
plot(apres, x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}

