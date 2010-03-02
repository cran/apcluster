\name{similarity}
\alias{negDistMat}
\alias{expSimMat}
\alias{linSimMat}
\alias{linKernel}
\title{Methods for computing similarity matrices}
\description{
   Compute similarity matrices from data set
}
\usage{
negDistMat(x, r=1, ...)
expSimMat(x, r=2, w=1, ...)
linSimMat(x, w=1, ...)
linKernel(x, normalize=FALSE)
}
\arguments{
  \item{x}{real-valued data matrix, every row is a sample, every
           column a feature/input dimension}
  \item{r}{exponent (see details below)}
  \item{w}{radius (see details below)}
  \item{normalize}{see details below}
  \item{...}{all other arguments are passed to
             \code{\link[stats:dist]{dist}} as they are;
             the default distance is \code{method="euclidean"},
             see \code{\link[stats:dist]{dist}} for other options.}
}
\details{\code{negDistMat} creates a square matrix of mutual
pairwise similarities of data vectors as negative distances. The
argument \code{r} (default is 1) is used to transform the resulting
distances by computing the r-th power (use \code{r=2} to obtain
negative squared distances as in Frey's and Dueck's demos), i.e.,
given a distance d, the resulting similarity is computed as
\eqn{s=-d^r}. Internally, the computation of distances is done using
\code{\link[stats:dist]{dist}}. All options of this function
except \code{diag} and \code{upper} can be used, especially
\code{method} that allows for selecting different distance measures.

\code{expSimMat} computes similarities in a way similar to
\code{negDistMat}, but the transformation of distances to similarities
is done in the following way:
\deqn{s=\exp\left(-\left(\frac{d}{w}\right)^r\right)}{s=exp(-(d/w)^r)}
As above, r is an exponent. The parameter w controls the speed of
descent. \code{r=2} in conjunction with Euclidean distances
corresponds to the well-known Gaussian/RBF kernel, 
whereas \code{r=1} corresponds to the Laplace kernel. Note that these
similarity measures can also be understood as fuzzy equality relations.

\code{linSimMat} provides another way of transforming distances
into similarities by applying the following transformation to a
distance d:
\deqn{s=\max\left(0,1-\frac{d}{w}\right)}{s=max(0,1-d/w)}
Here, w corresponds to a maximal radius of interest. Note that this
is a fuzzy equality relation with respect to the Lukasiewicz t-norm.

Unlike the above three functions, \code{linKernel} computes pairwise
similarities as scalar products of data vectors, i.e. it corresponds,
as the name suggests, to the \dQuote{linear kernel}. If
\code{normalize=TRUE}, the values are scaled to the unit sphere in the
following way (for two samples \code{x} and \code{y}:
\deqn{s=\frac{\vec{x}^T\vec{y}}{\|\vec{x}\| \|\vec{y}\|}}{s=(x^T y)/(|x| |y|)}
}
\value{
   all functions listed above return square matrices of similarities
}
\author{Ulrich Bodenhofer & Andreas Kothmeier
\email{apcluster@bioinf.jku.at}}
\references{\url{http://www.bioinf.jku.at/software/apcluster}

Frey, B. J. and Dueck, D. (2007) Clustering by passing messages
between data points. \emph{Science} \bold{315}, 972-976.

Micchelli, C. A. (1986) Interpolation of scattered data: distance
matrices and conditionally positive definite functions.
\emph{Constr. Approx.} \bold{2}, 11-20.

De Baets, B. and Mesiar, R. (1997) Pseudo-metrics and T-equivalences.
\emph{J. Fuzzy Math.} \bold{5}, 471-481.
}
\seealso{\code{\link[stats:dist]{dist}},
\code{\link{apcluster}}}
\examples{
## create two Gaussian clouds
cl1 <- cbind(rnorm(100,0.2,0.05),rnorm(100,0.8,0.06))
cl2 <- cbind(rnorm(50,0.7,0.08),rnorm(50,0.3,0.05))
x <- rbind(cl1,cl2)

## create negative distance matrix (default Euclidean)
sim1 <- negDistMat(x)

## compute similarities as squared negative distances
## (in accordance with Frey's and Dueck's demos)
sim2 <- negDistMat(x, r=2)

## compute RBF kernel
sim3 <- expSimMat(x, r=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}

