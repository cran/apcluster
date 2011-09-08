\name{preferenceRange}
\alias{preferenceRange}
\title{Determine Meaningful Ranges for Input Preferences}
\description{
   Determines meaningful ranges for affinity propagation input preference
}
\usage{
preferenceRange(s, exact=FALSE)
}
\arguments{
  \item{s}{an \eqn{l\times l}{lxl} similarity matrix}
  \item{exact}{flag indicating whether exact ranges should be computed,
    which is relatively slow; if bounds are sufficient,
    supply \code{FALSE} (default)}
}
\details{Affinity Propagation clustering relies on an appropriate choice
  of input preferences. This function helps in finding a good choice by
  determining meaningful lower and upper bounds.
}
\value{
  returns a vector with two entries, the first of which is the minimal
  input preference (which would lead to 1 or 2 clusters) and the second
  of which is the maximal input prefence (which would lead to as many
  clusters as data samples)
}
\author{Ulrich Bodenhofer & Andreas Kothmeier
\email{apcluster@bioinf.jku.at}}
\references{\url{http://www.bioinf.jku.at/software/apcluster}

Frey, B. J. and Dueck, D. (2007) Clustering by passing messages
between data points. \emph{Science} \bold{315}, 972-976.
DOI: \href{http://dx.doi.org/10.1126/science.1136800}{10.1126/science.1136800}.

Bodenhofer, U., Kothmeier, A., and Hochreiter, S. (2011)
APCluster: an R package for affinity propagation clustering.
\emph{Bioinformatics} \bold{27}, 2463-2464.
DOI: \href{http://dx.doi.org/10.1093/bioinformatics/btr406}{10.1093/bioinformatics/btr406}.
}
\seealso{\code{\link{apcluster}}}
\examples{
## create two Gaussian clouds
cl1 <- cbind(rnorm(100,0.2,0.05),rnorm(100,0.8,0.06))
cl2 <- cbind(rnorm(50,0.7,0.08),rnorm(50,0.3,0.05))
x <- rbind(cl1,cl2)

## create similarity matrix
sim <- negDistMat(x, r=2)

## determine bounds
preferenceRange(sim)

## determine exact range
preferenceRange(sim, exact=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}

