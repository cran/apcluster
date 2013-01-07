\name{heatmap}
\docType{methods}
\alias{heatmap}
\alias{heatmap-methods}
\alias{heatmap,APResult,missing-method}
\alias{heatmap,APResult,matrix-method}
\alias{heatmap,ExClust,missing-method}
\alias{heatmap,ExClust,matrix-method}
\alias{heatmap,AggExResult,missing-method}
\alias{heatmap,AggExResult,matrix-method}
\alias{heatmap,missing,matrix-method}
\alias{heatmap,matrix,missing-method}
\title{Plot Heatmap}
\description{
  Functions for Plotting of Heatmap
}
\usage{
\S4method{heatmap}{APResult,missing}(x, y, ...)
\S4method{heatmap}{APResult,matrix}(x, y, ...)
\S4method{heatmap}{ExClust,missing}(x, y, ...)
\S4method{heatmap}{ExClust,matrix}(x, y, ...)
\S4method{heatmap}{AggExResult,missing}(x, y, ...)
\S4method{heatmap}{AggExResult,matrix}(x, y, ...)
\S4method{heatmap}{matrix,missing}(x, y, ...)
\S4method{heatmap}{missing,matrix}(x, y, ...)
}
\arguments{
  \item{x}{a clustering result object of class
     \code{\linkS4class{APResult}}, \code{\linkS4class{ExClust}},
     or \code{\linkS4class{AggExResult}}; for compatibility, \code{x}
     may also be a similarity matrix (see details below).}
  \item{y}{a similarity matrix}
  \item{...}{all other arguments are passed to the plotting command 
    \code{\link[stats:heatmap]{heatmap}}.}
}
\details{
  The \code{heatmap} functions provide plotting of heatmaps from several
  different types of input object. The implementation uses the standard
  graphics function \code{\link[stats:heatmap]{heatmap}}. 
  Plotting heatmaps via the \code{\link{plot}} command as available in previous 
  versions of this package is still available for backward 
  compatibility.

  If \code{heatmap} is called for objects of classes
  \code{\linkS4class{APResult}} or \code{\linkS4class{ExClust}},
  a heatmap of the similarity matrix in slot \code{sim} of the parameter
  \code{x} is created with clusters grouped together and highlighted in
  different colors. The order of clusters is determined by running 
  \code{\link{aggExCluster}} on the clustering result \code{x}. This
  variant of \code{heatmap} returns an  invisible
  \code{\linkS4class{AggExResult}} object.

  If \code{heatmap} is called for an \code{\linkS4class{AggExResult}}
  object that contains all levels of clustering, the heatmap is
  displayed with the corresponding clustering dendrogram. If the
  \code{\linkS4class{AggExResult}} object is the result of running
  \code{\link{aggExCluster}} on a prior clustering result, the same heatmap
  plot is produced as if \code{heatmap} had been called on this
  prior clustering result, however, returning the cluster hierarchy's
  \code{\link{dendrogram}}.

  All variants described above only work if the input object \code{x}
  contains a slot \code{sim} with the similarity matrix (which is only
  the case if the preceding clustering method has been called with
  \code{includeSim=TRUE}). In case the slot \code{sim} of \code{x} does not
  contain the similarity matrix, the similarity matrix must be supplied
  as second argument \code{y}.

  The two variants with one of the two arguments being a matrix and one
  being missing are just wrappers around the standard
  \code{\link[stats:heatmap]{heatmap}} function with the aim to provide
  compatibility with this standard case.
}
\value{
  see details above
}
\author{Ulrich Bodenhofer, Andreas Kothmeier & Johannes Palme
\email{apcluster@bioinf.jku.at}}
\references{\url{http://www.bioinf.jku.at/software/apcluster}

Bodenhofer, U., Kothmeier, A., and Hochreiter, S. (2011)
APCluster: an R package for affinity propagation clustering.
\emph{Bioinformatics} \bold{27}, 2463-2464.
DOI: \href{http://dx.doi.org/10.1093/bioinformatics/btr406}{10.1093/bioinformatics/btr406}.

Frey, B. J. and Dueck, D. (2007) Clustering by passing messages
between data points. \emph{Science} \bold{315}, 972-976.
DOI: \href{http://dx.doi.org/10.1126/science.1136800}{10.1126/science.1136800}.
}
\seealso{\code{\linkS4class{APResult}},
  \code{\linkS4class{AggExResult}}, \code{\linkS4class{ExClust}},
  \code{\link{apcluster}}, \code{\link{apclusterL}}, 
  \code{\link{aggExCluster}}, \code{\link{cutree-methods}},
  \code{\link{plot-methods}}}
\examples{
## create two Gaussian clouds
cl1 <- cbind(rnorm(50,0.2,0.05),rnorm(50,0.8,0.06))
cl2 <- cbind(rnorm(50,0.7,0.08),rnorm(50,0.3,0.05))
x <- rbind(cl1,cl2)

## run affinity propagation using negative squared Euclidean
apres <- apcluster(negDistMat(r=2), x, p=-0.1)

## plot heatmap clustering run
heatmap(apres)

## rerun affinity propagation
## reuse similarity matrix from previous run
apres2 <- apcluster(s=apres@sim, q=0.6)

## plot heatmap of second run
heatmap(apres2, apres@sim)

## perform agglomerative clustering of affinity propagation clusters
aggres1 <- aggExCluster(apres@sim, apres2)

## plot heatmap
heatmap(cutree(aggres1, 2), apres@sim)

## perform agglomerative clustering of whole data set
aggres2 <- aggExCluster(negDistMat(r=2), x)

## show heatmap along with dendrogram
heatmap(aggres2)
}
\keyword{cluster}
\keyword{methods}
