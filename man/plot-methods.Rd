\name{plot-methods}
\docType{methods}
\alias{plot}
\alias{plot-methods}
\alias{plot,APResult,missing-method}
\alias{plot,APResult,matrix-method}
\alias{plot,ExClust,matrix-method}
\alias{plot,AggExResult,missing-method}
\alias{plot,AggExResult,matrix-method}
\title{Plot Clustering Results}
\description{
  Functions for Visualizing Clustering Results
}
\usage{
\S4method{plot}{APResult,missing}(x,
             type=c("netsim", "dpsim", "expref"),
             xlab="# Iterations",
             ylab="Similarity", ...)
\S4method{plot}{APResult,matrix}(x, y,
              connect=TRUE, xlab="", ylab="", ...)
\S4method{plot}{ExClust,matrix}(x, y,
              connect=TRUE, xlab="", ylab="", ...)
\S4method{plot}{AggExResult,missing}(x, main="Cluster dendrogram",
              xlab="", ylab="Balanced avg. similarity to exemplar",
              ticks=4, digits=2, ...)
\S4method{plot}{AggExResult,matrix}(x, y, k=NA, h=NA, ...)
}
\arguments{
  \item{x}{a clustering result object of class
     \code{\linkS4class{APResult}}, \code{\linkS4class{ExClust}},
     or \code{\linkS4class{AggExResult}}}
  \item{y}{a matrix (see details below)}
  \item{type}{a string or array of strings indicating which
              performance measures should be plotted; valid values are
              \code{"netsim"}, \code{"dpsim"}, and \code{"expref"}
              which can be used in any combination or order;
              all other strings are ignored (for the meaning see
              \code{\linkS4class{APResult}})}
  \item{xlab}{label for x axis of plot}
  \item{ylab}{label for y axis of plot}
  \item{connect}{used only if clustering is plotted on original data,
    ignored otherwise. If \code{connect} is \code{TRUE}, lines are
    drawn connecting exemplars with their cluster members.}
  \item{main}{title of plot}
  \item{ticks}{number of ticks used for the axis on the left side of the
    plot (applies to dendrogram plots only, see below)}
  \item{digits}{number of digits used for the axis tickmarks on the left
    side of the plot (applies to dendrogram plots only, see below)}
  \item{k}{level to be selected when plotting a single clustering
    level of cluster hierarchy (i.e. the number of clusters; see
    \code{\link{cutree-methods}})}
  \item{h}{cut-off to be used  when plotting a single clustering
    level of cluster hierarchy (see \code{\link{cutree-methods}})}
  \item{...}{all other arguments are passed to the plotting command that
    are used internally, \code{\link[graphics:plot]{plot}} or
    \code{\link{heatmap}}.}
}
\details{
  If \code{plot} is called for an \code{\linkS4class{APResult}} object
  without specifying the second argument \code{y},
  a plot is created that displays graphs of performance
  measures over execution time of the affinity propagation run.
  This only works if \code{\link{apcluster}} was called with
  \code{details=TRUE}, otherwise an error message is displayed.
  This variant of \code{plot} does not return any value.

  If \code{plot} is called for an \code{\linkS4class{APResult}} object
  along with a matrix as argument \code{y}, then the dimensions of
  the matrix determine the behavior of \code{plot}:
  \enumerate{
    \item{If the matrix \code{y} has two columns, \code{y} is
      interpreted as the original data set. Then a plot of
      the clustering result superimposed on the original
      data set is created. Each cluster is displayed in a
      different color. The exemplar of each cluster is highlighted
      by a black square. If \code{connect} is \code{TRUE}, lines
      connecting the cluster members to their exemplars are drawn.
      This variant of \code{plot} does not return any value.}
    \item{If the matrix \code{y} is quadratic, \code{y} is
      interpreted as the similarity matrix that was used for
      clustering. Then a heatmap of the similarity matrix is drawn, where
      the clusters are grouped together and highlighted in
      different colors. The order of clusters is determined by
      running \code{\link{aggExCluster}} on the clustering result
      \code{x}. This variant of \code{plot} returns an invisible
      \code{\linkS4class{AggExResult}} object.}
  }
  
  If \code{plot} is called for an \code{\linkS4class{ExClust}} object
  along with a matrix as argument \code{y}, then \code{plot} behaves
  exactly the same as described in the previous paragraph.
  
  If \code{plot} is called for an \code{\linkS4class{AggExResult}} object
  without specifying the second argument \code{y}, then a dendrogram
  plot is drawn. This variant returns an invisible
  \code{\link{dendrogram}} object.

  If \code{plot} is called for an \code{\linkS4class{AggExResult}} object
  along with a square matrix \code{y}, \code{y} is again interpreted
  as similarity matrix. If \code{x@maxNoClusters} is less than
  \code{x@l}, \code{plot} produces the same kind of heatmap as
  the variant for \code{\linkS4class{APResult}} and
  \code{\linkS4class{ExClust}} objects described above -
  with clusters grouped and highlighted
  in different colors. If \code{x@maxNoClusters} is equal to
  \code{x@l}, a heatmap is shown too, but instead of a color coding of
  clusters, the complete dendrogram of the cluster hierarchy is
  shown on the top and the left of the heatmap. This variant of
  \code{plot} does not return any value.

  If \code{plot} is called for an \code{\linkS4class{AggExResult}} object
  along with a two-column matrix \code{y}, \code{y} is again interpreted
  as original data set. If one of the two arguments \code{k} or
  \code{h} is present, a clustering is cut out from the cluster hierarchy using
  \code{\link{cutree}} and this clustering is displayed with the
  original data set as described above. This variant of
  \code{plot} returns an invisible \code{\linkS4class{ExClust}} object
  containing the extracted clustering.
}
\value{
  see details above
}
\author{Ulrich Bodenhofer & Andreas Kothmeier
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
  \code{\link{apcluster}}, \code{\link{aggExCluster}},
  \code{\link{cutree-methods}}}
\examples{
## create two Gaussian clouds
cl1 <- cbind(rnorm(50,0.2,0.05),rnorm(50,0.8,0.06))
cl2 <- cbind(rnorm(50,0.7,0.08),rnorm(50,0.3,0.05))
x <- rbind(cl1,cl2)

## compute similarity matrix (negative squared Euclidean)
sim <- negDistMat(x, r=2)

## run affinity propagation
apres <- apcluster(sim, q=0.6, details=TRUE)

## plot information about clustering run
plot(apres)

## plot clustering result
plot(apres, x)

## plot heatmap
plot(apres, sim)

## perform agglomerative clustering of affinity propagation clusters
aggres1 <- aggExCluster(sim, apres)

## show dendrogram
plot(aggres1)

## show heatmap
plot(aggres1, sim)

## show clustering result for 4 clusters
plot(aggres1, x, k=4)

## perform agglomerative clustering of whole data set
aggres2 <- aggExCluster(sim)

## show dendrogram
plot(aggres2)

## show heatmap along with dendrogram
plot(aggres2, sim)

## show clustering result for 2 clusters
plot(aggres2, x, k=2)
}
\keyword{cluster}
\keyword{methods}
