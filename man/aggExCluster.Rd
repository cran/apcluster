\name{aggExCluster}
\alias{aggExCluster}
\alias{aggexcluster}
\title{Exemplar-based Agglomerative Clustering}
\description{
  Runs exemplar-based agglomerative clustering for a given similarity matrix
}
\usage{
aggExCluster(s, cl)
}
\arguments{
  \item{s}{an \eqn{l\times l}{lxl} similarity matrix}
  \item{cl}{a prior clustering; if present, \code{cl} must be
    an object of class \code{\linkS4class{APResult}} or
    \code{\linkS4class{ExClust}}.}
}
\details{\code{aggExCluster} performs agglomerative clustering.
  Unlike other methods, e.g., the ones implemented in \code{\link{hclust}},
  \code{aggExCluster} is computing exemplars for each cluster and
  its merging objective is geared towards the identification of
  meaningful exemplars, too.

  For each pair of clusters, the merging objective is computed as
  follows:
  \enumerate{
    \item{An intermediate cluster is created as the union
      of the two clusters.}
    \item{The potential exemplar is selected from the intermediate
      cluster as the sample that has the largest average similarity
      to all other samples in the intermediate cluster.}
    \item{Then the average similarity of the exemplar with all
      samples in the first cluster and the average similarity with
      all samples in the second cluster is computed. These two values
      measure how well the joint exemplar describes the samples in the
      two clusters.}
    \item{The merging objective is finally computed as the
      average of the two measures above. Hence, we can consider the
      merging objective as some kind of \dQuote{balanced average
      similarity to the joint exemplar}.}
  }
  In each step, all pairs of clusters are considered and
  the pair with the largest merging objective is actually merged.
  The joint exemplar is then chosen as the exemplar of the merged
  cluster.

  \code{aggExCluster} can be used in two ways: Firstly, by calling it
  for a similarity matrix only (with \code{cl} argument missing).
  Then a full agglomeration run is performed that starts from \code{l}
  clusters (all samples in separate one-element clusters) and ends
  with one cluster (all samples in one single cluster). Secondly,
  \code{aggExCluster} can be called on a similarity matrix along with
  a prior clustering given by the \code{cl} argument. In this case,
  \code{aggExCluster} starts from the prior clustering. Consequently,
  a cluster hierarchy with numbers of clusters ranging from the
  number of clusters in \code{cl} down to 1 is created.

  The result is stored in an \code{\linkS4class{AggExResult}} object.
  If \code{aggExCluster} was called for a similarity matrix only
  (argument \code{cl} missing), \code{maxNoClusters} is equal to
  \code{l}, i.e. the result contains a complete
  cluster hierarchy from \code{l} clusters (every sample in a
  separate cluster) to one cluster (every sample in the same
  cluster). If \code{aggExCluster} was called for a similarity matrix
  along with a prior clustering (e.g. obtained from affinity
  propagation), then agglomeration starts from this clustering, hence,
  \code{maxNoClusters} is the number of clusters in the prior
  clustering. The slot \code{height} is filled with the merging
  objective of each of the \code{maxNoClusters-1} merges. The slot
  \code{order} contains a permutation of the samples/clusters for
  dendrogram plotting. The algorithm for computing this permutation
  is the same as the one used in \code{\link{hclust}}. If \code{aggExCluster}
  was called for a similarity matrix only, the slot \code{label}
  contains the names of the objects to be clustered (if available,
  otherwise the indices are used). If \code{aggExCluster} was called
  for a similarity matrix along with a prior clustering,
  then labels are set to \sQuote{Cluster 1}, \sQuote{Cluster 2}, etc.
}
\value{
  Upon successful completion, the function returns an
  \code{\linkS4class{AggExResult}} object.
}
\author{Ulrich Bodenhofer
\email{apcluster@bioinf.jku.at}}
\references{\url{http://www.bioinf.jku.at/software/apcluster}}
\seealso{\code{\linkS4class{AggExResult}},
  \code{\link{plot-methods}}, \code{\link{cutree-methods}}}
\examples{
## create two Gaussian clouds
cl1 <- cbind(rnorm(50,0.2,0.05),rnorm(50,0.8,0.06))
cl2 <- cbind(rnorm(50,0.7,0.08),rnorm(50,0.3,0.05))
x <- rbind(cl1,cl2)

## compute similarity matrix (negative squared Euclidean)
sim <- negDistMat(x, r=2)

## compute agglomerative clustering from scratch
aggres1 <- aggExCluster(sim)

## show results
show(aggres1)

## plot dendrogram
plot(aggres1)

## plot heatmap along with dendrogram
plot(aggres1, sim)

## plot level with two clusters
plot(aggres1, x, k=2)

## run affinity propagation
apres <- apcluster(sim, q=0.7)

## create hierarchy of clusters determined by affinity propagation
aggres2 <- aggExCluster(sim, apres)

## show results
show(aggres2)

## plot dendrogram
plot(aggres2)

## plot heatmap
plot(aggres2, sim)

## plot level with two clusters
plot(aggres2, x, k=2)
}
\keyword{cluster}

