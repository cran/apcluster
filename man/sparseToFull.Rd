\name{sparseToFull}
\alias{sparseToFull}
\title{Convert a Sparse Similarity Matrix to Standard Representation}
\description{
   Converts a similarity matrix in sparse format to a full l by l matrix
}
\usage{
sparseToFull(s, fill=-Inf)
}
\arguments{
  \item{s}{a similarity matrix in sparse format; s must have exactly
           three columns; the matrix
           interpreted as follows: assume that a given row contains an
           integer index i in the 1st column and an integer index j
           in the 2nd column. Then the value in the 3rd column is
           interpreted as the similarity of data point i and data
           point j.}
  \item{fill}{value to fill in for entries not specified in s}
}
\details{Frey's and Dueck's implementation can handle
  \eqn{l\times l}{l x l} similarity
  matrices and similarity matrices of the sparse format described above.
  Presently our implementation does not make any use of the sparsity of
  the similarity matrix. The functions \code{\link{apcluster}},
  \code{\link{apclusterK}}, and \code{\link{preferenceRange}},
  therefore, only accept square similarity matrices. This function
  is a utility for converting sparse similarity matrices - wherever
  they may come from. 

  Note that most of Frey's and Dueck's Matlab functions make an internal
  conversion. Guessing the right format from the matrix's dimension may
  be wrong for \eqn{3\times 3}{3x3} matrices.
  We avoid this problem by explicitly
  distinguishing between sparse and dense similarity matrices (although
  we admit that clustering of 3 values makes as little sense as
  considering similarity matrices with only three entries).
}
\value{returns a square similarity matrix}
\author{Ulrich Bodenhofer & Andreas Kothmeier
\email{apcluster@bioinf.jku.at}}
\references{\url{http://www.bioinf.jku.at/software/apcluster}

Frey, B. J. and Dueck, D. (2007) Clustering by passing messages
between data points. \emph{Science} \bold{315}, 972-976.
}
\examples{
## create similarity matrix in sparse format
sp <- matrix(c(1, 1, 0.5, 3, 1, 0.2, 5, 4, -0.2, 3, 3, 1.2), 4, 3, byrow=TRUE)
sp

## perform conversion
s <- sparseToFull(sp, fill=0)
s
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}

