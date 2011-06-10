\name{ch22Promoters}
\alias{ch22Promoters}
\title{ch22Promoters Data Set}
\description{Toy data set consisting of promoter regions of 150
  random genes from the human chromosome 22}

\usage{ch22Promoters}
\format{A character vector of length 150. Each string is
  a nucleotide sequence that corresponds to
  the promoter region of a gene from the human chromosome no. 22
  (according to the human genome assembly hg18). The sequences
  start 999 bases upstream of the transcription start site (TSS)
  and end with the TSS itself.
  The names attribute contains the RefSeq IDs of the genes.}
\examples{
## load data set
data(ch22Promoters)

## display names of first five samples
names(ch22Promoters)[1:5]

## display last 50 bases of the first five samples
substr(ch22Promoters[1:5], 951, 1000)
}
\keyword{datasets}

