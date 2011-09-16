\name{SMfootrule}
\alias{SMfootrule}
\title{Compute distance between two ranked lists.
}
\description{Compute the distances between the two input ranked lists.
}
\usage{
SMfootrule(R1, R2)
}
\arguments{
  \item{R1}{a ranked lists represented as a array, vector, data.frame, or matrix with one single column. 
}
  \item{R2}{with the same type of R1.
}
}
\details{Return a nonnegative number that represents the distances between the two input ranked lists using Spearman's algorithm. The two input ranked lists must be of the same length, which is used to measure the similarity of two ranked lists. In this package, this function can be used to compute distances between ranked lists which obtained for each gene expression profile by sorting the microarray probe-set identifiers according to the expression ratios (in decreasing order) with respect to the untreated hybridization.
}
\references{Diaconis, R.L Graham. (1977) \emph{ Spearman's footrule as a matter of disarray}
}
\seealso{\code{\link{FootruleMatrix}}
}
