\name{FootruleMatrix}
\alias{FootruleMatrix}
\title{Create a footrule-matrix
}
\description{Compute distances between any two ranked lists with the same length, and create a n-n matrix, where n is the length of the ranked lists.
}
\usage{
FootruleMatrix(Rankings, n)
}
\arguments{
  \item{Rankings}{a numeric matrix or a data.frame to be computed
}
  \item{n}{a number, a non zero value n means a normalized results matrix should be returned.
}
}
\details{This function uses \code{\link{SMfootrule}} to compute any two ranked lists in the first argument, a column represents a ranked list in the first argument. n is used to indicate whether the distance matrix is 
normalized. 
}
\value{This function returns a n-n numberic matrix as a distance matrix, where n is the number of column of the first argument. And the m-n elements in the result should be equal to the n-m elements in the result.
}

\seealso{\code{\link{SMfootrule}}, \code{\link{findclosestrank}}
}

