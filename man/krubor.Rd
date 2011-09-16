\name{krubor}
\alias{krubor}
\title{Aggregate all ranked lists into one list
}
\description{
Return a matrix with one column representing all the input ranked lists, get a single Prototype Ranked List(PRL)
}
\usage{
krubor(...)
}
\arguments{
  \item{\dots}{column vectors,matrices or data.frames. These can be given as named arguments. The mode of arguments must be numeric.
}

}

\details{This function is aim to aggregate all ranked lists with the same state into one single ranked list. First, remove the duplicate columns. If there are the same columns in combination, delete the same columns until only one of them left. Second, aggregate the lists with the same state using the Borda Merging Method until only one single list left.
 
The arguments can be a mix of matrices, vectors and data.frames. The length of the column of the matrices or data.frames and the length of the vectors must be equal.
}
\value{A matrix with one column as the aggregated list.
}

\seealso{\code{\link{aggregate}}which uses krubor to aggregate ranked lists according to the biological states.}
\examples{
## the inputs are in the same class
krubor(matrix(2,3,3),matrix(3,3,3))

## the inputs are mixed
krubor(matrix(2,3,3),as.data.frame(matrix(3,3,3)))
}