\name{findclosestrank}
\alias{findclosestrank}
\title{Find the closest ranks.
}
\description{Find the two closest ranks among ranks with the same state. Get the NO. of the two closest ranks.
}
\usage{
findclosestrank(SMDM)
}
\arguments{
  \item{SMDM}{a distance matrix or a data.frame represents the distances between any two ranked lists, which must be preprocessed before used (let the lower triangular part of the matrix is Inf).
}
}
\details{Get the distance matrix by using \code{\link{FootruleMatrix}} function. This function is used to find the two closest ranked lists to aggregate them into a new list.
}
\seealso{\code{\link{krubor}}, \code{\link{FootruleMatrix}}
}

