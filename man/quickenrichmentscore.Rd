\name{quickenrichmentscore}
\alias{quickenrichmentscore}
\title{Compute the Enrichment Score.
}
\description{Use Gene Set Enrichmentscore Analysis (GSEA) method to compute the Enrichment Score among PRLs .
}
\usage{
quickenrichmentscore(S, S1, List)
}
\arguments{
  \item{S}{optimal signature,the up-regulated genes of a PRL A
}
  \item{S1}{optimal signature, the down-regulated genes of a PRL A
}
  \item{List}{a PRL which distances from A (whose optimal signature is S and S1) will be computed.
}
}
\details{Once gene signature is given, we computed the GSEA enrichment score between ranked list A and B based on the gene signature of A, and vice versa. The average value of this two enrichment scores is used to quantify the distance between A and B. 
}
\references{Subramanian. (2005) \emph{Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles.}
}
\seealso{\code{\link{distances}}, \code{\link{integratePRL}}
}

