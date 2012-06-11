\name{Pihur.RankMerging}
\alias{Pihur.RankMerging}
\title{Merging the ranker lists with the same labels of the biological states into a single list with the Pihur's method.
}
\description{
Merging the assay data according to phenotypic data of the input ExpressionSet.  Each group of the ranked lists with the same phenotypic data is aggregated into a single list, return it as an ExpressionSet object.
}
\usage{
Pihur.RankMerging(exprSet, ...)
}
\arguments{
  \item{exprSet}{an ExpressionSet object, each column of assay data represents a ranked list obtained by preprocessing the corresponding gene expression profile, and phenotypic data represents the short description (characteristics of gene expression profile, such as the drug type, the disease state) about the assay data.
}
  \item{\dots}{additional arguments can be passed to the internal procedures
}
}
\details{
The function has the same function with Iorio.RankMerging,both this two functions use either the Spearman footrule distance or the Kendall's tau to measure the "closeness" of any two ordered lists. If choose Kendall as distance, the effectiveness of this function is certainly limited by the size of the  merging problem. 
}
\seealso{\code{\link{Iorio.RankMerging}}, \code{\link{SignatureDistance}}}
\examples{
#load the sample expressionSet
data(exampleSet)

## Merging each group of the ranked lists in the exampleSet with the same phenotypic data into a single PRL
MergingSet=Pihur.RankMerging(exampleSet[1:12])

}




