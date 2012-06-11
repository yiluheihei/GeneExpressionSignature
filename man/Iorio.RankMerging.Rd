\name{Iorio.RankMerging}
\alias{Iorio.RankMerging}
\title{
Merging the ranker lists with the same labels of the biological states into a single list with the Iorio's method.
}
\description{
Merging the assay data according to phenotypic data of the input ExpressionSet.  Each group of the ranked lists with the same phenotypic data is aggregated into a single list, return it as an ExpressionSet object.
}
\usage{
Iorio.RankMerging(exprSet, distance = c("Spearman", "Kendall"))
}
\arguments{
  \item{exprSet}{
an ExpressionSet object, each column of assay data represents a ranked list obtained by preprocessing the corresponding gene expression profile, and phenotypic data represents the short description (characteristics of gene expression profile, such as the drug type, the disease state) about the assay data.
}
  \item{distance}{
distance to be used which "measures" the similarity of ordered lists, the default is Spearman
}
}
\details{
The krubor function is used in the aggregating procedure. And the following methods are used in the implementation: a measure of the distance between two ranked lists (Spearman's Footrule), a method to merge two or more ranked lists the (Borda Merging Method), and a algorithm to obtain a single ranked list from a set of them in a hierarchical way (the Kruskal Algorithm). If choose Kendall as distance, the effectiveness of this function is certainly limited by the size of the  merging problem. 
}
\seealso{\code{\link{Pihur.RankMerging}}, \code{\link{SignatureDistance}}}
\examples{
#load the sample expressionSet
data(exampleSet)

## Merging each group of the ranked lists in the exampleSet with the same phenotypic data into a single PRL
MergingSet=Iorio.RankMerging(exampleSet,"Spearman")

}

