\name{SignatureDistance}
\alias{SignatureDistance}
\title{
Compute pairwise distances comprehensively.
}
\description{
 This function integrated the function for rank merging and distance scoring, we can do the rank merging and distance scoring simply with it.
}
\usage{
SignatureDistance(exprSet, MergingMethod = c("Iorio", "Pihur"), distance = c("Spearman", "Kendall"), ScoringMethod = c("GSEA", "PGSEA"), SignatureLength, DisMeasurement=c("avg", "max"), ...)
}
\arguments{
  \item{exprSet}{
an ExpressionSet object, each column of assay data represents a ranked list obtained by preprocessing the corresponding gene expression profile, and phenotypic data represents the short description (characteristics of gene expression profile, such as the drug type, the disease state) about the assay data.
}
  \item{MergingMethod}{
method to be used to perform rank merging,Iorio or Pihur
}
  \item{distance}{
distance to be used which "measures" the similarity of ordered lists, Spearman or Kendall
}
  \item{ScoringMethod}{
method to be used to perform distance scoring, GSEA or PGSEA
}
  \item{SignatureLength}{
the length of "gene signature". In order to compute pairwise distances among samples, genes lists are ranked according to the gene expression ratio (fold change). And the "gene signature" includes the most up-regulated genes (near the top of the list) and the most down-regulated genes (near the bottom of the list).     
}
  \item{DisMeasurement}{
the distance measurements between PRLs: the Average Enrichment Score Distance (avg), and the Maximum Enrichment Score Distance (max).    
}
  \item{\dots}{
additional arguments can be passed to the internal procedures
}
}
\seealso{
\code{\link{Iorio.RankMerging}}, \code{\link{Pihur.RankMerging}}, \code{\link{ScoreGSEA}}, \code{\link{ScorePGSEA}}
}
\examples{
#load the sample expressionSet
data(exampleSet)

#distance scoring
ds=SignatureDistance(exampleSet,MergingMethod="Iorio", distance="Spearman", ScoringMethod="GSEA", SignatureLength=250,DisMeasurement="avg")
}

