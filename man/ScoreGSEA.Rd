% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ScoreGSEA.R
\name{ScoreGSEA}
\alias{ScoreGSEA}
\title{Compute pairwise distances between samples with method in package GSEA}
\usage{
ScoreGSEA(
  MergingSet,
  SignatureLength,
  ScoringDistance = c("avg", "max"),
  p.value = FALSE
)
}
\arguments{
\item{MergingSet}{an \code{\link[Biobase:class.ExpressionSet]{Biobase::ExpressionSet}} object. The assay data
represents the PRLs of the samples, each column represents one PRL.
The number of sample of this argument must be greater than 1, otherwise,
this function is not meaningful.}

\item{SignatureLength}{the length of "gene signature". In order to compute
pairwise distances among samples, genes lists are ranked according to the
gene expression ratio (fold change). And the "gene signature" includes the
most up-regulated genes (near the top of the list) and the most
down-regulated genes (near the bottom of the list).}

\item{ScoringDistance}{the distance measurements between PRLs: the Average
Enrichment Score Distance (:avg"), and the Maximum Enrichment Score
Distance ("max").}

\item{p.value}{logical, if \code{TRUE} return a matrix of p.values of the
distance matrix, default \code{FALSE}.}
}
\value{
an distance-matrix, the max distance is more sensitive to weak
similarities, providing a lower precision but a larger recall.If \code{p.value}
is set to \code{TRUE}, then a list is returned that consists of the distance
matrix  as well as their p.values, otherwise, without p.values in the
result.
}
\description{
Compute pairwise distances between sample according to their (Prototype
Ranked List) PRL,  a \verb{N x N} distance matrix is generated by calling this
function, \code{N} is the length of PRL.
}
\details{
Once the PRL obtained for each sample, the distances between samples
are calculated base on gene signature, including the expression of genes
that seemed to consistently vary in response to the across different
experimental conditions (e.g., different cell lines and different dosages).
We take two distance measurements between PRLs: the Average
Enrichment-Score Distance \verb{Davg = (TES\{x,y\} + TES\{y,x\}) / 2}, and the Maximum
Enrichment-Score Distance \verb{Dmax = Min(TES\{x,y\}, TES\{y,x\}) / 2}.The avg is more
stringent than max, where max is more sensitive to weak similarities, with
lower precision but large recall.
}
\examples{
# load the sample expressionSet
data(exampleSet)
# Merging each group of the ranked lists in the exampleSet with the same
# phenotypic data into a single PRL
MergingSet <- RankMerging(exampleSet,"Spearman")
# get the distance matrix
ds <- ScoreGSEA(MergingSet, 250, "avg")
}
\seealso{
\code{\link[=ScorePGSEA]{ScorePGSEA()}},\code{\link[=SignatureDistance]{SignatureDistance()}}
}
