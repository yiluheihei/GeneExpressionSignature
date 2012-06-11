\name{GeneExpressionSignature}
\alias{GeneExpressionSignature}
\title{Introduction to the GeneExpressionSignature Package}
\description{

The \pkg{GeneExpressionSignature} add-on is an implementation of computing distances among preprocessed gene-expression profiles of samples for \R. 
The distances can be used to detect similarities among the signatures of drugs, diseases, and biological states of interest, and construct connectivity map. 
}
\details{
This package contains functions for the distances computation based on gene expression signature. First, list of genes is ranked according to their  expression ratios to produce the Prototype Ranked List (PRL). Second, all the PRLs with the same state are aggregated by \code{\link{Iorio.RankMerging}} and \code{\link{Pihur.RankMerging}} functions. Finally, all the ranked lists are made as one input of the \code{\link{ScorePGSEA}} and \code{\link{ScoreGSEA}} functions to compute the pairwise distances. 
}


