#' Compute pairwise distances comprehensively
#' 
#' This function integrated the function for rank merging and distance scoring, 
#' we can do the rank merging and distance scoring simply with it.
#' 
#' @param exprSet an `ExpressionSet` object, each column of assay data represents 
#'   a ranked list obtained by preprocessing the corresponding gene expression 
#'   profile, and phenotypic data represents the short description 
#'   (characteristics of gene expression profile, such as the drug type, the 
#'   disease state) about the assay data.
#' @param SignatureLength the length of "gene signature". In order to compute 
#'   pairwise distances among samples, genes lists are ranked according to the 
#'   gene expression ratio (fold change). And the "gene signature" includes the 
#'   most up-regulated genes (near the top of the list) and the most 
#'   down-regulated genes (near the bottom of the list). 
#' @param MergingDistance distance to be used which "measures" the similarity 
#'   of ordered lists, "Spearman" or "Kendall".
#' @param ScoringMethod method to be used to perform distance scoring, "GSEA" 
#'   or "PGSEA".
#' @param ScoringDistance the distance measurements between PRLs: the Average 
#'   Enrichment Score Distance ("avg"), or the Maximum Enrichment Score Distance 
#'   ("max"). 
#' @param weighted there are tow rank merging approaches for two cases: if 
#'   `weighted = FALSE`, all ranked list with the same biological state are 
#'   treated equally important, a simple but useful method average ranking 
#'   technique is selected; otherwise, `weighted = TRUE`, each individual 
#'   ranked lists has its own ranked weights, this takes the iterative 
#'   rank-aggregating algorithm, default is `TRUE`.
#' @param ... additional arguments can be passed to [`ScoreGSEA()`](while 
#'   `ScoringMethod = "GSEA"`) or to [`ScorePGSEA()`](while 
#'   `ScoringMethod = "PGSEA"`.
#' @seealso [`RankMerging()`],[`ScoreGSEA()`],[`ScorePGSEA()`]
#' @return the result from [`ScoreGSEA()`] or [`ScorePGSEA()`].
#' @examples 
#' #load the sample expressionSet
#' data(exampleSet)
#' 
#' # distance scoring
#' SignatureDistance(
#'   exampleSet,
#'   SignatureLength = 250,
#'   MergingDistance = "Spearman", 
#'   ScoringMethod = "GSEA",
#'   ScoringDistance = "avg",
#'   weighted = TRUE
#' )
#' @export
SignatureDistance <- function(exprSet,
                              SignatureLength, 
                              MergingDistance = c("Spearman", "Kendall"),
                              ScoringMethod = c("GSEA", "PGSEA"),
                              ScoringDistance = c("avg", "max"),
                              weighted = TRUE,
                              ...) {
    MergingDistance <- match.arg(MergingDistance, c("Spearman", "Kendall"))
    ScoringMethod <- match.arg(ScoringMethod, c("GSEA", "PGSEA"))
    ScoringDistance <- match.arg(ScoringDistance, c("avg", "max"))
    MergedRank <- RankMerging(exprSet, MergingDistance, weighted)
    if (ScoringMethod == "GSEA") {
        return(ScoreGSEA(MergedRank, SignatureLength, ScoringDistance, ...))
    } else {
        return(ScorePGSEA(MergedRank, SignatureLength, ScoringDistance, ...))
    }
}
