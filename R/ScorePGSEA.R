#' Compute pairwise distances between samples with method in package PGSEA
#' 
#' Compute pairwise distances between sample according to their 
#' (Prototype Ranked List) PRL, get a `N x N` distance matrix is generated by 
#' calling this function , `N` is the length of PRL.
#' 
#' @param MergingSet an [`Biobase::ExpressionSet`] object. The assay data 
#'   represents the PRLs of the samples, each column represents one PRL. The 
#'   number of sample must be greater than 1, otherwise, this function is not 
#'   meaningful.
#' @param SignatureLength the length of "gene signature". In order to compute 
#'   pairwise distances among samples, genes lists are ranked according to the 
#'   gene expression ratio (fold change). And the "gene signature" includes the 
#'   most up-regulated genes (near the top of the list) and the most 
#'   down-regulated genes (near the bottom of the list).   
#' @param ScoringDistance the distance measurements between PRLs: the Average 
#'   Enrichment Score Distance ("avg"), or the Maximum Enrichment Score 
#'   Distance ("max").  
#' @param p.value logical, if `TRUE` return a matrix of `p.values` of the 
#'   distance matrix, default `FALSE`.
#' @seealso [`ScoreGSEA()`],[`SignatureDistance()`]
#' @return an distance-matrix, the max distance is more sensitive to weak 
#'   similarities, providing a lower precision but a larger recall.If `p.value` 
#'   is set to `TRUE`, then a list is returned that consists of the distance 
#'   matrix  as well as their p.values, otherwise, without p.values in the 
#'   result. 
#' @examples
#' # load the sample expressionSet
#' data(exampleSet)
#' # Merging each group of the ranked lists in the exampleSet with the same
#' # phenotypic data into a single PRL
#' MergingSet <- RankMerging(exampleSet,"Spearman")
#' # get the distance matrix
#' ds <- ScorePGSEA(MergingSet,250, ScoringDistance="avg")
#' @importFrom Biobase exprs
#' @importFrom methods new
#' @importFrom stats pnorm
#' @export
ScorePGSEA <- function(MergingSet,
                       SignatureLength,
                       ScoringDistance = c("avg", "max"),
                       p.value = FALSE) {
    ScoringDistance <- match.arg(ScoringDistance, c("avg", "max"))
    PRLs <- exprs(MergingSet)
    n <- ncol(PRLs)
    m <- nrow(PRLs)
    pgscores <- matrix(0, 1, n)
    for (i in seq_len(n)) {
        prls <- as.matrix(PRLs[, i])
        up <- which(prls <= SignatureLength)
        down <- which(prls >= m - SignatureLength + 1)
        upsmc <- new("smc", ids = up)
        downsmc <- new("smc", ids = down)
        uppgscore <- PGSEA(MergingSet, list(upsmc), p.value = NA)
        downpgscore <- PGSEA(MergingSet, list(downsmc), p.value = NA)
        pgscore <- (uppgscore - downpgscore) / 2
        pgscores <- rbind(pgscores, pgscore)
    }
    pgscores <- as.matrix(pgscores)
    pgscores <- pgscores[-1, ]
    Mvalue <- max(abs(pgscores))
    pgscores <- pgscores / max(abs(pgscores))
    if (pgscores[1, 1] > 0) {
        pgscores <- pgscores
    }
    else {
        pgscores <- -pgscores
    }
    rownames(pgscores) <- colnames(pgscores)
    p.results <- matrix(0, n, n)
    if (ScoringDistance == "avg") {
        distances <- 1 - (pgscores + t(pgscores)) / 2
        for (i in seq_len(n)) {
            for (j in seq_len(n)) {
                if (distances[i, j] < 1) {
                    p.results[i, j] <- 2 * pnorm(distances[i, j],
                        mean = 1,
                        sd = 1 / (2 * Mvalue)
                    )
                } else {
                    p.results[i, j] <- 2 * (1 - pnorm(distances[i, j],
                        mean = 1,
                        sd = 1 / (2 * Mvalue)
                    ))
                }
            }
        }
    }
    else {
        distances <- pmin(1 - pgscores, t(1 - pgscores)) / 2
        p.results <- matrix(0, nrow(distances), ncol(distances))
        for (i in seq_len(n)) {
            for (j in seq_len(n)) {
                if (distances[i, j] < 1 / 2) {
                    p.results[i, j] <- 2 * pnorm(distances[i, j],
                        mean = 1 / 2,
                        sd = 1 / (8^(1 / 2) * Mvalue)
                    )
                } else {
                    p.results[i, j] <- 2 * (1 - pnorm(distances[i, j],
                        mean = 1 / 2,
                        sd = 1 / (8^(1 / 2) * Mvalue)
                    ))
                }
            }
        }
    }
    if (p.value == TRUE) {
        return(list(distances = distances, p.results = p.results))
    }
    return(distances)
}
