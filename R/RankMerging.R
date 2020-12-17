#' Merging the ranker lists with the same labels of the biological states into 
#' a single list with the Iorio's method
#' 
#' Merging the assay data according to phenotypic data of the input 
#' ExpressionSet. Each group of the ranked lists with the same phenotypic data 
#' is aggregated into a single list, return it as an ExpressionSet object.
#' 
#' @param exprSet an ExpressionSet object, each column of assay data represents 
#'   a ranked list obtained by preprocessing the corresponding gene expression 
#'   profile, and phenotypic data represents the short description 
#'   (characteristics of gene expression profile, such as the drug type, the 
#'   disease state) about the assay data.
#' @param MergingDistance distance to be used which "measures" the similarity 
#'   of ordered lists, the default is "Spearman"
#' @param weighted there are tow rank merging approaches for two cases: if 
#'   `weighted = FALSE`, all ranked list with the same biological state are 
#'   treated equally important, a simple but useful method average ranking 
#'   technique is selected; otherwise, `weighted = TRUE`, each individual ranked 
#'   lists has its own ranked weights, this takes the iterative 
#'   rank-aggregating algorithm, default is `TRUE`.
#'   
#' @details The krubor function is used in the aggregating procedure. And the 
#'   following methods are used in the implementation: a measure of the distance 
#'   between two ranked lists (Spearman's Footrule), a method to merge two or 
#'   more ranked lists the (Borda Merging Method), and a algorithm to obtain a 
#'   single ranked list from a set of them in a hierarchical way (the Kruskal 
#'   Algorithm). If choose Kendall as distance, the effectiveness of this 
#'   function is certainly limited by the size of the  merging problem. 
#' @seealso [`SignatureDistance()`]
#' @return a [`Biobase::ExpressionSet`] object.
#' @examples 
#' # load the sample expressionSet
#' data(exampleSet)
#' 
#' # Merging each group of the ranked lists in the exampleSet with the same
#' # phenotypic data into a single PRL
#' MergingSet <- RankMerging(exampleSet, "Spearman", weighted = TRUE)
#' 
#' @importFrom Biobase exprs phenoData
#' @importFrom methods new as
#' @export
RankMerging <- function(exprSet,
                        MergingDistance = c("Spearman", "Kendall"),
                        weighted = TRUE) {
    PRLs <- exprs(exprSet)
    MergingDistance <- match.arg(MergingDistance, c("Spearman", "Kendall"))

    # simple data preprocessing
    for (i in seq_len(ncol(PRLs))) {
        PRLs[, i] <- as.matrix(rank(PRLs[, i]))
    }

    phenodata <- as(as(phenoData(exprSet), "data.frame"), "matrix")
    if (ncol(PRLs) != length(phenodata)) {
        stop("the column of PRLS must be equal to the length of the phenodata")
    }
    FPRL <- matrix(0, nrow(PRLs))
    exp_names <- phenodata
    if (ncol(phenodata) > 1) {
        phenodata <- unique(phenodata, MARGIN = 2)
    } else {
        phenodata <- unique(phenodata, MARGIN = 1)
    }
    # phenodata=sort(phenodata)
    phenodata_num <- length(phenodata)
    gene_num <- nrow(PRLs)
    exp_num <- ncol(PRLs)
    tmp_indx <- matrix()
    for (n1 in seq_len(phenodata_num)) {
        tmp_indx <- matrix()
        diseasesI <- phenodata[n1]
        k <- 1
        for (n2 in seq_along(exp_names)) {
            if (diseasesI == exp_names[n2]) {
                tmp_indx[k] <- n2
                k <- k + 1
            }
        }
        R <- PRLs[, tmp_indx]
        R <- as.matrix(R)
        if (weighted) {
            R <- krubor(MergingDistance, R)
        } else {
            R <- rowMeans(R, na.rm = TRUE, dims = 1)
            R <- rank(R, ties.method = "first")
        }
        FPRL <- cbind(FPRL, R)
        FPRL <- as.matrix(FPRL)
    }
    FPRL <- as.matrix(FPRL[, -1])
    phenodata <- as.data.frame(phenodata)
    states <- NULL
    for (i in seq_len(nrow(phenodata))) {
        states <- cbind(states, paste(phenodata[i, 1]))
    }
    
    rownames(phenodata) <- states
    colnames(FPRL) <- rownames(phenodata)
    phenodata <- new("AnnotatedDataFrame", data = phenodata)
    return(new("ExpressionSet", exprs = FPRL, phenoData = phenodata))
}
