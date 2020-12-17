## The source code of this function is from the PGSEA.
## Since package PGSEA was deprecated in Bioconductor version 3.12 and removed
## from 3.13. To not depend on PGSEA, we create a function PGSEA by copying the
##  source code from package PGSEA.

setClass(
    "smc",
    representation(
        reference = "character",
        desc = "character",
        source = "character",
        design = "character",
        identifier = "character",
        species = "character",
        data = "character",
        private = "character",
        creator = "character",
        ids = "vector"
    )
)



#' @importFrom methods is
#' @importFrom stats pnorm sd t.test
#' @importFrom Biobase exprs
PGSEA <- function(exprs,
                  cl,
                  range = c(25, 500),
                  ref = NULL,
                  center = TRUE,
                  p.value = 0.005,
                  weighted = TRUE,
                  enforceRange = TRUE,
                  ...) {
    if (is(exprs, "ExpressionSet")) {
        exprs <- exprs(exprs)
    }

    if (!is.list(cl)) {
        stop("cl need to be a list")
    }
    if (!is.null(ref)) {
        if (!is.numeric(ref)) {
            stop("column index's required")
        }
    }
    if (!is.null(ref)) {
        if (options()$verbose) {
            cat("Creating ratios...", "\n")
        }
        ref_mean <- apply(exprs[, ref], 1, mean, na.rm = TRUE)
        exprs <- sweep(exprs, 1, ref_mean, "-")
    }
    if (center) {
        exprs <- scale(exprs, scale = FALSE)
    }
    results <- matrix(NA, length(cl), ncol(exprs))
    rownames(results) <- names(cl)
    colnames(results) <- colnames(exprs)
    mode(results) <- "numeric"
    if (is.logical(p.value)) {
        p.results <- results
    }

    for (i in seq_along(cl)) {
        if (is(cl[[i]], "smc")) {
            clids <- cl[[i]]@ids
        } else if (class(cl[[i]]) %in% c("GeneColorSet", "GeneSet")) {
            clids <- cl[[i]]@geneIds
        } else {
            clids <- cl[[i]]
        }
        if (options()$verbose) {
            cat("Testing region ", i, "\n")
        }
        ix <- match(clids, rownames(exprs))
        ix <- unique(ix[!is.na(ix)])
        present <- sum(!is.na(ix))
        if (present < range[1]) {
            if (options()$verbose) {
                cat("Skipping region ", i, " because too small-", present, ",\n")
            }
            next
        }
        if (present > range[2]) {
            if (options()$verbose) cat("Skipping region ", i, " because too large-", present, "\n")
            next
        }

        texprs <- exprs[ix, ]
        if (any(is.na(texprs))) {
            cat("Warning - 'NA' values within expression data, enrichment scores are estimates only.\n")
        }
        if (!is.matrix(texprs)) texprs <- as.matrix(texprs)

        if (!weighted) {
            stat <- try(apply(texprs, 2, t.test, ...))
        } else {
            try(mod <- (length(ix)^(1 / 2)) / apply(exprs, 2, sd, na.rm = TRUE))
            try(m <- apply(texprs, 2, mean, na.rm = TRUE) -
                apply(exprs, 2, mean, na.rm = TRUE))
            stat2 <- m * mod
            p.val <- 2 * pnorm(-abs(stat2))
            stat <- list()
            for (q in seq_along(stat2)) {
                stat[[q]] <- list(statistic = stat2[q], p.value = p.val[q])
            }
            names(stat) <- names(stat2)
        }
        if (is.list(stat)) {
            ps <- unlist(lapply(stat, function(x) x$p.value))
            stat <- unlist(lapply(stat, function(x) x$statistic))
            if (!is.na(p.value)) {
                if (is.numeric(p.value)) {
                    stat[ps > p.value] <- NA
                } else {
                    p.results[i, ] <- ps
                }
            }
        }
        results[i, ] <- as.numeric(stat)
        if (enforceRange) {
            for (w in seq_len(ncol(texprs))) {
                if (sum(!is.na(texprs[, w])) < range[1] | sum(!is.na(texprs[, w])) > range[2]) {
                    results[i, w] <- NA
                }
            }
        }
    }
    if (is.logical(p.value) & !is.na(p.value)) {
        return(list(results = results, p.results = p.results))
    }
    return(results)
}
