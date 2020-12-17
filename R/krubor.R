krubor <- function(MergingDistance = c("Spearman", "Kendall"), ...) {
    MergingDistance <- match.arg(MergingDistance, c("Spearman", "Kendall"))
    r <- data.frame(...)
    r <- as.matrix(r)
    r <- unique(r, MARGIN = 2)
    r <- as.matrix(r)
    nrank <- ncol(r)
    while (nrank != 0) {
        smdm <- FootruleMatrix(r, MergingDistance, 1)
        smdm <- as.matrix(smdm)
        smdm[lower.tri(smdm)] <- 0
        if (nrank == 1) {
            nrank <- 0
        } else {
            smdm[smdm == 0] <- "inf"
            mij <- find_closest_rank(smdm)
            i <- mij[[2]]
            j <- mij[[3]]
            r1 <- r[, i]
            r2 <- r[, j]
            r12 <- cbind(r1, r2)
            r12 <- as.matrix(r12)
            newrank <- bm_rank_merging(r12)
            newrank <- as.matrix(newrank)
            nrank <- nrank - 1
            r <- r[, c(-i, -j)]
            r <- cbind(r, newrank)
        }
    }
    r
}


# find the two closest ranks and merge them into a rank
find_closest_rank <- function(smdm) {
    if (length(smdm) != 0) {
        m <- apply(smdm, 2, min)
        i <- apply(smdm, 2, which.min)
        j <- which.min(m)
        m <- min(m)
        i <- i[j]
    }
    else {
        stop("zero-length smdm is illegal")
    }
    list(m, i, j)
}

# Borda Mering
bm_rank_merging <- function(rankings) {
    if (ncol(rankings) > 1) {
        majorities <- rowSums(rankings)
        tmp <- sort(majorities)
        sidxs <- order(majorities)
        tmp <- sort(sidxs)
        outrank <- order(sidxs)
    }
    else {
        outrank <- rankings
    }
    return(as.matrix(outrank))
}
