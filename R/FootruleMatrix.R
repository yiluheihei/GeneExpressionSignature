# footrule matrix
FootruleMatrix <- function(rankings,
                           MergingDistance = c("Spearman", "Kendall"),
                           n) {
    MergingDistance <- match.arg(MergingDistance, c("Spearman", "Kendall"))
    nrank <- ncol(rankings)
    if (length(rankings) != 0) {
        smdm <- matrix(0, nrow = nrank, ncol = nrank)
        if (nrank > 1) {
            for (i in seq_len(nrank - 1)) {
                for (j in (i + 1):nrank) {
                    # optimize
                    if (MergingDistance == "Spearman") {
                        t <- sm_footrule(rankings[, i], rankings[, j])
                    } else {
                        t <- kendall_footrule(rankings[, i], rankings[, j])
                    }
                    if (t > 0) {
                        smdm[i, j] <- t
                    } else {
                        smdm[i, j] <- 0
                    }
                }
            }
            smdm <- smdm + t(smdm)
            if (n != 0) {
                smdmrol <- nrow(smdm)
                smdmcol <- ncol(smdm)
                if (max(smdmrol, smdmcol) > 1) {
                    if (max(smdm[]) > 0) {
                        smdm <- smdm / max(smdm[])
                    }
                }
            }
        }
        else if (nrank == 1) {
            smdm <- 1
        }
        else {
            smdm <- 0
        }
    }
    return(smdm)
}


# Kendall distance,this source code is derived from R package RankAggreg
kendall_footrule <- function(x, y) {
    K <- 0
    n <- length(x)
    for (i in seq_len(n - 1)) {
        for (j in i:n) {
            if ((x[i] > x[j] & y[i] < y[j]) | (x[i] < x[j] & y[i] > y[j])) {
                K <- K + 1
            }
        }
    }
    K
}

# Spearman distance
sm_footrule <- function(r1, r2) {
    smd <- sum(abs(r1 - r2))
    smd
}
