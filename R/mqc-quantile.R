
mqc_quant_dist <- function(class_0, class_1, theta, quant_type, partition, var_lev) {

    n0 <- length(class_0)
    n1 <- length(class_1)
    part_len <- length(partition)

    if (identical(n0, 0L) || identical(n1, 0L)) {
        return (NULL)
    }

    if (isTRUE(all.equal(quant_type, "strict"))) {
        c0_quants <- class_0[ ceiling(theta * n0) ]
        c1_quants <- class_1[ ceiling(theta * n1) ]
    }
    else if (isTRUE(all.equal(quant_type, "interp"))) {
        c0_quants <- quantile(class_0, theta)
        c1_quants <- quantile(class_1, theta)
    }

    quant_dist <- abs(c0_quants - c1_quants)

    qlev <- vector("numeric", part_len)
    qdist <- vector("numeric", part_len)
    qval_c0 <- vector("numeric", part_len)
    qval_c1 <- vector("numeric", part_len)

    for (k in seq_len(part_len)) {

        curr_quant_idx <- dist_select_idx(quant_dist, partition[[k]])

        if (! is.null(curr_quant_idx)) {
            qlev[k] <- theta[curr_quant_idx]
            qdist[k] <- quant_dist[curr_quant_idx]
            qval_c0[k] <- c0_quants[curr_quant_idx]
            qval_c1[k] <- c1_quants[curr_quant_idx]
        }
        else {
            qlev[k] <- NA_real_
            qdist[k] <- NA_real_
            qval_c0[k] <- NA_real_
            qval_c1[k] <- NA_real_
        }
    }

    qweights <- dist_weights(qdist, var_lev)

    list(qlev     = qlev,
         qdist    = qdist,
         qweights = qweights,
         qval_c0  = qval_c0,
         qval_c1  = qval_c1)
}



dist_weights <- function(dist_vals, var_lev) {

    if (all(is.na(dist_vals) | (dist_vals == 0))) {
        return (rep(0, length(dist_vals)))
    }

    sort_idx <- order(dist_vals, decreasing=TRUE)
    var_cap <- var_lev * sum(dist_vals, na.rm=TRUE)

    tot_dist <- 0
    for (k in sort_idx) {
        if (tot_dist > var_cap) {
            dist_vals[k] <- 0
        }
        tot_dist <- tot_dist + dist_vals[k]
    }

    dist_vals / sum(dist_vals, na.rm=TRUE)
}




dist_select_idx <- function(distvals, idx) {

    x <- distvals[idx]
    if (all(x == 0)) {
        return (NULL)
    }

    # Find the set of indices with the maximum distance
    best_idx <- which(x == max(x))
    best_idx_len <- length(best_idx)
    local_best <- best_idx[ as.integer((best_idx_len + 1) / 2) ]

    idx[local_best]
}




mqc_quant_pred <- function(class_0, class_1, theta, quant_type, partitions,
                           var_lev, pred_rem_lev) {

    n0 <- length(class_0)
    n1 <- length(class_1)
    part_len <- length(partitions)

    if (identical(n0, 0L) || identical(n1, 0L)) {
        return (NULL)
    }

    if (isTRUE(all.equal(quant_type, "strict"))) {
        c0_quants <- class_0[ ceiling(theta * n0) ]
        c1_quants <- class_1[ ceiling(theta * n1) ]
    }
    else if (isTRUE(all.equal(quant_type, "interp"))) {
        c0_quants <- quantile(class_0, theta)
        c1_quants <- quantile(class_1, theta)
    }

    qlev <- vector("numeric", part_len)
    qpred <- vector("numeric", part_len)
    qval_c0 <- vector("numeric", part_len)
    qval_c1 <- vector("numeric", part_len)

    for (k in seq_len(part_len)) {

        rate <- classrate_noremove(class_0, class_1, theta[ partitions[[ k ]] ], quant_type)
        partit_idx <- select_idx(rate)
        best_idx <- partitions[[ k ]][ partit_idx ]

        # case: at least one quantile is different for the groups
        if (! is.null(partit_idx)) {
            qlev[k] <- theta[best_idx]
            qpred[k] <- rate$rate[partit_idx]
            qval_c0[k] <- c0_quants[best_idx]
            qval_c1[k] <- c1_quants[best_idx]
        }
        # case: no differences between groups
        else {
            qlev[k] <- NA_real_
            qpred[k] <- NA_real_
            qval_c0[k] <- NA_real_
            qval_c1[k] <- NA_real_
        }
    }

    qweights <- pred_weights(qpred, var_lev, pred_rem_lev)

    list(qlev = qlev,
         qpred = qpred,
         qweights = qweights,
         qval_c0 = qval_c0,
         qval_c1 = qval_c1)
}




pred_weights <- function(pred_vals, var_lev, pred_rem_lev) {

    EXP_VAL <- 8

    pred_vals[!is.na(pred_vals) & (pred_vals < pred_rem_lev)] <- 0

    rel_pred <- pred_vals / max(pred_vals, na.rm=TRUE)
    rel_exp <- rel_pred ^ EXP_VAL

    dist_weights(rel_exp, var_lev)
}
