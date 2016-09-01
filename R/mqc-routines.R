
check_loss <- function(x, q, tau) {
    ifelse(x > q, tau * (x - q), (1 - tau) * (q - x))
}




mqc_transform_x <- function(x, categ_prop_lev) {

    # Arbitrary choice for the minimum number of non-pm values needed to keep a
    # continuous variable around
    cont_rem_lev <- 0.05

    # Find which variables are numeric (as opposed to categorical)
    if (is.matrix(x)) {
        numer_col_bool <- rep(TRUE, p)
    } else {
        numer_col_bool <- sapply(x, is.numeric)
    }

    # Split x into categorical and continuous data
    x_numer_orig <- as.matrix( x[, numer_col_bool,  drop=FALSE] )
    x_categ <- x[, !numer_col_bool,  drop=FALSE]

    n <- NROW(x)
    p_numer <- NCOL(x_numer_orig)

    # Find vars with point masses
    pm_vals <- list()
    for (j in seq_len(p_numer)) {
        vals_count <- rle(sort( x_numer_orig[, j]  ))
        pm_vals[[ j ]] <- vals_count$values[
            vals_count$lengths > (categ_prop_lev * n)
        ]
    }

    # Track indices with point mass vars.  Note: which() returns integer(0) when
    # all vals are false
    pm_col_idx <- which( sapply(pm_vals, function(w) length(w) > 0) )

    # Convert point mass vars to categorical / continuous split
    convert_dat <- mqc_convert_pm(x_numer_orig, pm_col_idx, pm_vals)

    x_numer  <- convert_dat$x_numer
    x_derive <- convert_dat$x_derive
    pm_bool  <- convert_dat$pm_bool

    # Track numeric data without many continuous values (consequently to be
    # removed later)
    rem_cont_bool <- vector("logical", p_numer)
    for (j in pm_col_idx) {
        curr_ncont <- sum( pm_bool[[ j ]] )
        if (((curr_ncont / n) < cont_rem_lev) || (curr_ncont < 20)) {
            rem_cont_bool[j] <- TRUE
        }
    }
    rem_cont_idx <- which(rem_cont_bool)

    list(x_numer       = x_numer,
         x_categ       = x_categ,
         x_derive      = x_derive,
         numer_col_idx = which(numer_col_bool),
         pm_bool       = pm_bool,
         pm_col_idx    = pm_col_idx,
         pm_vals       = pm_vals,
         rem_cont_idx  = rem_cont_idx)
}




mqc_convert_pm <- function(x_numer, pm_col_idx, pm_vals,
                           rem_cont_idx=integer(0), low_pred_idx=integer(0)) {

    n <- NROW(x_numer)
    p_numer <- NCOL(x_numer)

    # Convert point mass vars to categorical / continuous split
    x_derive <- list()
    pm_bool <- list()
    for (j in seq_len(p_numer)) {

        # case: not a var with pm obs so nothing to do
        if (! (j %in% pm_col_idx)) {
            next
        }

        x_derive[[ j ]] <- rep("cont", n)
        for (currval in pm_vals[[ j ]]) {
            curr_pm_row_idx <- which( x_numer[, j] == currval )
            x_derive[[ j ]][curr_pm_row_idx] <- as.character(currval)
        }

        # Ignore vars with low number of cont. vals or low predictive power
        if ( !(j %in% rem_cont_idx) && !(j %in% low_pred_idx) ) {
            cont_row_bool <- (x_derive[[ j ]] == "cont")
            pm_bool[[ j ]] <- cont_row_bool
            x_numer[!cont_row_bool, j] <- 0
        }
    }

    # Combine the categorical variables derived from point mass vars into a
    # data.frame.  Note: data.frame call converts to factors.
    x_derive <- data.frame( Filter(function(w) !is.null(w), x_derive) )
    if (! NROW(x_derive)) {
        x_derive <- data.frame( matrix(0, n, 0) )
    }

    # Name data.frame columns
    if (! identical(NCOL(x_derive), 0L)) {
        if (is.null(colnames(x_numer))) {
            colnames(x_derive) <- paste0("col", pm_col_idx)
        } else {
            colnames(x_derive) <- colnames(x_numer)[pm_col_idx]
        }
    }

    list(x_numer  = x_numer,
         x_derive = x_derive,
         pm_bool  = pm_bool)
}




mqc_quantlev <- function(x_numer,
                         classbool,
                         pm_bool,
                         pm_col_idx,
                         rem_cont_idx,
                         pred_rem_lev,
                         theta,
                         quant_type,
                         simil_type,
                         npart,
                         var_lev) {

    p_numer <- NCOL(x_numer)

    quantlev  <- vector("numeric", p_numer)  # quantile level chosen for each var
    classrate <- vector("numeric", p_numer)  # leave-one-out classification rate
    c0_quant  <- vector("numeric", p_numer)  # quant vals at quantlevs, class 0
    c1_quant  <- vector("numeric", p_numer)  # quant vals at quantlevs, class 1

    low_pred_bool <- vector("logical", p_numer)
    quant_dat <- list()

    lo <- as.integer(0:(npart - 1L) * length(theta) / npart) + 1L
    hi <- as.integer(1:npart * length(theta) / npart)
    partitions <- mapply(function(a, b) a:b, lo, hi, SIMPLIFY=FALSE)

    # Each iteration calculates the best quantile level chosen for variable, the
    # leave-one-out classification rate for this variable, and the quantile
    # value at this level for classes 0 and 1
    for (j in seq_len(p_numer)) {

        # ignore variables with index in rem_cont_idx
        if (j %in% rem_cont_idx) {
            next
        }

        # case: pm var
        if (j %in% pm_col_idx) {

            # specifies the continuous set of values in var j
            curr_c0_bool <- pm_bool[[ j ]] & !classbool
            curr_c1_bool <- pm_bool[[ j ]] & classbool

            curr_c0 <- sort( x_numer[curr_c0_bool, j] )
            curr_c1 <- sort( x_numer[curr_c1_bool, j] )
        }
        # case: var is strictly continuous
        else {
            curr_c0 <- sort( x_numer[!classbool, j] )
            curr_c1 <- sort( x_numer[classbool, j] )
        }

        curr_n_c0 <- length(curr_c0)
        curr_n_c1 <- length(curr_c1)
        if (identical(curr_n_c0, 0L) || identical(curr_n_c1, 0L)) {
            next
        }

        if (isTRUE(all.equal(simil_type, "dist"))) {

            quant_dat[[ j ]] <- mqc_quant_dist(curr_c0,
                                               curr_c1,
                                               theta,
                                               quant_type,
                                               partitions,
                                               var_lev)
        }
        else if (isTRUE(all.equal(simil_type, "pred"))) {

            quant_dat[[ j ]] <- mqc_quant_pred(curr_c0,
                                               curr_c1,
                                               theta,
                                               quant_type,
                                               partitions,
                                               var_lev,
                                               pred_rem_lev)
        }
        else {
            stop("invalid arg for simil_type")
        }
    }

    quant_dat
}




mqc_fixed_quantlev <- function(x_numer, classbool, pm_bool, pm_col_idx, rem_cont_idx,
                               pred_rem_lev, quantlev, quant_type) {

    p_numer <- NCOL(x_numer)

    classrate <- vector("numeric", p_numer)  # leave-one-out classification rate
    c0_quant  <- vector("numeric", p_numer)  # quant vals at quantlevs, class 0
    c1_quant  <- vector("numeric", p_numer)  # quant vals at quantlevs, class 1

    low_pred_bool <- vector("logical", p_numer)

    # Each iteration calculates the best quantile level chosen for variable, the
    # leave-one-out classification rate for this variable, and the quantile
    # value at this level for classes 0 and 1
    for (j in 1:p_numer) {

        # ignore variables with index in rem_cont_idx
        if (j %in% rem_cont_idx) {
            next
        }

        # case: pm var
        if (j %in% pm_col_idx) {

            # specifies the continuous set of values in var j
            curr_c0_bool <- pm_bool[[ j ]] & !classbool
            curr_c1_bool <- pm_bool[[ j ]] & classbool

            curr_c0 <- sort( x_numer[curr_c0_bool, j] )
            curr_c1 <- sort( x_numer[curr_c1_bool, j] )
        }

        # case: var is strictly continuous
        else {
            curr_c0 <- sort( x_numer[!classbool, j] )
            curr_c1 <- sort( x_numer[classbool, j] )
        }

        curr_n_c0 <- length(curr_c0)
        curr_n_c1 <- length(curr_c1)

        # calculate leave-one-out classification rate.  best_idx is the index of
        # the best rate in the set, and low_pred_bool tracks whether the best
        # rate is below the threshold as determined by pred_rem_lev.
        if (identical(curr_n_c0, 0L) || identical(curr_n_c1, 0L)) {
            rate <- 0
            classrate[j] <- NA_real_
            c0_quant[j] <- NA_real_
            c1_quant[j] <- NA_real_
        }
        else {
            rate <- classrate_noremove(curr_c0, curr_c1, quantlev[j], quant_type)
            classrate[j] <- rate$rate

            if (isTRUE(all.equal(quant_type, "strict"))) {
                c0_quant[j] <- curr_c0[ ceiling(quantlev[j] * curr_n_c0) ]
                c1_quant[j] <- curr_c1[ ceiling(quantlev[j] * curr_n_c1) ]
            }
            else if (isTRUE(all.equal(quant_type, "interp"))) {
                c0_quant[j] <- quantile(curr_c0, quantlev[j])
                c1_quant[j] <- quantile(curr_c1, quantlev[j])
            }
        }

        if (rate$rate < pred_rem_lev) {
            low_pred_bool[j] <- TRUE
        }
    }

    list(quantlev     = quantlev,
         classrate    = classrate,
         c0_quant     = c0_quant,
         c1_quant     = c1_quant,
         low_pred_idx = which(low_pred_bool))
}




mqc_transform_dist <- function(x_numer,
                               quant_dat,
                               pm_bool,
                               pm_col_idx,
                               rem_cont_idx,
                               low_pred_idx,
                               std_parts) {

    n <- NROW(x_numer)
    p_numer <- NCOL(x_numer)
    K <- length( quant_dat[[ 1L ]]$qlev )

    all_keep_bool <- rep(TRUE, n)

    # Transform data
    trans0 <- matrix(0, n, p_numer)
    trans1 <- matrix(0, n, p_numer)

    for (j in seq_len(p_numer)) {

        # Ignore vars with low number of cont. vals or low predictive power
        if ((j %in% rem_cont_idx) || (j %in% low_pred_idx)) {
            next
        }

        # case: pm var, restrict transformation to continuous vals
        if (j %in% pm_col_idx) {
            curr_keep_bool <- pm_bool[[ j ]]
        }
        # case: var is strictly continuous, use all values
        else {
            curr_keep_bool <- all_keep_bool
        }

        n_curr <- sum(curr_keep_bool)
        t0_j <- rep(0, n_curr)
        t1_j <- rep(0, n_curr)

        x_numer_j <- x_numer[curr_keep_bool, j]

        for (k in seq_len(K)) {

            wt_k <- quant_dat[[ j ]]$qweights[k]
            if (is.na(wt_k) || (wt_k == 0)) {
                next
            }

            qlev_jk <- quant_dat[[ j ]]$qlev[k]
            qval_c0_jk <- quant_dat[[ j ]]$qval_c0[k]
            qval_c1_jk <- quant_dat[[ j ]]$qval_c1[k]

            t0_jk <- check_loss(x_numer_j, qval_c0_jk, qlev_jk)
            t1_jk <- check_loss(x_numer_j, qval_c1_jk, qlev_jk)

            if (isTRUE(all.equal(std_parts, "norm"))) {
                m0 <- mean(t0_jk)
                m1 <- mean(t1_jk)
                norm0 <- sqrt( sum( (t0_jk - m0)^2 ) / n )
                norm1 <- sqrt( sum( (t1_jk - m1)^2 ) / n )
                t0_jk <- t0_jk / norm0
                t1_jk <- t1_jk / norm1
            }
            else if (isTRUE(all.equal(std_parts, "none"))) {
                # noop
            }
            else {
                stop("illegal arg for std_parts")
            }

            t0_j <- t0_j + t0_jk
            t1_j <- t1_j + t1_jk

        }

        trans0[curr_keep_bool, j] <- t0_j
        trans1[curr_keep_bool, j] <- t1_j
    }

    # Return the difference of the distances
    keep_idx <- setdiff(seq_len(p_numer), union(rem_cont_idx, low_pred_idx))
    trans1[, keep_idx] - trans0[, keep_idx]
}




mqc_combine <- function(x_trans, x_categ, x_derive, x_numer, aug, keep_derive) {

    n <- NROW(x_trans)

    # Convert categorical to a design matrix
    if (identical(NCOL(x_derive), 0L)) {
        design_derive <- matrix(0, n, 0)
    } else {
        design_derive <- model.matrix( ~ ., x_derive)[, -1]
    }

    # Convert derived categorical to design matrix
    if (identical(NCOL(x_categ), 0L)) {
        design_categ <- matrix(0, n, 0)
    } else {
        design_categ <- model.matrix( ~ ., x_categ)[, -1]
    }

    # Combine data into a single matrix
    if (aug && keep_derive) {
        x_combine <- cbind(x_trans, design_derive, design_categ, x_numer)
    } else if (aug && !keep_derive) {
        x_combine <- cbind(x_trans, design_categ, x_numer)
    } else if (!aug && keep_derive) {
        x_combine <- cbind(x_trans, design_derive, design_categ)
    } else { # case: !aug && !keep_derive
        x_combine <- cbind(x_trans, design_categ)
    }

    x_combine
}
