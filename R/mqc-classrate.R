
# Find classification rate across theta.
#
# PRE: assumes that class_0 and class_1 are sorted

classrate_noremove <- function(class_0, class_1, theta, quant_type) {

    n0 <- length(class_0)
    n1 <- length(class_1)

    if (isTRUE(all.equal(quant_type, "strict"))) {
        c0_quants <- class_0[ ceiling(theta * n0) ]
        c1_quants <- class_1[ ceiling(theta * n1) ]
    }
    else if (isTRUE(all.equal(quant_type, "interp"))) {
        c0_quants <- quantile(class_0, theta)
        c1_quants <- quantile(class_1, theta)
    }

    equal_quant_bool <- (c0_quants == c1_quants)

    # Rows are theta, cols are data from c0
    c0_to_0 <- sapply(class_0, function(z) {
        ifelse(z > c0_quants, theta * (z - c0_quants), (1 - theta) * (c0_quants - z))
    })

    # Rows are theta, cols are data cols are data from c0
    c0_to_1 <- sapply(class_0, function(z) {
        ifelse(z > c1_quants, theta * (z - c1_quants), (1 - theta) * (c1_quants - z))
    })

    # Rows are theta, cols are data cols are data from c1
    c1_to_0 <- sapply(class_1, function(z) {
        ifelse(z > c0_quants, theta * (z - c0_quants), (1 - theta) * (c0_quants - z))
    })

    # Rows are theta, cols are data cols are data from c1
    c1_to_1 <- sapply(class_1, function(z) {
        ifelse(z > c1_quants, theta * (z - c1_quants), (1 - theta) * (c1_quants - z))
    })

    # Calculate rate correct for each theta
    if (identical(length(theta), 1L)) {
        rate <- (sum((c0_to_1 - c0_to_0) > 0) + sum((c1_to_0 - c1_to_1) > 0)) / (n0 + n1)
    } else {
        rate <- (rowSums((c0_to_1 - c0_to_0) > 0) + rowSums((c1_to_0 - c1_to_1) > 0)) / (n0 + n1)
    }

    # Return classification rates and theta after removing levels with ties
    list(rate  = rate[! equal_quant_bool],
         theta = theta[! equal_quant_bool])
}




# Choose the median index out of the set of indices that tie for the maximum
# value.  When there is an even number of such indices then the lower of the 2
# middle indices is chosen.

select_idx <- function(x) {

    # Level at which anything above is deemed too unstable to choose as a good
    # quantile level.  Arbitrarily chosen by looking at the graphs of
    # classification rates for real data.
    high_sd_lev <- 0.03

    # Find length of x.  If it is only 1 value then the following algorithm will
    # fail upon calling sd
    x_len <- length(x$rate)
    if (identical(x_len, 0L)) {
        return (NULL)
    } else if (identical(x_len, 1L)) {
        return (1L)
    }
    high_sd_bool <- vector("logical", x_len)

    # Arbitrary choice for default window width.  Once window width has been
    # decided, select left width and right width.
    win_width <- min(11L, x_len)
    lwidth <- as.integer(win_width / 2)
    rwidth <- win_width - lwidth - 1

    # Each iteration sets finds the window and checks the standard deviation of
    # the classification rate in that window;
    for (k in seq_len(x_len)) {

        # case: lower bound too low
        if ((k - lwidth) < 1) {
            local_idx <- 1:win_width
        }
        # case: upper bound too high
        else if ((k + rwidth) > x_len) {
            local_idx <- (x_len - win_width + 1) : x_len
        }
        # case: window fits
        else {
            local_idx <- (k - lwidth) : (k + rwidth)
        }

        # flag quantile level if the classification rate of the nearby vicinity
        # is high
        if (sd(x$rate[local_idx]) > high_sd_lev) {
            high_sd_bool[k] <- TRUE
        }
    }

    # Check that the entire range is not too variable.  In that case just return
    # the middle index
    if (all(high_sd_bool)) {
        return (NULL)
    }

    # Find the set of indices which achieve maximum classification rate
    best_idx <- which(x$rate == max( x$rate[! high_sd_bool]) )
    best_idx_len <- length(best_idx)

    best_idx[ as.integer((best_idx_len + 1) / 2) ]
}
