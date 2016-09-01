
aggregate <- function(x, y, theta=seq(0.01, 0.99, 0.01)) {

    # Center and scale predictors --------------------------

    # Should we standardize by groups??

    n <- nrow(x)
    p <- ncol(x)
    eps <- .Machine$double.eps

    # Center x
    meanx <- apply(x, 2, mean)
    x <- scale(x, meanx, FALSE)

    # Calculate variable norms
    normx <- sqrt(colSums(x^2))

    # Find variables with nontrivial amount of variation
    nonconst_idx <- which((normx / sqrt(n)) >= eps)
    if (! length(nonconst_idx)) {
        stop("All variables are constant", call.=FALSE)
    }

    p <- length(nonconst_idx)
    meanx <- meanx[nonconst_idx]
    normx <- normx[nonconst_idx]
    x <- scale(x[, nonconst_idx], FALSE, normx)


    # Choose variable quantile levels ----------------------

    classlev <- levels(y)
    classbool <- (y == classlev[2])

    c0 <- apply(x[!classbool, ], 2, sort)
    c1 <- apply(x[classbool, ], 2, sort)

    n0 <- nrow(c0)
    n1 <- nrow(c1)

    q0_map <- 1:n0 / n0
    q1_map <- 1:n1 / n1

    quantlev <- vector("numeric", p)
    classrate <- vector("numeric", p)
    c0_quant <- vector("numeric", p)
    c1_quant <- vector("numeric", p)

    for (j in 1:p) {
        rate <- classrate_noremove(c0[, j], c1[, j], theta, "interp")
        best_idx <- select_idx(rate)

        quantlev[j] <- theta[best_idx]
        classrate[j] <- rate[best_idx]
        c0_quant[j] <- c0[ head(which(q0_map > quantlev[j]), 1L), j ]
        c1_quant[j] <- c1[ head(which(q1_map > quantlev[j]), 1L), j ]
    }

    list(quantlev  = quantlev,
         classrate = classrate,
         c0_quant  = c0_quant,
         c1_quant  = c1_quant,
         keep_idx  = nonconst_idx,
         mean      = meanx,
         norm      = normx)
}




predict_aggr <- function(train, x, y) {

    x <- scale(x[, train$keep_idx], train$mean, train$norm)
    classlev <- levels(y)

    # predval <- apply(x, 1, function(z) {
    #     c0_diff <- train$c0_quant - z
    #     c1_diff <- train$c1_quant - z

    #     c0_dist <- ifelse(c0_diff > 0, train$quantlev * c0_diff, (train$quantlev - 1) * c0_diff)
    #     c1_dist <- ifelse(c1_diff > 0, train$quantlev * c1_diff, (train$quantlev - 1) * c1_diff)

    #     sum(c1_dist) > sum(c0_dist)
    # })

    predval <- vector("logical", nrow(x))
    dist_diff <- matrix(0, nrow=nrow(x), ncol=ncol(x))

    for (i in 1:nrow(x)) {
        c0_diff <- x[i, ] - train$c0_quant
        c1_diff <- x[i, ] - train$c1_quant

        c0_dist <- ifelse(c0_diff > 0, train$quantlev * c0_diff, (train$quantlev - 1) * c0_diff)
        c1_dist <- ifelse(c1_diff > 0, train$quantlev * c1_diff, (train$quantlev - 1) * c1_diff)

        dist_diff[i, ] <- c1_dist - c0_dist

        predval[i] <- ifelse((sum(c1_dist) < sum(c0_dist)),
                             classlev[2],
                             classlev[1])
    }

    # dist_diff
    mean(predval != y)
}




predict_single <- function(train, x, y, idx) {

    x <- scale(x, train$mean, train$norm)

    # predval <- apply(x, 1, function(z) {
    #     c0_diff <- train$c0_quant - z
    #     c1_diff <- train$c1_quant - z

    #     c0_dist <- ifelse(c0_diff > 0, train$quantlev * c0_diff, (train$quantlev - 1) * c0_diff)
    #     c1_dist <- ifelse(c1_diff > 0, train$quantlev * c1_diff, (train$quantlev - 1) * c1_diff)

    #     sum(c1_dist) > sum(c0_dist)
    # })

    predval <- vector("logical", nrow(x))

    for (i in 1:nrow(x)) {
        c0_diff <- train$c0_quant[idx] - x[i, idx]
        c1_diff <- train$c1_quant[idx] - x[i, idx]

        c0_dist <- ifelse(c0_diff > 0, train$quantlev * c0_diff, (train$quantlev - 1) * c0_diff)
        c1_dist <- ifelse(c1_diff > 0, train$quantlev * c1_diff, (train$quantlev - 1) * c1_diff)

        predval[i] <- (c1_dist < c0_dist)
    }
    #browser()
    mean(predval == y)
}




# Skewness transforms ----------------------------------------------------------


galton_transform <- function(x, y) {

  n <- nrow(x)
  p <- ncol(x)
  class_levels <- unique(y)

  skew <- matrix(0, 2L, p)
  for (i in 1:2) {
    skew[i,] <- apply(x[y == class_levels[i], , drop=FALSE], 2, galtonskew)
  }

  total_skew <- colSums(skew)
  total_skew <- ifelse(is.na(total_skew), 1, total_skew)
  x <- x * t( matrix(sign(total_skew), p, n) )


}


galtonskew <- function(x, y) {
  #
  #  Compute Galton's skewness measure for x
  #  NOTE: this procedure assumes no x values are missing
  #
  quarts <- as.numeric(quantile(x, probs = c(0.25, 0.5, 0.75)))
  num <- quarts[1L] + quarts[3L] - 2 * quarts[2L]
  denom <- quarts[3] - quarts[1]
  gskew <- num / denom
  gskew
}
