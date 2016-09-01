
library(MASS)


# exponentiated Gaussian r.v.'s ------------------------------------------------

sim_exp_gauss <- function(n,
                          p,
                          r,
                          prop_class_0 = 0.5,
                          mean_shift = 0.2,
                          rho = 0.8,
                          sigma = 1) {

    # number of observations in each class
    n0 <- n * prop_class_0
    n1 <- n - n0

    # number of relevant variables (i.e. the remaining vars are noise)
    nvars_relev <- r
    nvars_non <- p - nvars_relev

    # correlation matrix for underlying Gaussian data
    cov_relev <- ar1(nvars_relev, rho, sigma)

    c0_relev <- MASS::mvrnorm(n0, rep(0, nvars_relev), cov_relev)
    c1_relev <- MASS::mvrnorm(n1, rep(mean_shift, nvars_relev), cov_relev)

    meanvec <- rnorm(nvars_non, mean_shift / 2, mean_shift)
    x <- cbind(rbind(c0_relev, c1_relev),
               matrix(rnorm(n * nvars_non, mean_shift / 2, mean_shift / 4), n, nvars_non))

    y <- rep(c(0, 1), c(n0, n1))

    list(x = exp(x),
         y = y)
}


ar1 <- function(p, rho, sigma) {

    H <- abs(outer(1:p, 1:p, "-"))
    V <- sigma * rho^H
    diag(V) <- sigma * sigma
    V
}






