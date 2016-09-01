
library(ElemStatLearn)

library(quantileDA)
library(e1071)
library(MASS)


x_SAheart <- SAheart[, -10]
y_SAheart <- as.factor( SAheart[, 10] )
y_SAheart <- SAheart[, 10]


theta <- seq(0.01, 0.99, 0.01)
prop_train <- 0.5

train_idx <- sample(1:nrow(SAheart), trunc(prop_train * nrow(SAheart)))
test_idx <- setdiff(1:nrow(SAheart), train_idx)



out <- mqc(x = x_SAheart[train_idx, ],
           y = y_SAheart[train_idx],
           aug = FALSE,
           keep_derive = TRUE,
           split_prop = NULL,
           pred_rem_lev = 0.5,
           categ_prop_lev = 0.1,
           theta = seq(0.01, 0.99, 0.01),
           provide_quantlev = NULL,
           quant_type = "interp",
           cv_type = "class",
           simil_type = "dist",
           npart = 1,
           var_lev = 0.90,
           std_parts = "none")

1 - mean(predict_mqc(out,
                     x_SAheart[test_idx, ],
                     "probs") == y_SAheart[test_idx])

out <- compare(10,
               prop_train = 0.5,
               # args for mqc
               x = x_SAheart,
               y = y_SAheart,
               aug = 999,
               keep_derive = TRUE,
               split_prop = 999,
               pred_rem_lev = 0.5,
               categ_prop_lev = 0.1,
               theta = seq(0.01, 0.99, 0.01),
               provide_quantlev = NULL,
               quant_type = "interp",
               cv_type = "class",
               simil_type = 999,
               npart = 999,
               var_lev = 0.90,
               std_parts = 999)


err_mqc <- function(out) {

    mean_err <- vector("numeric", length(out))
    for (i in seq_along(out)) {
        # cat("aug = ",           out[[ i ]]$params$aug, "  ",
        #     "split_prop = ",    out[[ i ]]$params$split_prop, "  ",
        #     "simil_type = ",    out[[ i ]]$params$simil_type, "  ",
        #     "out$npart = ",     out[[ i ]]$params$npart, ",  ",
        #     "out$std_parts = ", out[[ i ]]$params$std_parts, sep="")
        mean_err[i] <- mean(out[[ i ]]$err)
        # cat("    err:  ", mean_err[i], "\n", sep="")
    }

    mean_err
}


print_summary <- function(out) {

    v1 <- rep(c("no", "yes"), each=16)
    v2 <- rep(rep(c("n/a", "0.5"), each=8, 2))
    v3 <- rep(rep(c("dist", "pred"), each=4, 4))
    v4 <- rep(rep(c("1", "10"), each=2, 8))
    v5 <- rep(c("n/a", "norm"), 16)

    tab <- data.frame(Incl_Orig    = v1,
                      Split_Prop   = v2,
                      Qlev_Choice  = v3,
                      Num_QPart    = v4,
                      Std_dist     = v5,
                      Class_Error  = err_mqc(out$mqc),
                      stringsAsFactors=FALSE)

    cat("agg:  ", mean(out$agg), "\n",
        "hen:  ", mean(out$hen), "\n",
        "plr:  ", mean(out$plr), "\n",
        "svm:  ", mean(out$svm), "\n",
        "lda:  ", mean(out$lda), "\n",
        "qda:  ", mean(out$qda), "\n",
        "\n", sep="")

    print(tab)
    cat("\n")

    for (j in 1:(length(tab) - 1)) {
        types_j <- unique(tab[, j])
        cat(format(names(tab)[j], width=15))
        for (k in seq_along(types_j)) {
            err_jk <- mean(tab[tab[, j] == types_j[k], NCOL(tab)])
            cat(format(types_j[k], width=10, justify="right"),
                ":  ", round(err_jk, 4), "    ", sep="")
        }
        cat("\n")
    }

}


















# ------------------------------------------------------------------------------


sim_dat_train <- sim_exp_gauss(50, 100, 20, 0.5, 1, 0.8, 1)
sim_dat_test <- sim_exp_gauss(1000, 100, 20, 0.5, 1, 0.8, 1)

train_x <- sim_dat_train$x
train_x_desmat <- sim_dat_train$x
train_y <- as.factor( sim_dat_train$y )

test_x <- sim_dat_test$x
test_x_desmat <- sim_dat_test$x
test_y <- as.factor( sim_dat_test$y )


# k-NN

NUM_K <- 9L
cv_correct <- vector("integer", length=NUM_K)
for (k in seq_len(NUM_K)) {
    cv_correct[k] <- sum( knn.cv(train_x_desmat, train_y, k) == train_y )
}
# which.max returns the smallest index that has value equal to the max
klev <- which.max(cv_correct)
mean(knn(train_x_desmat, test_x_desmat, train_y, klev) != test_y)


# Naive-Bayes

nai_fit <- naiveBayes(train_x, train_y)
err_nai <- mean(predict(nai_fit, test_x) != test_y)


# Penalized LDA.  Note: can have at most K-1 components of K classes
train_y_numeric <- (train_y == 1) + 1
test_y_numeric <- (test_y == 1) + 1
capture.output(out_cv <- PenalizedLDA.cv(train_x, train_y_numeric, K=1), file="/dev/null")
fit_lda <- PenalizedLDA(train_x, train_y_numeric, test_x, lambda=out_cv$bestlambda, K=1)
err_lda <- mean(as.factor(fit_lda$ypred - 1) != test_y)


# Nearest-shrunken centroid
shr_dat <- list(x=t(train_x), y=train_y)
capture.output(shr_tr <- pamr.train(shr_dat), file="/dev/null")
capture.output(shr_cv <- pamr.cv(shr_tr, shr_dat, nfold=10), file="/dev/null")
shr_threshold <- shr_cv$threshold[which.min(shr_cv$error)]
err_shr <- mean(pamr.predict(shr_tr, t(test_x), shr_threshold) != test_y)


# Decision trees
COMPLEXITY_VAL <- 0.2
dtr_dat <- data.frame(x=train_x, y=train_y)
dtr_fit <- rpart(y ~ ., dtr_dat)
err_dtr <- mean(predict(dtr_fit, data.frame(x=test_x), "class") != test_y)









# ------------------------------------------------------------------------------

x <- sim_dat$x
x_desmat <- x
y_fact <- factor(sim_dat$y)

train_idx <- sample(1:nrow(x), trunc(prop_train * nrow(x)))
test_idx <- setdiff(1:nrow(x), train_idx)

train_x <- x[train_idx, ]
train_x_desmat <- x_desmat[train_idx, ]
train_y <- y_fact[train_idx]

test_x <- x[test_idx, ]
test_desmat_x <- x_desmat[test_idx, ]
test_y <- y_fact[test_idx]

agg_fit    <- aggregate(train_x_desmat, train_y)
mean( predict_aggr(agg_fit, test_desmat_x, test_y) )

hen_fit    <- quantilecl(train_x_desmat, test_desmat_x, train_y, , test_y)
hen_fit$me.test

plr_fit    <- cv.glmnet(train_x_desmat, train_y, type.measure="class",
                        family = "binomial", parallel=TRUE)
1 - mean( predict(plr_fit,
                                test_desmat_x,
                                type = "class",
                                s="lambda.min") == test_y )



x <- compare(1,
             prop_train = prop_train,
             train_fcn = train_fcn,
             test_fcn = test_fcn,
             # args for mqc
             x = 999,
             y = 999,
             aug = 999,
             keep_derive = TRUE,
             split_prop = 999,
             pred_rem_lev = 0.5,
             categ_prop_lev = 0.1,
             theta = seq(0.01, 0.99, 0.01),
             provide_quantlev = NULL,
             quant_type = "interp",
             cv_type = "class",
             simil_type = 999,
             npart = 999,
             var_lev = 0.90,
             std_parts = 999)
