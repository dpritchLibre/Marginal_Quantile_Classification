
compare <- function(nrepl,
                    prop_train,
                    train_fcn,
                    test_fcn,
                    # args for mqc
                    x,
                    y,
                    aug,
                    keep_derive,
                    split_prop,
                    pred_rem_lev,
                    categ_prop_lev,
                    theta,
                    provide_quantlev,
                    quant_type,
                    cv_type,
                    simil_type,
                    npart,
                    var_lev,
                    std_parts) {

    eps = .Machine$double.eps

    # Check for constant cols
    if (! is.language(train_fcn)) {
        if (is.numeric(x)) {
            is_const_cols <- (apply(x, 2, var) < eps)
        }
        else {
            is_const_cols <- sapply(x, function(z) {
                if (is.numeric(z)) {
                    var(z) < eps
                }
                else {
                    (identical(length(unique(z)), 1L))
                }
            })
        }
        if (sum(is_const_cols) > 0) {
            cat("Some of the columns were removed:  \n")
            if (!is.null(colnames(x))) {
                colnames(x)[is_const_cols]
            } else {
                print(which(is_const_cols))
            }
            x <- x[, !is_const_cols]
        }

        if (! is.numeric(x)) {
            x_desmat <- model.matrix( ~ ., x)[, -1L]
        } else {
            x_desmat <- x
        }
        y_fact <- as.factor(y)
    }


    err_agg <- vector("numeric", nrepl)
    err_hen <- vector("numeric", nrepl)
    err_plr <- vector("numeric", nrepl)
    err_svm <- vector("numeric", nrepl)
    err_knn <- vector("numeric", nrepl)
    err_nai <- vector("numeric", nrepl)
    err_shr <- vector("numeric", nrepl)
    err_lda <- vector("numeric", nrepl)
    err_dtr <- vector("numeric", nrepl)
    err_mqc <- list()

    if (identical(Sys.info()["nodename"], c(nodename="dpritch-Satellite-C875"))) {
        use_parallel <- TRUE
        doMC::registerDoMC(cores=4)
    } else {
        use_parallel <- FALSE
    }


    for (s in 1:nrepl) {

        cat("repl  ", s, "\n", sep="")

        # case: real data. split into train and test sets
        if (is.null(train_fcn)) {
            train_idx <- sample(1:nrow(x), trunc(prop_train * nrow(x)))
            test_idx <- setdiff(1:nrow(x), train_idx)

            train_x <- x[train_idx, ]
            train_x_desmat <- x_desmat[train_idx, ]
            train_y <- y_fact[train_idx]

            test_x <- x[test_idx, ]
            test_x_desmat <- x_desmat[test_idx, ]
            test_y <- y_fact[test_idx]
        }
        # case: simulated data.  simulate new train and test sets
        else {
            train_data <- eval(train_fcn)
            train_x <- train_data$x
            train_x_desmat <- train_data$x
            train_y <- as.factor( train_data$y )

            test_data <- eval(test_fcn)
            test_x <- test_data$x
            test_x_desmat <- test_data$x
            test_y <- as.factor( test_data$y )
        }

        # Aggregate quantile classifier
        agg_fit    <- aggregate(train_x_desmat, train_y)
        err_agg[s] <- mean( predict_aggr(agg_fit, test_x_desmat, test_y) )

        # Hennig classifier
        hen_fit    <- quantilecl(train_x_desmat, test_x_desmat, train_y, , test_y)
        err_hen[s] <- hen_fit$me.test

        # plr
        plr_fit    <- cv.glmnet(train_x_desmat, train_y, type.measure="class",
                                grouped=FALSE, family = "binomial",
                                parallel=use_parallel)
        err_plr[s] <- mean(predict(plr_fit,
                                   test_x_desmat,
                                   type = "class",
                                   s="lambda.min") != test_y)

        # svm
        svm_cv <- tune.svm(train_x_desmat, train_y, gamma=c(0.001, 0.01, 0.1, 1, 2), cost=2^(0:4))
        svm_gam <- svm_cv$best.parameters$gamma
        svm_cost <- svm_cv$best.parameters$cost
        svm_fit    <- svm(train_x_desmat, train_y, gamma=svm_gam, cost=svm_cost)
        err_svm[s] <- mean( predict(svm_fit, test_x_desmat) != test_y )

        # k-NN
        NUM_K <- 9L
        cv_correct <- vector("integer", length=NUM_K)
        for (k in seq_len(NUM_K)) {
            cv_correct[k] <- sum( knn.cv(train_x_desmat, train_y, k) == train_y )
        }
        klev <- which.max(cv_correct)
        err_knn[s] <- mean(knn(train_x_desmat, test_x_desmat, train_y, klev) != test_y)

        # Naive-Bayes
        nai_fit <- naiveBayes(train_x, train_y)
        err_nai[s] <- mean(predict(nai_fit, test_x) != test_y)

        # Penalized LDA.  Note: can have at most K-1 components of K classes
        train_y_numeric <- (train_y == 1) + 1
        test_y_numeric <- (test_y == 1) + 1
        capture.output(out_cv <- PenalizedLDA.cv(train_x, train_y_numeric, K=1), file="/dev/null")
        fit_lda <- PenalizedLDA(train_x, train_y_numeric, test_x, lambda=out_cv$bestlambda, K=1)
        err_lda[s] <- mean(as.factor(fit_lda$ypred - 1) != test_y)

        # Nearest-shrunken centroid
        shr_dat <- list(x=t(train_x), y=train_y)
        capture.output(shr_tr <- pamr.train(shr_dat), file="/dev/null")
        capture.output(shr_cv <- pamr.cv(shr_tr, shr_dat, nfold=10), file="/dev/null")
        shr_threshold <- shr_cv$threshold[which.min(shr_cv$error)]
        err_shr[s] <- mean(pamr.predict(shr_tr, t(test_x), shr_threshold) != test_y)

        # Decision trees
        COMPLEXITY_VAL <- 0.2
        dtr_dat <- data.frame(x=train_x, y=train_y)
        dtr_fit <- rpart(y ~ ., dtr_dat)
        err_dtr[s] <- mean(predict(dtr_fit, data.frame(x=test_x), "class") != test_y)


        ctr <- 1L

        for (aug in c(FALSE, TRUE)) {
            for (split_prop in list(NULL, 0.5)) {
                for (simil_type in c("dist", "pred")) {
                    for (npart in c(1L, 10L)) {
                        for (std_parts in c("none", "norm")) {

                            cat(ctr, "  ", sep="")

                            mqcObj <- mqc(train_x,
                                          train_y,
                                          aug              = aug,
                                          keep_derive,
                                          split_prop       = split_prop,
                                          pred_rem_lev,
                                          categ_prop_lev,
                                          theta,
                                          provide_quantlev,
                                          quant_type,
                                          cv_type,
                                          simil_type       = simil_type,
                                          npart            = 10,
                                          var_lev = 0.8,
                                          std_parts        = std_parts,
                                          nboot = 21)

                            if (identical(s, 1L)) {
                                err_mqc[[ ctr ]]  <- list()
                                err_mqc[[ ctr ]]$err <- vector("numeric", nrepl)
                            }
                            err_mqc[[ ctr ]]$params <- list(aug        = aug,
                                                            split_prop = split_prop,
                                                            simil_type = simil_type,
                                                            npart      = npart,
                                                            std_parts  = std_parts)

                            pred_y <- predict_mqc(mqcObj, test_x, "probs")
                            err_mqc[[ ctr ]]$err[s] <- 1 - mean(test_y == pred_y)

                            ctr <- ctr + 1L
                        }
                    }
                }
            }
        }
        cat("\n")

    } # end test classifiers loop

    cat("\n\n")

    list(agg = err_agg,
         hen = err_hen,
         plr = err_plr,
         svm = err_svm,
         knn = err_knn,
         nai = err_nai,
         shr = err_shr,
         lda = err_lda,
         dtr = err_dtr,
         mqc = err_mqc,
         n   = NROW(train_x),
         p   = NCOL(train_x))
}
