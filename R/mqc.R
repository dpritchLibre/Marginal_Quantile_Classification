
mqc_ <- function(x,
                 y,
                 aug = TRUE,
                 keep_derive = TRUE,
                 split_prop = 0.5,
                 pred_rem_lev = 0.65,
                 categ_prop_lev = 0.2,
                 theta = seq(0.01, 0.99, 0.01),
                 provide_quantlev = NULL,
                 quant_type = "interp",
                 cv_type = "class",
                 simil_type = "dist",
                 npart = 10,
                 var_lev = 0.8,
                 std_parts = "norm") {

    cont_rem_lev <- 0.05

    n <- NROW(x)
    p <- NCOL(x)

    # Remove constant predictors ---------------------------

    # TODO


    # Transform non-continuous variables -------------------

    transform_dat <- mqc_transform_x(x, categ_prop_lev)

    x_numer       <- transform_dat$x_numer
    x_categ       <- transform_dat$x_categ
    x_derive      <- transform_dat$x_derive
    numer_col_idx <- transform_dat$numer_col_idx
    pm_bool       <- transform_dat$pm_bool
    pm_col_idx    <- transform_dat$pm_col_idx
    pm_vals       <- transform_dat$pm_vals
    rem_cont_idx  <- transform_dat$rem_cont_idx


    # Randomly split data ----------------------------------

    # Train / test split indices
    if (! is.null(split_prop)) {
        quant_idx <- sample.int(n, as.integer(split_prop * n))
        penal_idx <- setdiff(1:n, quant_idx)
    } else {
        quant_idx <- seq_len(n)
        penal_idx <- seq_len(n)
    }

    # classbool: TRUE / FALSE for y in class 1
    y <- as.factor(y)
    classlev <- levels(y)
    classbool <- (y == classlev[2])


    # Choose variable quantile levels ----------------------

    # Training data for quantile levels
    x_numer_quant <- x_numer[quant_idx, ]
    classbool_quant <- classbool[quant_idx]
    pm_var_cont_quant_bool <- lapply(pm_bool, function(w) w[quant_idx])

    # case: quantile levels not provided, have to select
    if (is.null(provide_quantlev)) {

        quant_dat <- mqc_quantlev(x_numer_quant,
                                  classbool_quant,
                                  pm_var_cont_quant_bool,
                                  pm_col_idx,
                                  rem_cont_idx,
                                  pred_rem_lev,
                                  theta,
                                  quant_type,
                                  simil_type,
                                  npart,
                                  var_lev)
    }
    # case: quantile levels provided
    else {

        stop("providing quantile levels not currently supported")
        # quant_dat <- mqc_fixed_quantlev(x_numer_quant,
        #                                 classbool_quant,
        #                                 pm_var_cont_quant_bool,
        #                                 pm_col_idx,
        #                                 rem_cont_idx,
        #                                 pred_rem_lev,
        #                                 provide_quantlev,
        #                                 quant_type,
        #                                 simil_type)
    }

    # TODO: what to do with this? (right now effectively a noop)
    low_pred_idx <- integer(0)


    # Transform data for model fitting ---------------------

    # Training data for plr penalty parameter

    x_numer_penal <- x_numer[penal_idx, , drop=FALSE]
    y_penal <- y[penal_idx]

    x_categ_penal <- x_categ[penal_idx, , drop=FALSE]
    x_derive_penal <- x_derive[penal_idx, , drop=FALSE]

    pm_penal_bool <- lapply(pm_bool, function(w) w[penal_idx])

    # Calculate transformed numeric data
    x_trans_penal <- mqc_transform_dist(x_numer_penal,
                                        quant_dat,
                                        pm_penal_bool,
                                        pm_col_idx,
                                        rem_cont_idx,
                                        low_pred_idx,
                                        std_parts)

    # Collect data into a single matrix
    x_combine_penal <- mqc_combine(x_trans_penal,
                                   x_categ_penal,
                                   x_derive_penal,
                                   x_numer_penal,
                                   aug,
                                   keep_derive)


    # Choose penalty parameter -----------------------------

    if (identical(Sys.info()["nodename"], c(nodename="dpritch-Satellite-C875"))) {
        use_parallel <- TRUE
        doMC::registerDoMC(cores=4)
    } else {
        use_parallel <- FALSE
    }

    # If data size < 25 the perform LOOCV, otherwise perform 10-fold CV
    n_penal <- length(penal_idx)
    nfolds <- ifelse(n_penal < 25, n_penal, 10)

    # Defaults to alpha=1, i.e. lasso
    fit <- glmnet::cv.glmnet(x            = x_combine_penal,
                             y            = y_penal,
                             type.measure = cv_type,
                             nfolds       = nfolds,
                             grouped      = FALSE,
                             parallel     = use_parallel,
                             family       = "binomial")


    # Return data ------------------------------------------

    list(quant_dat     = quant_dat,
         fit           = fit,
         aug           = aug,
         keep_derive   = keep_derive,
         pred_rem_lev  = pred_rem_lev,
         numer_col_idx = numer_col_idx,
         pm_col_idx    = pm_col_idx,
         pm_vals       = pm_vals,
         rem_cont_idx  = rem_cont_idx,
         low_pred_idx  = low_pred_idx,
         quant_type    = quant_type,
         cv_type       = cv_type,
         simil_type    = simil_type,
         std_parts     = std_parts,
         class_nm      = classlev)
}




predict_mqc_ <- function(mqcObj, newx, pred_type="class") {

    quant_dat     <- mqcObj$quant_dat
    fit           <- mqcObj$fit
    aug           <- mqcObj$aug
    keep_derive   <- mqcObj$keep_derive
    pred_rem_lev  <- mqcObj$pred_rem_lev
    numer_col_idx <- mqcObj$numer_col_idx
    pm_col_idx    <- mqcObj$pm_col_idx
    pm_vals       <- mqcObj$pm_vals
    rem_cont_idx  <- mqcObj$rem_cont_idx
    low_pred_idx  <- mqcObj$low_pred_idx
    quant_type    <- mqcObj$quant_type
    cv_type       <- mqcObj$cv_type
    simil_type    <- mqcObj$simil_type
    std_parts     <- mqcObj$std_parts

    # Split off categorical variables
    x_numer_orig <- as.matrix( newx[, numer_col_idx,  drop=FALSE] )
    x_categ <- newx[, setdiff(1:NCOL(newx), numer_col_idx),  drop=FALSE]

    n <- NROW(x_numer_orig)
    p_numer <- NCOL(x_numer_orig)

    # Convert point mass vars to categorical / continuous split
    convert_dat <- mqc_convert_pm(x_numer_orig, pm_col_idx, pm_vals, rem_cont_idx, low_pred_idx)

    x_numer  <- convert_dat$x_numer
    x_derive <- convert_dat$x_derive
    pm_bool  <- convert_dat$pm_bool

    # Transform numeric data to quantile distance difference
    x_trans <- mqc_transform_dist(x_numer, quant_dat, pm_bool, pm_col_idx,
                                  rem_cont_idx, low_pred_idx, std_parts)

    # Collect data into a single matrix
    x_combine <- mqc_combine(x_trans, x_categ, x_derive, x_numer, aug, keep_derive)

    # Predict new data observations
    predict(fit, x_combine, type=pred_type, s="lambda.min")
}




mqc <- function(x,
                y,
                aug = TRUE,
                keep_derive = TRUE,
                split_prop = 0.5,
                pred_rem_lev = 0.65,
                categ_prop_lev = 0.2,
                theta = seq(0.01, 0.99, 0.01),
                provide_quantlev = NULL,
                quant_type = "interp",
                cv_type = "class",
                simil_type = "dist",
                npart = 10,
                var_lev = 0.8,
                std_parts = "norm",
                nboot = 21) {

    if (is.null(split_prop)) {
        nboot <- 1L
    }

    out <- list()

    for (i in seq_len(nboot)) {
        out[[ i ]] <- mqc_(x, y, aug, keep_derive, split_prop, pred_rem_lev,
                           categ_prop_lev, theta, provide_quantlev, quant_type,
                           cv_type, simil_type, npart, var_lev, std_parts)
    }

    out
}




predict_mqc <- function(mqcObj, newx, classify_type="probs") {

    out_boot <- matrix(0, NROW(newx), length(mqcObj))

    if (isTRUE(all.equal(classify_type, "probs"))) {
        pred_type <- "response"
    }
    else if (isTRUE(all.equal(classify_type, "class"))) {
        pred_type <- "class"
    }

    for (i in seq_along(mqcObj)) {

        out_boot[, i] <- predict_mqc_(mqcObj[[ i ]], newx, pred_type)

        # quantlev      <- mqcObj[[ i ]]$quantlev
        # c0_quant      <- mqcObj[[ i ]]$c0_quant
        # c1_quant      <- mqcObj[[ i ]]$c1_quant
        # fit           <- mqcObj[[ i ]]$fit
        # aug           <- mqcObj[[ i ]]$aug
        # keep_derive   <- mqcObj[[ i ]]$keep_derive
        # pred_rem_lev  <- mqcObj[[ i ]]$pred_rem_lev
        # numer_col_idx <- mqcObj[[ i ]]$numer_col_idx
        # pm_col_idx    <- mqcObj[[ i ]]$pm_col_idx
        # pm_vals       <- mqcObj[[ i ]]$pm_vals
        # rem_cont_idx  <- mqcObj[[ i ]]$rem_cont_idx
        # low_pred_idx  <- mqcObj[[ i ]]$low_pred_idx
        # quant_type    <- mqcObj[[ i ]]$quant_type
        # cv_type       <- mqcObj[[ i ]]$cv_type

        # # Split off categorical variables
        # x_numer_orig <- as.matrix( x[, numer_col_idx,  drop=FALSE] )
        # x_categ <- x[, setdiff(1:NCOL(x), numer_col_idx),  drop=FALSE]

        # n <- NROW(x_numer_orig)
        # p_numer <- NCOL(x_numer_orig)

        # # Convert point mass vars to categorical / continuous split
        # convert_dat <- mqc_convert_pm(x_numer_orig, pm_col_idx, pm_vals,
        #                               rem_cont_idx, low_pred_idx)

        # x_numer  <- convert_dat$x_numer
        # x_derive <- convert_dat$x_derive
        # pm_bool  <- convert_dat$pm_bool

        # # Transform numeric data to quantile distance difference
        # x_trans <- mqc_transform_dist(x_numer, c0_quant, c1_quant, quantlev,
        #                               pm_bool, pm_col_idx, rem_cont_idx,
        #                               low_pred_idx, std_parts)

        # # Collect data into a single matrix
        # x_combine <- mqc_combine(x_trans, x_categ, x_derive, x_numer, aug, keep_derive)

        # if (isTRUE(all.equal(classify_type, "class"))) {
        #     out_boot[, i] <- predict(fit, x_combine, type="class", s="lambda.min")
        # }
        # else if (isTRUE(all.equal(classify_type, "probs"))) {
        #     out_boot[, i] <- predict(fit, x_combine, type="response", s="lambda.min")
        # }
        # else {
        #     stop("illegal arg")
        # }
    }


    class_nm <- mqcObj[[ 1L ]]$class_nm

    if (isTRUE(all.equal(classify_type, "class"))) {
        y_pred <- apply(out_boot, 1, function(w) {
            if (mean(w == class_nm[2L]) > 0.5) {
                class_nm[2L]
            } else {
                class_nm[1L]
            }
        })
    }
    else {
        y_pred <- rep(class_nm[1L], NROW(newx))
        y_pred[apply(out_boot, 1, mean) > 0.5] <- class_nm[2L]
    }

    y_pred
}

