
#library(devtools)

library(glmnet)        # plr
library(quantileDA)    # Hennig classifier
library(e1071)         # svm, naive bayes
library(penalizedLDA)  # penalized LDA
library(class)         # knn
library(pamr)          # nearest shrunken centroids
library(rpart)         # decision trees



if (! identical(Sys.info()["nodename"], c(nodename="dpritch-Satellite-C875"))) {
    setwd("~/MQC")
    source("R/aggregate.R")
    source("R/mqc-classrate.R")
    source("R/mqc-quantile.R")
    source("R/mqc.R")
    source("R/mqc-routines.R")
    #load_all()
}
source("tests/numerical/data_sims.R")
source("tests/numerical/compare.R")

NREPL <- 25
NTEST <- 1000
r <- 50

out_exp_gauss <- list()

for (n in c(500, 250, 100, 50)) {
    for (p in c(50, 100, 250, 500)) {

        train_fcn <- quote( sim_exp_gauss(n, p, r, 0.5, 1, 0.8, 1) )
        test_fcn <- quote( sim_exp_gauss(NTEST, p, r, 0.5, 1, 0.8, 1) )

        nm <- paste0("n", n, "_p", p, "_r", r)

        out_exp_gauss[[ nm ]]  <- compare(NREPL,
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
    }
}


save(out_exp_gauss, file="tests/numerical/results/exp_gauss.RData")
