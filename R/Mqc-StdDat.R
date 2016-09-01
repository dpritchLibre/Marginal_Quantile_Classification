
mqc_stdDat <- function(splitDat, sd_min, standardize=TRUE) {

  if (!standardize) {
    splitDat$useVar <- rep(TRUE, ncol(splitDat$train$obs))
    return (splitDat)
  }
  
  trn_obs <- splitDat$train$obs
  tes_obs <- splitDat$test$obs
  trn_nr <- nrow(trn_obs)

  meanVals <- colMeans(trn_obs)

  # tr_obs_center: column-centered training data
  trn_obs_center <- sweep(trn_obs, 2, meanVals)
  tes_obs_center <- sweep(tes_obs, 2, meanVals)

  # st_dev: training data column standard deviations
  st_dev <- sqrt( apply(trn_obs_center, 2, crossprod) / (trn_nr - 1) )

  # If standard deviation is too low then we don't use variable
  useVar <- (st_dev > sd_min)

  # Change small standard deviation values to ensure that later we are not
  # dividing by 0 when standardizing the data.  This changes more than just
  # those values, but since we are not using those variables anyway this not a
  # problem.
  st_dev[!useVar] <- 1

  # Create train and test data with observations standardized by the training
  # data
  stdDat <- list()
  stdDat$train <- list(
    obs      = sweep(trn_obs_center, 2, st_dev, "/"),
    class_id = splitDat$train$class_id
  )
  stdDat$test <- list(
    obs      = sweep(tes_obs_center, 2, st_dev, "/"),
    class_id = splitDat$test$class_id
  )
  stdDat$useVar <- useVar

  return (stdDat)
}


# mqc_stdDat <- function(splitDat) {
#
#   if ( identical(length(split), 1L) ) {
#     return ( mqc_getStdDat_both(splitDat) )
#   }
#   else {
#     return ( mqc_getStdDat_sep(splitDat) )
#   }
# }
#
#
#
#
#
# mqc_stdDat_both <- function(splitDat) {
#   obs <- splitDat$obs
#   grp <- splitDat$grp
#   split <- splitDat$split
#   remIdx <- splitDat$remIdx
#
#   # idx: the column indexes in obs which contain data
#   idx <- setdiff(seq_len( NCOL(obs) ), remIdx)
#   p <- length(idx)
#
#   # Find the number of groups
#   grp_uniq <- unique(grp[split])
#   grp_nm <- paste0("grp_", grp_uniq)
#   nGrp <- length(grp_nm)
#   if ((nGrp <= 1L) || identical(nGrp, sum(split))) {
#     stop(paste("The number of groups in the training set must be",
#                "greater than 1 and less than the number of obserations\n"))
#   }
#
#   # Calculate training data means / sd
#   if (is.matrix(obs)) {
#     trn_mn <- apply(obs[split, idx], 2, mean)
#     trn_sd <- apply(obs[split, idx], 2, sd)
#   }
#   else {
#     trn_mn <- sapply(obs[split, idx], mean)
#     trn_sd <- sapply(obs[split, idx], sd)
#   }
#
#   # Calculate indices for group in training set
#   grpIdx <- lapply(grp_uniq, function(x) which((grp == x) & split))
#
#   # Create training data separated by group (and centered and scaled)
#   trn <- vector("list", nGrp)
#   for (i in seq_len(nGrp)) {
#     trn[[i]] <- scale(obs[grpIdx[[i]], idx], center=trn_mn, scale=trn_sd)
#   }
#   names(trn) <- grp_nm
#
#   # Create centered and scaled test data
#   tes <- scale(obs[!split, idx], center=trn_mn, scale=trn_sd)
#
#   list( trn = trn,
#         tes = tes )
# }
#
#
#
#
# mqc_stdDat_sep <- function(splitDat) {
#
#   # trnDatIdx: the column indexes in splitDat$trn$obs which contain data
#   trnDatIdx <- setdiff(seq_len( NCOL(splitDat$trn$obs) ), splitDat$trn$remIdx)
#
#   # tesDatIdx: the column indexes in splitDat$tes$obs which contain data
#   tesDatIdx <- setdiff(seq_len( NCOL(splitDat$tes$obs) ), splitDat$tes$remIdx)
#
#   # Find the dimension of the observations
#   p <- length(trnDatIdx)
#   q <- length(tesDatIdx)
#   if ( !identical(p, q) ) {
#     stop("Training and test data must have the same number of columns of data\n")
#   }
#
#   # Find the number of groups
#   grp_uniq <- unique(splitDat$trn$grp)
#   grp_nm <- paste0("grp_", grp_uniq)
#   nGrp <- length(grp_nm)
#   if (nGrp <= 1L) {
#     stop("The number of groups in the training set must be greater than 1\n")
#   }
#
#   # Calculate training data means / sd
#   if (is.matrix(splitDat$trn$obs)) {
#     trn_mn <- apply(splitDat$trn$obs, 2, mean)[trnDatIdx]
#     trn_sd <- apply(splitDat$trn$obs, 2, sd)[trnDatIdx]
#   }
#   else {
#     trn_mn <- sapply(splitDat$trn$obs, mean)[trnDatIdx]
#     trn_sd <- sapply(splitDat$trn$obs, sd)[trnDatIdx]
#   }
#
#   # Calculate indices for group in training set
#   grpIdx <- lapply(grp_uniq, function(x) which(splitDat$trn$grp == grp_uniq))
#
#   # Create training data separated by group (and centered and scaled)
#   trn <- vector("list", nGrp)
#   for (i in seq_len(nGrp)) {
#     trn[[i]] <- (splitDat$trn$obs[grpIdx[[i]], trnDatIdx] - trn_mn) / trn_sd
#   }
#   names(trn) <- grp_nm
#
#   # Create centered and scaled test data
#   tes <- (splitDat$tes$obs - trn_mn) / trn_sd
#   tes <- tes[tesDatIdx]
#
#   list( trn = trn,
#         tes = tes )
# }







