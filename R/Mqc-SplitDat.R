
mqc_splitDat <- function(obs, class_id, split) {

  if ( is.matrix(obs) || is.data.frame(obs) ) {
    return ( mqc_splitDat_single(obs, class_id, split) )
  }
  else {
    return ( list( trn = mqc_splitDat_single(obs$train, class_id, split),
                   tes = mqc_splitDat_single(obs$test, class_id, split) ) )
  }
}




mqc_splitDat_single <- function(obs, class_id, split) {

  # Construct var specifying class information
  class_info <- mqc_class_info(obs, class_id)

  # case: obs is exactly one of train or test data, and hence does not need to
  # be split into seperate datasets
  if ( is.null(split) ) {
    splitObj <- list(
      obs      = obs[, -class_info$classIdx],
      class_id = class_info$class_vec
    )
  }

  # case: obs contains both train and data which need to be split into seperate
  # datasets.
  else {
    # construct var specifying split information
    split_info <- mqc_split_info(obs, split)

    # remIdx: a length-0, 1, or 2 vector with the indices of the columns in obs
    # containing the class info and/or split info (if any)
    remIdx <- c(class_info$classIdx, split_info$splitIdx)

    splitObj <- list()

    splitObj$train <- list(
      obs      = obs[split_info$split_vec, setdiff(1:NCOL(obs), remIdx)],
      class_id = class_info$class_id_vec[split_info$split_vec]
    )

    splitObj$test <- list(
      obs      = obs[!split_info$split_vec, setdiff(1:NCOL(obs), remIdx)],
      class_id = class_info$class_id_vec[!split_info$split_vec]
    )
  }

  return (splitObj)
}




mqc_split_info <- function(obs, split) {

  # case: split is length-1 and tells us which column the class info is located in
  if ( identical(length(split), 1L) ) {

    # case: split is length-1 numeric vector
    if (is.numeric(split)) {
      splitIdx <- as.integer(split)
      if ((splitIdx < 1L) || (splitIdx > NCOL(obs))) {
        stop(paste("If the argument for split is a single numeric value then it ",
                   "cannot be smaller than one or larger than the number of ",
                   "columns in the argument for obs"), call.=FALSE)
      }
    } # end split is a length-1 numeric vector case

    # case: split is a length-1 character vector
    else {
      splitIdx <- grep(split, colnames(obs))
      if ( !identical(length(splitIdx), 1L) ) {
        stop(paste("If the argument for split is a single character string, ",
                   "then it must have exactly 1 partial match in the column ",
                   "names for the argument for obs"), call.=FALSE)
      }
    } # end split is a length-1 character vector case

    split_vec <- obs[, splitIdx]
  } # end group length is 1 case

  # case: split is a vector of the same length as obs
  else {
    split_vec <- split
    splitIdx <- integer(0)
  }

  list( split_vec = split_vec,
        splitIdx  = splitIdx )
}













# Returns a vector containing the class information with length equal to the
# number of observations in obs

mqc_class_info <- function(obs, class_id) {

  # case: class_id is length-1 and tells us which column the class info is located in
  if ( identical(length(class_id), 1L) ) {

    # case: class_id is length-1 numeric vector
    if (is.numeric(class_id)) {
      classIdx <- as.integer(class_id)
      if ((classIdx < 1L) || (classIdx > obs_ncol)) {
        stop(paste("If the argument for class_id is a single numeric value ",
                   "then it cannot be smaller than one or larger than the ",
                   "number of columns in the argument for obs"), call.=FALSE)
      }
    } # end class_id is a length-1 numeric vector case

    # case: class_id is a length-1 character vector
    else {
      classIdx <- grep(class_id, colnames(obs))
      if ( !identical(length(classIdx), 1L) ) {
        stop(paste("If the argument for class_id is a single character string ",
                   "then it must have exactly 1 partial match in the column ",
                   "names for the argument for obs"), call.=FALSE)
      }
    } # end class_id is a length-1 character vector case
    class_id_vec <- obs[, classIdx]
  } # end group length is 1 case

  # case: class_id has length equal to the number of observations and provides
  # the class info data
  else {
    class_id_vec <- as.character(class_id)
    classIdx <- integer(0)
  }

  list( class_id_vec = class_id_vec,
        classIdx     = classIdx )
}




