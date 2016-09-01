
# email spam example -----------------------------------------------------------

library(ElemStatLearn)
data(spam)

# Create a list of vectors of length len, each have the proportion of TRUEs
# given by an element of prop.  The generation of the vectors is random.

sample_splits <- function(len, prop) {
  lapply(prop, function(x) {
    n_tru <- floor(x * len)
    n_fal <- len - n_tru
    tf_vals <- c(rep(TRUE, n_tru), rep(FALSE, n_fal))
    tf_vals[ sample.int(len, len) ]
  })
}

# split_list <- sample_splits(nrow(spam), c(0.05, 1:8 / 10))


# obs <- spam
# class_id <- "spam"
# split <- split_list[[6]]
# sd_min <- 0.05
# epsilon <- 0.25
# quants <- seq(0.1, 0.9, 0.1)



# Regular FAN approach ---------------------------------------------------------

split_rep_50 <- sample_splits(nrow(spam), rep(0.5, 100))

x05 <- sapply(sample_splits(nrow(spam), rep(0.05, 100)), function(x)
  mqc(spam, "spam", x, sd_min=0.05, epsilon=0.25, method="fan"))

x10 <- sapply(sample_splits(nrow(spam), rep(0.10, 100)), function(x)
  mqc(spam, "spam", x, sd_min=0.05, epsilon=0.25, method="fan"))

x30 <- sapply(sample_splits(nrow(spam), rep(0.30, 100)), function(x)
  mqc(spam, "spam", x, sd_min=0.05, epsilon=0.25, method="fan"))

x40 <- sapply(sample_splits(nrow(spam), rep(0.40, 100)), function(x)
  mqc(spam, "spam", x, sd_min=0.05, epsilon=0.25, method="fan"))

x50 <- sapply(sample_splits(nrow(spam), rep(0.50, 100)), function(x)
  mqc(spam, "spam", x, sd_min=0.05, epsilon=0.25, method="fan"))

x70 <- sapply(sample_splits(nrow(spam), rep(0.70, 100)), function(x)
    mqc(spam, "spam", x, sd_min=0.05, epsilon=0.25, method="fan"))

x80 <- sapply(sample_splits(nrow(spam), rep(0.80, 100)), function(x)
  mqc(spam, "spam", x, sd_min=0.05, epsilon=0.25, method="fan"))






# FAN2 approach ----------------------------------------------------------------

x05_augm <- sapply(sample_splits(nrow(spam), rep(0.05, 50)), function(x)
  mqc(spam, "spam", x, sd_min=0.05, epsilon=0.25, method="fan", augm=TRUE))

x10_augm <- sapply(sample_splits(nrow(spam), rep(0.10, 50)), function(x)
  mqc(spam, "spam", x, sd_min=0.05, epsilon=0.25, method="fan", augm=TRUE))

x30_augm <- sapply(sample_splits(nrow(spam), rep(0.30, 50)), function(x)
  mqc(spam, "spam", x, sd_min=0.05, epsilon=0.25, method="fan", augm=TRUE))

x50_augm <- sapply(sample_splits(nrow(spam), rep(0.50, 50)), function(x)
  mqc(spam, "spam", x, sd_min=0.05, epsilon=0.25, method="fan", augm=TRUE))

x80_augm <- sapply(sample_splits(nrow(spam), rep(0.80, 50)), function(x)
  mqc(spam, "spam", x, sd_min=0.05, epsilon=0.25, method="fan", augm=TRUE))




# Hennig approach --------------------------------------------------------------

x50_hen <- sapply(sample_splits(nrow(spam), rep(0.50, 50)), function(trBool) {
  quantfit <- quantileDA::quantilecl(train   = spam[trBool, -58],
                                     test    = spam[!trBool, -58],
                                     cl      = spam[trBool, 58],
                                     cl.test = spam[!trBool, 58])
  quantfit$me.test
})

# A single run
split <- unlist( sample_splits(nrow(spam), 0.5) )
test_hennig <- quantileDA::quantilecl(train=spam[split, -58],
                                      test=spam[!split, -58],
                                      cl=spam[split, 58],
                                      cl.test=spam[!split, 58])
test_hennig$me.test




# Hennig extension approach ----------------------------------------------------

obs <- spam
class_id <- "spam"
split <- unlist( sample_splits(nrow(spam), 0.5) )
quants <- seq(0.1, 0.9, 0.1)
sd_min  <- -1
epsilon <- 0.25
method <- "hall"

test <- mqc(obs, class_id, split, sd_min, epsilon, quants, method)
