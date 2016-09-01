
# prostate example -------------------------------------------------------------

# pg 49 ESL
library(ElemStatLearn)
data(prostate)

# gleason and pgg45 are categorical vars
prostate <- prostate[, setdiff(colnames(prostate),  c("gleason", "pgg45"))]

obs <- prostate[, setdiff(colnames(prostate), c("train"))]
class_id <- "svi"
split <- prostate[, "train"]
sd_min <- 0.05
epsilon <- 0.25
quants <- seq(0.1, 0.9, 0.1)
method <- "fan"
augm <- TRUE

prostate_mqc <- mqc(prostate[, setdiff(colnames(prostate), c("train"))],
                    "svi", prostate[, "train"], 0.25, 0.25, augm=TRUE)

## Using qda

library(MASS)

prostate_qda <- qda(svi ~ lcavol + lweight + age + lbph + lcp + lpsa,
                    data=prostate[prostate$train, ])
prostate_qda_pred <- predict(prostate_qda, prostate[!prostate$train, ])

table(prostate_qda_pred$class, prostate[!prostate$train, "svi"])


train <- list(x = subset(prostate, train, c("lcavol", "lweight", "age", "lbph", "lcp", "lpsa")),
              y = with(prostate, svi[train]))
test <- list(x  = subset(prostate, !train, c("lcavol", "lweight", "age", "lbph", "lcp", "lpsa")),
             y  = with(prostate, svi[!train]))






# bone example -----------------------------------------------------------------

library(ElemStatLearn)

data(bone, package="ElemStatLearn")

data(orange4.train, package="ElemStatLearn")
data(orange4.test, package="ElemStatLearn")
data(orange10.train, package="ElemStatLearn")













