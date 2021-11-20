## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- echo = FALSE------------------------------------------------------------
library(formatR)

## ---- eval=FALSE--------------------------------------------------------------
#  install.packages("SIS", repos = "http://cran.us.r-project.org")

## -----------------------------------------------------------------------------
library(SIS)

## ---- tidy=TRUE, tidy.opts=list(width.cutoff=70)------------------------------
set.seed(0)
n = 400; p = 50; rho = 0.5
corrmat = diag(rep(1-rho, p)) + matrix(rho, p, p)
corrmat[,4] = sqrt(rho)
corrmat[4, ] = sqrt(rho)
corrmat[4,4] = 1
corrmat[,5] = 0
corrmat[5, ] = 0
corrmat[5,5] = 1
cholmat = chol(corrmat)
x = matrix(rnorm(n*p, mean=0, sd=1), n, p)
x = x%*%cholmat

# gaussian response 
set.seed(1)
b = c(4,4,4,-6*sqrt(2),4/3)
y=x[, 1:5]%*%b + rnorm(n)

## ---- tidy=TRUE, tidy.opts=list(width.cutoff=70)------------------------------
# SIS without regularization
model10 = SIS(x, y, family='gaussian', iter = FALSE)

# Getting the final selected variables after regularization step
model10$ix

# Getting the ranked list of variables from the screening step
model10$sis.ix0

# The top 10 ranked variables from the screening step
model10$ix0[1:10]

## ---- tidy=TRUE, tidy.opts=list(width.cutoff=70)------------------------------
# SIS with regularization
model11 = SIS(x, y, family='gaussian', penalty = 'SCAD', iter = TRUE)

# Getting the final selected variables
model10$ix

# The top 10 ranked variables for the final screening step
model11$ix0[1:10]

# The top 10 ranked variables for each screening step
lapply(model11$ix_list,f<-function(x){x[1:10]})

## ---- tidy=TRUE, tidy.opts=list(width.cutoff=70)------------------------------
set.seed(2)
feta <- x[, 1:5] %*% b
fprob <- exp(feta) / (1 + exp(feta))
y <- rbinom(n, 1, fprob)
model21 <- SIS(x, y, family = "binomial", tune = "bic")

# Getting the final selected variables
model21$ix

# The top 10 ranked variables for the final screening step
model11$ix0[1:10]

# The top 10 ranked variables for each screening step
lapply(model11$ix_list,f<-function(x){x[1:10]})

## ---- tidy=TRUE, tidy.opts=list(width.cutoff=70)------------------------------
set.seed(4)
b <- c(4, 4, 4, -6 * sqrt(2), 4 / 3)
myrates <- exp(x[, 1:5] %*% b)
Sur <- rexp(n, myrates)
CT <- rexp(n, 0.1)
Z <- pmin(Sur, CT)
ind <- as.numeric(Sur <= CT)
y <- survival::Surv(Z, ind)
model41 <- SIS(x, y,
  family = "cox", penalty = "lasso", tune = "bic",
  varISIS = "aggr", seed = 41
)
model42 <- SIS(x, y,
  family = "cox", penalty = "lasso", tune = "bic",
  varISIS = "cons", seed = 41
)
model41$ix
model42$ix

## ---- tidy=TRUE, tidy.opts=list(width.cutoff=70)------------------------------
y <- as.factor(iris$Species)
noise <- matrix(rnorm(nrow(iris)*200),nrow(iris),200)
x <- cbind(as.matrix(iris[,-5]),noise)

model21 <- SIS(x, y, family = "multinom", penalty = 'lasso')

# Getting the final selected variables
model21$ix

# The top 10 ranked variables for the final screening step
model21$ix0[1:10]

# The top 10 ranked variables for each screening step
lapply(model21$ix_list,f<-function(x){x[1:10]})

