library(glmnet)
library(survival)
# hack for testing
surv_clin$status[surv_clin$status == 2] <- 1


cvgfit <- cv.glmnet(
  x = as.matrix(surv_clin[, -c(1, 2)]),
  y =  Surv(surv_clin$time, surv_clin$status),
  family = "cox",
  #foldid = fwelnet::stratified_cv_folds(surv_geno)$fold
)

plot(cvgfit)

surv_x <- model.matrix( ~ -1 + ., surv_clin[, -c(1, 2)])
cv.glmnet(
  x = surv_x,
  y =  Surv(surv_clin$time, surv_clin$status),
  family = "cox"
)

survival::veteran

xx <- model.matrix( ~ -1 + ., survival::veteran[, (!names(survival::veteran) %in% c("time", "status"))])
glmnet::cv.glmnet(
  x = xx,
  y = Surv(survival::veteran$time, survival::veteran$status),
  family = "cox"
)


xt <- data.frame(
  x2 = sample(c("A", "B"), 10, TRUE),
  x3 = sample(c("E", "F", "G"), 10, TRUE),
  y = rnorm(10)
)

model.matrix( ~ -1 + ., data = xt)
model.matrix( ~ ., data = xt)
