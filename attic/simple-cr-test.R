# if not using renv, update to latest fwelnet fork
# remotes::install_github("jemus42/fwelnet")
#

library(survival)
library(data.table)

xdat <- survival::pbc[, -1]
xdat <- as.data.table(xdat)
xdat[, sex := as.integer(sex == "m")]
dim(xdat)
xdat <- na.omit(xdat)
dim(xdat)

str(xdat)
table(xdat$status)

# Note: input data needs to have (time, status) columns, status in (0,1,2)
# predictors need to be numeric/binary, dummy-encoding needs to be done beforehand because we run into
# issues with Z otherwise
# Also: no missing values
fit <- fwelnet::fwelnet_mt_cox(
  xdat,
  alpha = 1,
  mt_max_iter = 3, # how many multi-task iterations, keep low since no convergence anyway
  a = 0.5,
  t = 100, # start learning rate for theta optim
  thresh = 1e-7, # threshold for theta optim, the lower the longer, better ods of getting non-shit theta
  stratify_by_status = TRUE, # internal CV stratified by "status" variable for low-prevalence tasks
  nfolds = 3, # number of folds in internal CV, default 10 but might cause issues with small case counts
  include_mt_beta_history = TRUE # store history of betas per cause in $beta1, $beta2
)

# Extract coefficients, first col is initial cause-specific glmnet fit
glmnet_beta1 <- fit$beta1[, 1]
glmnet_beta2 <- fit$beta2[, 1]
# Same for fwelnet estimates, last col are final coefs
fwel_beta1 <- fit$beta1[, ncol(fit$beta1)]
fwel_beta2 <- fit$beta2[, ncol(fit$beta2)]

# use predict.cv.fwelnet method directly, gives linear predictor, see fwelnet:::predict.cv.fwelnet calling fwelnet:::predict.fwelnet
predict(fit$fwfit1, xnew = xdat[, -c(1,2)], s = "lambda.min")
