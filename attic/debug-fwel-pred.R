# Trying to get predictions for evaluation out of fwelnet_mt object
library(riskRegression)
library(survival)
library(data.table)
library(fwelnet) # remotes::install_github("jemus42/fwelnet")

# Low-dimensional train-/test data --------------------------------------------------------------------------------
# set.seed(8)
# dtrain <- sampleData(1000, outcome = "competing.risk")
# dtest <- sampleData(400, outcome = "competing.risk")
#
# # Only keep time, event + predictors, rename event to status for fwelnet WIP/NYI reasons
# dtrain[, `:=`(eventtime1 = NULL, eventtime2 = NULL, censtime = NULL, status = event, event = NULL)]
# dtest[, `:=`(eventtime1 = NULL, eventtime2 = NULL, censtime = NULL, status = event, event = NULL)]
#
# # There's something fishy going on with factor encoding in fwelnet_mt which I'll have to debug later
# # So for now we integerize the factor variables to (0, 1)
# dtrain[, `:=`(X1 = as.integer(X1) - 1, X2 = as.integer(X2) - 1, X3 = as.integer(X3) - 1,
#               X4 = as.integer(X4) - 1, X5 = as.integer(X5) - 1)]
# dtest[, `:=`(X1 = as.integer(X1) - 1, X2 = as.integer(X2) - 1, X3 = as.integer(X3) - 1,
#              X4 = as.integer(X4) - 1, X5 = as.integer(X5) - 1)]
#
# # Take all predictors, since `y ~ .` doesn't work we have to get weird about it, sorry
# fpreds <- names(dtrain)[startsWith(names(dtrain), "X")]
# fformula <- as.formula(sprintf("Hist(time, status) ~ %s", paste0(fpreds, collapse = " + ")))
# # Formula checks out
# fformula
#
# # Fit regular CSC as sacrificial model
# m_csc = CSC(formula = fformula, data = dtrain, singular.ok = TRUE, iter.max = 1)
#
# # Initiate models for fwelnet and glmnet
# m_fwelnet <- m_csc
# m_glmet <- m_csc
#
# # Original coefficients
# m_csc$models[[1]]$coefficients
# m_csc$models[[2]]$coefficients
#
# # Fit fwelnet_mt model
# mt_max_iter <- 2 # How many multi-task iterations to fit - "step 1" is always original glmnet solution
# m_fwelmt <- fwelnet::fwelnet_mt_cox(
#   dtrain, mt_max_iter = mt_max_iter, alpha = 1, t = 100, thresh = 1e-3,
#   include_mt_beta_history = TRUE
#   )
#
# # fwelnet_mt model fit includes "final" coefficients per model for fwelnet and underlying original cs-glmnet solution
# m_glmet$models[[1]]$coefficients <- m_fwelmt$beta1[, 1]
# m_glmet$models[[2]]$coefficients <- m_fwelmt$beta2[, 1]
#
# m_fwelnet$models[[1]]$coefficients <- m_fwelmt$beta1[, mt_max_iter + 1]
# m_fwelnet$models[[2]]$coefficients <- m_fwelmt$beta2[, mt_max_iter + 1]
#
# mod_scores <- Score(
#   list(
#     cs_glmnet = m_glmet,
#     cs_fwelnet = m_fwelnet,
#     csc = m_csc
#   ),
#   formula = Hist(time, status) ~ 1, data = dtest
# )
#
# summary(mod_scores)
#

# If CSC is fit on subset of available predictors -----------------------------------------------------------------
# Since we're looking at high-dim data with 5000 predictors, fitting regular CSC on it won't work,
# but fitting fwelnet_mt will. So can we fudge more coefficients into the sacrifical CSC model than that itself
# was fit on?


set.seed(8)
dtrain <- sampleData(1000, outcome = "competing.risk")
dtest <- sampleData(400, outcome = "competing.risk")

# Only keep time, event + predictors, rename event to status for fwelnet WIP/NYI reasons
dtrain[, `:=`(eventtime1 = NULL, eventtime2 = NULL, censtime = NULL, status = event, event = NULL)]
dtest[, `:=`(eventtime1 = NULL, eventtime2 = NULL, censtime = NULL, status = event, event = NULL)]

# There's something fishy going on with factor encoding in fwelnet_mt which I'll have to debug later
# So for now we integerize the factor variables to (0, 1)
dtrain[, `:=`(X1 = as.integer(X1) - 1, X2 = as.integer(X2) - 1, X3 = as.integer(X3) - 1,
              X4 = as.integer(X4) - 1, X5 = as.integer(X5) - 1)]
dtest[, `:=`(X1 = as.integer(X1) - 1, X2 = as.integer(X2) - 1, X3 = as.integer(X3) - 1,
             X4 = as.integer(X4) - 1, X5 = as.integer(X5) - 1)]

# Fit regular CSC as sacrificial model, fit on only 3 predictors
m_csc = CSC(formula = Hist(time, status) ~ X1 + X2 + X9, data = dtrain, singular.ok = TRUE)

m_glmet <- m_csc

dummy_model_matrix <- model.matrix(prodlim::Hist(time, status) ~ . -1, data = dtrain)

# assign object is created from formula, relevant especially if factors are dummy-coded
# survival::coxph uses attrassign internally to construct that same object
dummy_assign <- survival::attrassign(dummy_model_matrix, terms(Surv(time, status) ~ ., data = dtrain))
m_glmet$models[[1]]$assign <- dummy_assign
m_glmet$models[[2]]$assign <- dummy_assign

# Have to remove that attribute apparently, but it has to be present for attrassign beforehand
# attr(dummy_model_matrix, "assign") <- NULL

m_glmet$models[[1]]$x <- dummy_model_matrix
m_glmet$models[[2]]$x <- dummy_model_matrix

m_glmet$models[[1]]$means <- colMeans(m_glmet$models[[1]]$x)
m_glmet$models[[2]]$means <- colMeans(m_glmet$models[[2]]$x)

# Do we have to fusge a covariance matrix of the correct length?
# p <- length(m_glmet$models[[1]]$means)
#
# m_glmet$models[[1]]$var <- matrix(numeric(p * p), ncol = p)
# m_glmet$models[[2]]$var <- matrix(numeric(p * p), ncol = p)


# terms object can be created from formula fairly easily but has to be consistent with actual data/coefs
m_glmet$models[[1]]$terms <- terms(Surv(time, status) ~ ., data = dtrain)
m_glmet$models[[2]]$terms <- terms(Surv(time, status) ~ ., data = dtrain)

# copy fudged CSC-franken-glmnet object to hold equally fudged fwelnet object
m_fwelnet <- m_glmet

# Original coefficients for checking
# m_csc$models[[1]]$coefficients
# m_csc$models[[2]]$coefficients

# Fit fwelnet_mt model, using all predictors in dtrain (10)
mt_max_iter <- 5 # How many multi-task iterations to fit - "step 1" is always original glmnet solution
m_fwelmt <- fwelnet::fwelnet_mt_cox(
  dtrain, mt_max_iter = mt_max_iter,
  alpha = 1, t = 100, thresh = 1e-7, include_mt_beta_history = TRUE
)

# fwelnet_mt model fit includes "final" coefficients per model for fwelnet and underlying original cs-glmnet solution
m_glmet$models[[1]]$coefficients <- m_fwelmt$beta1[, 1]
m_glmet$models[[2]]$coefficients <- m_fwelmt$beta2[, 1]

m_fwelnet$models[[1]]$coefficients <- m_fwelmt$beta1[, mt_max_iter + 1]
m_fwelnet$models[[2]]$coefficients <- m_fwelmt$beta2[, mt_max_iter + 1]

mod_scores <- Score(
  list(
    cs_glmnet = m_glmet,
    cs_fwelnet = m_fwelnet
  ),

  formula = Hist(time, status) ~ 1, data = dtest,
  times = quantile(dtest$time, probs = seq(.1, .9, .1))
)

summary(mod_scores)
mod_scores$Brier$score

# debugonce(Score.list)
# -> debugonce(riskRegression:::getPerformanceData)

# Direct predictions?
predict(m_glmet, times = 5, newdata = dtest[1:5,])
predict(m_fwelnet, times = 5, newdata = dtest[1:5,])
predict(m_csc, times = 5, newdata = dtest[1:5,])

predictRisk(m_csc, newdata = dtest[1:10,], times = c(1:5))
predictRisk(m_glmet, newdata = dtest[1:10,], times = c(1:5))

library(ggplot2)
mod_scores$Brier$score |>
  ggplot(aes(x = times, y = Brier, color = model, fill = model)) +
  facet_wrap(~model, nrow = 1) +
  #geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .1) +
  geom_path(size = 1.5) +
  scale_color_brewer(palette = "Dark2", aesthetics = c("color", "fill")) +
  theme_minimal() +
  theme(legend.position = "top")
