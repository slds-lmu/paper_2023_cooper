convert_models_csc <- function(train_data, fwel_fit, dummy_csc) {
  # Take csc model and fwelnet fit
  # create glmet/fwelnet "CSC" models
  dummy_glmnet <- dummy_csc

  # Recreate a suitable model matrix with all predictors
  dummy_model_matrix <- model.matrix(prodlim::Hist(time, status) ~ . -1, data = train_data)

  # assign object is created from formula, relevant especially if factors are dummy-coded
  # survival::coxph uses attrassign internally to construct that same object
  dummy_assign <- survival::attrassign(dummy_model_matrix, terms(Surv(time, status) ~ ., data = train_data))
  dummy_glmnet$models[[1]]$assign <- dummy_assign
  dummy_glmnet$models[[2]]$assign <- dummy_assign

  dummy_glmnet$models[[1]]$x <- dummy_model_matrix
  dummy_glmnet$models[[2]]$x <- dummy_model_matrix

  dummy_glmnet$models[[1]]$means <- colMeans(dummy_glmnet$models[[1]]$x)
  dummy_glmnet$models[[2]]$means <- colMeans(dummy_glmnet$models[[2]]$x)

  # terms object can be created from formula fairly easily but has to be consistent with actual data/coefs
  dummy_glmnet$models[[1]]$terms <- terms(Surv(time, status) ~ ., data = train_data)
  dummy_glmnet$models[[2]]$terms <- terms(Surv(time, status) ~ ., data = train_data)

  # copy fudged CSC-franken-glmnet object to hold equally fudged fwelnet object
  dummy_fwelnet <- dummy_glmnet

  # fwelnet_mt model fit includes "final" coefficients per model for fwelnet and
  # underlying original cs-glmnet solution
  dummy_glmnet$models[[1]]$coefficients <- fwel_fit$beta1[, 1]
  dummy_glmnet$models[[2]]$coefficients <- fwel_fit$beta2[, 1]

  dummy_fwelnet$models[[1]]$coefficients <- fwel_fit$beta1[, ncol(fwel_fit$beta1)]
  dummy_fwelnet$models[[2]]$coefficients <- fwel_fit$beta2[, ncol(fwel_fit$beta2)]

  # Return list of glmnet and fwelnet "CSC" models
  list(
    csc_glmnet = dummy_glmnet,
    csc_fwelnet = dummy_fwelnet
  )
}
