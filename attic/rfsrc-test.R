library(randomForestSRC)
library(ggplot2)

instance <- sim_surv_binder(n_train = 400, n_test = 200, p = 5000, ce = 0.5, lambda = 0.1, lambda_c = 0.1)


## Analysis 1
## modified Gray's weighted log-rank splitting
## (equivalent to cause=c(1,1) and splitrule="logrankCR")
o1 <- rfsrc(Surv(time, status) ~ ., instance[["data"]], ntree = 1000)

## Analysis 2
## log-rank cause-1 (death) specific splitting and targeted VIMP
o2 <- rfsrc(Surv(time, status) ~ ., instance[["data"]],
            splitrule = "logrank", cause = c(1,0), importance = TRUE)

## Analysis 3
## log-rank cause-2 (transplant) specific splitting and targeted VIMP
o3 <- rfsrc(Surv(time, status) ~ ., instance[["data"]],
            splitrule = "logrank", cause = c(0,1), importance = TRUE)

## extract VIMP from the log-rank forests: event-specific
## extract minimal depth from the Gray log-rank forest: non-event specific
vimpOut <- data.frame(md = max.subtree(o1)$order[, 1],
                      vimp.c1 = 100 * o2$importance[ ,1],
                      vimp.c2 = 100 * o3$importance[ ,2])
saveRDS(vimpOut, "attic/vimpOut-testfit.rds")
# vimpOut <- readRDS("attic/vimpOut-testfit.rds")

print(vimpOut[order(vimpOut$md), ], digits = 2)

any(vimpOut$md <= 0)
any(vimpOut$vimp.c1 <= 0)
any(vimpOut$vimp.c2 <= 0)

sum(vimpOut$vimp.c1 <= 0)
sum(vimpOut$vimp.c2 <= 0)

library(ggplot2)

vimpOut |>
  tidyr::gather() |>
  ggplot(aes(x = value, fill = key)) +
  facet_wrap(vars(key), ncol = 1, scales = "free_x") +
  geom_histogram(position = "dodge", color = "darkgray") +
  geom_vline(xintercept = 0, color = "red") +
  theme_minimal()


vimpOut |>
  tidyr::gather() |>
  dplyr::filter(key == "md") |>
  ggplot(aes(x = value, fill = key)) +
  facet_wrap(vars(key), ncol = 1, scales = "free_x") +
  geom_density() +
  #geom_histogram(position = "dodge", color = "darkgray") +
  #geom_vline(xintercept = 0, color = "red") +
  theme_minimal()

sum(mx$order[, 1] <= mx$threshold)

instance$covar_true_effect |> sapply(length) |> sum()



rfsrc_importance <- function(data, importance = "anti", cause = c(1, 1)) {
  # importances
  # default TRUE => importance = "anti"
  checkmate::assert_subset(importance, c("random", "permute", "anti", "anti.join", "permute.joint", "random.joint"))


  o <- rfsrc(
    Surv(time, status) ~ ., data = data,
    splitrule = "logrank",
    cause = cause,
    importance = importance,
    perf.type = "none" # no crank, for speedup
  )

  vimp <- o$importance[, which(cause == 1), drop = FALSE]
  colnames(vimp) <- paste0(importance, ".", colnames(vimp))
  vimp
}


# Trying out all the importances ----------------------------------------------------------------------------------
# Except for the `logrankCR` Gray one which is weird

rf_c1 <- rfsrc(
  Surv(time, status) ~ ., data = instance$data,
  splitrule = "logrank",
  cause = c(1, 0),
  importance = FALSE,
  perf.type = "none" # no crank, for speedup
)

rf_c2 <- rfsrc(
  Surv(time, status) ~ ., data = instance$data,
  splitrule = "logrank",
  cause = c(0, 1),
  importance = FALSE,
  perf.type = "none" # no crank, for speedup
)

vimps_all <- vector(mode = "list")

for (importance in  c("random", "permute", "anti")) {
   vtemp1 <- vimp(rf_c1, importance = importance, joint = FALSE)[["importance"]][, 1, drop = FALSE]
   vtemp2 <- vimp(rf_c2, importance = importance, joint = FALSE)[["importance"]][, 2, drop = FALSE]
   vtemp <- cbind(vtemp1, vtemp2)
   vtemp <- data.table::as.data.table(vtemp, keep.rownames = "variable")
   vtemp$importance <- importance

   vimps_all[[importance]] <- vtemp
}

vimps_all <- data.table::rbindlist(vimps_all)
saveRDS(vimps_all, "attic/vimps_all.rds")
vimps_all <- readRDS("attic/vimps_all.rds")

vimps_all |>
  data.table::melt(measure.vars = c("event.1", "event.2"), variable.name = "cause", value.name = "vi") |>
  ggplot(aes(x = vi * 100, fill = importance)) +
  facet_wrap(cause~importance, ncol = 3, scales = "free") +
  geom_histogram() +
  geom_vline(xintercept = 0, color = "red") +
  theme_minimal()


vimps_all$event.1[which(vimps_all$event.1 <= abs(min(vimps_all$event.1)))] |> length()
sum(vimps_all$event.1 <= 0)

vimps_all[, event.1r := ifelse(event.1 <= abs(min(event.1)), 0, event.1), by = importance]
vimps_all[, event.2r := ifelse(event.2 <= abs(min(event.2)), 0, event.2), by = importance]

vimps_all |>
  data.table::melt(measure.vars = c("event.1", "event.2"), variable.name = "cause", value.name = "vi") |>
  ggplot(aes(x = vi * 100, fill = importance)) +
  facet_wrap(cause~importance, ncol = 3, scales = "free") +
  geom_histogram() +
  geom_vline(xintercept = 0, color = "red") +
  theme_minimal()

vimps_all |>
  data.table::melt(measure.vars = c("event.1r", "event.2r"), variable.name = "cause", value.name = "vi") |>
  ggplot(aes(x = vi * 100, fill = importance)) +
  facet_wrap(cause~importance, ncol = 3, scales = "free") +
  geom_histogram() +
  geom_vline(xintercept = 0, color = "red") +
  theme_minimal()



# Tuning test -----------------------------------------------------------------------------------------------------

instance <- sim_surv_binder(n_train = 400, p = 5000, ce = 0.5, lambda = 0.1, lambda_c = 0.1)

# ceiling((ncol(instance$data) - 2)/3)
# Winning combs
# nodesize: 55
# mtry: 2083

tictoc::tic()
rftuned <- tune(
  formula = Surv(time, status) ~ .,
  data = instance$data,
  # Default nodesize for survival is 15, probably shouldn't go lower
  nodesizeTry = 15,
  # due to high-dim setting recommended > p /5, possibly p/2
  mtryStart = ceiling((ncol(instance$data) - 2)/3),
  # Fixing this to 500 (default) for now
  ntreeTry = 500,
  doBest = TRUE,
  # Passed to rfsrc
  splitrule = "logrank",
  cause = c(1, 0),
  importance = "random"
)
tictoc::toc()

rftuned$optimal

# Benchmark to check which fitting procedure is faster ------------------------------------------------------------

bench::mark(
  aio = {
    set.seed(2)
    rf_c1 <- rfsrc(
      Surv(time, status) ~ ., data = instance$data,
      splitrule = "logrank",
      cause = c(1, 0),
      importance = "random"
    )
    rf_c1[["importance"]][, 1, drop = FALSE]

  },
  sequential = {
    set.seed(2)
    rf_c1 <- rfsrc(
      Surv(time, status) ~ ., data = instance$data,
      splitrule = "logrank",
      cause = c(1, 0),
      importance = FALSE, perf.type = "none" # no crank, for speedup
    )
    set.seed(2)
    vimp(rf_c1, importance = "random", joint = FALSE)[["importance"]][, 1, drop = FALSE]
  },
  check = FALSE,
  iterations = 3, min_iterations = 3
) -> bm


# Different cause-settings ----------------------------------------------------------------------------------------

library(randomForestSRC)

instance <- sim_surv_binder(n_train = 400, n_test = 200, p = 400, ce = 0.5, lambda = 0.1, lambda_c = 0.1)

instance[["data"]]$status

## Analysis 1
## modified Gray's weighted log-rank splitting
## (equivalent to cause=c(1,1) and splitrule="logrankCR")

set.seed(2)
rf_c1 <- rfsrc(
  Surv(time, status) ~ ., data = instance$data,
  splitrule = "logrank",
  cause = c(1, 0),
  importance = FALSE,
  perf.type = "none"
)

# vimp(rf_c1, importance = "random", joint = FALSE)[["importance"]][, 1, drop = FALSE]

# recode to censor event 2
instance$data$status2 <- instance$data$status
instance$data$status2[instance$data$status2 == 2] <- 0

set.seed(2)
rf_c1cens <- rfsrc(
  Surv(time, status2) ~ ., data = instance$data,
  splitrule = "logrank",
  #cause = c(1, 0),
  importance = FALSE,
  perf.type = "none"
)

plot.competing.risk(rf_c1)
plot.survival(rf_c1cens)

# get.brier.survival(rf_c1, cens.mode = "km")$brier.score # only standard setting
get.brier.survival(rf_c1cens, cens.mode = "km")$brier.score


## cr splitrule
rf_cr <- rfsrc(
  Surv(time, status) ~ ., data = instance$data,
  splitrule = "logrankCR",
  cause = c(1, 1),
  importance = FALSE, perf.type = "none" # no crank, for speedup
)

plot.competing.risk(rf_c1)
plot.competing.risk(rf_cr)

cr_maxsub <- max.subtree(rf_cr)
cr_maxsub$order |> head()
cr_maxsub$order[, 1] # cause 1 only
