#renv::install("bioc::TCGAbiolinks")
#renv::install("bioc::SummarizedExperiment")

library(TCGAbiolinks)
library(SummarizedExperiment)
library(data.table)
library(readxl)
library(fwelnet)
library(riskRegression)
source("4-varsel-prediction/algorithms.R")

# Get RNAseq data
query <- GDCquery(project = "TCGA-BLCA",
                  data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  experimental.strategy = "RNA-Seq",
                  platform = "Illumina HiSeq",
                  file.type = "results",
                  legacy = TRUE)
GDCdownload(query)
rnaseq <- GDCprepare(query)

# Transform to matrix and data table
mat <- assay(rnaseq, "raw_count")
gdt <- as.data.table(t(mat))

# Just use the first 5000 genes for now
gdt <- gdt[, 1:5000]
colnames(gdt) <- gsub("\\|", "_", colnames(gdt))

# Add ID
gdt[, id := substr(colnames(mat), 1, 12)]

# Get survival data
# Paper here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6066282/
#bb <- as.data.table(read_excel("https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6066282/bin/NIHMS978596-supplement-1.xlsx", sheet = 3))

tcga_excel <- here::here("data-raw", "NIHMS978596-supplement-1.xlsx")
if (!file.exists(tcga_excel)) {
  download.file("https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6066282/bin/NIHMS978596-supplement-1.xlsx", tcga_excel)
}
bb <- as.data.table(read_excel(here::here("data-raw", "NIHMS978596-supplement-1.xlsx"), sheet = 3))

# check different types/endpoints
table(bb$type, bb$DSS_cr) |>
  as.data.frame() |>
  arrange(desc(`2`))
  rowSums()

bb |>
  count(type, DSS_cr) |>
  tidyr::pivot_wider(names_from = DSS_cr, values_from = n) |>
  arrange(desc(`2`))

cc <- bb[, .(bcr_patient_barcode, type, DSS.time.cr, DSS_cr)]
#cc[, table(DSS_cr, type)]
dd <- cc[type == "BLCA", ] # check differen types
colnames(dd) <- c("id", "type", "time", "status")

# Merge RNAseq and survival data
dt <- merge(dd, gdt, by = "id")
dt[, type := NULL]
dt[, id := NULL]
#dt[status == 2, status := 0] # Only used to convert to non-cr survival (testing)
dt <- dt[!is.na(time) & !is.na(status), ]
dt <- dt[time > 0, ]

dim(dt)

# De-select vars with zero sd, i.e. constant
dt <- dt[, which(sapply(dt, sd) > 0), with = FALSE]

# Train/test split
split <- 2/3
dt[, row_id := .I]
# FIXME only works with fixed seed b/c otherwise we get an error due
# to too many censored obs in glmnet CV - not sure how to fix
set.seed(25)
train <- dt[ , .SD[sample(.N, size = round(split * .N), replace = FALSE)], by = status]
test <- dt[setdiff(row_id, train$row_id), ]
train[, row_id := NULL]
test[, row_id := NULL]

# get column-wise means and sd of genetic data only (1, 2 are time, status)
gen_means <- apply(train[, -c(1,2)], 2, mean)
gen_sds <- apply(train[, -c(1,2)], 2, sd)
# apply to only genetic data, cbind with time, status. There's probably a more elegant solution, sorry.
train <- cbind(train[, 1:2], scale(train[, -c(1,2)], center = gen_means, scale = gen_sds))
test <- cbind(test[, 1:2], scale(test[, -c(1,2)], center = gen_means, scale = gen_sds))

#no missings
stopifnot(!anyNA(train))
sds <- apply(train, 2, sd)
stopifnot(!any(sds == 0))
table(train$status)
#train <- as.data.table(train)

# lapply(train, \(x) {
#   data.frame(
#     n_distinct = length(unique(x)),
#     min_freq = min(table(x)),
#     max_freq = max(table(x))
#   )
# }) |> data.table::rbindlist(idcol = TRUE) |>
#   dplyr::filter(n_distinct < 50) |>
#   dplyr::arrange(n_distinct, min_freq) |>
#   View()

# debugonce(fwelnet)
# Feature selection
fit <- fwelnet::fwelnet_mt_cox(
  train,
  mt_max_iter = 3,
  z_method = "original",
  alpha = 1,
  t = 100,
  a = 0.5,
  thresh = 1e-7, standardize = FALSE,
  include_mt_beta_history = TRUE, verbose = FALSE
)

p <- ncol(train)

glmnet_beta1 <- fit$beta1[, 1]
glmnet_beta2 <- fit$beta2[, 1]

# Same for fwelnet estimates
fwel_beta1 <- fit$beta1[, ncol(fit$beta1)]
fwel_beta2 <- fit$beta2[, ncol(fit$beta1)]

fw_coefs <- list(
  fwelnet = list(
    cause1 = fwel_beta1[fwel_beta1 != 0],
    cause2 = fwel_beta2[fwel_beta2 != 0]
  ),
  glmnet = list(
    cause1 = glmnet_beta1[glmnet_beta1 != 0],
    cause2 = glmnet_beta2[glmnet_beta2 != 0]
  )
)

# Performance with CSC
# sFIXME: beta2 is all zeros!
scores <- data.table::rbindlist(list(
  fit_csc(train, test, model = "fwelnet", coefs = fw_coefs$fwelnet$cause1, cause = 1),
  #fit_csc(train, test, model = "fwelnet", coefs = fw_coefs$fwelnet$cause2, cause = 2),
  fit_csc(train, test, model = "glmnet", coefs = fw_coefs$glmnet$cause1, cause = 1)
  #fit_csc(train, test, model = "glmnet", coefs = fw_coefs$glmnet$cause2, cause = 2)
))

list(scores = scores, coefs = fw_coefs)

scores |>
  filter(score > 0 & score < 1) |>
  ggplot(aes(x = times, y = score, color = model)) +
  facet_wrap(facets = vars(metric), scales = "free") +
  geom_line()
