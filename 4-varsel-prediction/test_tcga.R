
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
bb <- as.data.table(read_excel("~/sandbox/NIHMS978596-supplement-1.xlsx", sheet = 3))
cc <- bb[, .(bcr_patient_barcode, type, DSS.time.cr, DSS_cr)]
#cc[, table(DSS_cr, type)]
dd <- cc[type == "BLCA", ]
colnames(dd) <- c("id", "type", "time", "status")

# Merge RNAseq and survival data
dt <- merge(dd, gdt, by = "id")
dt[, type := NULL]
dt[, id := NULL]
#dt[status == 2, status := 0] # Only used to convert to non-cr survival (testing)
dt <- dt[!is.na(time) & !is.na(status), ]
dt <- dt[time > 0, ]

# Train/test split
split <- 2/3
dt[, row_id := .I]
train <- dt[ , .SD[sample(.N, size = round(split * .N), replace = FALSE)], by = status]
test <- dt[setdiff(row_id, train$row_id), ]
train[, row_id := NULL]
test[, row_id := NULL]

# Feature selection
fit <- fwelnet::fwelnet_mt_cox(
  train,
  mt_max_iter = 3,
  z_method = "original",
  alpha = 1,
  t = 100,
  a = 0.5,
  thresh = 1e-7,
  include_mt_beta_history = TRUE
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
scores <- data.table::rbindlist(list(
  fit_csc(train, test, model = "fwelnet", coefs = fw_coefs$fwelnet$cause1, cause = 1),
  fit_csc(train, test, model = "fwelnet", coefs = fw_coefs$fwelnet$cause2, cause = 2),
  fit_csc(train, test, model = "glmnet", coefs = fw_coefs$glmnet$cause1, cause = 1),
  fit_csc(train, test, model = "glmnet", coefs = fw_coefs$glmnet$cause2, cause = 2)
))

list(scores = scores, coefs = fw_coefs)

