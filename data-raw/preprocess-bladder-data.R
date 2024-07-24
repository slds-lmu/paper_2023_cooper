# Recreate bladder cancer data as described by and used in:
# - Dyrskjot et al. (2005.) "A Molecular Signature in Superficial Bladder Carcinoma Predicts Clinical Outcome." Clinical Cancer Research: An Official Journal of the American Association for Cancer Research 11 (11): 4029–36. https://doi.org/10.1158/1078-0432.CCR-04-2095.
# - Binder et al. (2009). "Boosting for High-Dimensional Time-to-Event Data with Competing Risks." Bioinformatics 25 (7): 890–96. https://doi.org/10.1093/bioinformatics/btp088.

# Output files are:
# - data/bladder-binder-geno.rds
# - data/bladder-binder-clinical.rds
# - data/bladder-binder-clinical_geno.rds

dat1url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE5nnn/GSE5479/suppl/GSE5479_Final_processed_data_1.txt"
dat2url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE5nnn/GSE5479/suppl/GSE5479_Final_processed_data_2.txt"

dat1file <- here::here("data-raw", "GSE5479_Final_processed_data_1.txt")
dat2file <- here::here("data-raw", "GSE5479_Final_processed_data_2.txt")

if (!file.exists(dat1file)) download.file(dat1url, destfile = dat1file)
if (!file.exists(dat2file)) download.file(dat2url, destfile = dat2file)

exprdat1 <- read.delim(dat1file)
exprdat2 <- read.delim(dat2file)

exprmat <- t(cbind(as.matrix(exprdat1[, 2:ncol(exprdat1)]), as.matrix(exprdat2[, 2:ncol(exprdat2)])))
colnames(exprmat) <- as.character(exprdat1[, 1])
rownames(exprmat)[substr(rownames(exprmat), 1, 1) == "X"] <- substr(rownames(exprmat), 2, 1000)[substr(rownames(exprmat), 1, 1) == "X"]
rownames(exprmat) <- sub("[.]", "-", rownames(exprmat))
rownames(exprmat)[rownames(exprmat) == "1082-1_DK"] <- "1082-1"

# Supp table 1 (xlsx) from
# https://aacrjournals.org/clincancerres/article/13/12/3545/13137/Gene-Expression-Signatures-Predict-Outcome-in-Non
# downloads 10780432ccr062940-sup-supplemental_file_1.xls
rawclinical <- readxl::read_excel(here::here("data-raw", "10780432ccr062940-sup-supplemental_file_1.xls"), .name_repair = make.names)

rawclinical$Sample.ID <- as.character(rawclinical$Sample.ID)
rawclinical$Sample.ID[rawclinical$Sample.ID == "20421_S (91?)"] <- "20421_S"

#   There are two samples with the same ID
#   We assume that the order is the same as in the file with the microarray measurements
rawclinical$Sample.ID[309] <- "692-1.1"

#   check whether the entries match
which(is.na(match(rownames(exprmat), rawclinical$Sample.ID)))

exprmat <- exprmat[match(rawclinical$Sample.ID, rownames(exprmat)), ]
stopifnot("Order is not ok, clinical mismatch genetic" = all.equal(rawclinical$Sample.ID, rownames(exprmat)))

# SEQ1384 appears twice, like many others:
grep("^SEQ1384", colnames(exprmat), value = TRUE)

length(unique(colnames(exprmat))) == ncol(exprmat)
anyDuplicated(colnames(exprmat))

duplicate_name_ids <- which(duplicated(colnames(exprmat)))
duplicate_names <- colnames(exprmat)[duplicate_name_ids]

# Rename by appending "a" as simple workaround
colnames(exprmat)[duplicate_name_ids] <- paste0(duplicate_names, "a")
stopifnot("Duplicate colnames for genetic data" = anyDuplicated(colnames(exprmat)) == 0)

# As per Binder et al, refusing to touch this _at all_
clinical <- data.frame(
  deathcause = rawclinical$Death.cause..0.alive..1.dead.from.bladder.cancer..2.other.specific.cause..3.unknown.cause.,
  progsurvival = rawclinical$Progression.free.survival,
  progression = rawclinical$Progression..0.no.progression..1.progression.to.T1..2.progression.to.T2..,
  progmodel = as.numeric(!is.na(rawclinical$Samples.used.for.training.progression.classifier) | !is.na(rawclinical$Samples.used.for.validating.progression.classifier)),
  clinicalrisk = rawclinical$clinical.risk..1.high.risk..0.low.risk.,
  followup = rawclinical$Follow.up..months.from.tumor.to.last.visit.to.the.clinic..or.to.cystectomy.,
  age = rawclinical$AGE, female = ifelse(rawclinical$SEX == "F" | rawclinical$SEX == " F ", 1, ifelse(rawclinical$SEX == "M" | rawclinical$SEX == " M ", 0, NA)),
  treatment = ifelse(rawclinical$BCG.MMC.treatment == "", 0, 1),
  cystectomy = ifelse(rawclinical$Cystectomy == "cystectomy", 1, 0),
  T1 = ifelse(substr(as.character(rawclinical$Reevaluated.pathological.disease.stage...no.reevaluation), 1, 3) == "pTa", 0, 1),
  gradehigh = ifelse(substr(as.character(rawclinical$Reevaluated.WHO.grade...no.reevaluation), 1, 4) == "HIGH", 1,
    ifelse(rawclinical$Reevaluated.WHO.grade...no.reevaluation == "" | rawclinical$Reevaluated.WHO.grade...no.reevaluation == "-", NA, 0)
  )
)

clinical$progstatus <- ifelse(clinical$progression > 0, 1, ifelse(clinical$deathcause > 1, 2, 0))

if (FALSE) {
  library(survival)
  summary(coxph(Surv(progsurvival,ifelse(progstatus != 0,1,0)) ~ age + female + T1 + gradehigh + treatment,
                data=clinical[clinical$progmodel == 1 & !is.na(clinical$age) & !is.na(clinical$female) & !is.na(clinical$gradehigh),]))
}

table(clinical$progstatus)
table(clinical[clinical$progmodel == 1, ]$progstatus)

# expected 194; 74; 33, sum 301
final_status_table <- table(clinical[clinical$progmodel == 1 & !is.na(clinical$age) & !is.na(clinical$female), "progstatus"])
final_status_table

# triple sanity check
stopifnot("Status counts don't match" = {
  identical(as.integer(final_status_table), c(194L, 74L, 33L))
})

sapply(clinical, \(x) sum(is.na(x)))

# Postprocess: Subset to only covariates used in paper and to non-missing obs with progmodel == 1
non_missing_idx <- which(clinical$progmodel == 1 & !is.na(clinical$age) & !is.na(clinical$female))

clinical <- clinical[non_missing_idx, ]
exprmat <- exprmat[non_missing_idx, ] # have same order as per previous match()

# status = progstatus
# time = progsurvival
# predictors:
# - age
# - sex --> female
# - stage (pTa vs pT1) --> T1
# - grade (PUNLMP/low versus high) --> gradehigh
clinical <- clinical[, c("progsurvival", "progstatus", "age", "female", "T1", "gradehigh")]

# rename for convenience
names(clinical) <- c("time", "status", "age", "female", "T1", "gradehigh")

stopifnot("NAs left in clinical data" = !any(sapply(clinical, \(x) sum(is.na(x))) > 0))
stopifnot("NAs left in genetic data" = !any(sapply(exprmat, \(x) sum(is.na(x))) > 0))

# Adjust event time == 0 to be very small instead, affects 1 obs only
clinical$time[which(clinical$time == 0)] <- 0.001

# Clinical only
saveRDS(clinical, file = here::here("data", "bladder-binder-clinical.rds"))
# saveRDS(rawclinical, file = here::here("data", "bladder-binder-rawclinical.rds"))

# Clinical + genetic vars
clinical_geno <- cbind(clinical, exprmat)
saveRDS(clinical_geno, file = here::here("data", "bladder-binder-clinical_geno.rds"))

# Genetic only but with time/status
geno <- cbind(time = clinical$time, status = clinical$status, exprmat)
saveRDS(exprmat, file = here::here("data", "bladder-binder-geno.rds"))

# dim(geno) == (dim(exprmat) + c(0, 2))
# dim(clinical_geno) == c(nrow(clinical), ncol(clinical) + ncol(exprmat))
