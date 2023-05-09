load(here::here("data-raw/tumor.Rdata"))
load(here::here("data-raw/geno_data.Rdata"))

# Data prep via bimj2384-sup-0001-suppmat/simulation/case_study.R

tumor$stage <- tumor$Reevaluated.pathological.disease.stage...no.reevaluation
tumor$stage[grep("T1", tumor$Reevaluated.pathological.disease.stage...no.reevaluation, ignore.case = T)] <- "PT1"
tumor$stage[grep("T[2-4]", tumor$Reevaluated.pathological.disease.stage...no.reevaluation, perl = T, ignore.case = T)] <- "PT2 or Higher"
tumor$stage[grep("T[ai]", tumor$Reevaluated.pathological.disease.stage...no.reevaluation, ignore.case = T)] <- "PTa"
# table(tumor$stage, exclude = NULL)
tumor$grade <- tumor$Reevaluated.WHO.grade...no.reevaluation
tumor$grade[grep("HIGH", tumor$Reevaluated.WHO.grade...no.reevaluation, ignore.case = T)] <- "HIGH"
tumor$grade[grep("PUNLMP|LOW", tumor$Reevaluated.WHO.grade...no.reevaluation, ignore.case = T)] <- "PUNLMPorLOW"
tumor$grade[tumor$Reevaluated.WHO.grade...no.reevaluation %in% c("-", "")] <- "PUNLMPorLOW"

tumor$tumor <- tumor$Reevaluated.tumor.configuration...no.reevaluation
tumor$tumor[grep("?", tumor$Reevaluated.tumor.configuration...no.reevaluation)] <- NA
tumor$tumor[grep("pap|pip", tumor$Reevaluated.tumor.configuration...no.reevaluation)] <- "papillary"
tumor$tumor[grep("squa|Squa", tumor$Reevaluated.tumor.configuration...no.reevaluation)] <- "papillary"
tumor$tumor[grep("flat|solid", tumor$Reevaluated.tumor.configuration...no.reevaluation)] <- "solid_or_flat"
tumor$CIS <- tumor$CIS.diagnosis..CIS.cis.in.more.than.one.random.biopsy.in.the.disease.course.before.muscle.invasive.tumor..or.CIS.in.diagnostic.tumor.sectio.for.analyzed.tumor.
sex <- gsub("\\s+", "", tumor$SEX)
# table(sex, exclude=NULL)
tumor$SEX <- sex
tumor$trt <- tumor$Cystectomy
tumor$trt[tumor$BCG.MMC.treatment != "" & tumor$Cystectomy != ""] <- paste0(tumor$BCG.MMC.treatment, tumor$Cystectomy)[tumor$BCG.MMC.treatment != "" & tumor$Cystectomy != ""]
tumor$trt[tumor$BCG.MMC.treatment != "" & tumor$Cystectomy == ""] <- tumor$BCG.MMC.treatment[tumor$BCG.MMC.treatment != "" & tumor$Cystectomy == ""]
# table(tumor$trt, exclude=NULL)
tumor$service <- "No"
tumor$service[tumor$trt != ""] <- "Yes"
# table(tumor$service, exclude = NULL)


rmsp <- which(tumor$SEX == "" | is.na(tumor$AGE) | is.na(tumor$tumor))
# length(rmsp)
clin <- tumor[-rmsp, ]
dim(clin) # 385  40
table(clin$Death.cause..0.alive..1.dead.from.bladder.cancer..2.other.specific.cause..3.unknown.cause.)
clin$cens <- clin$Death.cause..0.alive..1.dead.from.bladder.cancer..2.other.specific.cause..3.unknown.cause.
clin$cens[clin$cens > 1] <- 2

clin$time <- clin$Progression.free.survival

dim(clin)
table(clin$cens)
range(clin$time)

# Avoid time == 0
set.seed(24)
clin$time <- clin$time + stats::runif(nrow(clin), min = 0.001, max = .01)
length(unique(clin$time))
range(clin$time)

# clin$Cystectomy[clin$Cystectomy=='']<-'No'
# clin$Cystectomy[clin$Cystectomy=='No']<-'_No'
# cdata<-clin[,c('time','cens','AGE','SEX','stage','grade','service','tumor','CIS','Country')]

#sum(tumor$Reevaluated.WHO.grade...no.reevaluation %in% c('','-')|is.na(tumor$Reevaluated.WHO.grade...no.reevaluation))

stopifnot(all(rownames(geno1) == tumor$Sample.ID))
dim(geno1)

geno1 <- geno1[-rmsp,]

stopifnot(all(rownames(geno1) == clin$Sample.ID))

dim(geno1)
surv_geno <- data.frame(time = clin$time, status = clin$cens, geno1)

if (!file.exists(here::here("data"))) dir.create(here::here("data"))
saveRDS(surv_geno, here::here("data/bladder_surv_geno.rds"))

# Version with genetic + clinical data
#
surv_clin_geno <- data.frame(
  time = clin$time, status = clin$cens,
  sex = factor(clin$SEX), age = clin$AGE,
  grade = factor(clin$grade), factor(clin$stage),
  geno1
)

# Some sanity checking just in case
if (FALSE) {
  names(clin)

  table(clin$AGE)
  table(clin$SEX)
  table(clin$grade)
  table(clin$stage)
  table(clin$cens)

  clin |>
    dplyr::select(AGE, SEX, grade, stage) |>
    sapply(anyNA)

}

