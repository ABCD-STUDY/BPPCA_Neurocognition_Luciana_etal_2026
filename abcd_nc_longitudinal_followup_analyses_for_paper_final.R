#######################################
#######################################
## ABCD Bayesian Longitudinal PCA    ##
## Wes Thompson (wes.stat@gmail.com) ##
## Feb 16, 2026                      ##
#######################################
#######################################

rm(list = ls()); gc()

############################
## 1) Libraries (minimal) ##
############################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(mgcv)   # gamm()
  library(nlme)   # lme()
})

############################
## 2) Helpers / utilities ##
############################

## Simple spaghetti plot (kept for identical output logic)
spaghetti <- function(y, times, ids,
                      ylab = "y", xlab = "times",
                      xlim = c(min(times, na.rm = TRUE), max(times, na.rm = TRUE)),
                      ylim = c(min(y, na.rm = TRUE), max(y, na.rm = TRUE)),
                      cols = rep(1, length(unique(ids))), lwd = 1) {
  plot(times, y, type = "n", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab)
  uid <- unique(ids)
  for (i in seq_along(uid)) {
    ii <- ids == uid[i]
    lines(times[ii], y[ii], col = cols[i], lwd = lwd)
  }
}

## Convert "special missing" codes to NA across numeric columns
na_from_codes <- function(df, codes = c(777, 999)) {
  df %>%
    mutate(across(where(is.numeric), ~ replace(., . %in% codes, NA_real_)))
}

## Standardize using externally supplied mean/sd (for applying fixed PCA loadings)
z_with <- function(x, mu, sig) (x - mu) / sig

## Make a contiguous integer ID per group (faster + clearer than manual loops)
dense_int <- function(x) as.integer(factor(x))

####################################
## 3) Data locations + constants  ##
####################################

## These are the original paths. Update if needed.
DATA_DIR_ABCD <- "../../../../../../../../aaa_mac/research/brains/abcd/data/RDS_files/abcd_5.1/rds_file_v5.1"

RDS_IMAGING <- file.path(DATA_DIR_ABCD, "all_imaging_data_5.1.rds")
RDS_NONIMG  <- file.path(DATA_DIR_ABCD, "nonimaging_data_excluding_nt_5.1.rds")

KEEP_EVENTS <- c("baseline_year_1_arm_1",
                 "2_year_follow_up_y_arm_1",
                 "4_year_follow_up_y_arm_1")

## Neuropsych measures used in PCA (5 NIH Toolbox variables)
NP_VARS <- c("nihtbx_picvocab_uncorrected",
             "nihtbx_flanker_uncorrected",
             "nihtbx_pattern_uncorrected",
             "nihtbx_picture_uncorrected",
             "nihtbx_reading_uncorrected")

## Outcome variables (CBCL) - same set as original (ADHD excluded later)
CBCL_VARS <- c("cbcl_scr_syn_internal_r",
               "cbcl_scr_syn_external_r",
               "cbcl_scr_dsm5_adhd_r",
               "cbcl_scr_07_stress_r",
               "cbcl_scr_syn_totprob_r")

## Demographics: new = old
DEMOG_MAP <- c(
  age     = "interview_age",
  sex     = "demo_sex_v2.x",
  race    = "race_ethnicity",
  educ    = "demo_prnt_ed_v2",
  married = "demo_prnt_marital_v2",
  income  = "demo_comb_income_v2"
)

## Nesting: new = old
NEST_MAP <- c(
  id     = "src_subject_id",
  site   = "site_id_l",
  fam    = "rel_family_id.x",
  school = "school_id"
)
  
## fMRI ROI/contrast variables (exact same list as original script)
FMRI_VARS <- c(
  "tfmri_nback_all_545","tfmri_nback_all_579","tfmri_nback_all_568","tfmri_nback_all_602",
  "tfmri_nback_all_544","tfmri_nback_all_578","tfmri_nback_all_567","tfmri_nback_all_601",
  "tfmri_nback_all_557","tfmri_nback_all_591","tfmri_nback_all_133","tfmri_nback_all_147"
)
BEHAV_VARS <- c("tfmri_nb_all_beh_c0b_rate", "tfmri_nb_all_beh_c2b_rate")


####################################
## 4) Load and filter ABCD data   ##
####################################

## Imaging
abcd_imaging <- readRDS(RDS_IMAGING) %>%
  ## Match original filtering:
  ##  iqc_nback_ok_ser.x > 0 AND tfmri_nback_beh_performflag == 1 and not NA
  filter(iqc_nback_ok_ser.x > 0,
         tfmri_nback_beh_performflag == 1,
         !is.na(tfmri_nback_beh_performflag)) %>%
  filter(eventname %in% KEEP_EVENTS) %>%
  ## Keep identifiers + requested fMRI and behavior vars
  select(
    id = src_subject_id, eventname,
    # keep a small stable set of early columns (original kept 1:6, then vars)
    any_of(names(.)[1:6]),
    any_of(FMRI_VARS),
    any_of(BEHAV_VARS)
  )

## Non-imaging
abcd <- readRDS(RDS_NONIMG) # %>%
#  dplyr::rename(
#    id = dplyr::any_of(c("src_subject_id", "src_subject_id.x")),
#    site = dplyr::any_of(c("site_id_l")),
#    fam = dplyr::any_of(c("rel_family_id.x")),
#	school = dplyr::any_of(c("school_id"))
#	sex = dplyr::any_of(c("demo_sex_v2.x"))
#  )
  
#########################################
## 5) Select + process analysis dataset ##
#########################################

data <- abcd %>%
  filter(eventname %in% KEEP_EVENTS) %>%
  select(
    eventname,
    any_of(unname(NEST_MAP)),     # select old names
    any_of(unname(DEMOG_MAP)),    # select old names
    any_of(NP_VARS),
    any_of(CBCL_VARS)
  ) %>%
  rename(!!!NEST_MAP, !!!DEMOG_MAP) %>%
  ## Match original encoding: sex == 1 -> 1 (male), else 0
  mutate(sex = as.numeric(sex == 1)) %>%
  ## Join imaging on (id, eventname) with left join (all.x=TRUE in original)
  left_join(
    abcd_imaging %>% select(id, eventname, any_of(FMRI_VARS), any_of(BEHAV_VARS)),
    by = c("id", "eventname")
  ) %>%
  ## Recode 777/999 -> NA everywhere EXCEPT family & school
  mutate(
    across(
      where(is.numeric),
      ~ if (cur_column() %in% c("fam", "school")) .x
        else replace(.x, .x %in% c(777, 999), NA_real_)
    )
  )
  
rm(abcd_imaging, abcd); gc()

## 5a) Carry baseline demographics/nesting across visits (replaces explicit per-ID loops)
## Original behavior: for each id, if baseline exists, overwrite sex/income/married/educ/race/fam/school
## at all visits with baseline values.
baseline_info <- data %>%
  filter(eventname == "baseline_year_1_arm_1") %>%
  select(id, sex, income, married, educ, race, fam, school)

data <- data %>%
  left_join(baseline_info, by = "id", suffix = c("", ".bl")) %>%
  mutate(
    sex    = coalesce(sex.bl, sex),
    income = coalesce(income.bl, income),
    married= coalesce(married.bl, married),
    educ   = coalesce(educ.bl, educ),
    race   = coalesce(race.bl, race),
    fam    = coalesce(fam.bl, fam),
    school = coalesce(school.bl, school)
  ) %>%
  select(-ends_with(".bl"))

## 5b) Require complete NIH toolbox measures (original used complete.cases on cols 12:16)
data <- data %>%
  filter(if_all(all_of(NP_VARS), ~ !is.na(.)))

## 5c) Keep only subjects with all 3 visits
data <- data %>%
  group_by(id) %>%
  mutate(visit_num = n()) %>%
  ungroup() %>%
  filter(visit_num == 3)

## 5d) Compute family size (number of unique IDs in the family), keep as in original
data <- data %>%
  group_by(fam) %>%
  mutate(fam_num = n_distinct(id)) %>%
  ungroup()

## 5e) Recode ids/fam/school/site to contiguous integers (same intent as manual loops)
## Also reproduce the site22 -> site07 mapping before integer coding.
data <- data %>%
  mutate(
    site = as.character(site),
    site = if_else(site == "site22", "site07", site),
    site = factor(site),
    id     = dense_int(id),
    fam    = dense_int(fam),
    school = dense_int(school),
    site   = dense_int(site)
  ) %>%
  arrange(id, age)

## data2: families with exactly 2 subjects, and then re-code ids/fams to be contiguous within this subset
data2 <- data %>%
  filter(fam_num == 2) %>%
  mutate(
    id  = dense_int(id),
    fam = dense_int(fam)
  ) %>%
  arrange(id, age)

############################################
## 6) Fixed PCA loadings + compute scores ##
############################################

## Loadings from Kelly (same as original)
pca_loadings <- matrix(
  c(0.778, 0.205, 0.142,
    0.205, 0.931, 0.099,
    0.333, 0.879, 0.154,
    0.165, 0.147, 0.986,
    0.824, 0.256, 0.099),
  nrow = 5, byrow = TRUE,
  dimnames = list(c("PicVocab","Flanker","Pattern","Picture","Reading"),
                  c("PC1","PC2","PC3"))
)

## Compute means/sds from data (same as original script)
np_means <- data %>% summarise(across(all_of(NP_VARS), mean)) %>% unlist(use.names = FALSE)
np_sds   <- data %>% summarise(across(all_of(NP_VARS), sd))   %>% unlist(use.names = FALSE)

print(pca_loadings)
print(np_means)
print(np_sds)

## Standardize and compute PCA scores in one shot
np_to_pcs <- function(df) {
  X <- df %>%
    select(all_of(NP_VARS)) %>%
    as.matrix()

  ## z-score using baseline-fitted (global) means/sds
  for (p in seq_len(ncol(X))) X[, p] <- z_with(X[, p], np_means[p], np_sds[p])

  pcs <- X %*% pca_loadings
  df %>%
    mutate(
      pca1 = as.numeric(pcs[, "PC1"]),
      pca2 = as.numeric(pcs[, "PC2"]),
      pca3 = as.numeric(pcs[, "PC3"])
    )
}

data  <- np_to_pcs(data)
data2 <- np_to_pcs(data2)

## Save workspace for downstream section (kept to preserve original flow)
save.image()

##############################################################
################### 7) PCA Regression Models #################
##############################################################

## Clear and reload (matches original script structure)
rm(list = ls()); gc()
load(".RData")

## Create centered sex, visit factor (1:3 within subject), and male factor
data <- data %>%
  mutate(
    male = factor(sex == 1),
    sex  = sex - mean(sex, na.rm = TRUE)
  ) %>%
  group_by(id) %>%
  mutate(visit = factor(row_number())) %>%
  ungroup()

## Linear mixed effects
pc1.lme <- lme(pca1 ~ age + sex + visit, random = ~ 1 | site/fam/id, data = data)
pc2.lme <- lme(pca2 ~ age + sex + visit, random = ~ 1 | site/fam/id, data = data)
pc3.lme <- lme(pca3 ~ age + sex + visit, random = ~ 1 | site/fam/id, data = data)

print(summary(pc1.lme))
print(summary(pc2.lme))
print(summary(pc3.lme))

## GAMMs (nonlinear age) and LRT vs linear-age GAMM (same as original)
pc1.gamm <- gamm(pca1 ~ s(age) + sex + visit, random = list(site = ~1, fam = ~1, id = ~1), data = data)
pc2.gamm <- gamm(pca2 ~ s(age) + sex + visit, random = list(site = ~1, fam = ~1, id = ~1), data = data)
pc3.gamm <- gamm(pca3 ~ s(age) + sex + visit, random = list(site = ~1, fam = ~1, id = ~1), data = data)

pc1.linear.gamm <- gamm(pca1 ~ age + sex + visit, random = list(site = ~1, fam = ~1, id = ~1), data = data)
pc2.linear.gamm <- gamm(pca2 ~ age + sex + visit, random = list(site = ~1, fam = ~1, id = ~1), data = data)
pc3.linear.gamm <- gamm(pca3 ~ age + sex + visit, random = list(site = ~1, fam = ~1, id = ~1), data = data)

print(summary(pc1.gamm$gam)); plot(pc1.gamm$gam)
print(summary(pc2.gamm$gam)); plot(pc2.gamm$gam)
print(summary(pc3.gamm$gam)); plot(pc3.gamm$gam)

print(summary(pc1.linear.gamm$gam)); print(anova(pc1.linear.gamm$lme, pc1.gamm$lme))
print(summary(pc2.linear.gamm$gam)); print(anova(pc2.linear.gamm$lme, pc2.gamm$lme))
print(summary(pc3.linear.gamm$gam)); print(anova(pc3.linear.gamm$lme, pc3.gamm$lme))

## Variance explained proxy (same computation as original)
print(var(predict(pc1.gamm$lme, level = 4)) / var(data$pca1))
print(var(predict(pc2.gamm$lme, level = 4)) / var(data$pca2))
print(var(predict(pc3.gamm$lme, level = 4)) / var(data$pca3))

############################
## 8) Spaghetti plots PCs  ##
############################

age_grid <- seq(min(data$age), max(data$age), by = 0.01)
newdata  <- data.frame(sex = 0, age = age_grid, visit = factor(1))

jpeg("pca-spaghetti-plots.jpeg")
par(mfrow = c(1, 3))

nkeep <- 3 * 400
spaghetti(y = data$pca1[1:nkeep], times = data$age[1:nkeep], ids = data$id[1:nkeep],
          ylab = "pc1 scores", xlab = "age (months)", ylim = c(-6.5, 6.5))
lines(age_grid, predict(pc1.gamm$gam, newdata = newdata), col = 2, lwd = 5); title(main = "PC1")

spaghetti(y = data$pca2[1:nkeep], times = data$age[1:nkeep], ids = data$id[1:nkeep],
          ylab = "pc2 scores", xlab = "age (months)", ylim = c(-6.5, 6.5))
lines(age_grid, predict(pc2.gamm$gam, newdata = newdata), col = 2, lwd = 5); title(main = "PC2")

spaghetti(y = data$pca3[1:nkeep], times = data$age[1:nkeep], ids = data$id[1:nkeep],
          ylab = "pc3 scores", xlab = "age (months)", ylim = c(-6.5, 6.5))
lines(age_grid, predict(pc3.gamm$gam, newdata = newdata), col = 2, lwd = 5); title(main = "PC3")

dev.off()

##############################################################
################### 9) CBCL analyses #########################
##############################################################

## Match original CBCL indexing logic, and remove ADHD from the set
ind_cbcl <- grep("cbcl", names(data))
ind_cbcl <- ind_cbcl[-3]  # remove ADHD, matching original comment

jpeg("cbcl-spaghetti-plots.jpeg")
par(mfrow = c(2, 4))
hist(data[, ind_cbcl[1]][[1]], prob = TRUE, xlab = "score", main = "Internalizing")
hist(data[, ind_cbcl[2]][[1]], prob = TRUE, xlab = "score", main = "Externalizing")
hist(data[, ind_cbcl[3]][[1]], prob = TRUE, xlab = "score", main = "Stress Reactivity")
hist(data[, ind_cbcl[4]][[1]], prob = TRUE, xlab = "score", main = "Total Problems")

nkeep <- 3 * 250
spaghetti(y = data[, ind_cbcl[1]][[1]][1:nkeep], times = data$age[1:nkeep], ids = data$id[1:nkeep],
          ylab = "cbcl internalizing", xlab = "age (months)",
          ylim = c(min(data[, ind_cbcl[1]], na.rm = TRUE), 40))
title(main = "Internalizing")

spaghetti(y = data[, ind_cbcl[2]][[1]][1:nkeep], times = data$age[1:nkeep], ids = data$id[1:nkeep],
          ylab = "cbcl externalizing", xlab = "age (months)",
          ylim = c(min(data[, ind_cbcl[2]], na.rm = TRUE), 40))
title(main = "Externalizing")

spaghetti(y = data[, ind_cbcl[3]][[1]][1:nkeep], times = data$age[1:nkeep], ids = data$id[1:nkeep],
          ylab = "cbcl stress reactivity", xlab = "age (months)",
          ylim = c(min(data[, ind_cbcl[3]], na.rm = TRUE), 40))
title(main = "Stress Reactivity")

spaghetti(y = data[, ind_cbcl[4]][[1]][1:nkeep], times = data$age[1:nkeep], ids = data$id[1:nkeep],
          ylab = "cbcl total problems", xlab = "age (months)",
          ylim = c(min(data[, ind_cbcl[4]], na.rm = TRUE), 100))
title(main = "Total Problems")
dev.off()

## Run Gamma GAMMs for CBCL outcomes (same models as original)
fit_cbcl_gamma <- function(outcome_var) {
  data %>%
    mutate(y = .data[[outcome_var]] + 1) %>%
    { gamm(y ~ age + sex * (pca1 + pca2 + pca3),
           random = list(site = ~1, fam = ~1, id = ~1),
           data = .,
           family = "Gamma") }
}

cbcl_internalizing.gamm <- fit_cbcl_gamma("cbcl_scr_syn_internal_r")
cbcl_externalizing.gamm <- fit_cbcl_gamma("cbcl_scr_syn_external_r")
cbcl_stress.gamm        <- fit_cbcl_gamma("cbcl_scr_07_stress_r")
cbcl_totprob.gamm       <- fit_cbcl_gamma("cbcl_scr_syn_totprob_r")

print(summary(cbcl_internalizing.gamm$gam))
print(summary(cbcl_externalizing.gamm$gam))
print(summary(cbcl_stress.gamm$gam))
print(summary(cbcl_totprob.gamm$gam))

##############################################################
################### 10) fMRI analyses ########################
##############################################################

## Keep complete cases for all fMRI + behavioral vars, then z-scale those predictors
data_fmri <- data %>%
  filter(if_all(all_of(c(FMRI_VARS, BEHAV_VARS)), ~ !is.na(.))) %>%
  mutate(across(all_of(c(FMRI_VARS, BEHAV_VARS)), ~ as.numeric(scale(.))))

## Helper: run one lme per predictor and extract (estimate, p) for that predictor term
lme_scan_one <- function(outcome, predictor) {
  form <- as.formula(paste0(outcome, " ~ ", predictor, " + sex + age"))
  fit  <- lme(fixed = form, random = ~ 1 | site/fam/id, data = data_fmri)
  tt   <- summary(fit)$tTable
  tibble(
    predictor = predictor,
    estimate  = tt[2, "Value"],
    p_value   = tt[2, "p-value"]
  )
}

## Vectorized scans using purrr (replaces repeated loops + manual tables)
scan_predictors <- function(outcome, predictors) {
  map_dfr(predictors, ~ lme_scan_one(outcome, .x))
}

pca1_tfmri_table <- scan_predictors("pca1", FMRI_VARS) %>%
  mutate(across(c(estimate, p_value), as.numeric))
pca2_tfmri_table <- scan_predictors("pca2", FMRI_VARS) %>%
  mutate(across(c(estimate, p_value), as.numeric))
pca3_tfmri_table <- scan_predictors("pca3", FMRI_VARS) %>%
  mutate(across(c(estimate, p_value), as.numeric))

pca1_behav_table <- scan_predictors("pca1", BEHAV_VARS) %>%
  mutate(across(c(estimate, p_value), as.numeric))
pca2_behav_table <- scan_predictors("pca2", BEHAV_VARS) %>%
  mutate(across(c(estimate, p_value), as.numeric))
pca3_behav_table <- scan_predictors("pca3", BEHAV_VARS) %>%
  mutate(across(c(estimate, p_value), as.numeric))

print(pca1_tfmri_table)
print(pca2_tfmri_table)
print(pca3_tfmri_table)

print(pca1_behav_table)
print(pca2_behav_table)
print(pca3_behav_table)

## If you want to exactly mimic the original "array with rownames":
as_array_table <- function(tbl) {
  m <- as.matrix(tbl %>% select(estimate, p_value))
  colnames(m) <- c("value", "p-value")
  rownames(m) <- tbl$predictor
  m
}

pca1_tfmri_array  <- as_array_table(pca1_tfmri_table)
pca2_tfmri_array  <- as_array_table(pca2_tfmri_table)
pca3_tfmri_array  <- as_array_table(pca3_tfmri_table)
pca1_behav_array  <- as_array_table(pca1_behav_table)
pca2_behav_array  <- as_array_table(pca2_behav_table)
pca3_behav_array  <- as_array_table(pca3_behav_table)

pca1_tfmri_array
pca2_tfmri_array
pca3_tfmri_array
pca1_behav_array
pca2_behav_array
pca3_behav_array


