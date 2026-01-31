# ==============================================================================
# Program: create_adsl.R
# Purpose: Create ADaM ADSL (Subject Level Analysis Dataset) from SDTM domains
# Input:   pharmaversesdtm::dm, vs, ex, ds, ae
# Output:  adsl.rds, adsl.csv
# Author:  Chunyi Wang
# Date:    2025-01-30
# 
# Note: Set working directory to question_2_adam/ before running this script
# ==============================================================================

# Load required packages
library(admiral)
library(pharmaversesdtm)
library(dplyr)
library(lubridate)
library(stringr)

# Read source data
dm <- pharmaversesdtm::dm
vs <- pharmaversesdtm::vs
ex <- pharmaversesdtm::ex
ds <- pharmaversesdtm::ds
ae <- pharmaversesdtm::ae

# Create ADSL base dataset
adsl <- dm %>%
  select(-DOMAIN)

# Derive TRTSDTM and TRTSTMF
# Process EX data to derive datetime variables
ex_ext <- ex %>%
  derive_vars_dtm(
    dtc = EXSTDTC,
    new_vars_prefix = "EXST",
    time_imputation = "first"
  ) %>%
  derive_vars_dtm(
    dtc = EXENDTC,
    new_vars_prefix = "EXEN",
    time_imputation = "last"
  )

# Merge first exposure datetime to ADSL
adsl <- adsl %>%
  derive_vars_merged(
    dataset_add = ex_ext,
    filter_add = (EXDOSE > 0 | (EXDOSE == 0 & str_detect(EXTRT, "PLACEBO"))) & !is.na(EXSTDTM),
    new_vars = exprs(TRTSDTM = EXSTDTM, TRTSTMF = EXSTTMF),
    order = exprs(EXSTDTM, EXSEQ),
    mode = "first",
    by_vars = exprs(STUDYID, USUBJID)
  )

# Merge last exposure datetime to ADSL
adsl <- adsl %>%
  derive_vars_merged(
    dataset_add = ex_ext,
    filter_add = (EXDOSE > 0 | (EXDOSE == 0 & str_detect(EXTRT, "PLACEBO"))) & !is.na(EXENDTM),
    new_vars = exprs(TRTEDTM = EXENDTM),
    order = exprs(EXENDTM, EXSEQ),
    mode = "last",
    by_vars = exprs(STUDYID, USUBJID)
  )

# Convert datetime to date
adsl <- adsl %>%
  derive_vars_dtm_to_dt(source_vars = exprs(TRTSDTM, TRTEDTM))

# Derive AGEGR9 and AGEGR9N
adsl <- adsl %>%
  mutate(
    AGEGR9 = case_when(
      AGE < 18 ~ "<18",
      AGE >= 18 & AGE <= 50 ~ "18 - 50",
      AGE > 50 ~ ">50",
      TRUE ~ NA_character_
    ),
    AGEGR9N = case_when(
      AGE < 18 ~ 1,
      AGE >= 18 & AGE <= 50 ~ 2,
      AGE > 50 ~ 3,
      TRUE ~ NA_real_
    )
  )

# Derive ITTFL
adsl <- adsl %>%
  mutate(
    ITTFL = ifelse(!is.na(ARM), "Y", "N")
  )

# Derive LSTAVLDT
adsl <- adsl %>%
  derive_vars_extreme_event(
    by_vars = exprs(STUDYID, USUBJID),
    events = list(
      # (1) VS: valid test result and datepart not missing
      event(
        dataset_name = "vs",
        order = exprs(VSDTC, VSSEQ),
        condition = !is.na(VSDTC) & !(is.na(VSSTRESN) & is.na(VSSTRESC)),
        set_values_to = exprs(
          LSTAVLDT = convert_dtc_to_dt(VSDTC, highest_imputation = "M")
        )
      ),
      # (2) AE: onset date not missing
      event(
        dataset_name = "ae",
        order = exprs(AESTDTC, AESEQ),
        condition = !is.na(AESTDTC),
        set_values_to = exprs(
          LSTAVLDT = convert_dtc_to_dt(AESTDTC, highest_imputation = "M")
        )
      ),
      # (3) DS: disposition date not missing
      event(
        dataset_name = "ds",
        order = exprs(DSSTDTC, DSSEQ),
        condition = !is.na(DSSTDTC),
        set_values_to = exprs(
          LSTAVLDT = convert_dtc_to_dt(DSSTDTC, highest_imputation = "M")
        )
      ),
      # (4) ADSL: treatment end date
      event(
        dataset_name = "adsl",
        condition = !is.na(TRTEDT),
        set_values_to = exprs(LSTAVLDT = TRTEDT)
      )
    ),
    source_datasets = list(ae = ae, vs = vs, ds = ds, adsl = adsl),
    tmp_event_nr_var = event_nr,
    order = exprs(LSTAVLDT, event_nr),
    mode = "last",
    new_vars = exprs(LSTAVLDT)
  )

# Select and order final variables
adsl_final <- adsl %>%
  select(
    STUDYID, USUBJID, SUBJID, SITEID,
    AGE, AGEU, AGEGR9, AGEGR9N, SEX, RACE, ETHNIC,
    ARM, ARMCD, ACTARM, ACTARMCD, COUNTRY,
    TRTSDTM, TRTSTMF, TRTSDT, TRTEDTM, TRTEDT,
    ITTFL, LSTAVLDT
  )

# Save output
saveRDS(adsl_final, "adsl.rds")
write.csv(adsl_final, "adsl.csv", row.names = FALSE, na = "")

# Print summary for verification
cat("\n=== ADSL Summary ===\n")
cat("Dimensions:", nrow(adsl_final), "obs x", ncol(adsl_final), "vars\n")
cat("Variables:", paste(names(adsl_final), collapse = ", "), "\n")

cat("\nAGEGR9 distribution:\n")
print(table(adsl_final$AGEGR9, useNA = "ifany"))

cat("\nITTFL distribution:\n")
print(table(adsl_final$ITTFL, useNA = "ifany"))

cat("\nMissing values:\n")
print(colSums(is.na(adsl_final)))