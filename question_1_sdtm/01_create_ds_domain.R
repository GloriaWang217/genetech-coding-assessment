# ==============================================================================
# Program: 01_create_ds_domain.R
# Purpose: Create SDTM DS (Disposition) domain from raw clinical trial data
# Input:   pharmaverseraw::ds_raw, pharmaverseraw::ec_raw, sdtm_ct.csv
# Output:  ds.rds, ds.csv
# Author:  Chunyi
# Date:    2025-01-30
# 
# Note: Set working directory to question_1_sdtm/ before running this script
# ==============================================================================

# Load required packages -------------------------------------------------------
library(sdtm.oak)
library(pharmaverseraw)
library(dplyr)
library(readr)
library(lubridate)

# Read source data -------------------------------------------------------------
ds_raw <- pharmaverseraw::ds_raw
ec_raw <- pharmaverseraw::ec_raw
study_ct <- read_csv("sdtm_ct.csv", show_col_types = FALSE)

# Generate oak_id_vars and fix data quality issues -----------------------------
ds_raw <- ds_raw %>%
  generate_oak_id_vars(
    pat_var = "PATNUM",
    raw_src = "ds_raw"
  ) %>%
  mutate(
    # Fix spelling to match controlled terminology
    INSTANCE = ifelse(INSTANCE == "Ambul Ecg Removal", "Ambul ECG Removal", INSTANCE),
    # Fill IT.DSTERM and IT.DSDECOD with OTHERSP when missing
    IT.DSTERM = ifelse(is.na(IT.DSTERM), toupper(OTHERSP), IT.DSTERM),
    IT.DSDECOD = ifelse(is.na(IT.DSDECOD), toupper(OTHERSP), IT.DSDECOD)
  )

# Create DM domain for derive_study_day ----------------------------------------
dm <- ec_raw %>%
  mutate(ECSTDT = dmy(IT.ECSTDAT)) %>%
  group_by(PATNUM) %>%
  summarise(RFSTDTC = min(ECSTDT, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    USUBJID = paste0(unique(ds_raw$STUDY), "-", PATNUM),
    RFXSTDTC = as.character(RFSTDTC)
  ) %>%
  select(USUBJID, RFXSTDTC)

# Create base DS dataset -------------------------------------------------------
ds <- ds_raw %>%
  select(oak_id, raw_source, patient_number)

# Map variables using sdtm.oak functions ---------------------------------------
ds <- ds %>%
  # DSTERM: verbatim disposition term
  assign_no_ct(
    raw_dat = ds_raw,
    raw_var = "IT.DSTERM",
    tgt_var = "DSTERM",
    id_vars = oak_id_vars()
  ) %>%
  # DSDECOD: standardized disposition term (C66727 codelist)
  assign_ct(
    raw_dat = ds_raw,
    raw_var = "IT.DSDECOD",
    tgt_var = "DSDECOD",
    ct_spec = study_ct,
    ct_clst = "C66727",
    id_vars = oak_id_vars()
  ) %>%
  # VISIT: standardized visit name
  assign_ct(
    raw_dat = ds_raw,
    raw_var = "INSTANCE",
    tgt_var = "VISIT",
    ct_spec = study_ct,
    ct_clst = "VISIT",
    id_vars = oak_id_vars()
  ) %>%
  # VISITNUM: numeric visit identifier
  assign_ct(
    raw_dat = ds_raw,
    raw_var = "INSTANCE",
    tgt_var = "VISITNUM",
    ct_spec = study_ct,
    ct_clst = "VISITNUM",
    id_vars = oak_id_vars()
  ) %>%
  # DSSTDTC: disposition event start date (ISO8601)
  assign_datetime(
    raw_dat = ds_raw,
    raw_var = "IT.DSSTDAT",
    tgt_var = "DSSTDTC",
    raw_fmt = c("m-d-y"),
    id_vars = oak_id_vars()
  ) %>%
  # DSDTC: disposition date (ISO8601)
  assign_datetime(
    raw_dat = ds_raw,
    raw_var = "DSDTCOL",
    tgt_var = "DSDTC",
    raw_fmt = c("m-d-y"),
    id_vars = oak_id_vars()
  )

# Derive SDTM standard variables -----------------------------------------------
ds <- ds %>%
  mutate(
    STUDYID = unique(ds_raw$STUDY),
    DOMAIN = "DS",
    USUBJID = paste0("01-", patient_number),
    DSTERM = toupper(DSTERM),
    # DSCAT: categorize based on DSDECOD
    DSCAT = case_when(
      DSDECOD %in% c("FINAL LAB VISIT", "FINAL RETRIEVAL VISIT") ~ "OTHER EVENT",
      DSDECOD == "RANDOMIZED" ~ "PROTOCOL MILESTONE",
      TRUE ~ "DISPOSITION EVENT"
    ),
    VISITNUM = as.numeric(VISITNUM),
    VISITNUM = ifelse(is.na(VISITNUM) & grepl("UNSCHEDULED", VISIT),
                      as.numeric(gsub("UNSCHEDULED ", "", VISIT)),
                      VISITNUM)
  ) %>%
  # DSSEQ: sequence number within subject, sorted by date
  derive_seq(
    tgt_var = "DSSEQ",
    rec_vars = c("STUDYID", "USUBJID", "DSSTDTC")
  ) %>%
  # DSSTDY: study day relative to first exposure
  derive_study_day(
    sdtm_in = .,
    dm_domain = dm,
    tgdt = "DSSTDTC",
    refdt = "RFXSTDTC",
    study_day_var = "DSSTDY"
  )

# Select and order final variables ---------------------------------------------
ds_final <- ds %>%
  select(
    STUDYID, DOMAIN, USUBJID, DSSEQ, DSTERM, DSDECOD,
    DSCAT, VISITNUM, VISIT, DSDTC, DSSTDTC, DSSTDY
  ) %>%
  arrange(STUDYID, USUBJID, DSSEQ)

# Save output ------------------------------------------------------------------
saveRDS(ds_final, "ds.rds")
write.csv(ds_final, "ds.csv", row.names = FALSE, na = "")

# Print summary for verification -----------------------------------------------
cat("\n=== DS Domain Summary ===\n")
cat("Dimensions:", nrow(ds_final), "obs x", ncol(ds_final), "vars\n")
cat("Variables:", paste(names(ds_final), collapse = ", "), "\n")
cat("Unique subjects:", n_distinct(ds_final$USUBJID), "\n")
cat("\nDSCAT distribution:\n")
print(table(ds_final$DSCAT))
cat("\nMissing values:\n")
print(colSums(is.na(ds_final)))