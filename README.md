# Genentech Analytical Data Science Programmer Coding Assessment

**Author:** Chunyi Wang  
**Date:** January 2025

---

## Repository Structure
```
├── README.md
├── question_1_sdtm/
│   ├── 01_create_ds_domain.R
│   ├── sdtm_ct.csv
│   ├── ds.rds, ds.csv
│   └── q1_log.txt
├── question_2_adam/
│   ├── 02_create_adsl.R
│   ├── adsl.rds, adsl.csv
│   └── q2_log.txt
├── question_3_tlg/
│   ├── 01_create_ae_summary_table.R
│   ├── 02_create_visualizations.R
│   ├── ae_summary_table.html
│   ├── plot1_severity.png, plot2_top10_ae.png
│   └── q3_01_log.txt, q3_01_log.txt
└── question_4_python/
    ├── genai_assistant.ipynb
    └── adae.csv
```

---

## Question 1: SDTM DS Domain Creation

Create SDTM Disposition (DS) domain using `{sdtm.oak}` from `pharmaverseraw::ds_raw`.

**Output:** 850 obs × 12 vars (STUDYID, DOMAIN, USUBJID, DSSEQ, DSTERM, DSDECOD, DSCAT, VISITNUM, VISIT, DSDTC, DSSTDTC, DSSTDY)

---

## Question 2: ADaM ADSL Dataset Creation

Create ADSL dataset using `{admiral}` with derived variables: AGEGR9, AGEGR9N, TRTSDTM, TRTSTMF, ITTFL, LSTAVLDT.

**Output:** 306 obs × 23 vars

---

## Question 3: TLG - Adverse Events Reporting

- **Task 1:** TEAE summary table using `{gtsummary}` → `ae_summary_table.html`
- **Task 2:** Two visualizations using `{ggplot2}`:
  - `plot1_severity.png` - AE severity distribution by treatment
  - `plot2_top10_ae.png` - Top 10 AEs with 95% CI

---

## Question 4: GenAI Clinical Data Assistant (Bonus)

Python-based AI assistant using OpenAI API to translate natural language questions into structured Pandas queries.

**Features:**
- Parses natural language into structured JSON (`target_column`, `filter_value`)
- Returns count of unique subjects (USUBJID) and list of matching IDs
- Interactive Gradio web interface with visualizations

**Example:**
```
Question: "Find patients with cardiac disorders"
Parsed: {"target_column": "AESOC", "filter_value": "CARDIAC DISORDERS"}
Result: 44 unique subjects
```

---

## How to Run

**R (Q1-Q3):** Set working directory to question folder, then `source()` the R script.


**Python (Q4):** Open `genai_assistant.ipynb` in Jupyter/Colab, add your OpenAI API key, run all cells.

## Video Explanation

Here is a brief screen-share walkthrough of my repository. https://www.loom.com/share/cbbf66753e794592adb3c276cd70668d
