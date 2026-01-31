# ==============================================================================
# Program: 02_create_visualizations.R
# Purpose: Create AE visualizations - severity distribution and top 10 AEs
# Input:   pharmaverseadam::adae, pharmaverseadam::adsl
# Output:  plot1_severity.png, plot2_top10_ae.png
# Author:  Chunyi Wang
# Date:    2025-01-30
# 
# Note: Set working directory to question_3_tlg/ before running this script
# ==============================================================================

# Load required packages
library(ggplot2)
library(dplyr)
library(pharmaverseadam)

# Read source data
adsl <- pharmaverseadam::adsl
adae <- pharmaverseadam::adae

# Filter TEAEs
adae_teae <- adae %>%
  filter(TRTEMFL == "Y")

# ==============================================================================
# Plot 1: AE Severity Distribution by Treatment
# ==============================================================================
p1 <- ggplot(adae_teae, aes(x = ACTARM, fill = AESEV)) +
  geom_bar() +
  labs(
    title = "AE severity distribution by treatment",
    x = "Treatment Arm",
    y = "Count of AEs",
    fill = "Severity/Intensity"
  ) +
  theme_gray() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

# Save Plot 1
ggsave("plot1_severity.png", p1, width = 8, height = 6, dpi = 300)

# ==============================================================================
# Plot 2: Top 10 Most Frequent AEs with 95% CI
# ==============================================================================

# Calculate total number of subjects with TEAEs
n_total <- n_distinct(adae_teae$USUBJID)

# Calculate subject counts per AETERM
ae_counts <- adae_teae %>%
  group_by(AETERM) %>%
  summarise(n = n_distinct(USUBJID), .groups = "drop") %>%
  arrange(desc(n)) %>%
  head(10)

# Calculate percentage and 95% Clopper-Pearson CI
ae_counts <- ae_counts %>%
  rowwise() %>%
  mutate(
    pct = n / n_total * 100,
    ci_low = binom.test(n, n_total)$conf.int[1] * 100,
    ci_high = binom.test(n, n_total)$conf.int[2] * 100
  ) %>%
  ungroup()

# Reorder AETERM by frequency for plotting
ae_counts <- ae_counts %>%
  mutate(AETERM = reorder(AETERM, pct))

# Create forest plot
p2 <- ggplot(ae_counts, aes(x = pct, y = AETERM)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.2) +
  labs(
    title = "Top 10 Most Frequent Adverse Events",
    subtitle = paste0("n = ", n_total, " subjects; 95% Clopper-Pearson CIs"),
    x = "Percentage of Patients (%)",
    y = NULL
  ) +
  theme_gray() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

# Save Plot 2
ggsave("plot2_top10_ae.png", p2, width = 8, height = 6, dpi = 300)

# Print summary for verification
cat("\n=== Visualization Summary ===\n")
cat("TEAEs included:", nrow(adae_teae), "\n")
cat("Unique subjects with TEAEs:", n_total, "\n")
cat("\nPlot 1 saved: plot1_severity.png\n")
cat("Plot 2 saved: plot2_top10_ae.png\n")