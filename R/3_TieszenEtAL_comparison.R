# ==========================================================
# Tieszen Comparison and One-Sample Tests
# ==========================================================
rm(list = ls())

library(tidyverse)

setwd("~/Documents/paleodietThesisData")

if (!dir.exists("output")) dir.create("output")

diet_data <- read.csv("data/dietTSZN.csv")

target_sites <- c("Majagora", "Hoya Fria", "Uchova")

this_means <- diet_data %>%
  filter(source == "This",
         Settlement %in% target_sites) %>%
  group_by(Settlement) %>%
  summarise(
    Mean_d13C_This = mean(d13Ccol, na.rm = TRUE),
    Mean_d15N_This = mean(d15N,    na.rm = TRUE),
    n_samples      = n(),
    .groups        = "drop"
  )

tieszen_values <- diet_data %>%
  filter(source == "Tieszen",
         Settlement %in% target_sites) %>%
  select(
    Settlement,
    d13C_Tieszen = d13Ccol,
    d15N_Tieszen = d15N
  )

comparison_final <- inner_join(
  this_means,
  tieszen_values,
  by = "Settlement"
) %>%
  mutate(
    Diff_C = Mean_d13C_This - d13C_Tieszen,
    Diff_N = Mean_d15N_This - d15N_Tieszen
  )

write.csv(comparison_final,
          "output/Tieszen_Mean_Differences.csv",
          row.names = FALSE)

# ----------------------------------------------------------
# One-sample t-tests
# ----------------------------------------------------------
t_test_results <- diet_data %>%
  filter(source == "This",
         Settlement %in% target_sites) %>%
  group_by(Settlement) %>%
  summarise(
    p_val_C13 = t.test(
      d13Ccol,
      mu = comparison_final$d13C_Tieszen[
        comparison_final$Settlement == first(Settlement)
      ]
    )$p.value,
    p_val_N15 = t.test(
      d15N,
      mu = comparison_final$d15N_Tieszen[
        comparison_final$Settlement == first(Settlement)
      ]
    )$p.value,
    .groups = "drop"
  )

write.csv(t_test_results,
          "output/Tieszen_TTest_Results.csv",
          row.names = FALSE)
