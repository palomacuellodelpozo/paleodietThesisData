# ==========================================================
# Collagen-Apatite Comparisons
# ==========================================================
rm(list = ls())

library(MixSIAR)
library(tidyverse)
library(patchwork)

setwd("~/Documents/paleodietThesisData")

if (!dir.exists("output")) dir.create("output")

ALLmydata <- read.csv("data/ALLmydata.csv")

human_isotopes_clean <- ALLmydata %>%
  filter(
    Category == "Human",
    Site != "Hoya Brunco",
    Site != "Tenefe"
  ) %>%
  select(Site, SampleID, d13Cap, d13Ccol, d15N, dates)

write.csv(human_isotopes_clean,
          "output/human_isotopes_clean_filtered.csv",
          row.names = FALSE)


human_isotopes_clean <- human_isotopes_clean %>%
  mutate(
    spacing = d13Cap - d13Ccol,
    diet_type = case_when(
      spacing < 5.5 ~ "High carbohydrate",
      spacing < 7   ~ "Mixed",
      TRUE          ~ "High protein"
    )
  )

spacing_summary <- human_isotopes_clean %>%
  group_by(Site) %>%
  summarise(
    n            = n(),
    mean_spacing = mean(spacing, na.rm = TRUE),
    sd_spacing   = sd(spacing,   na.rm = TRUE),
    .groups      = "drop"
  )

write.csv(spacing_summary,
          "output/Collagen_Apatite_Spacing_Summary.csv",
          row.names = FALSE)

d13C_C3_apatite     <- -27
d13C_marine_apatite <- -12

human_isotopes_clean <- human_isotopes_clean %>%
  mutate(
    marine_pct_apatite = pmin(pmax(
      (d13Cap - d13C_C3_apatite) /
        (d13C_marine_apatite - d13C_C3_apatite),
      0), 1) * 100
  )

marine_summary <- human_isotopes_clean %>%
  group_by(Site) %>%
  summarise(
    mean_marine_apatite = mean(marine_pct_apatite, na.rm = TRUE),
    .groups             = "drop"
  )

write.csv(marine_summary,
          "output/Marine_Contribution_from_Apatite.csv",
          row.names = FALSE)

load("output/Scenario_3_Informed_WithCon.RData")

mix <- load_mix_data(
  filename     = "data/thesis_CN.csv",
  iso_names    = c("d13C", "d15N"),
  factors      = c("Site"),
  fac_random   = c(FALSE),
  fac_nested   = c(FALSE),
  cont_effects = NULL
)

attach.jags(jags_3)

sites <- mix$FAC[[1]]$labels

collagen_estimates <- map_dfr(seq_along(sites), function(i) {
  data.frame(
    Site          = sites[i],
    Fish_collagen = mean(p.fac1[, i, 2]) * 100
  )
})

detach.jags()

# ----------------------------------------------------------
# Compare collagen vs apatite
# ----------------------------------------------------------
comparison <- left_join(collagen_estimates,
                        marine_summary,
                        by = "Site")

write.csv(comparison,
          "output/Collagen_vs_Apatite_Comparison.csv",
          row.names = FALSE)

# ----------------------------------------------------------
# Temporal categories — only individuals WITH dates
# ----------------------------------------------------------
human_dated <- human_isotopes_clean %>%
  filter(!is.na(dates)) %>%
  mutate(
    Period = case_when(
      dates >= 533 & dates <= 877 ~ "Early",
      dates <= 1221               ~ "Middle",
      dates <= 1565               ~ "Late",
      TRUE                        ~ NA_character_
    ),
    Period = factor(Period,
                    levels  = c("Early", "Middle", "Late"),
                    ordered = TRUE)
  )

cat("Period distribution among dated individuals:\n")
print(table(human_dated$Period, useNA = "always"))

# ----------------------------------------------------------
# Period ANOVA (only if at least 2 periods present)
# ----------------------------------------------------------
n_periods <- length(unique(na.omit(human_dated$Period)))

if (n_periods >= 2) {
  anova_result <- aov(spacing ~ Period, data = human_dated)
  summary(anova_result)
  
  anova_summary <- data.frame(
    Effect  = "Period",
    F_value = summary(anova_result)[[1]]$`F value`[1],
    p_value = summary(anova_result)[[1]]$`Pr(>F)`[1]
  )
}
  write.csv(anova_summary,
            "output/Spacing_ANOVA_by_Period.csv",
            row.names = FALSE)
  