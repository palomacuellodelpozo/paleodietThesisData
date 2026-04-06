# ==========================================================
# Collagen–Apatite Comparisons
# ==========================================================

rm(list = ls())

library(tidyverse)
library(ggplot2)
library(patchwork)

data_file <- "data/ALLmydata.csv"

ALLmydata <- read.csv(data_file)

# ----------------------------------------------------------
# Filter humans
# ----------------------------------------------------------

human_isotopes_clean <- ALLmydata %>%
  filter(Category=="Human",
         Site!="Hoya Brunco",
         Site!="Tenefe") %>%
  select(Site, SampleID, d13Cap, d13Ccol, d15N, dates)

write.csv(human_isotopes_clean,
          "output/human_isotopes_clean_filtered.csv",
          row.names=FALSE)

# ----------------------------------------------------------
# Spacing calculation
# ----------------------------------------------------------

human_isotopes_clean <- human_isotopes_clean %>%
  mutate(
    spacing = d13Cap - d13Ccol,
    diet_type = case_when(
      spacing < 5.5 ~ "High carbohydrate",
      spacing < 7   ~ "Mixed",
      TRUE ~ "High protein"
    )
  )

spacing_summary <- human_isotopes_clean %>%
  group_by(Site) %>%
  summarise(
    n=n(),
    mean_spacing=mean(spacing,na.rm=TRUE),
    sd_spacing=sd(spacing,na.rm=TRUE)
  )

write.csv(spacing_summary,
          "output/Collagen_Apatite_Spacing_Summary.csv",
          row.names=FALSE)

# ----------------------------------------------------------
# Marine contribution from apatite
# ----------------------------------------------------------

d13C_C3_apatite <- -27
d13C_marine_apatite <- -12

human_isotopes_clean <- human_isotopes_clean %>%
  mutate(
    marine_pct_apatite =
      pmin(pmax(
        (d13Cap - d13C_C3_apatite) /
          (d13C_marine_apatite - d13C_C3_apatite),
        0),1)*100
  )

marine_summary <- human_isotopes_clean %>%
  group_by(Site) %>%
  summarise(
    mean_marine_apatite =
      mean(marine_pct_apatite,na.rm=TRUE)
  )

write.csv(marine_summary,
          "output/Marine_Contribution_from_Apatite.csv",
          row.names=FALSE)

# ----------------------------------------------------------
# Load collagen (Scenario 3)
# ----------------------------------------------------------

load("output/Scenario_3_Informed_WithCon.RData")

attach.jags(jags_3)
sites <- mix$FAC[[1]]$labels

collagen_estimates <- map_dfr(seq_along(sites), function(i){
  data.frame(
    Site=sites[i],
    Fish_collagen=mean(p.fac1[,i,2])*100
  )
})

detach.jags()

comparison <- left_join(collagen_estimates,
                        marine_summary,
                        by="Site")

write.csv(comparison,
          "output/Collagen_vs_Apatite_Comparison.csv",
          row.names=FALSE)

# ----------------------------------------------------------
# Temporal categories
# ----------------------------------------------------------

human_isotopes_clean <- human_isotopes_clean %>%
  mutate(
    Period = case_when(
      dates >= 533 & dates <= 877 ~ "Early",
      dates <= 1221 ~ "Middle",
      dates <= 1565 ~ "Late",
      TRUE ~ NA_character_
    )
  )

human_isotopes_clean$Period <-
  factor(human_isotopes_clean$Period,
         levels=c("Early","Middle","Late"),
         ordered=TRUE)

# Period ANOVA example
summary(aov(spacing ~ Period,
            data=human_isotopes_clean))