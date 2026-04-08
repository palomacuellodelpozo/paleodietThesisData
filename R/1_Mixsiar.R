rm(list = ls())

library(MixSIAR)
library(tidyverse)
library(R2jags)
library(rjags)

setwd("~/Documents/paleodietThesisData")

if (!dir.exists("output")) dir.create("output")

mix <- load_mix_data(
  filename     = "data/data/thesis_CN.csv",
  iso_names    = c("d13C", "d15N"),
  factors      = c("Site"),
  fac_random   = c(FALSE),
  fac_nested   = c(FALSE),
  cont_effects = NULL
)

discr <- load_discr_data("data/data_discrimination.csv", mix)

source_s1 <- load_source_data(
  filename       = "data/diet_sourcesConc.csv",
  source_factors = NULL,
  conc_dep       = FALSE,
  data_type      = "means",
  mix            = mix
)

source_s2 <- load_source_data(
  filename       = "data/diet_sourcesConc.csv",
  source_factors = NULL,
  conc_dep       = TRUE,
  data_type      = "means",
  mix            = mix
)

source_s3 <- source_s2

load("output/Scenario_1_Uninformed_NoCon.RData")
load("output/Scenario_2_Uninformed_WithCon.RData")
load("output/Scenario_3_Informed_WithCon.RData")

# ----------------------------------------------------------
# Extract full posterior draws per site per source
# ----------------------------------------------------------
extract_posterior <- function(jags_model, scenario_name, mix) {
  post_draws    <- jags_model$BUGSoutput$sims.list$p.fac1
  n_sites       <- dim(post_draws)[2]
  n_sources     <- dim(post_draws)[3]
  site_labels   <- mix$FAC[[1]]$labels
  source_labels <- c("C3 Plants", "Fish", "Shellfish", "Terrestrial Fauna")
  
  map_dfr(1:n_sites, function(i) {
    map_dfr(1:n_sources, function(j) {
      data.frame(
        Scenario   = scenario_name,
        Site       = site_labels[i],
        Source     = source_labels[j],
        Proportion = post_draws[, i, j] * 100
      )
    })
  })
}

# ----------------------------------------------------------
# Combine all three scenarios
# ----------------------------------------------------------
all_posterior <- bind_rows(
  extract_posterior(jags_1_verylong, "S1: Uninformed, No Concentration", mix),
  extract_posterior(jags_2,          "S2: Uninformed, Concentration",    mix),
  extract_posterior(jags_3,          "S3: Informed, Concentration",      mix)
) %>%
  mutate(
    Source = factor(Source,
                    levels = c("C3 Plants", "Fish",
                               "Shellfish", "Terrestrial Fauna")),
    Scenario = factor(Scenario,
                      levels = c("S1: Uninformed, No Concentration",
                                 "S2: Uninformed, Concentration",
                                 "S3: Informed, Concentration")),
    Region = case_when(
      Site %in% c("Hoya Fria", "Bco Santos", "Maspalomas") ~ "Lowland",
      Site %in% c("Pico Yeje", "Majagora", "Uchova")       ~ "Highland",
      Site %in% c("Cruz Animas", "Mt Guerra", "Risco Perro") ~ "Midland",
      TRUE ~ "Other"
    ),
    Region = factor(Region, levels = c("Highland", "Midland", "Lowland")),
    Site = factor(Site,
                  levels = c(
                    "Pico Yeje", "Majagora", "Uchova",
                    "Cruz Animas", "Mt Guerra", "Risco Perro",
                    "Hoya Fria", "Bco Santos", "Maspalomas"
                  ))
  )

# ----------------------------------------------------------
# Figure 4 — Posterior boxplots, all sites, all scenarios
# ----------------------------------------------------------
p2 <- ggplot(all_posterior,
             aes(x = Source, y = Proportion, fill = Scenario)) +
  geom_boxplot(
    outlier.shape = NA,
    width         = 0.6,
    position      = position_dodge(width = 0.7),
    alpha         = 0.85
  ) +
  facet_wrap(~ Site, ncol = 3) +
  scale_fill_manual(values = c(
    "S1: Uninformed, No Concentration" = "grey70",
    "S2: Uninformed, Concentration"    = "steelblue",
    "S3: Informed, Concentration"      = "coral"
  )) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 25),
    labels = function(x) paste0(x, "%")
  ) +
  labs(
    x        = NULL,
    y        = "Diet Proportion (%)",
    fill     = "Model Scenario",
    title    = "Posterior Dietary Proportions by Site",
    subtitle = "Left column: Highland | Middle column: Midland | Right column: Lowland"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 0,
                                    hjust = 0.5,
                                    size  = 12),
    strip.text       = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "grey92"),
    legend.position  = "bottom",
    legend.title     = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.spacing    = unit(0.8, "lines")
  )

print(p2)

ggsave("output/Posterior_Boxplots_AllSites.pdf",
       plot = p2, width = 14, height = 14, dpi = 300)
ggsave("output/Posterior_Boxplots_AllSites.png",
       plot = p2, width = 14, height = 14, dpi = 300)