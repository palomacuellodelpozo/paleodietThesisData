# ==========================================================
# BSIMMS / MixSIAR
# Three Scenarios:
# S1: Uninformed – No Concentration
# S2: Uninformed – With Concentration
# S3: Informed – With Concentration
# ==========================================================

rm(list = ls())

library(MixSIAR)
library(tidyverse)
library(R2jags)
library(rjags)

# ----------------------------------------------------------
# File paths (relative for GitHub)
# ----------------------------------------------------------

mix_file    <- "data/thesis_CN.csv"
source_file <- "data/diet_sourcesConc.csv"
discr_file  <- "data/data_discrimination.csv"

# ----------------------------------------------------------
# Load mixture data
# ----------------------------------------------------------

mix <- load_mix_data(
  filename = mix_file,
  iso_names = c("d13C","d15N"),
  factors = c("Site"),
  fac_random = c(FALSE),
  fac_nested = c(FALSE),
  cont_effects = NULL
)

discr <- load_discr_data(discr_file, mix)

# ==========================================================
# SCENARIO 1 — Uninformed, NO concentration
# ==========================================================

source_s1 <- load_source_data(
  filename = source_file,
  source_factors = NULL,
  conc_dep = FALSE,
  data_type = "means",
  mix = mix
)

model_s1 <- "output/Model_1_Uninformed_NoCon.txt"

write_JAGS_model(model_s1, resid_err=TRUE, process_err=FALSE, mix, source_s1)

jags_1_verylong <- run_model(
  run="very long",
  mix, source_s1, discr, model_s1,
  alpha.prior=c(1,1,1,1)
)

save(jags_1_verylong, file="output/Scenario_1_Uninformed_NoCon.RData")

# ==========================================================
# SCENARIO 2 — Uninformed, WITH concentration
# ==========================================================

source_s2 <- load_source_data(
  filename = source_file,
  source_factors=NULL,
  conc_dep=TRUE,
  data_type="means",
  mix
)

model_s2 <- "output/Model_2_Uninformed_WithCon.txt"

write_JAGS_model(model_s2, resid_err=TRUE, process_err=FALSE, mix, source_s2)

jags_2 <- run_model(
  run="very long",
  mix, source_s2, discr, model_s2,
  alpha.prior=c(1,1,1,1)
)

save(jags_2, file="output/Scenario_2_Uninformed_WithCon.RData")

# ==========================================================
# SCENARIO 3 — Informed priors, WITH concentration
# ==========================================================

informed_prior <- c(3,1,0.5,2)

source_s3 <- source_s2
model_s3  <- "output/Model_3_Informed_WithCon.txt"

write_JAGS_model(model_s3, resid_err=TRUE, process_err=FALSE, mix, source_s3)

jags_3 <- run_model(
  run="very long",
  mix, source_s3, discr, model_s3,
  alpha.prior=informed_prior
)

save(jags_3, file="output/Scenario_3_Informed_WithCon.RData")

# ==========================================================
# Extract comparison across scenarios
# ==========================================================

extract_summary <- function(jags_model, scenario_name){
  
  post_draws <- jags_model$BUGSoutput$sims.list$p.fac1
  sites <- mix$FAC[[1]]$labels
  
  map_dfr(seq_along(sites), function(i){
    data.frame(
      Scenario = scenario_name,
      Site = sites[i],
      C3 = mean(post_draws[,i,1])*100,
      Fish = mean(post_draws[,i,2])*100,
      Shellfish = mean(post_draws[,i,3])*100,
      Terrestrial = mean(post_draws[,i,4])*100
    )
  })
}

all_scenarios <- bind_rows(
  extract_summary(jags_1_verylong, "S1: No Concentration"),
  extract_summary(jags_2, "S2: With Concentration"),
  extract_summary(jags_3, "S3: Informed Priors")
)

write.csv(all_scenarios,
          "output/All_Three_Scenarios_Complete.csv",
          row.names=FALSE)

# Comparison figure

all_long <- all_scenarios %>%
  pivot_longer(C3:Terrestrial,
               names_to="Source",
               values_to="Proportion")

p <- ggplot(all_long,
            aes(Source, Proportion, fill=Scenario)) +
  geom_col(position="dodge") +
  facet_wrap(~Site, ncol=3) +
  theme_bw()

ggsave("output/Final_All_Scenarios_Comparison.pdf",
       p, width=16, height=12)