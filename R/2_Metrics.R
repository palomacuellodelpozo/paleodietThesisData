# ==========================================================
# Isospace Metrics and Enhanced Biplot
# ==========================================================
rm(list = ls())

library(MixSIAR)
library(tidyverse)
library(geometry)
library(sp)

setwd("~/Documents/paleodietThesisData")

if (!dir.exists("output")) dir.create("output")

mix <- load_mix_data(
  filename     = "data/thesis_CN.csv",
  iso_names    = c("d13C", "d15N"),
  factors      = c("Site"),
  fac_random   = c(FALSE),
  fac_nested   = c(FALSE),
  cont_effects = NULL
)

source <- load_source_data(
  filename       = "data/diet_sourcesConc.csv",
  source_factors = NULL,
  conc_dep       = TRUE,
  data_type      = "means",
  mix            = mix
)

discr <- load_discr_data("data/data_discrimination.csv", mix)

# ----------------------------------------------------------
# Consumer isotope summary by site
# ----------------------------------------------------------
consumer_summary <- mix$data %>%
  group_by(Site) %>%
  summarise(
    n          = n(),
    mean_d13C  = mean(d13C),
    sd_d13C    = sd(d13C),
    min_d13C   = min(d13C),
    max_d13C   = max(d13C),
    mean_d15N  = mean(d15N),
    sd_d15N    = sd(d15N),
    min_d15N   = min(d15N),
    max_d15N   = max(d15N),
    range_d13C = max(d13C) - min(d13C),
    range_d15N = max(d15N) - min(d15N),
    .groups    = "drop"
  )

write.csv(consumer_summary,
          "output/Consumer_Isotope_Summary_by_Site.csv",
          row.names = FALSE)

overall_consumer <- consumer_summary %>%
  summarise(
    n          = sum(n),
    mean_d13C  = mean(mean_d13C),
    mean_d15N  = mean(mean_d15N),
    range_d13C = max(max_d13C) - min(min_d13C),
    range_d15N = max(max_d15N) - min(min_d15N)
  )

# ----------------------------------------------------------
# Source values + TDF correction
# ----------------------------------------------------------
source_data <- data.frame(
  Source    = source$source_names,
  mean_d13C = source$S_MU[, 1],
  sd_d13C   = source$S_SIG[, 1],
  mean_d15N = source$S_MU[, 2],
  sd_d15N   = source$S_SIG[, 2]
)

TDF_d13C <- 1.0
TDF_d15N <- 4.0

source_corrected <- source_data %>%
  mutate(
    corrected_d13C = mean_d13C + TDF_d13C,
    corrected_d15N = mean_d15N + TDF_d15N
  )

write.csv(source_corrected,
          "output/Source_Values_TDF_Corrected.csv",
          row.names = FALSE)

# ----------------------------------------------------------
# Convex hull metrics
# ----------------------------------------------------------
hull_area        <- calc_area(source = source, mix = mix, discr = discr)
consumer_coords  <- cbind(mix$data$d13C, mix$data$d15N)
consumer_hull    <- convhulln(consumer_coords, options = "FA")

# ----------------------------------------------------------
# Mixing polygon inclusion
# ----------------------------------------------------------
mixing_polygon <- source_corrected %>%
  select(corrected_d13C, corrected_d15N)

consumers_in_polygon <- sum(
  point.in.polygon(
    mix$data$d13C,
    mix$data$d15N,
    mixing_polygon$corrected_d13C,
    mixing_polygon$corrected_d15N
  ) > 0
)

# ----------------------------------------------------------
# Pairwise source distances
# ----------------------------------------------------------
source_distances <- as.matrix(
  dist(source_corrected[, c("corrected_d13C", "corrected_d15N")])
)
rownames(source_distances) <- source$source_names
colnames(source_distances) <- source$source_names

# ----------------------------------------------------------
# Trophic position
# ----------------------------------------------------------
baseline_d15N <- min(source_corrected$corrected_d15N)

consumer_trophic <- mix$data %>%
  mutate(
    trophic_position = 1 + ((d15N - baseline_d15N) / TDF_d15N)
  )

trophic_summary <- consumer_trophic %>%
  group_by(Site) %>%
  summarise(
    mean_trophic_position = mean(trophic_position),
    sd_trophic_position   = sd(trophic_position),
    .groups               = "drop"
  )

write.csv(trophic_summary,
          "output/Trophic_Position_by_Site.csv",
          row.names = FALSE)

# ----------------------------------------------------------
# Summary table
# ----------------------------------------------------------
isospace_summary <- data.frame(
  Metric = c(
    "Number of consumers",
    "Consumer δ13C range",
    "Consumer δ15N range",
    "Source convex hull area",
    "Consumer convex hull area",
    "Consumers in mixing polygon (%)",
    "Mean trophic position"
  ),
  Value = c(
    nrow(mix$data),
    round(overall_consumer$range_d13C, 2),
    round(overall_consumer$range_d15N, 2),
    round(hull_area, 2),
    round(consumer_hull$vol, 2),
    round(consumers_in_polygon / nrow(mix$data) * 100, 1),
    round(mean(consumer_trophic$trophic_position), 2)
  )
)

write.csv(isospace_summary,
          "output/Isospace_Summary_Metrics.csv",
          row.names = FALSE)

# ----------------------------------------------------------
# Enhanced isospace plot
# ----------------------------------------------------------
p_isospace <- ggplot() +
  geom_point(
    data  = mix$data,
    aes(d13C, d15N, color = Site),
    size  = 3,
    alpha = 0.6
  ) +
  geom_point(
    data  = source_corrected,
    aes(corrected_d13C, corrected_d15N),
    size  = 5,
    shape = 17,
    color = "black"
  ) +
  geom_polygon(
    data  = source_corrected,
    aes(corrected_d13C, corrected_d15N),
    alpha = 0.1,
    fill  = "grey60"
  ) +
  geom_text(
    data  = source_corrected,
    aes(corrected_d13C, corrected_d15N, label = Source),
    vjust = -1,
    size  = 3.5
  ) +
  theme_bw() +
  labs(
    title = "Bone Collagen Isotope Biplot",
    x     = "δ13C (‰, VPDB)",
    y     = "δ15N (‰, AIR)"
  )

ggsave("output/Isospace_Plot_Enhanced.pdf",
       p_isospace, width = 10, height = 8)

ggsave("output/Isospace_Plot_Enhanced.png",
       p_isospace, width = 10, height = 8, dpi = 300)