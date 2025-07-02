# Purpose:
# This script analyzes root trait data across multiple years.
# It extracts root trait subsets for specific years, checks correlations among traits,
# and visualizes the correlation matrices using heatmaps.

# ---- Load Packages ----
library(tidyverse)
library(ggcorrplot)
library(patchwork)

# ---- Read and Prepare Data ----
Lys_data <- read_csv("data/RGG_Lys_Gradient.csv") %>%
  unite("Plot_ID", c("Lys_no", "date"), remove = FALSE)

# Function to filter root traits by year
get_traits_by_year <- function(data, yr) {
  data %>%
    filter(year == yr) %>%
    select(SRL_m_g, RD_mm, RTD_g_cm3, hyphae)
}

traits_2023 <- get_traits_by_year(Lys_data, 2023)
traits_2024 <- get_traits_by_year(Lys_data, 2024)

# Helper function to rename columns, compute correlation, and return the plot
make_corr_plot <- function(df, title = NULL) {
  df_renamed <- df %>%
    rename(
      "Specific root length" = SRL_m_g,
      "Root diameter" = RD_mm,
      "Root tissue density" = RTD_g_cm3,
      "Roots colonized with AMF" = hyphae
    )
  corr_mat <- round(cor(df_renamed, use = "pairwise.complete.obs", method = "pearson"), 2)
  ggcorrplot(
    corr_mat,
    hc.order = FALSE, type = "lower",
    lab = TRUE, lab_size = 4,
    colors = c("red", "white", "blue"),
    legend.title = "Pearson r"
  ) +
    ggtitle(title)
}

# ---- Create the three plots ----
p_all <- make_corr_plot(
  Lys_data %>% select(SRL_m_g, RD_mm, RTD_g_cm3, hyphae),
  title = "Root Traits Correlation (2023 and 2024)"
)
p_2023 <- make_corr_plot(traits_2023, title = "Root Traits Correlation (2023)")
p_2024 <- make_corr_plot(traits_2024, title = "Root Traits Correlation (2024)")

# Combine plots 

combined_plot <- p_all + p_2023 + p_2024 +
  plot_layout(
  design = "
  AA
  BC", 
  guides = "collect") &
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold', size=20))

print(combined_plot)

ggsave("figures/correlation_plot.png", combined_plot, width = 12, height = 12, dpi = 150)
