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

traits_all <- Lys_data%>%
  filter(!year == 2025) %>%
  select(SRL_m_g, RD_mm, RTD_g_cm3, hyphae)

traits_2023 <- get_traits_by_year(Lys_data, 2023)
traits_2024 <- get_traits_by_year(Lys_data, 2024)
traits_2025 <- Lys_data %>% filter(year == 2025) %>%
  select(SRL_m_g, RD_mm, RTD_g_cm3, leafN_mgg, rootN_mgg)

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


make_corr_plot2025 <- function(df, title = NULL) {
  df_renamed <- df %>%
    rename(
      "Specific root length" = SRL_m_g,
      "Root diameter" = RD_mm,
      "Root tissue density" = RTD_g_cm3,
      "Leaf N" =leafN_mgg,
      "Root N" =rootN_mgg
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

# ---- Create plots ----
p_all <- make_corr_plot(traits_2023, title = "Overall Root Trait Correlation")
p_2023 <- make_corr_plot(traits_2023, title = "First Sampling Year")
p_2024 <- make_corr_plot(traits_2024, title = "Second Sampling Year")
# p_2025 <- make_corr_plot2025(traits_2025, title = "2025")

# Combine plots 

combined_plot <- 
  p_all + p_2023 + p_2024  +
  plot_layout(nrow = 1, ncol = 3, guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold', size=15),
        plot.tag.position  = c(0.3, 1.07),
        plot.title = element_text(hjust = 0.5, face = 'bold', size=11))

print(combined_plot)

ggsave("figures/correlation_plot.png", combined_plot, width = 12, height = 12, dpi = 150)
