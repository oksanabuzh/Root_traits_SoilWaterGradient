# Purpose: Fit ordination (RDA) for root traits and their predictors for 2023 and 2024

# ---- Load packages ----
library(tidyverse)
library(vegan)
library(ggrepel)
library(tibble)

# ---- Read and Prepare Data ----
Lys_data <- read_csv("data/RGG_Lys_Gradient.csv") %>%
  unite("Plot_ID", c("Lys_no", "date"), remove = TRUE)

# ---- Utility function to extract and subset data by year ----
get_traits_by_year <- function(data, yr) {
  data %>%
    filter(year == yr) %>%
    select(Plot_ID, GW_level_cm, month, SRL_m_g, RD_mm, RTD_g_cm3, hyphae)
}

# ---- Prepare predictors and trait matrices for both years ----
prep_predictors <- function(data, yr) {
  get_traits_by_year(data, yr) %>%
    select(Plot_ID, GW_level_cm, month)
}

prep_traits <- function(data, yr) {
  get_traits_by_year(data, yr) %>%
    select(Plot_ID, SRL_m_g, RD_mm, RTD_g_cm3, hyphae) %>%
    column_to_rownames("Plot_ID")
}

predictors_2023 <- prep_predictors(Lys_data, 2023)
predictors_2024 <- prep_predictors(Lys_data, 2024)
traits_2023 <- prep_traits(Lys_data, 2023)
traits_2024 <- prep_traits(Lys_data, 2024)

# ---- Check DCA axis length (to justify linear RDA) ----
decorana(traits_2023)
decorana(traits_2024)

# ---- Run Redundancy Analysis for both years ----
ordin_rda <- function(traits, predictors) {
  rda(traits ~ GW_level_cm + month, data = predictors, scale = TRUE)
}

ordin2023 <- ordin_rda(traits_2023, predictors_2023)
ordin2024 <- ordin_rda(traits_2024, predictors_2024)

# ---- Model evaluation ----
vif.cca(ordin2023)
vif.cca(ordin2024)
summary(eigenvals(ordin2023))
summary(eigenvals(ordin2024))
RsquareAdj(ordin2023)
RsquareAdj(ordin2024)

# ---- Permutation tests ----
set.seed(1)
perm_test <- function(model) {
  perm_stat <- permustats(anova(model))
  list(
    summary = summary(perm_stat),
    densityplot = densityplot(perm_stat),
    anova = anova(model),
    anova_margin = anova(model, by= "margin"),
    anova_axis = anova(model, by= "axis")
  )
}

perm_2023 <- perm_test(ordin2023)
perm_2024 <- perm_test(ordin2024)

# ---- Goodness of fit for species (traits) ----
goodness(ordin2023, display = "species", summarize = FALSE)
goodness(ordin2024, display = "species", summarize = FALSE)

# ---- Plot ordination for both years ----
plot(ordin2023, scaling = "species")
plot(ordin2024, scaling = "species")

# ---- Extract RDA scores for plotting ----
extract_env_scrs <- function(ordin) {
  scores(ordin, display = "bp") %>%
    as_tibble(rownames="predictors") %>%
    filter(predictors=="GW_level_cm") %>%
    mutate(predictors_names="Groundwater level")
}

extract_centroids <- function(ordin) {
  scores(ordin, display="cn", scaling="species") %>%
    as_tibble(rownames = "month") %>%
    mutate(month=stringr::str_sub(month, 6))
}

extract_trait_scrs <- function(ordin) {
  scores(ordin, display = "species", scaling = "species") %>%
    as_tibble(rownames = "traits") %>%
    mutate(traits_names=case_when(
      traits=="SRL_m_g" ~ "Specific root length",
      traits=="RD_mm" ~ "Root diameter",
      traits=="RTD_g_cm3" ~ "Root tissue density",
      traits=="hyphae" ~"Roots AMF"
    ))
}

extract_observ_scrs <- function(ordin, predictors, centroids) {
  scores(ordin, display = "sites", scaling = "species") %>%
    as_tibble(rownames="Plot_ID") %>%
    left_join(predictors, by="Plot_ID") %>%
    left_join(centroids %>% rename(RDA1_centr=RDA1, RDA2_centr=RDA2), by="month")
}

# --- For 2023 ---
env.scrs_2023 <- extract_env_scrs(ordin2023)
centroids_2023 <- extract_centroids(ordin2023)
trait.scrs_2023 <- extract_trait_scrs(ordin2023)
observ.scrs_2023 <- extract_observ_scrs(ordin2023, predictors_2023, centroids_2023)

# --- For 2024 ---
env.scrs_2024 <- extract_env_scrs(ordin2024)
centroids_2024 <- extract_centroids(ordin2024)
trait.scrs_2024 <- extract_trait_scrs(ordin2024)
observ.scrs_2024 <- extract_observ_scrs(ordin2024, predictors_2024, centroids_2024)


# ---- Plotting ordination results: 2023 and 2024 ----

# ---- 2023: Plot all elements ----
ggplot(observ.scrs_2023) +
  geom_hline(yintercept = 0, color = "grey", lty = 1) +
  geom_vline(xintercept = 0, color = "grey", lty = 1) +
  # stat_ellipse(aes(RDA1, RDA2, fill=month), alpha=.1,type='t', linewidth =1, geom="polygon") +
  eom_segment(aes(x = RDA1, y = RDA2, xend = RDA1_centr, yend = RDA2_centr, color = month),
              alpha = .4, linewidth = 0.2, linetype = 5) +
  geom_point(data = centroids_2023, aes(x = RDA1, y = RDA2, color = month), size = 4, alpha = .5, stroke = 2) +
  geom_point(aes(RDA1, RDA2, color = month), size = 3, alpha = .5) +
  geom_text(data = trait.scrs_2023, aes(RDA1, RDA2, label = traits_names),
            color = "black", vjust = c(-1.1, 1.2, -0.6, 0.2), hjust = c(0.8, 0, 0, -0.1)) +
  geom_segment(data = trait.scrs_2023, aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black", linewidth = 0.9) +
  geom_text(data = env.scrs_2023, aes(RDA1, RDA2, label = predictors_names),
            color = "blue", size = 5, vjust = -0.5, hjust = 0.7) +
  geom_segment(data = env.scrs_2023, aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.3, "cm")), color = "blue", linewidth = 1) +
  labs(color = "month", x = "RDA1 (32.7 %)", y = "RDA2 (4.8 %)") +
  theme_bw() +
  scale_color_manual(values = c("#3CB22D", "#FF8000")) +
  scale_fill_manual(values = c("#3CB22D", "#FF8000"))

# ---- 2023: Plot only traits and predictors ----
ggplot() +
  geom_hline(yintercept = 0, color = "grey", lty = 1) +
  geom_vline(xintercept = 0, color = "grey", lty = 1) +
  geom_point(data = centroids_2023, aes(x = RDA1, y = RDA2, color = month), size = 4, alpha = .5, stroke = 2) +
  geom_text(data = centroids_2023, aes(RDA1, RDA2, label = month, color = month), size = 5, vjust = c(-1, -1), hjust = c(0, 1)) +
  geom_text(data = trait.scrs_2023, aes(RDA1, RDA2, label = traits_names),
            color = "black", vjust = c(-1.1, 1.2, -0.6, 0.2), hjust = c(0.5, 0, 0.5, -0.1)) +
  geom_segment(data = trait.scrs_2023, aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black", linewidth = 0.9) +
  geom_text(data = env.scrs_2023, aes(RDA1, RDA2, label = predictors_names),
            color = "blue", size = 5, vjust = -0.5, hjust = 0.7) +
  geom_segment(data = env.scrs_2023, aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.3, "cm")), color = "blue", linewidth = 1) +
  labs(color = "month", x = "RDA1 (32.7 %)", y = "RDA2 (4.8 %)") +
  theme_bw() +
  scale_color_manual(values = c("#3CB22D", "#FF8000")) +
  scale_fill_manual(values = c("#3CB22D", "#FF8000")) +
  ylim(-0.5, 0.9) +
  xlim(-1.2, 1.2)

# ---- 2024: Plot all elements ----
ggplot(observ.scrs_2024) +
  geom_hline(yintercept = 0, color = "grey", lty = 1) +
  geom_vline(xintercept = 0, color = "grey", lty = 1) +
  # stat_ellipse(aes(RDA1, RDA2, fill=month), alpha=.1,type='t', linewidth =1, geom="polygon") +
  geom_segment(aes(x = RDA1, y = RDA2, xend = RDA1_centr, yend = RDA2_centr, color = month),
               alpha = .4, linewidth = 0.2, linetype = 5) +
  geom_point(data = centroids_2024, aes(x = RDA1, y = RDA2, color = month), size = 4, alpha = .5, stroke = 2) +
  geom_point(aes(RDA1, RDA2, color = month), size = 3, alpha = .5) +
  geom_text(data = trait.scrs_2024, aes(RDA1, RDA2, label = traits_names),
            color = "black", vjust = c(-1.1, 1.2, 1.2, 1.2), hjust = c(0.7, 0.7, 0.5, 0.7)) +
  geom_segment(data = trait.scrs_2024, aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black", linewidth = 0.9) +
  geom_text(data = env.scrs_2024, aes(RDA1, RDA2, label = predictors_names),
            color = "blue", size = 5, vjust = -0.5, hjust = 0.7) +
  geom_segment(data = env.scrs_2024, aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.3, "cm")), color = "blue", linewidth = 1) +
  labs(color = "month", x = "RDA1 (32.1 %)", y = "RDA2 (13 %)") +
  theme_bw() +
  scale_color_manual(values = c("#3CB22D", "#FF8000")) +
  scale_fill_manual(values = c("#3CB22D", "#FF8000"))+
  xlim(-1.5, 2.2)

# ---- 2024: Plot only traits and predictors ----
ggplot() +
  geom_hline(yintercept = 0, color = "grey", lty = 1) +
  geom_vline(xintercept = 0, color = "grey", lty = 1) +
  geom_point(data = centroids_2024, aes(x = RDA1, y = RDA2, color = month), size = 4, alpha = .5, stroke = 2) +
  geom_text(data = centroids_2024, aes(RDA1, RDA2, label = month, color = month), ,show_guide = F, size = 5, vjust = c(-1, -1), hjust = c(0, 1)) +
  geom_text(data = trait.scrs_2024, aes(RDA1, RDA2, label = traits_names),
            color = "black", vjust = c(-1.1, 1.2, 1.2, 1.2), hjust = c(0.7, 0.7, 0.5, 0.7)) +
  geom_segment(data = trait.scrs_2024, aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black", linewidth = 0.9) +
  geom_text(data = env.scrs_2024, aes(RDA1, RDA2, label = predictors_names),
            color = "blue", size = 5, vjust = -0.5, hjust = 0.7) +
  geom_segment(data = env.scrs_2024, aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.3, "cm")), color = "blue", linewidth = 1) +
  labs(color = "month", x = "RDA1 (32.1 %)", y = "RDA2 (13 %)") +
  theme_bw() +
  scale_color_manual(values = c("#3CB22D", "#FF8000")) +
  scale_fill_manual(values = c("#3CB22D", "#FF8000")) +
  ylim(-1, 1) +
  xlim(-1.4, 1.6)
