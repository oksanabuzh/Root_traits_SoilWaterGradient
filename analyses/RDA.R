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

# ---- Permutation tests (model fit) ----
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
perm_2023

perm_2024 <- perm_test(ordin2024)
perm_2024

# ---- Goodness of fit for traits ----
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
plot_2023_v1 <- ggplot(observ.scrs_2023) +
  geom_hline(yintercept = 0, color = "grey", lty = 1) +
  geom_vline(xintercept = 0, color = "grey", lty = 1) +
  # stat_ellipse(aes(RDA1, RDA2, fill=month), alpha=.1,type='t', linewidth =1, geom="polygon") +
  geom_segment(aes(x = RDA1, y = RDA2, xend = RDA1_centr, yend = RDA2_centr, color = month),
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
  labs(color = "month", x = "RDA1 (32.7 %)", y = "RDA2 (4.8 %)", title = "2023") +
  theme_bw() +
  scale_color_manual(values = c("#3CB22D", "#FF8000")) +
  scale_fill_manual(values = c("#3CB22D", "#FF8000"))

plot_2023_v1
  
# ---- 2023: Plot only traits and predictors ----
plot_2023_v2 <- ggplot() +
  geom_hline(yintercept = 0, color = "grey", lty = 1) +
  geom_vline(xintercept = 0, color = "grey", lty = 1) +
  geom_point(data = centroids_2023, aes(x = RDA1, y = RDA2, color = month), size = 4, alpha = .5, stroke = 2) +
  geom_text(data = centroids_2023, aes(RDA1, RDA2, label = month, color = month), size = 5, show.legend = F, vjust = c(-1, -1), hjust = c(0, 1)) +
  geom_text(data = trait.scrs_2023, aes(RDA1, RDA2, label = traits_names),
            color = "black", vjust = c(-1.1, 1.2, -0.6, 0.2), hjust = c(0.5, 0, 0.5, -0.1)) +
  geom_segment(data = trait.scrs_2023, aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black", linewidth = 0.9) +
  geom_text(data = env.scrs_2023, aes(RDA1, RDA2, label = predictors_names),
            color = "blue", size = 5, vjust = -0.5, hjust = 0.7) +
  geom_segment(data = env.scrs_2023, aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.3, "cm")), color = "blue", linewidth = 1) +
  labs(color = "month", x = "RDA1 (32.7 %)", y = "RDA2 (4.8 %)", title = "2023") +
  theme_bw() +
  scale_color_manual(values = c("#3CB22D", "#FF8000")) +
  scale_fill_manual(values = c("#3CB22D", "#FF8000")) +
  ylim(-0.5, 0.9) +
  xlim(-1.2, 1.2)

plot_2023_v2

# ---- 2024: Plot all elements ----
plot_2024_v1 <- ggplot(observ.scrs_2024) +
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
  labs(color = "month", x = "RDA1 (32.1 %)", y = "RDA2 (13 %)", title = "2024") +
  theme_bw() +
  scale_color_manual(values = c("#3CB22D", "#FF8000")) +
  scale_fill_manual(values = c("#3CB22D", "#FF8000"))+
  xlim(-1.5, 2.2)

plot_2024_v1
# ---- 2024: Plot only traits and predictors ----
plot_2024_v2 <- ggplot() +
  geom_hline(yintercept = 0, color = "grey", lty = 1) +
  geom_vline(xintercept = 0, color = "grey", lty = 1) +
  geom_point(data = centroids_2024, aes(x = RDA1, y = RDA2, color = month), size = 4, alpha = .5, stroke = 2) +
  geom_text(data = centroids_2024, aes(RDA1, RDA2, label = month, color = month), show_guide = F, size = 5, vjust = c(-1, -1), hjust = c(0, 1)) +
  geom_text(data = trait.scrs_2024, aes(RDA1, RDA2, label = traits_names),
            color = "black", vjust = c(-1.1, 1.2, 1.2, 1.2), hjust = c(0.7, 0.7, 0.5, 0.7)) +
  geom_segment(data = trait.scrs_2024, aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black", linewidth = 0.9) +
  geom_text(data = env.scrs_2024, aes(RDA1, RDA2, label = predictors_names),
            color = "blue", size = 5, vjust = -0.5, hjust = 0.7) +
  geom_segment(data = env.scrs_2024, aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.3, "cm")), color = "blue", linewidth = 1) +
  labs(color = "month", x = "RDA1 (32.1 %)", y = "RDA2 (13 %)", title = "2024") +
  theme_bw() +
  scale_color_manual(values = c("#3CB22D", "#FF8000")) +
  ylim(-1, 1) +
  xlim(-1.4, 1.6)

plot_2024_v2

# ---- Combine plots ----
library(patchwork)

combined_plot_v1 <- plot_2023_v1 + plot_2024_v1 +
  plot_layout(nrow = 1, ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold', size=20),
        plot.tag.position  = c(0.1, 1.01),
        plot.margin = margin(t = 20, r = 5, b = 5, l = 5),
        plot.title = element_text(hjust = 0.5, face = 'bold', size=15))

print(combined_plot_v1)


ggsave("figures/RDA_plot_v1.png", combined_plot_v1, width = 13, height = 6, dpi = 150)



combined_plot_v2 <- plot_2023_v2 + plot_2024_v2 +
  plot_layout(nrow = 1, ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold', size=20),
        plot.tag.position  = c(0.1, 1.01),
        plot.margin = margin(t = 20, r = 5, b = 5, l = 5),
        plot.title = element_text(hjust = 0.5, face = 'bold', size=15))

print(combined_plot_v2)

ggsave("figures/RDA_plot_v2.png", combined_plot_v2, width = 13, height = 6, dpi = 150)


# ---- Extract and present permutation test results from perm_2023 and perm_2024 as publication-ready tables ----

library(dplyr)
library(knitr)
library(tibble)

# Vector of years and types 
years <- c(2023, 2024)
types <- c("anova", "anova_margin", "anova_axis")

# Create a list of all results with labels 
anova_list <- lapply(years, function(yr) {
  perm_obj <- get(paste0("perm_", yr))
  lapply(types, function(tp) {
    as.data.frame(perm_obj[[tp]]) %>%
      tibble::rownames_to_column("Effect") %>%
      mutate(Test = tp, Year = yr)
  }) %>% bind_rows()
}) %>% bind_rows() %>% 
  mutate(Variance=round(Variance, 2),
         F=round(F, 2)) %>% 
  mutate(Type=case_when(Test=="anova" ~ "model fit",
                        Test=="anova_margin" ~ "effects of predictors",
                        Test=="anova_axis" ~ "fit for RDA axes")) %>% 
 mutate(Effect=case_when(Effect=="GW_level_cm" ~ "Groundwater level",
                        Effect=="month" ~ "Months",
                        .default = Effect))

# Save combined table
kable(anova_list, caption = "Permutation ANOVA Results for 2023 and 2024")

library(knitr)
library(kableExtra)

# Show as a nice table for presentation
kable(anova_list, caption = "Permutation Results for 2023 and 2024") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                full_width = FALSE,
                position = "center") %>%
  column_spec(1, bold = TRUE) # Make 'Effect' column bold for clarity

# To split and print by year and test type, using split/filter:
for (yr in years) {
  for (tp in types) {
    tab <- dplyr::filter(anova_list, Year == yr, Test == tp)
    kable(tab, caption = paste(yr, toupper(gsub("anova", "Permutation Test", tp))))
  }
}

# Export all as one clean CSV for supplementary material
write.csv(anova_list, "tables/RDA_permutation_results.csv", row.names = FALSE)
