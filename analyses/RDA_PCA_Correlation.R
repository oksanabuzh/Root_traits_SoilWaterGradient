# Purpose: Fit ordination (RDA) for root traits and their predictors for 2023 and 2024

# ---- Load packages ----
library(tidyverse)
library(vegan)
library(ggrepel)
library(tibble)
library(lattice)
library(patchwork)


# ---- Read and Prepare Data ----
Lys_data <- read_csv("data/RGG_Lys_Gradient.csv") %>%
  unite("Plot_ID", c("Lys_no", "date"), remove = TRUE)

names(Lys_data)

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

names(perm_2024)


p1 <- ggplotify::as.ggplot(perm_2023$densityplot)+
  ggtitle("First Sampling Year") +
  theme(plot.title = element_text(hjust = 0.5))

p2 <- ggplotify::as.ggplot(perm_2024$densityplot) +
ggtitle("Second Sampling Year") +
  theme(plot.title = element_text(hjust = 0.5))


perm_density_combined <- p1 + p2 +  
 # plot_layout(nrow = 1, ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold', size=20),
        plot.tag.position  = c(0.2, 1.04),
        plot.margin = margin(t = 20, r = 5, b = 5, l = 5),
        plot.title = element_text(hjust = 0.5, face = 'bold', size=15))

print(perm_density_combined)

ggsave("figures/perm_density_plot.png", perm_density_combined, width = 10, height = 5, dpi = 150)


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
      traits=="hyphae" ~"Roots AMF",
      traits=="leafN_mgg" ~"Leaf N",
      traits=="rootN_mgg" ~"Root N"
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
  stat_ellipse(aes(RDA1, RDA2, fill=month), alpha=.1,type='t', linewidth =1, geom="polygon") +
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
  labs(color = "Month", fill = "Month", x = "RDA1 (32.7 %)", y = "RDA2 (4.8 %)", 
       title = "First year") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold")) +
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
  labs(color = "Month", fill = "Month", x = "RDA1 (32.7 %)", y = "RDA2 (4.8 %)", 
       title = "First year") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold")) +
  scale_color_manual(values = c("#3CB22D", "#FF8000")) +
  ylim(-0.5, 0.9) +
  xlim(-1.2, 1.2)

plot_2023_v2

# ---- 2024: Plot all elements ----
plot_2024_v1 <- ggplot(observ.scrs_2024) +
  geom_hline(yintercept = 0, color = "grey", lty = 1) +
  geom_vline(xintercept = 0, color = "grey", lty = 1) +
  stat_ellipse(aes(RDA1, RDA2, fill=month), alpha=.1,type='t', linewidth =1, geom="polygon") +
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
  labs(color = "Month", fill = "Month", x = "RDA1 (32.1 %)", y = "RDA2 (13 %)", 
       title = "Second year") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold")) +
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
  labs(color = "Month", fill = "Month", x = "RDA1 (32.1 %)", y = "RDA2 (13 %)", 
       title = "Second year") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold")) +
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


ggsave("figures/RDA_plot_v1.png", combined_plot_v1, width = 13, height = 7, dpi = 150)



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

library(knitr)
library(kableExtra)
kable(anova_list, caption = "Permutation ANOVA Results for 2023-2025")


# Show as a nice table for presentation
kable(anova_list, caption = "Permutation Results for 2023-2025") %>%
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
# write.csv(anova_list, "tables/RDA_permutation_results.csv", row.names = FALSE)



# ---- Run PCA ----

PCA2023 <- rda(traits_2023,  scale = TRUE)
PCA2024 <- rda(traits_2024, scale = TRUE)
PCA_all <- rda(Lys_data %>% 
                 filter(!year==2025) %>% 
                 select(Plot_ID, SRL_m_g, RD_mm, RTD_g_cm3, hyphae) %>%
                 column_to_rownames("Plot_ID"), scale = TRUE)

summary(eigenvals(PCA2023))
summary(eigenvals(PCA2024))
summary(eigenvals(PCA_all))

# Extract scores ----------
# --- For 2023 ---
PCA_trait.scrs_2023 <- extract_trait_scrs(PCA2023)
PCA_env.scrs_2023 <- scores(PCA2023, display = "sites", scaling = "species") %>%
  as_tibble(rownames="Plot_ID") %>%
  left_join(predictors_2023, by="Plot_ID") 

# --- For 2024 ---
PCAtrait.scrs_2024 <- extract_trait_scrs(PCA2024)
PCA_env.scrs_2024 <- scores(PCA2024, display = "sites", scaling = "species") %>%
  as_tibble(rownames="Plot_ID") %>%
  left_join(predictors_2024, by="Plot_ID") 


# --- For all ---
PCAtrait.scrs_all <- extract_trait_scrs(PCA_all)
PCA_env.scrs_all <- scores(PCA_all, display = "sites", scaling = "species") %>%
  as_tibble(rownames="Plot_ID") %>%
  left_join(Lys_data %>% 
              filter(!year==2025) %>%
              select(Plot_ID, year, month, GW_level_cm), 
            by="Plot_ID") 


# ---- Plotting ordination results ----

# ---- 2023 ----
plot_2023_PCA <- ggplot(PCA_trait.scrs_2023) +
  geom_hline(yintercept = 0, color = "grey", lty = 1) +
  geom_vline(xintercept = 0, color = "grey", lty = 1) +
  geom_point(data=PCA_env.scrs_2023, aes(PC1, PC2, color = month), size = 3, alpha = .5) +
    geom_text(data = PCA_trait.scrs_2023, aes(PC1, PC2, label = traits_names),
            color = "black", vjust = c(-1.1, 1.2, -0.6, 0.2), hjust = c(0.8, 0, 0, -0.1)) +
  geom_segment(data = PCA_trait.scrs_2023, aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black", linewidth = 0.9) +
 labs(x = "PC1 (72.1 %)", y = "PC2 (1.5 %)", 
       title = "First year") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold")) +
  scale_color_manual(values = c("#3CB22D", "#FF8000")) +
   # ylim(-0.5, 0.9) +
  xlim(-2, 2)

plot_2023_PCA

# ---- 2024 ----
plot_2024_PCA <- ggplot(PCAtrait.scrs_2024) +
  geom_hline(yintercept = 0, color = "grey", lty = 1) +
  geom_vline(xintercept = 0, color = "grey", lty = 1) +
  geom_point(data=PCA_env.scrs_2024, aes(PC1, PC2, color = month), size = 3, alpha = .5) +
  geom_text(data = PCAtrait.scrs_2024, aes(PC1, PC2, label = traits_names),
            color = "black", vjust = c(1.8, 1.2, 1.1, 0.2), 
                             hjust = c(0.8, -0.2, 0, -0.1)) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black", linewidth = 0.9) +
  labs(x = "PC1 (41.4 %)", y = "PC2 (33.8 %)", 
       title = "Second year") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold")) +
  scale_color_manual(values = c("#3CB22D", "#FF8000")) +
  # ylim(-0.5, 0.9) +
  xlim(-2, 2)

plot_2024_PCA


# ---- all  ----
plot_all_PCA <- ggplot(PCAtrait.scrs_all) +
  geom_hline(yintercept = 0, color = "grey", lty = 1) +
  geom_vline(xintercept = 0, color = "grey", lty = 1) +
  geom_point(data=PCA_env.scrs_all, aes(PC1, PC2, color = factor(year)), size = 3, alpha = .5) +
  geom_text(data = PCAtrait.scrs_all, aes(PC1, PC2, label = traits_names),
            color = "black", vjust = c(1.8, 1.2, 1.1, 0.2), 
            hjust = c(0.5, -0.2, 0.7, -0.1)) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black", linewidth = 0.9) +
  labs(# title = "Overall dataset"
    x = "PC1 (46.6 %)", y = "PC2 (34.8 %)",
    color="Year") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold")) +
  scale_color_manual(values = c("#E03C53", "#4C72FE"),
                                labels=c("2023" = "first year",
                                         "2024" = "second year")) +
  xlim(-2, 2)

plot_all_PCA
ggsave("figures/PCA_plot_allData.png", plot_all_PCA, width = 10, height = 7.5, dpi = 150)


# ---- Combine plots for 2023 and 2024 ----
library(patchwork)

combined_PCA_plots <- plot_2023_PCA + plot_2024_PCA +  
  plot_layout(nrow = 1, ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold', size=20),
        plot.tag.position  = c(0.1, 1.01),
        plot.margin = margin(t = 20, r = 5, b = 5, l = 5),
        plot.title = element_text(hjust = 0.5, face = 'bold', size=15))

print(combined_PCA_plots)


ggsave("figures/PCA_plot_2023_2024.png", combined_PCA_plots, width = 13, height = 7, dpi = 150)



# Correlation analysis ----------------

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
   # colors = c("red", "white", "blue"),
    colors = c("#6D9EC1", "white", "#E46726"),
        legend.title = "Pearson r") +
    ggtitle(title)
}

?ggcorrplot
# ---- Create correlation plots ----

p_2023 <- make_corr_plot(traits_2023, title = "First Sampling Year")
p_2024 <- make_corr_plot(traits_2024, title = "Second Sampling Year")
p_all <- make_corr_plot(Lys_data %>% 
                          filter(!year==2025) %>% 
                          select(Plot_ID, SRL_m_g, RD_mm, RTD_g_cm3, hyphae) %>%
                          column_to_rownames("Plot_ID"), title = "Overall Root Trait Correlation")


# Combine plots 

library(patchwork)

A <- plot_all_PCA +
  theme(plot.tag.position = c(0.1, 1.05))

B <- p_all +
  theme(plot.tag.position = c(0.4, 0.91))   

C <- p_2023 +
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        plot.tag.position = c(0.03, 0.91))

D <- p_2024 +
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        plot.tag.position = c(0.03, 0.91))


top <- (plot_spacer() + A + plot_spacer()) +
  plot_layout(ncol = 3, widths = c(1, 6, 1), guides = "collect")&
  theme(plot.margin = margin(t = 20, r = 0, b = 0, l = 0))


bottom <- (B + C + D) +
  plot_layout(ncol = 3, widths = c(1, 1, 1), 
              guides = "collect") &
  theme(plot.margin = margin(t = 10, r = 5, b = 5, l = 20))

combined <- top / bottom +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold', size = 15),
        # plot.tag.position  = c(0.18, 1.05),
        plot.title = element_text(hjust = 0.6, vjust=-18, #face = 'bold',
                                  size = 11))

print(combined)

ggsave("figures/correlation_plot.png", combined, width = 12, height = 12, dpi = 150)

