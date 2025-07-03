# Purpose: Fit (G)LMM models to root trait data and visualize results

# ----------------- Load Packages -----------------
library(tidyverse)    # Data manipulation and plotting
library(performance)  # Model checks (collinearity, overdispersion)
library(ggsignif)     # Significance annotations on ggplots
library(sjPlot)       # Model visualization
library(ggeffects)    # For ggpredict() predictions
library(car)          # For model fitresults

# ----------------- Data Preparation -----------------
# Read in the root traits data and convert year to factor
Lys_data <- read_csv("data/RGG_Lys_Gradient.csv") %>%
  mutate(year = factor(year),
         hyphae_fract = no_counts / no_intersects)
str(Lys_data) # Check structure

# ----------------- Helper Functions -----------------
# Extract significance stars AND p-value for a predictor from ANOVA table
get_pval_and_sign <- function(anova_tab, pred) {
  # Find the p-value column name (works for lm, glm, etc.)
  pcol <- names(anova_tab)[grepl("^Pr\\(", names(anova_tab))]
  anova_tab %>%
    as_tibble(rownames = "predictor") %>%
    filter(predictor == pred) %>%
    mutate(
      p.value = !!sym(pcol),
      sign = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01 ~ "**",
        p.value < 0.05 ~ "*",
        TRUE ~ ""
      )
    ) %>%
    dplyr::select(predictor, p.value, sign)
}

# ----------------- (1) Specific Root Length (SRL) -----------------

# Fit linear model with interactions
SRLran1 <- lm(SRL_m_g ~ year * GW_level_cm + GW_level_cm * month, data = Lys_data)

# check residuals and collinearity
par(mfrow = c(2,2)); plot(SRLran1); par(mfrow = c(1,1))
check_collinearity(SRLran1)
car::Anova(SRLran1)

# Remove non-significant interactions
SRLran2 <- lm(SRL_m_g ~ year + GW_level_cm + month, data = Lys_data)
car::Anova(SRLran2)

# Likelihood ratio test
anova(SRLran1, SRLran2)

# final model
SRLran2

## Plots----

### Plot 1 (faceted by year and month)----
# Generate predicted values for plot 
SRLran2_pred <- ggpredict(SRLran2, terms = c("GW_level_cm [all]", "year", "month")) %>%
  as.data.frame() %>% rename(month = facet, year = group)

# Plot predicted SRL for each year/month, overlaying raw data
Plot_SRL_1 <- ggplot(SRLran2_pred, aes(x = x, y = predicted, color = month, fill = month)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1, color = NA) +
  geom_line(size = 0.8) +
  geom_point(data = Lys_data, aes(x = GW_level_cm, y = SRL_m_g, color = month), size = 2) +
  facet_wrap(~year) +
  labs(x = "Groundwater Level [cm]", y = "Specific root length [m/g]", color = "Month", fill = "Month") +
  scale_color_manual(values = c("#3CB22D", "#FF8000")) +
  scale_fill_manual(values = c("#3CB22D", "#FF8000")) +
  theme_bw()
Plot_SRL_1

### Plot 2 (groundwater level)----
SRLran2_pred2 <- ggpredict(SRLran2, terms = c("GW_level_cm [all]", "month")) %>%
  as.data.frame() %>% rename(month = group)

Plot_SRL_2 <- ggplot(SRLran2_pred2, aes(x = x, y = predicted, color = month, fill = month)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1, color = NA) +
  geom_line(size = 0.8) +
  geom_point(data = Lys_data, aes(x = GW_level_cm, y = SRL_m_g, color = month, pch = year), size = 2) +
  labs(x = "Groundwater Level [cm]", y = "Specific root length [m/g]") +
  scale_color_manual(values = c("#3CB22D", "#FF8000")) +
  theme_bw()
Plot_SRL_2

### Plot 3 (month, violin plot) ----

SRL_month_stats <- get_pval_and_sign(car::Anova(SRLran2), "month")
SRL_month_annot <- paste0(round(signif(SRL_month_stats$p.value, 3), 3), SRL_month_stats$sign)

Plot_SRL_3 <- ggplot(Lys_data, aes(x = month, y = SRL_m_g, fill = month)) +
  geom_violin(alpha = 0.7) +
  geom_signif(comparisons = list(c("June", "October")),
              y_position = max(Lys_data$SRL_m_g, na.rm = TRUE) * 1.01,
              annotations = c(SRL_month_annot)) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  scale_fill_manual(values = c("#3CB22D", "#FF8000")) +
  labs(x = "Seasonality", y = "Specific root length [m/g]") +
  theme_minimal() +
  ylim(155, 440)
Plot_SRL_3

### Plot 4 (year, violin plot)----

SRL_year_stats <- get_pval_and_sign(car::Anova(SRLran2), "year")
SRL_year_annot <- paste0(round(signif(SRL_year_stats$p.value, 3), 2), SRL_year_stats$sign)

Plot_SRL_4 <- ggplot(Lys_data, aes(x = year, y = SRL_m_g, fill = year)) +
  geom_violin(alpha = 0.7) +
  geom_signif(comparisons = list(c("2023", "2024")),
              y_position = max(Lys_data$SRL_m_g, na.rm = TRUE) * 1.01,
              annotations = c(SRL_year_annot)) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  labs(x = "Year", y = "Specific root length [m/g]") +
  theme_minimal() +
  ylim(155, 440)
Plot_SRL_4

# ----------------- (2) Root Tissue Density (RTD) -----------------
# Fit initial and log-transformed models, check diagnostics

RTDran1 <- lm(RTD_g_cm3 ~ year * GW_level_cm + GW_level_cm * month, data = Lys_data)
par(mfrow = c(2,2)); plot(RTDran1); par(mfrow = c(1,1))

RTDran1b <- lm(log(RTD_g_cm3) ~ year * GW_level_cm + GW_level_cm * month, data = Lys_data)
par(mfrow = c(2,2)); plot(RTDran1b); par(mfrow = c(1,1))

# check collinearity 
check_collinearity(RTDran1b)
car::Anova(RTDran1b)

# Remove non-significant interaction
RTDran2 <- lm(log(RTD_g_cm3) ~ year * GW_level_cm + month, data = Lys_data)

# model fit
car::Anova(RTDran2)

# check collinearity 
check_collinearity(RTDran2)

# Likelihood-ratio test
anova(RTDran1b, RTDran2)

# final model
RTDran2

## Plots----

### Plot 1 (faceted by year and month)----

RTDran2_pred <- ggpredict(RTDran2, terms = c("GW_level_cm [all]", "year", "month")) %>%
  as.data.frame() %>% rename(month = facet, year = group)

Plot_RTD_1 <- ggplot(RTDran2_pred, aes(x = x, y = predicted, color = month, fill = month)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1, color = NA) +
  geom_line(size = 0.8) +
  geom_point(data = Lys_data, aes(x = GW_level_cm, y = RTD_g_cm3, color = month), size = 2) +
  facet_wrap(~year) +
  labs(x = "Groundwater Level [cm]", y = "Root Tissue Density", color = "Month", fill = "Month") +
  scale_color_manual(values = c("#3CB22D", "#FF8000")) +
  scale_fill_manual(values = c("#3CB22D", "#FF8000")) +
  theme_bw()
Plot_RTD_1

### Plot 2 (groundwater level)----

RTDran2_pred2 <- ggpredict(RTDran2, terms = c("GW_level_cm [all]", "year")) %>%
  as.data.frame() %>% rename(year = group)

Plot_RTD_2 <- ggplot(RTDran2_pred2, aes(x = x, y = predicted, color = year, fill = year)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1, color = NA) +
  geom_line(size = 0.8) +
  geom_point(data = Lys_data, aes(x = GW_level_cm, y = RTD_g_cm3, color = year, pch = year), size = 2) +
  labs(x = "Groundwater Level [cm]", y = "Root Tissue Density") +
  theme_bw()

Plot_RTD_2

### Plot 3 (year, violin plot) ----
RTD_month_stats <- get_pval_and_sign(car::Anova(RTDran2), "month")
RTD_month_annot <- paste0("<0.001", RTD_month_stats$sign)

Plot_RTD_3 <- ggplot(Lys_data, aes(x = month, y = RTD_g_cm3, fill = month)) +
  geom_violin(alpha = 0.7) +
  geom_signif(comparisons = list(c("June", "October")),
              y_position = max(Lys_data$RTD_g_cm3, na.rm = TRUE) * 1.05,
              annotations = c(RTD_month_annot)) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  scale_fill_manual(values = c("#3CB22D", "#FF8000")) +
  labs(x = "Seasonality", y = "Root Tissue Density") +
  theme_minimal() + ylim(0.1, 0.32)

Plot_RTD_3

### Plot 4 (year, violin plot)----

RTD_year_stats <- get_pval_and_sign(car::Anova(RTDran2), "year")
RTD_year_annot <- paste0(">0.001", RTD_year_stats$sign)

Plot_RTD_4 <- ggplot(Lys_data, aes(x = year, y = RTD_g_cm3, fill = year)) +
  geom_violin(alpha = 0.7) +
  geom_signif(comparisons = list(c("2023", "2024")),
              y_position = max(Lys_data$RTD_g_cm3, na.rm = TRUE) * 1.05,
              annotations = c(RTD_year_annot)) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  labs(x = "Year", y = "Root Tissue Density") +
  theme_minimal() + ylim(0.1, 0.32)
Plot_RTD_4


# ----------------- (3) Average Root Diameter (RD) -----------------

# Fit linear model with interactions
RDran1 <- lm(RD_mm ~ year * GW_level_cm + GW_level_cm * month, data = Lys_data)

# check residuals and collinearity
par(mfrow = c(2,2)); plot(RDran1); par(mfrow = c(1,1))

# model fit
car::Anova(RDran1)

# Remove non-significant interaction
RDran2 <- lm(RD_mm ~ year + GW_level_cm + month, data = Lys_data)
car::Anova(RDran2)

# Likelihood-ratio test
anova(RDran1, RDran2)

# final model
RDran2

## Plots----

### Plot 1 (faceted by year and month)----

RDran2_pred <- ggpredict(RDran2, terms = c("GW_level_cm [all]", "year", "month")) %>%
  as.data.frame() %>% rename(month = facet, year = group)

Plot_RD_1 <- ggplot(RDran2_pred, aes(x = x, y = predicted, color = month, fill = month)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1, color = NA) +
  geom_line(size = 0.8) +
  geom_point(data = Lys_data, aes(x = GW_level_cm, y = RD_mm, color = month), size = 2) +
  facet_wrap(~year) +
  labs(x = "Groundwater Level [cm]", y = "Average Root Diameter [mm]", color = "Month", fill = "Month") +
  scale_color_manual(values = c("#3CB22D", "#FF8000")) +
  scale_fill_manual(values = c("#3CB22D", "#FF8000")) +
  theme_bw()

Plot_RD_1

### Plot 2 (groundwater level)----

RDran2_pred2 <- ggpredict(RDran2, terms = c("GW_level_cm [all]", "year")) %>%
  as.data.frame() %>% rename(year = group)

Plot_RD_2 <- ggplot(RDran2_pred2, aes(x = x, y = predicted, fill = year, color = year)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1, color = NA) +
  geom_line(size = 0.8) +
  geom_point(data = Lys_data, aes(x = GW_level_cm, y = RD_mm, color = year, pch = year), size = 2) +
  labs(x = "Groundwater Level [cm]", y = "Average Root Diameter [mm]") +
  theme_bw()

Plot_RD_2

### Plot 3 (month, violin plot) ----

RD_month_stats <- get_pval_and_sign(car::Anova(RDran2), "month")
RD_month_annot <- paste0(round(signif(RD_month_stats$p.value, 3), 3), RD_month_stats$sign)

Plot_RD_3 <- ggplot(Lys_data, aes(x = month, y = RD_mm, fill = month)) +
  geom_violin(alpha = 0.7) +
  geom_signif(comparisons = list(c("June", "October")),
              y_position = max(Lys_data$RD_mm, na.rm = TRUE) * 1.01,
              annotations = c(RD_month_annot)) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  scale_fill_manual(values = c("#3CB22D", "#FF8000")) +
  labs(x = "Seasonality", y = "Average Root Diameter [mm]") +
  theme_minimal()

Plot_RD_3

### Plot 4 (year, violin plot)----

RD_year_stats <- get_pval_and_sign(car::Anova(RDran2), "year")
RD_year_annot <- paste0("<0.001", RD_year_stats$sign)

Plot_RD_4 <- ggplot(Lys_data, aes(x = year, y = RD_mm, fill = year)) +
  geom_violin(alpha = 0.7) +
  geom_signif(comparisons = list(c("2023", "2024")),
              y_position = max(Lys_data$RD_mm, na.rm = TRUE) * 1.01,
              annotations = c(RD_year_annot)) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  labs(x = "Year", y = "Average Root Diameter [mm]") +
  theme_minimal() + ylim(0.13, 0.22)

Plot_RD_4

# ----------------- (4) Roots colonized with AMF -----------------

# Fit initial model (binomial), check for overdispersion
Hyphran1 <- glm(hyphae_fract ~ year * GW_level_cm + GW_level_cm * month,
                data = Lys_data, family = binomial, weights = no_intersects)
summary(Hyphran1)
check_overdispersion(Hyphran1)

# Overdispersed -> use quasibinomial family
Hyphran2 <- glm(hyphae_fract ~ year * GW_level_cm + GW_level_cm * month,
                data = Lys_data, family = quasibinomial, weights = no_intersects)
summary(Hyphran2)
car::Anova(Hyphran2)
drop1(Hyphran2, test = "Chisq")

# Remove non-significant interaction
Hyphran3 <- glm(hyphae_fract ~ year + GW_level_cm + month,
                data = Lys_data, family = quasibinomial, weights = no_intersects)

# Likelihood-ratio test
anova(Hyphran2, Hyphran3, test="F")

# final model
Hyphran3

# model fit
car::Anova(Hyphran3)

## Plots----

### Plot 1 (faceted by year and month)----

Hyphran3_pred <- ggpredict(Hyphran3, terms = c("GW_level_cm [all]", "year", "month")) %>%
  as.data.frame() %>% rename(month = facet, year = group)

Plot_Hyph_1 <- ggplot(Hyphran3_pred, aes(x = x, y = predicted, color = month, fill = month)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1, color = NA) +
  geom_line(size = 0.8) +
  geom_point(data = Lys_data, aes(x = GW_level_cm, y = hyphae, color = month), size = 2) +
  facet_wrap(~year) +
  labs(x = "Groundwater Level [cm]", y = "Roots colonized with AMF", color = "Month", fill = "Month") +
  scale_color_manual(values = c("#3CB22D", "#FF8000")) +
  scale_fill_manual(values = c("#3CB22D", "#FF8000")) +
  theme_bw()+
  scale_y_continuous(labels = scales::label_percent())


Plot_Hyph_1

### Plot 2 (groundwater level)----

Hyphran3_pred2 <- ggpredict(Hyphran3, terms = c("GW_level_cm [all]", "year")) %>%
  as.data.frame() %>% rename(year = group)

Plot_Hyph_2 <- ggplot(Hyphran3_pred2, aes(x = x, y = predicted, color = year, fill = year)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1, color = NA) +
  geom_line(size = 0.8) +
  geom_point(data = Lys_data, aes(x = GW_level_cm, y = hyphae, color = year), size = 2) +
  labs(x = "Groundwater Level [cm]", y = "Roots colonized with AMF", color = "Year", fill = "Year") +
  theme_bw()+
  scale_y_continuous(labels = scales::label_percent())


Plot_Hyph_2

### Plot 3 (month, violin plot) ----

AMF_month_stats <- get_pval_and_sign(car::Anova(Hyphran3), "month")

AMF_month_annot <- paste0(round(signif(AMF_month_stats$p.value, 3), 3), AMF_month_stats$sign)

Plot_Hyph_3 <- ggplot(Lys_data, aes(x = month, y = hyphae, fill = month)) +
  geom_violin(alpha = 0.7) +
  geom_signif(comparisons = list(c("June", "October")),
              y_position = max(Lys_data$hyphae, na.rm = TRUE) * 1.01,
              annotations = c(AMF_month_annot)) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  scale_fill_manual(values = c("#3CB22D", "#FF8000")) +
  labs(x = "Seasonality", y = "Roots colonized with AMF") +
  theme_minimal()+
  scale_y_continuous(labels = scales::label_percent())

Plot_Hyph_3

### Plot 4 (year, violin plot)----

AMF_year_stats <- get_pval_and_sign(car::Anova(Hyphran3), "year")
AMF_year_annot <- paste0(round(signif(AMF_year_stats$p.value, 3), 3), AMF_year_stats$sign)

Plot_Hyph_4 <- ggplot(Lys_data, aes(x = year, y = hyphae, fill = year)) +
  geom_violin(alpha = 0.7) +
  geom_signif(comparisons = list(c("2023", "2024")),
              y_position = max(Lys_data$hyphae, na.rm = TRUE) * 1.01,
              annotations = c(AMF_year_annot)) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  labs(x = "Year", y = "Roots colonized with AMF") +
  theme_minimal() + 
  scale_y_continuous(labels = scales::label_percent(), limits = c(0, 1))

Plot_Hyph_4


# ---- Combine plots ----
library(patchwork)

combined_plot_1 <- (Plot_SRL_1 + Plot_RTD_1) / 
  plot_spacer() / 
  (Plot_RD_1 + Plot_Hyph_1) +
  plot_layout(heights = c(10, 0.5, 10), guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold', size=18),
        plot.tag.position  = c(0.11, 1.05),
        plot.margin = margin(t = 20, r = 5, b = 5, l = 5),
        plot.title = element_text(hjust = 0.5, face = 'bold', size=15))

print(combined_plot_1)

ggsave("figures/interaction_effects.png", combined_plot_1, width = 10, height = 7.5, dpi = 150)

combined_plot_3 <- (Plot_SRL_3 + Plot_RTD_3) / 
  plot_spacer() / 
  (Plot_RD_3 + Plot_Hyph_3) +
  plot_layout(heights = c(10, 1, 10), guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold', size=20),
        plot.tag.position  = c(0.2, 1.07),
        plot.margin = margin(t = 20, r = 5, b = 5, l = 5),
        plot.title = element_text(hjust = 0.5, face = 'bold', size=15))

print(combined_plot_3)

ggsave("figures/month_effcet.png", combined_plot_3, width = 7, height = 9, dpi = 150)


combined_plot_4 <- (Plot_SRL_4 + Plot_RTD_4) / 
  plot_spacer() / 
  (Plot_RD_4 + Plot_Hyph_4) +
  plot_layout(heights = c(10, 1, 10), guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold', size=20),
        plot.tag.position  = c(0.2, 1.07),
        plot.margin = margin(t = 20, r = 5, b = 5, l = 5),
        plot.title = element_text(hjust = 0.5, face = 'bold', size=15))

print(combined_plot_4)

ggsave("figures/year_effcet.png", combined_plot_4, width = 7, height = 9, dpi = 150)



# --- Extract model tables  ----

library(MuMIn)
library(tibble)
library(dplyr)
library(car)

# --- Helper: Get model formula as character string ---
pretty_formula <- function(model) {
  f <- deparse(formula(model))
  f <- gsub("\\s+", " ", f)
  f <- gsub(" +", " ", f)
  f
}

get_family_chr <- function(model) {
  if ("family" %in% names(model)) {
    fam <- model$family$family
    link <- model$family$link
    paste0("GLM (", fam, ", link = ", link, ")")
  } else if ("lm" %in% class(model)) {
    "Linear Model (LM)"
  } else {
    class(model)[1]
  }
}

# --- Helper:  R2 (MuMIn) ---
get_r2_theoretical <- function(model) {
  # For LM and GLM, MuMIn::r.squaredGLMM works and gives R2m and R2c
  # We'll use only the marginal R2 (R2m) as we have no random effcets
  r2 <- tryCatch(MuMIn::r.squaredGLMM(model), error = function(e) NULL)
  if (!is.null(r2)) {
    # Handles both LM and GLMM/GLM
    r2m <- as.numeric(r2[1])
    return(sprintf("%.3f", r2m))
  }
  # fallback for plain LM
  if (inherits(model, "lm")) {
    r2 <- summary(model)$r.squared
    return(sprintf("%.3f", r2))
  }
  return(NA)
}

# --- 1. Final Model Formulas Table ---
final_models_list <- list(
  SRLran2 = SRLran2,
  RTDran2 = RTDran2,
  RDran2 = RDran2,
  Hyphran3 = Hyphran3
)

final_model_formulas <- tibble::tibble(
  Trait = c(
    "Specific Root Length (SRL)",
    "Root Tissue Density (RTD)",
    "Root Diameter (RD)",
    "AMF Colonization (AMF)"
  ),
  Model_Object = names(final_models_list),
  Model_Formula = vapply(final_models_list, pretty_formula, character(1)),
  Model_Type = vapply(final_models_list, get_family_chr, character(1)),
  R2 = vapply(final_models_list, get_r2_theoretical, character(1))
)
write.csv(final_model_formulas, "tables/final_model_information.csv", row.names = FALSE)

# --- 2. Likelihood Ratio Test Results ---
get_lrt_table <- function(mod1, mod2, name1, name2, test_type = "F") {
  aov_tab <- anova(mod1, mod2, test = test_type)
  tibble::tibble(
    Comparison = paste(name1, "vs", name2),
    Model_1_Formula = pretty_formula(mod1),
    Model_2_Formula = pretty_formula(mod2),
    Test_Statistic = round(aov_tab[2, "F"], 3),
    Df = paste(aov_tab[2, "Df"], collapse = "/"),
    p_value = signif(aov_tab[2, "Pr(>F)"], 3)
  )
}
lrt_SRL <- get_lrt_table(SRLran1, SRLran2, "SRLran1", "SRLran2")
lrt_RTD <- get_lrt_table(RTDran1b, RTDran2, "RTDran1b", "RTDran2")
lrt_RD <- get_lrt_table(RDran1, RDran2, "RDran1", "RDran2")
lrt_Hyph <- {
  aov_tab <- anova(Hyphran2, Hyphran3, test = "F")
  tibble::tibble(
    Comparison = "Hyphran2 vs Hyphran3",
    Model_1_Formula = pretty_formula(Hyphran2),
    Model_2_Formula = pretty_formula(Hyphran3),
    Test_Statistic = round(aov_tab[2, "F"], 3),
    Df = paste(aov_tab[2, "Df"], collapse = "/"),
    p_value = signif(aov_tab[2, "Pr(>F)"], 3)
  )
}
lrt_table <- dplyr::bind_rows(lrt_SRL, lrt_RTD, lrt_RD, lrt_Hyph)
write.csv(lrt_table, "tables/likelihood_ratio_tests.csv", row.names = FALSE)

# --- 3. car::Anova Results for Each Final Model ---
get_anova_table <- function(model, model_name) {
  a <- car::Anova(model)
  predictors <- rownames(a)
  results <- tibble::tibble(
    Model = model_name,
    Model_Formula = pretty_formula(model),
    Predictor = predictors,
    Statistic = round(a[, 1], 3),
    Df = a[, "Df"],
    p_value = signif(a[, ncol(a)], 3),
    Significance = case_when(
      a[, ncol(a)] < 0.001 ~ "***",
      a[, ncol(a)] < 0.01 ~ "**",
      a[, ncol(a)] < 0.05 ~ "*",
      TRUE ~ ""
    )
  )
  results
}
anova_SRL <- get_anova_table(SRLran2, "SRLran2")
anova_RTD <- get_anova_table(RTDran2, "RTDran2")
anova_RD <- get_anova_table(RDran2, "RDran2")
anova_Hyph <- get_anova_table(Hyphran3, "Hyphran3")
anova_all <- dplyr::bind_rows(anova_SRL, anova_RTD, anova_RD, anova_Hyph)
write.csv(anova_all, "tables/Anova_final_models.csv", row.names = FALSE)

