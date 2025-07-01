# Purpose: to fit (G)LMM models

# packeges ----
library(tidyverse)
library(readxl)
library(lmerTest) # for GLMs
library(performance) # for checking overdispersion
library(emmeans) # for posthoc test
library(multcomp) # for posthoc tests wit the letters
library(ggsignif)

library(sjPlot) # to plot predictions from the model

# -----------------------------------------------------------------------------#
# Data preparation -----
# -----------------------------------------------------------------------------#

# read data ----
Lys_data <- read_excel("data/RGG_Lys_Gradient.xlsx") %>% 
  mutate(year=factor(year))
str(Lys_data) # check the structure of the data


# -----------------------------------------------------------------------------#
# Analysis -----
# -----------------------------------------------------------------------------#

# (1) Specific root length------

Lys_data %>% 
  ggplot(aes(GW_level_cm, SRL_m_g))+
           geom_point()
         

SRLran1 <- lm(SRL_m_g ~ year * GW_level_cm + GW_level_cm * month, 
                data = Lys_data)

# residuals:
par(mfrow=c(2,2))
plot(SRLran1)
par(mfrow=c(1,1))

# check singularity of predictors
check_collinearity(SRLran1)

# model fit
car::Anova(SRLran1)

# remove non significant interaction: 
SRLran2 <- lm(SRL_m_g ~ year + GW_level_cm + month, #rewett_d,  
              data = Lys_data)


car::Anova(SRLran2)



# Plots:

## faceted plot ---- 
library(ggeffects) # for ggpredict()
SRLran2_pred <-  ggpredict(SRLran2,
                           terms = c("GW_level_cm [all]", "year", "month")) %>% 
  as.data.frame() %>% 
  # ggpredict returns a data frame with columns:
  # x (GW_level_cm), predicted (on original scale), conf.low, conf.high, group (month), facet (year)
  rename(month=facet, year=group)

# Plot
library(ggplot2)

ggplot(SRLran2_pred, aes(x = x, y = predicted, color = month, fill = month)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.1, color = NA) +
  geom_line(size = 0.8) +  
  geom_point(
    data = Lys_data,
    aes(x = GW_level_cm, y = SRL_m_g, color = month),
    size = 2
  ) +
  facet_wrap(~year) +
  labs(
    x = "Groundwater Level [cm]",
    y = "Specific root length [m/g]",
    color = "month",
    fill = "month"
  ) +
  scale_color_manual(values = c("#3CB22D", "#FF8000")) +
  scale_fill_manual(values = c("#3CB22D", "#FF8000")) +
  theme_bw() 

## water levels effects ------
library(ggeffects) # for ggpredict()
SRLran2_pred2 <-  ggpredict(SRLran2,
                           terms = c("GW_level_cm [all]", "month")) %>% 
  as.data.frame() %>% 
  # ggpredict returns a data frame with columns:
  # x (GW_level_cm), predicted (on original scale), conf.low, conf.high, group (month), facet (year)
  rename(month=group)

# Plot
ggplot(SRLran2_pred2, aes(x = x, y = predicted, color = month, fill = month)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.1, color = NA) +
  geom_line(size = 0.8) +  
  geom_point(
    data = Lys_data,
    aes(x = GW_level_cm, y = SRL_m_g, color=month, pch=year),
    size = 2
  ) +
 # facet_wrap(~year) +
  labs(
    x = "Groundwater Level [cm]",
    y = "Specific root length [m/g]") +
  scale_color_manual(values = c("#3CB22D", "#FF8000")) +
#  scale_fill_manual(values = c("#3CB22D", "#FF8000")) +
  theme_bw() 

#------------------------------------------------------------------------#
##  month effects -----

# Extract p-value for June vs October
p_value <- car::Anova(SRLran2) %>% 
  as_tibble(rownames = "predictor") %>% 
  filter(predictor=="month") %>% 
  mutate(sign = ifelse( `Pr(>F)` < 0.001, "***", 
                        ifelse(`Pr(>F)` < 0.01, "**", 
                               ifelse( `Pr(>F)`< 0.05, "*", "ns"))))


P1 <- ggplot(Lys_data, aes(x = month, y = SRL_m_g, fill = month)) +
  geom_violin(alpha = 0.7) +
  geom_signif(
    comparisons = list(c("June", "October")),
    y_position = c(max(Lys_data$SRL_m_g, na.rm = TRUE) * 1.01), # Define y position for the significance bar (adjust as needed)
    annotations = c(p_value$sign)
    ) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  scale_fill_manual(values = c("#3CB22D", "#FF8000")) +
  labs(x = "Seasonality", y = "Specific root length [m/g]") +
  theme_minimal()+
  ylim(155, 440)

P1
ggsave(filename = "images/SRLseason.png", plot = P1, width = 18, height = 12, units = "cm")

## year effect -----
# Extract p-value for June vs October
p_value2 <- car::Anova(SRLran2) %>% 
  as_tibble(rownames = "predictor") %>% 
  filter(predictor=="year") %>% 
  mutate(sign = ifelse( `Pr(>F)` < 0.001, "***", 
                        ifelse(`Pr(>F)` < 0.01, "**", 
                               ifelse( `Pr(>F)`< 0.05, "*", "ns"))))


ggplot(Lys_data, aes(x = year, y = SRL_m_g, fill = year)) +
  geom_violin(alpha = 0.7) +
  geom_signif(
    comparisons = list(c("2023", "2024")),
    y_position = c(max(Lys_data$SRL_m_g, na.rm = TRUE) * 1.01), # Define y position for the significance bar (adjust as needed)
    annotations = c(p_value2$sign)
  ) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
 # scale_fill_manual(values = c("#3CB22D", "#FF8000")) +
  labs(x = "Year", y = "Specific root length [m/g]") +
  theme_minimal()+
  ylim(155, 440)

# -----------------------------------------------------------------------------#

# (2) Root tissue density ------

Lys_data %>% 
  ggplot(aes(GW_level_cm, RTD_g_cm3))+
  geom_point()


RTDran1 <- lm(RTD_g_cm3 ~ year * GW_level_cm + GW_level_cm*month, 
              data = Lys_data)

# residuals:
par(mfrow=c(2,2))
plot(RTDran1)
par(mfrow=c(1,1))

# transform data to improve normal distribution of residuals:
RTDran1b <- lm(log(RTD_g_cm3) ~ year * GW_level_cm + GW_level_cm*month, 
              data = Lys_data)

# residuals:
par(mfrow=c(2,2))
plot(RTDran1b)
par(mfrow=c(1,1))

# check singularity of predictors
check_collinearity(RTDran1b)

car::Anova(RTDran1b)

# remove non significant interaction 
RTDran2 <- lm(log(RTD_g_cm3) ~ year * GW_level_cm +  month,   
              data = Lys_data)

car::Anova(RTDran2)
# check singularity of predictors
check_collinearity(RTDran2)


# Create the plot
## faceted plot -----
library(ggeffects) # for ggpredict()
RTDran2_pred <-  ggpredict(RTDran2,
  terms = c("GW_level_cm [all]", "year", "month")) %>% 
  as.data.frame() %>% 
   # ggpredict returns a data frame with columns:
   # x (GW_level_cm), predicted (on original scale), conf.low, conf.high, group (month), facet (year)
 rename(month=facet, year=group)

# Plot
library(ggplot2)

P2 <- ggplot(RTDran2_pred, aes(x = x, y = predicted, color = month, fill = month)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.1, color = NA) +
  geom_line(size = 0.8) +  
  geom_point(
    data = Lys_data,
    aes(x = GW_level_cm, y = RTD_g_cm3, color = month),
    size = 2
  ) +
  facet_wrap(~year) +
  labs(
    x = "Groundwater Level [cm]",
    y = "Root Tissue Density",
    color = "month",
    fill = "month"
  ) +
  scale_color_manual(values = c("#3CB22D", "#FF8000")) +
  scale_fill_manual(values = c("#3CB22D", "#FF8000")) +
  theme_bw() 


P2
ggsave(filename = "images/RTDran.png", plot = P2, width = 18, height = 12, units = "cm")

## water levels effects ------

RTDran2_pred2 <-  ggpredict(RTDran2,
                           terms = c("GW_level_cm [all]", "year")) %>% 
  as.data.frame() %>% 
  # ggpredict returns a data frame with columns:
  # x (GW_level_cm), predicted (on original scale), conf.low, conf.high, group (month), facet (year)
  rename(year=group)

# Plot
ggplot(RTDran2_pred2, aes(x = x, y = predicted, color = year, fill = year)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.1, color = NA) +
  geom_line(size = 0.8) +  
  geom_point(
    data = Lys_data,
    aes(x = GW_level_cm, y = RTD_g_cm3, color=year, pch=year),
    size = 2
  ) +
  # facet_wrap(~year) +
  labs(
    x = "Groundwater Level [cm]",
    y = "Root Tissue Density") +
 # scale_color_manual(values = c("#3CB22D", "#FF8000")) +
  #  scale_fill_manual(values = c("#3CB22D", "#FF8000")) +
  theme_bw() 

#------------------------------------------------------------------------#
## month effects -----

# Extract p-value for June vs October
p_value <- car::Anova(RTDran2) %>% 
  as_tibble(rownames = "predictor") %>% 
  filter(predictor=="month") %>% 
  mutate(sign = ifelse( `Pr(>F)` < 0.001, "***", 
                        ifelse(`Pr(>F)` < 0.01, "**", 
                               ifelse( `Pr(>F)`< 0.05, "*", "ns"))))


ggplot(Lys_data, aes(x = month, y = RTD_g_cm3, fill = month)) +
  geom_violin(alpha = 0.7) +
  geom_signif(
    comparisons = list(c("June", "October")),
    y_position = c(max(Lys_data$RTD_g_cm3, na.rm = TRUE) * 1.05), # Define y position for the significance bar (adjust as needed)
    annotations = c(p_value$sign)
  ) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  scale_fill_manual(values = c("#3CB22D", "#FF8000")) +
  labs(x = "Seasonality", y = "Root Tissue Density") +
  theme_minimal()+
  ylim(0.1,0.32)


## year effect -----------
# Extract p-value for June vs October
p_value2 <- car::Anova(RTDran2) %>% 
  as_tibble(rownames = "predictor") %>% 
  filter(predictor=="year") %>% 
  mutate(sign = ifelse( `Pr(>F)` < 0.001, "***", 
                        ifelse(`Pr(>F)` < 0.01, "**", 
                               ifelse( `Pr(>F)`< 0.05, "*", "ns"))))


ggplot(Lys_data, aes(x = year, y = RTD_g_cm3, fill = year)) +
  geom_violin(alpha = 0.7) +
  geom_signif(
    comparisons = list(c("2023", "2024")),
    y_position = c(max(Lys_data$RTD_g_cm3, na.rm = TRUE) * 1.05), # Define y position for the significance bar (adjust as needed)
    annotations = c(p_value2$sign)
  ) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
 # scale_fill_manual(values = c("#3CB22D", "#FF8000")) +
  labs(x = "Year", y = "Root Tissue Density") +
  theme_minimal()+
  ylim(0.1,0.32)
# -----------------------------------------------------------------------------#

# (3) Average root diameter------


Lys_data %>% 
  ggplot(aes(GW_level_cm, RD_mm))+
  geom_point()


RDran1 <- lm(RD_mm ~ year * GW_level_cm +  GW_level_cm*month, 
             data = Lys_data)

# residuals:
par(mfrow=c(2,2))
plot(RDran1)
par(mfrow=c(1,1))

car::Anova(RDran1)

# remove nonsignificant interaction 
RDran2 <- lm(RD_mm ~ year + GW_level_cm + month,  
             data = Lys_data)

car::Anova(RDran2)


# Create the plot

library(ggeffects) # for ggpredict()
RDran2_pred <-  ggpredict(RDran2,
                             terms = c("GW_level_cm [all]", "year", "month")) %>% 
  as.data.frame() %>% 
  # ggpredict returns a data frame with columns:
  # x (GW_level_cm), predicted (on original scale), conf.low, conf.high, group (month), facet (year)
  rename(month=facet, year=group)

# Plot
## faceted plot ----

P3 <- ggplot(RDran2_pred, aes(x = x, y = predicted, color = month, fill = month)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.1, color = NA) +
  geom_line(size = 0.8) +  
  geom_point(
    data = Lys_data,
    aes(x = GW_level_cm, y = RD_mm, color = month),
    size = 2
  ) +
  facet_wrap(~year) +
  labs(
    x = "Groundwater Level [cm]",
    y = "Average Root Diameter [mm]",
    color = "month",
    fill = "month"
  ) +
  scale_color_manual(values = c("#3CB22D", "#FF8000")) +
  scale_fill_manual(values = c("#3CB22D", "#FF8000")) +
  theme_bw() 

P3


ggsave(filename = "images/RDran.png", plot = P3, width = 18, height = 12, units = "cm")

#---------------------------------------------------------------------------#
## water levels effects ------

RDran2_pred2 <-  ggpredict(RDran2,
                           terms = c("GW_level_cm [all]", "year")) %>% 
  as.data.frame()  %>% 
# ggpredict returns a data frame with columns:
# x (GW_level_cm), predicted (on original scale), conf.low, conf.high, group (month), facet (year)
 rename(year=group)

# Plot
ggplot(RDran2_pred2, aes(x = x, y = predicted, fill=year, color=year)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.1, color = NA) +
  geom_line(data=RDran2_pred2, size = 0.8) +  
#  geom_ribbon(data=RDran2_pred3, aes(ymin = conf.low, ymax = conf.high,
 #                                    fill=year), 
  #            alpha = 0.1, color = NA) +
  geom_line(data=RDran2_pred3, size = 0.8, aes(color=year)) +  
  geom_point(
    data = Lys_data,
    aes(x = GW_level_cm, y = RD_mm, color = year, pch=year),
    size = 2
  ) +
 # facet_wrap(~year) +
  labs(
    x = "Groundwater Level [cm]",
    y = "Average Root Diameter [mm]"
  ) +
#  scale_color_manual(values = c("#3CB22D", "#FF8000")) +
 # scale_fill_manual(values = c("#3CB22D", "#FF8000")) +
  theme_bw() 


#------------------------------------------------------------------------#
## month effects -----

# Extract p-value for June vs October
p_value <- car::Anova(RDran2) %>% 
  as_tibble(rownames = "predictor") %>% 
  filter(predictor=="month") %>% 
  mutate(sign = ifelse( `Pr(>F)` < 0.001, "***", 
                        ifelse(`Pr(>F)` < 0.01, "**", 
                               ifelse( `Pr(>F)`< 0.05, "*", "ns"))))


ggplot(Lys_data, aes(x = month, y = RD_mm, fill = month)) +
  geom_violin(alpha = 0.7) +
  geom_signif(
    comparisons = list(c("June", "October")),
    y_position = c(max(Lys_data$RD_mm, na.rm = TRUE) * 1.01), # Define y position for the significance bar (adjust as needed)
    annotations = c(p_value$sign)
  ) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  scale_fill_manual(values = c("#3CB22D", "#FF8000")) +
  labs(x = "Seasonality", y = "Average Root Diameter [mm]") +
  theme_minimal()


## year effects -------

# Extract p-value for June vs October
p_value2 <- car::Anova(RDran2) %>% 
  as_tibble(rownames = "predictor") %>% 
  filter(predictor=="year") %>% 
  mutate(sign = ifelse( `Pr(>F)` < 0.001, "***", 
                        ifelse(`Pr(>F)` < 0.01, "**", 
                               ifelse( `Pr(>F)`< 0.05, "*", "ns"))))


ggplot(Lys_data, aes(x = year, y = RD_mm, fill = year)) +
  geom_violin(alpha = 0.7) +
  geom_signif(
    comparisons = list(c("2023", "2024")),
    y_position = c(max(Lys_data$RD_mm, na.rm = TRUE) * 1.01), # Define y position for the significance bar (adjust as needed)
    annotations = c(p_value2$sign)
  ) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
 # scale_fill_manual(values = c("#3CB22D", "#FF8000")) +
  labs(x = "Year", y = "Average Root Diameter [mm]") +
  theme_minimal()+
  ylim(0.13, 0.22
       )
# -----------------------------------------------------------------------------#

# (4) hyphae------

Hyphran1 <- glm(hyphae ~  year + GW_level_cm + month + GW_level_cm:month, 
            data = Lys_data,
            family = binomial, 
            weights = no_intersects) # total number of intersects

# check how the model look
summary(Hyphran1)
str(Lys_data)
# check overdispersion:
# Residual deviance/degrees of freedom
1036.44/35 
check_overdispersion(Hyphran1)
# dispersion ratio should be less than 1.5, othervise you have overdispersion
# we have strong overdispersion, thus we have to use quasibinomial family

# use quasibinomial
Hyphran2 <- glm(hyphae ~  year * GW_level_cm + GW_level_cm*month, 
                data = Lys_data,
                family = quasibinomial, 
                weights = no_intersects) # total number of intersects


summary(Hyphran2)

# see the effects 
car::Anova(Hyphran2)
drop1(Hyphran2, test = "Chisq")

Hyphran3 <- glm(hyphae ~  year + GW_level_cm + month, 
                data = Lys_data,
                family = quasibinomial, 
                weights = no_intersects) # total number of intersects


car::Anova(Hyphran3)


# Plots:

## faceted plot ---- 
library(ggeffects) # for ggpredict()
Hyphran3_pred <-  ggpredict(Hyphran3,
                           terms = c("GW_level_cm [all]", "year", "month")) %>% 
  as.data.frame() %>% 
  # ggpredict returns a data frame with columns:
  # x (GW_level_cm), predicted (on original scale), conf.low, conf.high, group (month), facet (year)
  rename(month=facet, year=group)

# Plot
ggplot(Hyphran3_pred, aes(x = x, y = predicted, color = month, fill = month)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.1, color = NA) +
  geom_line(size = 0.8) +  
  geom_point(
    data = Lys_data,
    aes(x = GW_level_cm, y = hyphae, color = month),
    size = 2
  ) +
  facet_wrap(~year) +
  labs(
    x = "Groundwater Level [cm]",
    y = "Roots colonized with AMF [%]",
    color = "month",
    fill = "month"
  ) +
  scale_color_manual(values = c("#3CB22D", "#FF8000")) +
  scale_fill_manual(values = c("#3CB22D", "#FF8000")) +
  theme_bw() 


#---------------------------------------------------------------------------#
## water levels effects ------

Hyphran3_pred2 <-  ggpredict(Hyphran3,
                            terms = c("GW_level_cm [all]", "year")) %>% 
  as.data.frame() %>% 
  # ggpredict returns a data frame with columns:
  # x (GW_level_cm), predicted (on original scale), conf.low, conf.high, group (month), facet (year)
  rename(year=group)

# Plot
ggplot(Hyphran3_pred2, aes(x = x, y = predicted, color = year, fill = year)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.1, color = NA) +
  geom_line(size = 0.8) +  
  geom_point(
    data = Lys_data,
    aes(x = GW_level_cm, y = hyphae, color = year),
    size = 2
  ) +
 # facet_wrap(~year) +
  labs(
    x = "Groundwater Level [cm]",
    y = "Roots colonized with AMF [%]",
    color = "month",
    fill = "month"
  ) +
#  scale_color_manual(values = c("#3CB22D", "#FF8000")) +
 # scale_fill_manual(values = c("#3CB22D", "#FF8000")) +
  theme_bw() 

#-----------------------------------------------------------------------#
## month effcet -----

# Extract p-value for June vs October
p_value2 <- car::Anova(Hyphran3) %>% 
  as_tibble(rownames = "predictor") %>% 
  filter(predictor=="year") %>% 
  mutate(sign = ifelse( `Pr(>Chisq)` < 0.001, "***", 
                        ifelse(`Pr(>Chisq)` < 0.01, "**", 
                               ifelse( `Pr(>Chisq)`< 0.05, "*", "ns"))))


ggplot(Lys_data, aes(x = year, y = hyphae, fill = year)) +
  geom_violin(alpha = 0.7) +
  geom_signif(
    comparisons = list(c("2023", "2024")),
    y_position = c(max(Lys_data$hyphae, na.rm = TRUE) * 1.01), # Define y position for the significance bar (adjust as needed)
    annotations = c(p_value2$sign)
  ) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
 # scale_fill_manual(values = c("#3CB22D", "#FF8000")) +
  labs(x = "Year", y = "Roots colonized with AMF [%]") +
  theme_minimal()+
  ylim(0, 1)
