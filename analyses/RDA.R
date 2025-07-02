# Purpose: fit ordination for the root traits and their predictors

# packages ----
library(tidyverse)
library(vegan)
library(ggrepel)

# ---- Read and Prepare Data ----
Lys_data <- read_csv("data/RGG_Lys_Gradient.csv") %>%
  unite("Plot_ID", c("Lys_no", "date"), remove = TRUE)

# Filter root traits by year
get_traits_by_year <- function(data, yr) {
  data %>%
    filter(year == yr) %>%
    select(Plot_ID, GW_level_cm, month, SRL_m_g, RD_mm, RTD_g_cm3, hyphae)
}

predictors_2023 <- get_traits_by_year(Lys_data, 2023) %>% 
  select(Plot_ID, GW_level_cm, month) 

predictors_2024 <- get_traits_by_year(Lys_data, 2024) %>% 
  select(Plot_ID, GW_level_cm, month) 

traits_2023 <- get_traits_by_year(Lys_data, 2023) %>% 
  select(Plot_ID, SRL_m_g, RD_mm, RTD_g_cm3, hyphae)%>% 
  column_to_rownames("Plot_ID")

traits_2024 <- get_traits_by_year(Lys_data, 2024) %>% 
  select(Plot_ID, SRL_m_g, RD_mm, RTD_g_cm3, hyphae)%>% 
  column_to_rownames("Plot_ID")

# check gradient length of first DCA axis (only optional)
decorana(traits_2023) #  < 3 SD = linear methods are applicable
# Axis lengths         0.007819

# PCA on scaled data
ordin2023 <- rda(traits_2023 ~ GW_level_cm + month, data = predictors_2023,
                scale = TRUE) # scale data to have the same units
ordin2023

vif.cca(ordin2023)

# Check eigenvalues of each PC
summary(eigenvals(ordin2023)) # cumulative proportion of PC1 + PC2 = 34.6%


# R2 Adjusted for ordination models 
RsquareAdj(ordin2023)

# Permutation tests ------
set.seed(1)

perm_stat <- permustats(anova(ordin2023))
summary(perm_stat)
densityplot(perm_stat)

anova(ordin2023)
anova(ordin2023, by= "margin") # test for conditional effects (does not depend on the order)
anova(ordin2023, by= "axis")  # significance of individual constrained axes


# extracts goodness of fit for individual sites or species 
goodness(ordin2023, display = "species", 
         summarize = F) # summarize = TRUE show only the accumulated total across all axes

goodness(ordin2023, display = "species", 
         summarize = TRUE) # show only the accumulated total across all axes



plot(ordin2023, scaling = "species")


# extract scores  --------------------------------------



# vector for groundwater level effects
env.scrs <- scores(ordin2023, display = "bp") %>% 
  as_tibble(rownames="predictors") %>% 
  filter(predictors=="GW_level_cm") %>% 
  mutate(predictors_names=c("Groundwater level"))

env.scrs

# centroids for months effects
centroids <-scores(ordin2023, 
                   display="cn",  
                   scaling="species") %>%   
  as_tibble(rownames = "month") %>%  
  mutate(month=stringr::str_sub(month, 6)) 

centroids

# extract scores for each trait
trait.scrs <- scores(ordin2023, display = "species",
                     scaling = "species") %>% 
  as_tibble(rownames = "traits") %>% 
  mutate(traits_names=case_when(traits=="SRL_m_g" ~ "Specific root length",
                                traits== "RD_mm" ~ "Root diameter",
                                traits== "RTD_g_cm3" ~ "Root tissue density",
                                traits== "hyphae" ~"Roots AMF"))

trait.scrs

# extract scores for each observation
observ.scrs <- scores(ordin2023, display = "sites",
                      scaling = "species") %>% 
  as_tibble(rownames="Plot_ID") %>% 
  left_join(predictors_2023, by="Plot_ID") %>% 
  left_join(centroids %>% 
              rename(RDA1_centr=RDA1,
                     RDA2_centr=RDA2), by="month")

observ.scrs

# plot  --------------------------------------

observ.scrs %>% 
  ggplot()+
  geom_hline(yintercept = 0, color="grey", lty =1) +
  geom_vline(xintercept = 0, color="grey", lty =1) +
  # centroids:
  # stat_ellipse(aes(RDA1, RDA2, fill=month), 
  #              alpha=.1,type='t', linewidth =1, geom="polygon")+
  geom_segment(data=observ.scrs, 
              mapping = aes(x = RDA1, y = RDA2, 
                           xend = RDA1_centr, yend = RDA2_centr,
                             color=month), alpha=.4, linewidth=0.2,
               linetype=5) + 
    geom_point(data=centroids, aes(x=RDA1, y=RDA2, color=month), 
                 size=4, alpha=.5, stroke=2) +
    # observations:
  geom_point(aes(RDA1, RDA2, color=month), # shape=factor(plant_richness)), 
             size = 3, alpha=.5)+
  # traits:
  geom_text(data=trait.scrs, aes(RDA1, RDA2,
                                 label=traits_names),
            color="black",
            vjust=c(-1.1, 1.2, -0.6, 0.2),
            hjust=c(0.8, 0, 0, -0.1)) +
   geom_segment(data=trait.scrs, 
               aes(x=0, y=0, xend=RDA1, yend=RDA2), 
               arrow=arrow(length=unit(0.2,"cm")), 
               color="black", linewidth=0.9) +
  # water level:
  geom_text(data=env.scrs, aes(RDA1, RDA2,
                                 label=predictors_names),
            color="blue",size=5,
            vjust=c(-0.5),
            hjust=c(0.7)) +
  geom_segment(data=env.scrs, 
               aes(x=0, y=0, xend=RDA1, yend=RDA2), 
               arrow=arrow(length=unit(0.3,"cm")), 
               color="blue", linewidth=1) +
  labs(color="month", 
       x="RDA1 (32.7%)", y="RDA2 (4.8%)") +
  theme_bw()+
  scale_color_manual(values = c("#3CB22D", "#FF8000")) +
  scale_fill_manual(values = c("#3CB22D", "#FF8000")) 



# only traits and predictors
observ.scrs %>% 
  ggplot()+
  geom_hline(yintercept = 0, color="grey", lty =1) +
  geom_vline(xintercept = 0, color="grey", lty =1) +
  # centroids:
  # stat_ellipse(aes(RDA1, RDA2, fill=month), 
  #              alpha=.1,type='t', linewidth =1, geom="polygon")+
  #  geom_segment(data=observ.scrs, 
  #               mapping = aes(x = RDA1, y = RDA2, 
  #                             xend = RDA1_centr, yend = RDA2_centr,
  #                             color=month), alpha=.4, linewidth=0.2,
  #               linetype=5) + 
  geom_point(data=centroids, aes(x=RDA1, y=RDA2, color=month), 
             size=4, alpha=.5, stroke=2) +
  geom_text(data=centroids, aes(RDA1, RDA2,
                               label=month, color=month),
            size=5,show_guide = F,
            vjust=c(-1, -1),
            hjust=c(0, 1)) +
  # observations:
  #  geom_point(aes(RDA1, RDA2, color=month), # shape=factor(plant_richness)), 
  #             size = 3, alpha=.5)+
  # traits:
  geom_text(data=trait.scrs, aes(RDA1, RDA2,
                                 label=traits_names),
            color="black",
            vjust=c(-1.1, 1.2, -0.6, 0.2),
            hjust=c(0.5, 0, 0.5, -0.1)) +
  geom_segment(data=trait.scrs, 
               aes(x=0, y=0, xend=RDA1, yend=RDA2), 
               arrow=arrow(length=unit(0.2,"cm")), 
               color="black", linewidth=0.9) +
  # water level:
  geom_text(data=env.scrs, aes(RDA1, RDA2,
                               label=predictors_names),
            color="blue",size=5,
            vjust=c(-0.5),
            hjust=c(0.7)) +
  geom_segment(data=env.scrs, 
               aes(x=0, y=0, xend=RDA1, yend=RDA2), 
               arrow=arrow(length=unit(0.3,"cm")), 
               color="blue", linewidth=1) +
  labs(color="month", 
       x="RDA1 (32.7%)", y="RDA2 (4.8%)") +
  theme_bw()+
  scale_color_manual(values = c("#3CB22D", "#FF8000")) +
  scale_fill_manual(values = c("#3CB22D", "#FF8000")) +
  ylim(-0.5, 0.9)+
  xlim(-1.2, 1.2)
