# Purpose: fit ordination for the root traits and their predictors

# packages ----
library(tidyverse)
library(vegan)

# read data 
Lys_data <- read_csv("data/RGG_Lys_Gradient.csv") %>% 
  unite("Plot_ID", c("Lys_no", "date"))

traits_2023 <- Lys_data %>% 
  filter(year==2023) %>% 
  dplyr::select(SRL_m_g, RD_mm, RTD_g_cm3, hyphae)

traits_2024 <- Lys_data %>% 
  filter(year==2024) %>% 
  dplyr::select(SRL_m_g, RD_mm, RTD_g_cm3, hyphae)

## Check correlation 
library(ggcorrplot)

round(cor(Lys_data %>% 
                     dplyr::select(SRL_m_g, RD_mm, RTD_g_cm3, hyphae) %>% 
                     rename("Specific root length" = SRL_m_g,
                            "Root diameter" = RD_mm,
                            "Root tissue density" = RTD_g_cm3,
                            "Roots colonized with AMF" = hyphae), 
                   use="pairwise.complete.obs", method = c("pearson")),2) %>% 
ggcorrplot(hc.order = F, type = "lower",
           lab = TRUE, lab_size = 4, 
           colors = c("red", "white", "blue"))

# for 2023
round(cor(traits_2023 %>% 
                     rename("Specific root length" = SRL_m_g,
                            "Root diameter" = RD_mm,
                            "Root tissue density" = RTD_g_cm3,
                            "Roots colonized with AMF" = hyphae), 
                   use="pairwise.complete.obs", method = c("pearson")),2) %>% 
ggcorrplot(hc.order = F, type = "lower",
           lab = TRUE, lab_size = 4, 
           colors = c("red", "white", "blue"))

# for 2024
round(cor(traits_2024 %>% 
            rename("Specific root length" = SRL_m_g,
                   "Root diameter" = RD_mm,
                   "Root tissue density" = RTD_g_cm3,
                   "Roots colonized with AMF" = hyphae), 
          use="pairwise.complete.obs", method = c("pearson")),2) %>% 
  ggcorrplot(hc.order = F, type = "lower",
             lab = TRUE, lab_size = 4, 
             colors = c("red", "white", "blue"))

names(Lys_data)


# check gradient length of first DCA axis (only optional)
decorana(traits) #  < 3 SD = linear methods are applicable

# PCA on on Hellinger-transformed data
pca_spec <- rda(traits ~ year + GW_level_cm + month, data = Lys_data,
                scale = TRUE) # scale data to have the same units
pca_spec

anova(pca_spec)
# Check eigenvalues of each PC

summary(eigenvals(pca_spec)) # cumulative proportion of PC1 + PC2 = 34.6%


plot(pca_spec, scaling = "species")


set.seed(1835)
fit <- envfit(pca_spec ~ year + GW_level_cm + month, data = Lys_data,
              choices = 1:2, # for which ordination axes
              scaling = "species", # symmetric scaling (as we want to see both sizes and species and their relationship)
              permutations = 1000) # can increase permutations if needed
fit

# extract only significant p-values
env.scrs <- scores(fit, display = "vectors") %>% 
  as_tibble(rownames="predictors") %>% 
  # remove "_" from predictor names
  mutate(predictors = gsub("_", " ", predictors))

plot(pca_spec, scaling = "species", display="species")
plot(fit)


### Plotting PCA results using the ggplot --------------------------------------

spec.scrs <- scores(pca_spec, display = "species",
                    scaling = "species") %>% 
  as_tibble(rownames = "Plot_Code") 
spec.scrs



# extract only significant p-values
env.scrs <- scores(pca_spec, display = "bp") %>% 
  as_tibble(rownames="predictors") #%>% 
  # remove "_" from predictor names
 # mutate(predictors = gsub("_", " ", predictors))

site.scrs <- scores(pca_spec, display = "sites",
                    scaling = "sites") %>% 
  as_tibble(rownames = "Plot_Code") %>% 
  left_join(envar, by="Plot_Code") # add site column from plant treatment data to PCA site scores, needed for coloring points by plant species richness

site.scrs

centroids <-scores(pca_spec, # from envfit() 
                   display="cn",  # centroid
                   scaling="species") %>%   # scaling methods
  as_tibble(rownames = "Centroids") #%>%  
 # mutate(Centroids=stringr::str_sub(Centroids, 13)) %>%  # remove letters before 3 letter
 # separate(Centroids, into = c("Plant_SR", "year")) %>% 
#  mutate(Plant_SR=as.numeric(Plant_SR))


ggplot(spec.scrs)+
  geom_hline(yintercept = 0, color="grey", lty =1) +
  geom_vline(xintercept = 0, color="grey", lty =1) +
  geom_point(aes(RDA1, RDA2), # shape=factor(plant_richness)), 
             size = 3,
             alpha=0.5)+
  # add centroid:
  geom_text_repel(data=centroids, 
                  aes(x=RDA1, y=RDA2), 
                  size=5, fontface="bold", show_guide = F) +
  #  geom_point(data=centroids, 
  #             aes(x=RDA1, y=RDA2, shape=factor(Plant_SR), color =year), 
  #             size=5, alpha=1, stroke=2) +
  
  scale_shape_manual(values = c(3,4,5,9, 7, 8))+ 
  labs(color="Year", shape="Plant species richness",
       x="RDA1 (24.7%)", y="RDA2 (9.8%)") +
  theme_bw()
