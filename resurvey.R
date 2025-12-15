rm(list=ls())     #remove list: cleans workspace
gc()              #garbage collector

setwd("C:/Users/chris/Documents/Projekte/PhenObs/sensitivity_trends")

library(readxl)
library(dplyr)
library(reshape2)    # for casting matrices
library(data.table)  # for aggregating data
library(ggplot2)     # for graphics
library(DescTools)   # Gini coefficient and Lorenz curve
library(BSDA)        # sign test
library(Hmisc)       # error bars in Fig. 4
library(vegan)       # Shannon diversity  
library(matrixStats) # to calculate rowRanks
library(lmerTest)    # for mixed effects models (lmer)
library(parameters)  # for credible intervals for lmer

slope <- read.csv("Overview_Slopes_German.csv", sep = ",")

cultivated_or_nonDE <- c(
  "Amorpha canescens",
  "Arnica longifolia",
  "Bergenia purpurascens",
  "Camassia cusickii",
  "Clematis integrifolia",
  "Crocus speciosus",
  "Cyclamen coum",
  "Hyacinthoides hispanica",
  "Lysichiton americanus",
  "Lysichiton camtschatcensis",
  "Paeonia delavayi",
  "Paeonia peregrina",
  "Penstemon barbatus",
  "Phlomis russeliana",
  "Platycodon grandiflorus",
  "Podophyllum peltatum",
  "Silphium integrifolium",
  "Solidago rigida",
  "Trillium sessile",
  "Tulipa saxatilis",
  "Anemone hupehensis",
  "Gunnera manicata",
  "Helleborus orientalis",
  "Hypericum olympicum",
  "Hibiscus moschatus",
  "Lavandula angustifolia",
  "Iberis sempervirens",
  "Fuchsia magellanica",
  "Asphodelus albus",
  "Primula denticulata",
  "Veratrum nigrum",
  "Althaea officinalis",
  "Brunnera macrophylla",
  "Centranthus ruber",
  "Corydalis solida",
  "Dictamnus albus",
  "Eranthis hyemalis",
  "Galanthus nivalis",
  "Salvia nemorosa",
  "Salvia officinalis",
  "Galega officinalis",
  "Glycyrrhiza glabra",
  "Helleborus niger",
  "Hemerocallis fulva",
  "Leonurus cardiaca",
  "Leucojum aestivum",
  "Narcissus pseudonarcissus",
  "Paeonia officinalis",
  "Scopolia carniolica",
  "Tulipa sylvestris"
)


slope_sp <- slope_sp %>%
  filter(!Species %in% cultivated_or_nonDE)

slope_sp <- slope_sp %>%
  mutate( Species = case_when(
      Species == "Atropa belladonna"        ~ "Atropa bella-donna",
      Species == "Stachys officinalis"      ~ "Betonica officinalis",
      Species == "Ficaria verna"            ~ "Ranunculus ficaria",
      Species == "Comarum palustre"         ~ "Potentilla palustris",
      Species == "Lotus maritimus"          ~ "Tetragonolobus maritimus",
      Species == "Anemone hepatica"         ~ "Hepatica nobilis",
      Species == "Arum italicum"            ~ "Arum maculatum agg.",
      Species == "Alchemilla xanthochlora"  ~ "Alchemilla vulgaris agg.",
      Species == "Primula acaulis"          ~ "Primula vulgaris",
      Species == "Viscaria vulgaris"        ~ "Silene viscaria",
      TRUE                                          ~ Species))




str(slope)

slope_sp <- slope_sp %>%
  select(Species, mean_slope)
slope_sp$Species <- as.factor(slope_sp$Species)

str(slope_sp)

###preparing to join
library(stringr)

clean_species <- function(x) {
  x %>%
    tolower() %>%                       # alles klein
    str_replace_all("[,;]", " ") %>%    # Komma, Semikolon entfernen
    str_replace_all("\\bagg\\.?\\b", "") %>%  # "agg." oder "agg" entfernen
    str_replace_all("\\bsp\\.?\\b", "") %>%   # optional: "sp." entfernen
    str_replace_all("\\bspp\\.?\\b", "") %>%  # optional: "spp." entfernen
    str_replace_all("[^a-z0-9äöüß \\-]", " ") %>%  # Sonderzeichen raus
    str_squish()                        # Mehrfach-Leerzeichen reduzieren
}

slope_sp <- slope_sp %>%
  mutate(species_clean = clean_species(Species))


change<-read_excel('resurvey_species_change.xlsx')%>%
  rename(Species = species) %>%
  mutate(species_clean = clean_species(Species))

summary(change)


# Calculate the probability to increase in cover for the 
# 161 species with a significantly negative or positive change
# (according to a binomial test at p<0.05, with Holm correction) 
# and at least 100 change observations.
# Holm adjustment of p values
#change$p.adjust <- p.adjust(change$p.values.binom, method="holm")
str(change)
change2b <- change %>%
  filter(p.adjust < 0.05, n >= 100)
str(change2b) #161 obs. of  11 variables:
change2b <- change2b[order(change2b$est.binom),]
dim(change2b) # 161 11 Fig. 4 and Fig. 6

# merging into one
combo <- merge(slope_sp, change2b, by = "species_clean")
str(combo)



# linear model
m1 <-lm(combo$mean.absolute.change ~ combo$mean_slope)
summary(m1) 
ggplot(combo, aes(x=mean_slope, y=mean.absolute.change)) +  
  geom_point(alpha = 0.4,  
             size = 2,
             shape = 19) +
  geom_smooth(method = "lm",  
              linetype = "dashed",
              linewidth = 1.2) 


######Alternative: probbaility with all species (significant or not in change)

combo_alternative <- merge(slope_sp, change, by = "species_clean")
str(combo_alternative)

hist(combo_alternative$n,50)
combo_alternative2<-combo_alternative%>%
  filter( n >= 10)
# linear model
m1 <-lm(combo_alternative2$est.binom  ~ combo_alternative2$mean_slope)
summary(m1) 
  ggplot(combo_alternative2, aes(x=mean_slope, y=est.binom)) +  
  geom_point(alpha = 0.4,  
             size = 1.7,
             shape = 19) +
    geom_vline(xintercept = 0, color = "gray", size = 0.4) +
    geom_hline(yintercept = 0.5, color = "gray", size = 0.4) +
  geom_smooth( method = "lm", 
               linetype = "dashed",
               linewidth = 1.2) +
  labs(x="Phenological sensistivity",y="Probability of increasing in cover")+
  theme_bw()
  ggsave("resurvey2.pdf", width = 18, height = 10, units = "cm")
  


###new tabel CHange list

# from:  group_by(Species) %>%
#summarise(relative.change = mean(relative.change, na.rm = TRUE), slope = mean(slope, na.rm = TRUE),.groups = 'drop')

change_list<-read_excel('resurvey_species_change_list.xlsx')%>%
  mutate(species_clean = clean_species(Species))%>%
  filter( n_records >= 10)



# merging into one
combo2 <- merge(slope_sp, change_list, by ="species_clean")
slope_only_species <- slope_sp %>%
  anti_join(change_list, by = "species_clean")


# linear model of relative change
m2 <- lm(relative.change ~ mean_slope, data = combo2)
summary(m2)
#plot(m2)
ggplot(combo2, aes(x=mean_slope, y=relative.change)) +  
  geom_point(alpha = 0.4,  
             size = 1.7,
             shape = 19) +
  geom_vline(xintercept = 0, color = "gray", size = 0.4) +
  geom_hline(yintercept = 0, color = "gray", size = 0.4) +
  geom_smooth( method = "lm", 
              linetype = "dashed",
              linewidth = 1.2) +
  labs(x="Phenological sensistivity",y="Mean relative change across resurveys")+
  theme_bw()

ggsave("resurvey.pdf", width = 18, height = 10, units = "cm")


#
setDT(species.change.list2)

# getting a mean slope of species.chang.list2
species_summary <- species.change.list2 %>%
  filter(!is.na(slope)) %>%
  group_by(Species) %>%
  summarise(
    mean_slope_change = mean(slope, na.rm = TRUE),
    frequency = n(),  
    .groups = 'drop' 
  ) %>%
  filter(frequency > 0) %>%
  mutate(species_clean = clean_species(Species))

# merging
comb_3 <- merge(slope_sp, species_summary, by = "species_clean")


# linear model
m3 <- lm(mean_slope_change ~ mean_slope, data = comb_3)
summary(m3)
ggplot(comb_3, aes(x=mean_slope, y=mean_slope_change)) +  
  geom_point(alpha = 0.4,  
             size = 1.7,
             shape = 19) +
  geom_vline(xintercept = 0, color = "gray", size = 0.4) +
  geom_hline(yintercept = 0, color = "gray", size = 0.4) +
  geom_smooth( method = "lm", 
               linetype = "dashed",
               linewidth = 1.2) +
  labs(x="Phenological sensistivity",y="Mean slope across resurveys")+
  theme_bw()

ggsave("resurvey_slopes.pdf", width = 18, height = 10, units = "cm")

# for relative rank change
species_summary_rel <- species.change.list2 %>%
  filter(!is.na(relative.rank.change)) %>%
  group_by(species) %>%
  summarise(
    mean_rel_rank = mean(relative.rank.change, na.rm = TRUE),
    frequency = n(),  
    .groups = 'drop'  
  ) %>%
  filter(frequency > 0) %>%
  rename(Species = species) %>%
  mutate(species_clean = clean_species(Species))

# merging
comb_3_rel <- merge(slope_sp, species_summary_rel, by = "species_clean")

# linear model
m3_rel <- lm(mean_rel_rank ~ mean_slope, data = comb_3_rel)
summary(m3_rel)
ggplot(comb_3_rel, aes(x=mean_slope, y=mean_rel_rank)) +  
  geom_point(alpha = 0.4,  
             size = 1.7,
             shape = 19) +
  geom_vline(xintercept = 0, color = "gray", size = 0.4) +
  geom_hline(yintercept = 0, color = "gray", size = 0.4) +
  geom_smooth( method = "lm", 
               linetype = "dashed",
               linewidth = 1.2) +
  labs(x="Phenological sensistivity",y="Mean relative ranks across resurveys")+
  theme_bw()

ggsave("resurvey_slopes_ranks.pdf", width = 18, height = 10, units = "cm")

# absolute change
species_summary_ab <- species.change.list2 %>%
  filter(!is.na(absolute.change)) %>%
  group_by(species) %>%
  summarise(
    mean_ab_chg = mean(absolute.change, na.rm = TRUE),
    frequency = n(),  
    .groups = 'drop'  
  ) %>%
  filter(frequency > 0) %>%
  rename(Species = species)%>%
  mutate(species_clean = clean_species(Species))

# merging
comb_3_ab <- merge(slope_sp, species_summary_ab, by = "species_clean")

# linear model
m3_ab <- lm(mean_ab_chg ~ mean_slope, data = comb_3_ab)
summary(m3_ab)
