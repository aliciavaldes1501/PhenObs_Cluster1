####Load libraries ####
library(readr)
library(here)
library(rWCVP)
library(tidyverse)
library(purrr)
library(sf)
library(glmmTMB)
library(ggeffects)

# Try rWCVP or GIFT R packages, which allow to access data from the World 
# Checklist of Vascular Plants and the Global Inventory of Floras and Traits, 
# respectively.
# Trt also the rgbif package which allows retrieving data from GBIF. 

#### Try package rWCVP ####

species <- read_csv(here("data", "raw", "list_species_PhenObs.csv"))
species

# taxonomy data
names <- rWCVPdata::wcvp_names

# distribution data
distributions <- rWCVPdata::wcvp_distributions

# Filtering the WCVP to generate lists or summaries of vascular plant species 
# in particular areas. 
# Two arguments:
# taxon: the name of a valid taxon with a taxonomic rank of species or higher 
# (e.g. the species “Myrcia almasensis”, the genus “Myrcia”, or the family “Myrtaceae”).
# area: a vector of WGSRPD level 3 codes for the regions you want to focus on.
# These arguments can be combined in the wcvp_checklist, wcvp_occ_mat, 
# and wcvp_summary to generate outputs for focal taxa in a desired area. 

matches <- wcvp_match_names(species, name_col = "species", fuzzy = TRUE,
                            progress_bar = FALSE)

# Generate occurency matrix

# Extract the genus
genera <- species %>%
  mutate(genus = word(species, 1)) %>%
  distinct(genus) %>%
  pull(genus)

# Get country codes for all countries in PhenObs
get_wgsrpd3_codes("Canada") # More than one
get_wgsrpd3_codes("Germany")
get_wgsrpd3_codes("Austria")
get_wgsrpd3_codes("Switzerland")
get_wgsrpd3_codes("Czechoslovakia")
get_wgsrpd3_codes("Russia") # More than one
get_wgsrpd3_codes("Norway") # More than one
get_wgsrpd3_codes("India") # More than one
get_wgsrpd3_codes("Spain") # "SPA"
get_wgsrpd3_codes("Iran")

# Points for gardens where I have doubts about the code
points_df <- data.frame(
  lon = c(-64.3683299719477, 34.38191913554825, 10.45239702332544, 74.87900495304142),   
  lat = c(45.08747478246653, 61.84258691264972, 63.448299228918934, 34.09344415438499),
  name = c("Acadia", "Petrozavodsk", "Trondheim", "Srinagar")
)

points_sf <- st_as_sf(points_df, coords = c("lon", "lat"), crs = 4326)

# Map
ggplot() +
  geom_sf(data = wgsrpd3 %>%
            filter(LEVEL3_COD %in% 
                     (wgsrpd_mapping %>%
                        filter(COUNTRY %in% c("Canada", "Germany", "Austria", 
                                              "Switzerland", "Czechoslovakia",
                                              "Russia", "Norway", "India",
                                              "Spain", "Iran")) %>% 
                        distinct(LEVEL3_COD) %>% pull(LEVEL3_COD))), 
          fill = "lightgray", color = "white") +
  geom_sf_text(data = wgsrpd3 %>%
                 filter(LEVEL3_COD %in% 
                          (wgsrpd_mapping %>%
                             filter(COUNTRY %in% c("Canada", "Germany", "Austria", 
                                                   "Switzerland", "Czechoslovakia",
                                                   "Russia", "Norway", "India",
                                                   "Spain", "Iran")) %>% 
                             distinct(LEVEL3_COD) %>% pull(LEVEL3_COD))), 
               aes(label = LEVEL3_COD), size = 3, color = "black") +
  geom_sf(data = points_sf, color = "red", size = 3) +
  theme_minimal()

# Map of all gardens and countries

# Apply wcvp_occ_mat to each genus
safe_occ_mat <- possibly(
  function(genus) {
    wcvp_occ_mat(taxon = genus, taxon_rank = "genus",
                 area_codes = c(get_wgsrpd3_codes("Canada"), # More than one,
                                get_wgsrpd3_codes("Germany"),
                                get_wgsrpd3_codes("Austria"),
                                get_wgsrpd3_codes("Switzerland"),
                                get_wgsrpd3_codes("Czechoslovakia"),
                                get_wgsrpd3_codes("Russia"), # More than one
                                get_wgsrpd3_codes("Norway"), # More than one
                                get_wgsrpd3_codes("India"), # More than one
                                get_wgsrpd3_codes("Spain"), # More than one
                                get_wgsrpd3_codes("Iran"))
                 ) %>%
      mutate(input_genus = genus)
  },
  otherwise = NULL
)

# Apply to all genera
occurrence_matrix_genera <- map_df(genera, safe_occ_mat)

occurrence_matrix_genera

# Filter for species that are on the PhenObs list
occurrence_matrix_species <- occurrence_matrix_genera %>%
  filter(taxon_name %in% species$species)
nrow(occurrence_matrix_species) # 161 species of 177 in PhenObs list

# Species from PhenObs list that were not matched
species %>% filter(!species %in% occurrence_matrix_species$taxon_name)

# Anemone hepatica = Hepatica nobilis
wcvp_occ_mat(taxon = "hupehensis", taxon_rank = "species",
             area_codes = c("GER", "AUT", "SWI", "CZE", "RUN", "NOR", "WHM",
                            "SPA", "IRN")) %>%
  filter(taxon_name == "Hepatica nobilis")

# Species not appearing in any of the area_codes from PhenObs
wcvp_distribution(taxon ="Camassia cusickii", taxon_rank = "species") %>%
  wcvp_distribution_map() + ggtitle("Camassia cusickii")
wcvp_distribution(taxon ="Fuchsia magellanica", taxon_rank = "species") %>%
  wcvp_distribution_map() + ggtitle("Fuchsia magellanica")
wcvp_distribution(taxon ="Gunnera manicata", taxon_rank = "species") %>%
  wcvp_distribution_map() + ggtitle("Gunnera manicata")
wcvp_distribution(taxon ="Paeonia delavayi", taxon_rank = "species") %>%
  wcvp_distribution_map() + ggtitle("Paeonia delavayi")
wcvp_distribution(taxon ="Paeonia peregrina", taxon_rank = "species") %>%
  wcvp_distribution_map() + ggtitle("Paeonia peregrina")
wcvp_distribution(taxon = "Trillium sessile", taxon_rank = "species") %>%
  wcvp_distribution_map() + ggtitle("Trillium sessile")
wcvp_distribution(taxon = "Tulipa saxatilis", taxon_rank = "species") %>%
  wcvp_distribution_map() + ggtitle("Tulipa saxatilis")

# Add synonyms
species <- species %>%
  mutate(synonym = case_when(
    species == "Anemone hepatica" ~ "Hepatica nobilis",
    species == "Anemone hupehensis" ~ "Eriocapitella hupehensis",
    species == "Atropa belladonna" ~ "Atropa bella-donna",
    species == "Centranthus ruber" ~ "Valeriana rubra",
    species == "Ficaria verna" ~ "Ranunculus ficaria",
    species == "Hibiscus moschatus" ~ "Abelmoschus moschatus",
    species == "Primula acaulis" ~ "Primula vulgaris",
    species == "Securigera varia" ~ "Coronilla varia",
    species == "Stachys officinalis" ~ "Betonica officinalis",
    TRUE ~ NA_character_
  ))

# Redo

genera <- species %>%
  mutate(genus = ifelse(
    !is.na(synonym), word(synonym, 1),
    word(species, 1)
    )
    ) %>%
  distinct(genus) %>%
  pull(genus)

# Apply to all genera
occurrence_matrix_genera <- map_df(genera, safe_occ_mat)

occurrence_matrix_genera

# Filter for species that are on the PhenObs list
species <- species %>%
  mutate(species_corrected = ifelse(
    !is.na(synonym), synonym,
    species
  ))
occurrence_matrix_species <- occurrence_matrix_genera %>%
  filter(taxon_name %in% species$species_corrected)
nrow(occurrence_matrix_species) # 170 species of 177 in PhenObs list

# Species from PhenObs list that were not matched
# According to WCVP, these species are not appearing in any of the countries
# (Canada, Germany, Austria, Switzerland, Czechia, Russia, Norway, India, Spain,
# Iran)
species_absent <- species %>% filter(!species_corrected %in% occurrence_matrix_species$taxon_name)
species_absent
species_present <-species %>% filter(species_corrected %in% occurrence_matrix_species$taxon_name)
species_present

# Occurence matrix for countries
occurrence_matrix_species <- occurrence_matrix_species %>%
  select(-input_genus) %>%
  pivot_longer(cols = 3:48, names_to = "code") %>% 
  mutate(country = case_when(
    code %in% c("ABT", "BRC", "LAB", "MAN", "NWT", "NBR", "NFL", "NSC", "NUN",
                "ONT", "QUE", "SAS", "YUK") ~ "Canada",
    code == "GER" ~ "Germany",
    code == "AUT" ~ "Austria",
    code == "SWI" ~ "Switzerland",
    code == "CZE" ~ "Czech Republic",
    code %in% c("ALT", "AMU", "RUC", "RUN", "BRY", "RUS", "RUW", "CTA", "IRK",
                "KAM", "KHA", "KRA", "KUR", "MAG", "NCS", "PRM", "RUE", "SAK",
                "TVA", "WSB", "YAK") ~ "Russia",
    code %in% c("NOR", "SVA") ~ "Norway",
    code %in% c("ASS", "IND", "WHM") ~ "India",
    code %in% c("SPA", "BAL") ~ "Spain",
    code == "IRN" ~ "Iran",
    TRUE ~ NA_character_
  )) %>%
  group_by(taxon_name, country) %>%
  summarise(occurrence = max(value), .groups = "drop")

occurrence_matrix_species_wide <- occurrence_matrix_species %>%
  pivot_wider(names_from = country, values_from = occurrence)
occurrence_matrix_species_wide

#### Try package GIFT ####

#### Try package rgbif ####

### MOTIVATE Red List data ####

redlist <- read_tsv(here("data", "raw", "redlist_MOTIVATE_subset.csv"))
redlist

countries <- c("Canada", "Germany", "Austria", "Switzerland", "Czech Republic",
               "Russia", "Norway", "India", "Spain")

# Canada and India not in Red List data

redlist_countries <- redlist %>%
  # Filter only the countries for which we have phenological sensitivity
  filter(Country %in% countries) 

redlist_phenobs <-bind_rows(
  species_present %>%
    left_join(redlist_countries, 
              by = c("species" = "parsed_species_name")), 
  species_present %>%
    left_join(redlist_countries, 
              by = c("species_corrected" = "parsed_species_name"))
  ) %>%
  distinct() %>%
  # Keep columns that I think I will use
  select(species, species_corrected, Redlist_cat, Threat, Country, Publ_year, 
         redlist_source, Continent)
redlist_phenobs

# Species with no info in Red List data

redlist_phenobs %>% filter(is.na(Threat)) %>% distinct(species_corrected)
sp# We suppose these are not threatened in any country, but should double check

# We remove them from the red list
redlist_phenobs <- redlist_phenobs %>% filter(!is.na(Threat))

# Cases where there is info for more than one year for a species and country

redlist_phenobs_dupl_diff <- redlist_phenobs %>%
  group_by(species_corrected, Country) %>%
  filter(n() > 1, n_distinct(Threat) > 1) %>%
  ungroup()

redlist_phenobs_dupl_diff

# Keep threat info for the latest publication year

# In Russia, there is no info about publ year

# Keep Russia separated
redlist_phenobs_Russia <- redlist_phenobs %>%
  filter(Country == "Russia")

# Remove Russia
redlist_phenobs_noRussia <- redlist_phenobs %>%
  filter(Country != "Russia")

# For other countries apart from Russia,
# if there are several rows for same publ year, 
# pick the most severe red list category

severity_order <- c("CR", "EN", "VU", "NT", "LC")

redlist_phenobs_noRussia_latest <- redlist_phenobs_noRussia %>%
  group_by(species_corrected, Country) %>%
  filter(if (all(is.na(Publ_year))) TRUE else Publ_year == max(Publ_year, na.rm = TRUE)) %>%
  slice(which.min(match(Redlist_cat, severity_order))) %>%
  ungroup()

redlist_phenobs_noRussia_latest

# Check that there is only one entry for each species and country
redlist_phenobs_noRussia_latest %>% group_by(species_corrected, Country) %>%
  count() %>% filter(n > 1)

# For Russia

severity_threatened <- c("CR", "EN", "VU")
severity_not_threatened <- c("NT", "LC")

redlist_phenobs_Russia <- redlist_phenobs_Russia %>%
  group_by(species_corrected) %>%
  # Apply custom logic to pick ONE row but keep all columns
  slice({
    threatened_count <- sum(Threat == "threatened")
    not_threatened_count <- sum(Threat == "not_threatened")
    
    if (threatened_count > not_threatened_count) {
      which.max(Threat == "threatened") # pick first threatened row
    } else if (not_threatened_count > threatened_count) {
      which.max(Threat == "not_threatened") # pick first not_threatened row
    } else {
      # Tie: break by severity groups
      threatened_severity_count <- sum(Redlist_cat %in% severity_threatened)
      not_threatened_severity_count <- sum(Redlist_cat %in% severity_not_threatened)
      
      if (threatened_severity_count >= not_threatened_severity_count) {
        which.max(Redlist_cat %in% severity_threatened)
      } else {
        which.max(Redlist_cat %in% severity_not_threatened)
      }
    }
  }) %>%
  ungroup()

# Check that there is only one entry for each species
redlist_phenobs_Russia %>% group_by(species_corrected) %>%
  count() %>% filter(n > 1)

# Join noRussia and Russia
redlist_phenobs_final <- bind_rows(redlist_phenobs_noRussia_latest,
                                   redlist_phenobs_Russia)

# Note that species that are not threatened in any country 
# are absent from redlist_phenobs_final

#### Add occurrence data to Red List data ####

countries_redlist <- c(#"Canada", 
  "Germany", "Austria", "Switzerland", "Czech Republic", "Russia", "Norway", 
  #"India", 
  "Spain")

redlist_phenobs_final_occ <- 
  # Here, all 170 present species are included
  occurrence_matrix_species %>% filter(country %in% countries_redlist) %>%
  # Here, only the species included in red lists in ANY country are included
  left_join(redlist_phenobs_final, by = c("taxon_name" = "species_corrected",
                                          "country" = "Country")) %>%
  # Remove unneded columns
  select(-species, -Redlist_cat, -Publ_year, -redlist_source, - Continent) %>%
  # Rename
  rename(species = taxon_name) %>%
  # Order columns
  select(species, country, occurrence, Threat) %>%
  # Set those where Threat is NA or not_evaluated_dd as not_threatened
  mutate(Threat = case_when(
    Threat == "not_evaluated_dd" ~ "not_threatened",
    is.na(Threat) ~ "not_threatened",
    TRUE ~ Threat
  ))

#### Calculate % endangerment ####

perc_endangerment <- redlist_phenobs_final_occ %>%
  filter(occurrence > 0) %>%  # only countries where species is present
  group_by(species) %>%
  summarise(
    total_countries = n(),
    threatened_countries = sum(Threat == "threatened"),
    perc_endangerment = (threatened_countries / total_countries) * 100,
    .groups = "drop"
  ) %>%
  arrange(species)

perc_endangerment

#### Sensitivity data (slopes from Robert) ####
sensitivity_data <- read_csv(here("data", "raw", "Overview_Slopes_all.csv")) %>%
  select(-...1) %>%
  # Modify some species names to match those in red list data
  mutate(Species = case_when(
    Species == "Anemone hepatica" ~ "Hepatica nobilis",
    Species == "Anemone hupehensis" ~ "Eriocapitella hupehensis",
    Species == "Atropa belladonna" ~ "Atropa bella-donna",
    Species == "Centranthus ruber" ~ "Valeriana rubra",
    Species == "Ficaria verna" ~ "Ranunculus ficaria",
    Species == "Hibiscus moschatus" ~ "Abelmoschus moschatus",
    Species == "Primula acaulis" ~ "Primula vulgaris",
    Species == "Securigera varia" ~ "Coronilla varia",
    Species == "Stachys officinalis" ~ "Betonica officinalis",
    TRUE ~ Species
  )) %>%
  # Add variable for "sens" vs "not-sens"
  mutate(sens_01 = factor(ifelse(ci_lower > 0 | ci_upper < 0, 1, 0)))
sensitivity_data

#### Models % endangerment vs sensitivity ####

# Join data

sens_per_endanger <- perc_endangerment %>% 
  left_join(sensitivity_data %>% 
              select(Species, mean_slope, sens_01) %>% 
              rename(species = Species))

sens_per_endanger %>% filter(is.na(perc_endangerment))
sens_per_endanger %>% filter(is.na(mean_slope))

# Histograms
ggplot(sens_per_endanger, aes(x = perc_endangerment)) +
  geom_histogram(color = "black", fill = "lightblue") +
  theme_bw() +
  xlab("% endangerment") + ylab("N species")

ggplot(sens_per_endanger, aes(x = mean_slope)) +
  geom_histogram(color = "black", fill = "lightblue") +
  theme_bw() +
  xlab("Mean slope") + ylab("N species")

ggplot(sens_per_endanger, aes(sens_01, perc_endangerment)) +
  geom_boxplot()

sens_per_endanger %>% group_by(sens_01) %>%
  summarise(mean_perc = mean(perc_endangerment))

# Fit a binomial model

# With 0/1
sens_per_endanger <- sens_per_endanger %>% 
  mutate(endangered = ifelse(perc_endangerment == 0, 0, 1))
model_01_endangerment <- glmmTMB(endangered ~
    mean_slope, family = binomial, data = sens_per_endanger)
summary(model_01_endangerment)
plot(ggpredict(model_01_endangerment, terms="mean_slope [all]"))

# With proportions

# perc_endangerment is in percentage, so use counts:
# successes = threatened_countries
# failures = total_countries - threatened_countries

model_perc_endangerment <- glmmTMB(
  cbind(threatened_countries, total_countries - threatened_countries) ~
    mean_slope, family = binomial, data = sens_per_endanger %>% 
    filter(mean_slope <= 0  & mean_slope > -17))
model_perc_endangerment_quad <- glmmTMB(
  cbind(threatened_countries, total_countries - threatened_countries) ~
    mean_slope + I(mean_slope^2), family = binomial, 
  data = sens_per_endanger %>% filter(mean_slope <= 0 & mean_slope > -17))
# Does it make sense to fit a quadratic model?
summary(model_perc_endangerment)
summary(model_perc_endangerment_quad)

# Get predictions
pred <- ggpredict(model_perc_endangerment, terms = "mean_slope [all]")
pred_quad <- ggpredict(model_perc_endangerment_quad, terms = "mean_slope [all]")

# Plot predictions + raw data
ggplot() +
  geom_line(data = pred, aes(x = x, y = predicted),
            color = "blue", linewidth = 1) +
  geom_ribbon(data = pred, 
              aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high),
              fill = "blue", alpha = 0.1) +
  # geom_line(data = pred_quad, aes(x = x, y = predicted),
  #           color = "red", linewidth = 1) +
  # geom_ribbon(data = pred_quad, 
  #             aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high),
  #             fill = "red", alpha = 0.1) +
  geom_point(data = sens_per_endanger %>% filter(mean_slope <= 0),
             aes(x = mean_slope, y = threatened_countries / total_countries),
             color = "black", alpha = 0.6) +
  labs(x = "Mean slope", y = "Predicted probability of endangerment") +
  theme_minimal()


model_perc_endangerment_quad_subset <- glmmTMB(
  cbind(threatened_countries, total_countries - threatened_countries) ~
    mean_slope + I(mean_slope^2), family = binomial, 
  data = sens_per_endanger %>% filter(mean_slope> -12 & mean_slope < 5))
summary(model_perc_endangerment_quad_subset)
pred_quad_subset <- ggpredict(model_perc_endangerment_quad_subset,
                              terms = "mean_slope [all]")

ggplot() +
  geom_line(data = pred_quad_subset, aes(x = x, y = predicted),
            color = "red", linewidth = 1) +
  geom_ribbon(data = pred_quad_subset, 
              aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high),
              fill = "red", alpha = 0.1) +
  geom_point(data = sens_per_endanger  %>% 
               filter(mean_slope> -12 & mean_slope < 5),
             aes(x = mean_slope, y = threatened_countries / total_countries),
             color = "black", alpha = 0.6) +
  labs(x = "Mean slope", y = "Predicted probability of endangerment") +
  theme_minimal()

# Zero inflation
sens_per_endanger <- sens_per_endanger %>%
  mutate(prop_endangerment = perc_endangerment/100)

model_perc_endangerment_zi <- glmmTMB(prop_endangerment ~ mean_slope,
                                      family = "beta_family", ziformula=~., 
                                      data = sens_per_endanger %>% 
                                        filter(mean_slope <= 0 & mean_slope > -17))
model_perc_endangerment_quad_zi <- glmmTMB(prop_endangerment ~ mean_slope +
                                             I(mean_slope^2), 
                                           family = "beta_family", ziformula=~.,
                                           data = sens_per_endanger %>% 
                                             filter(mean_slope <= 0 & mean_slope > -17))
# Does it make sense to fit a quadratic model?
summary(model_perc_endangerment_zi)
summary(model_perc_endangerment_quad_zi)
