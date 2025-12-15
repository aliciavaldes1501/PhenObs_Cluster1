rm(list = ls())
gc()

setwd("C:/Users/chris/Documents/Projekte/PhenObs/sensitivity_trends")

library(readr)
library(taxize)
library(dplyr)
library(rgbif)
library(V.PhyloMaker)
library(ggtree)
library(ggtreeExtra)
library(ggstar)
library(phytools)

df = read_csv('Overview_Slopes_all.csv')

# Only using unique species
unique_species <- unique(df$Species)

# Due to gna_verifier only being able to handle 50 names at once subsetting
first <- unique_species[1:50]
second <- unique_species[51:100]
third <- unique_species[101:150]
fourth <- unique_species[151:177]

# Unconventional way for getting current canonical names
resolved_1 <- gna_verifier(first)
resolved_2 <- gna_verifier(second)
resolved_3 <- gna_verifier(third)
resolved_4 <- gna_verifier(fourth)

# Combining all of the names
all_resolved <- rbind(resolved_1, resolved_2, resolved_3, resolved_4)

# Add to df
df_tax <- df %>%
  cbind(., all_resolved$currentCanonicalFull) %>%
  rename(name_res = `all_resolved$currentCanonicalFull`)

# Getting backbone information
gbif_id <- get_gbifid_(df_tax$name_res)

backbone <- lapply(gbif_id, function(df) df[,c("genus", "family", "order" ) ])

# Bit cheated
backbone_combined <- bind_rows(lapply(backbone, function(df) df[1, ]))

# Fusing with df_tax
df_tax_bb <- df_tax %>%
  select(Species, mean_slope, name_res) %>%
  cbind(., backbone_combined) %>%
  rename(name_old = Species,
         Species = name_res,
         Genus = genus,
         Family = family,
         Order = order)

# Prep for Tree
df_tax_tree <- df_tax_bb %>%
  select(Species, Genus, Family, Order, mean_slope)


## Creating Trees

phyloResult = phylo.maker(df_tax_tree, scenarios = "S3",output.tree=T)

# Getting clean names again
phyloResult$scenario.3[["tip.label"]] <- gsub("_", " ", phyloResult$scenario.3[["tip.label"]])

#
phyloTree_sp = phyloResult$scenario.3
phyloTree_all = phyloResult$tree.scenario.3

# phyloTree_q = chronos(phyloTree_sp, lambda = 0, model = "relaxed")

dd = data.frame(species = phyloResult$species.list[["species"]],
                mean_slope = phyloResult$species.list[["mean_slope"]],
                family = phyloResult$species.list[["family"]])


# checking if all species match in scenario.3 and species.list
dd$species <- trimws(dd$species)
phyloResult$scenario.3[["tip.label"]] <- trimws(phyloResult$scenario.3[["tip.label"]])

# Create a clean matching dataset
matching_species <- intersect(dd$species, phyloResult$scenario.3[["tip.label"]])
cat("Number of matching species:", length(matching_species), "\n")

# Filter data to only include matching species
dd_filtered <- dd[dd$species %in% matching_species, ]
library(ggplot2)
# Creating tree visualization
p <- ggtree(phyloResult$scenario.3, layout='fan', open.angle = 10)

# Adding tip labels
p <- p + geom_tiplab(size = 0.8, aes(angle = angle))

# Adding barplot
p2 <- p + ggnewscale::new_scale_fill() +
  geom_fruit(
    data=dd_filtered,
    geom=geom_col,
    mapping=aes(y=species, x=mean_slope),
    offset = 0.75,
    pwidth=0.4,
    axis.params=list(
      axis='x',
      text="mean slope",
      text.size=1,# add axis text of the layer.
      text.angle=-90, # the text size of axis.
      hjust=0  # adjust the horizontal position of text of axis.
    ),
    grid.params=list(),
    # add the grid line of the external bar plot.
  ) +
  theme( #when adding legend
    legend.background=element_rect(fill=NA),
    legend.title=element_text(size=7),
    legend.text=element_text(size=6),
    legend.spacing.y = unit(0.02, "cm") 
  )

print(p2)
ggsave("tree.pdf", width = 18, height = 18, units = "cm")

######## Testing

index <- match(phyloResult$scenario.3[["tip.label"]], dd_filtered$species)
tab1 <- as.data.frame(cbind(phyloResult$scenario.3[["tip.label"]], dd_filtered$Species[index]))
head(tab1)

data2<-dd_filtered[index,] #now the traits dataframe has the species in the same order from the tree file
dimnames(data2)[[1]] <- phyloResult$scenario.3[["tip.label"]] #  obs. of    variables

#Results
slope <- as.matrix(data2[,2])
dimnames(slope)[[1]] <- dimnames(data2)[[1]]
slope

#calculate Blomberg's K and Pagels Lambda
phylosig(phyloResult$scenario.3,slope,method="K",test=TRUE, nsim=9999)
#Phylogenetic signal K : 0.0943208
#P-value (based on 999 randomizations) : 0.757758
phylosig(phyloResult$scenario.3 ,slope,method="lambda",test=TRUE, nsim=9999)
#Phylogenetic signal lambda : 7.34769e-05
#logL(lambda) : -524.568
#LR(lambda=0) : -0.00220029
#P-value (based on LR test) : 1
