
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





## Ancillary Function
# The records of species cover values were made by layer
# (tree layer, shrub layer, herb layer, moss layer)
# The percentage cover values of taxa occurring in different 
# layers of the same plot will be  merged, assuming a 
# random overlap of their, cover values and making sure 
# that the combined cover values cannot exceed 100. 
combine.cover <- function(x){
  # x= COV_PERC
  while (length(x)>1){
    x[2] <- x[1]+(100-x[1])*x[2]/100
    x <- x[-1]
  }
  return(x)
}

###  Load data
DT <- read.csv(file="ReSurveyGermany.csv")
DT <- data.table(DT)
str(DT) #610313 obs. of  6 variables:
# DT holds 
# PROJECT_ID:     ID of the project
# RELEVE_NR:      ID for plot records within the project
# PROJECT_ID_RELEVE_NR: ID across all projects, links DT to header
# LAYER:          code for the layer
# TaxonName:      Harmonized taxonomic name
# Cover_Perc:     Per cent cover of the species in the plot

# Load header data
header <- read.csv("Header_ReSurveyGermany.csv")
str(header) # 23641 obs. of  47 variables:

# Link DT to header
index1 <- match(DT$PROJECT_ID_RELEVE_NR, header$PROJECT_ID_RELEVE_NR)

header$RS_PROJECT_PLOT <- paste(header$PROJECT_ID, header$RS_PLOT, sep="_")
length(unique(header$RS_PROJECT_PLOT))
# In total, there are 7738 time series

# Combine cover of species across layers
DT2 <- DT[,list(Cover_Perc=combine.cover(Cover_Perc)), by=list(PROJECT_ID_RELEVE_NR,TaxonName)]
str(DT2)
# 583395 obs. of  3 variables

# Split RELEVE_NR_comb in DT into project-specific RELEVE_NR and project ID
DT2$PROJECT_ID   <- as.numeric(tstrsplit(DT2$PROJECT_ID_RELEVE_NR,"_")[[1]])
DT2$RELEVE_NR <- as.numeric(tstrsplit(DT2$PROJECT_ID_RELEVE_NR,"_")[[2]])

length(unique(DT2$PROJECT_ID))
# In total, there are 92 projects

##Create cross-link table between project names and sequential id
Projects <- header[!duplicated(header[,c("PROJECT_ID", "RS_PROJECT")]), c("PROJECT_ID", "RS_PROJECT")]

### Step. 2: Calculation of change 
# both by resurvey ID x species x time interval combinations 
# (species.change.list) and by
# resurvey ID x time interval combinations

# The analysis is done project-wise, which allows to analyse community-level changes 
# together with plot-specific changes (i.e. those taken on (semi-)permanent plots)
species.change.list <- data.frame(RS_PROJECT=NULL, RS_PLOT=NULL,
                                  from.n=NULL, to.n=NULL,from=NULL, to=NULL,species=NULL,
                                  absolute.change=NULL, relative.change=NULL, relative.rank.change=NULL,
                                  log.repsonse.change=NULL, absolute.change.colonizer=NULL,
                                  absolute.change.extinct=NULL)
# richness.change.list holds the different metrics for the
# resurvey ID x time interval combinations
richness.change.list <- data.frame(RS_PROJECT=NULL, RS_PLOT=NULL,
                                   from.n=NULL, to.n=NULL,from=NULL, to=NULL,log.richness.change=NULL,
                                   log.shannon.change=NULL, log.evenness.change=NULL, curve.diff=NULL,
                                   log.rank.change=NULL,log.rank.change.sum=NULL,
                                   mean.cover.change=NULL, median.cover.change=NULL,
                                   diff.gains.losses=NULL)

# loop for all 92 projects
for (i in c(1:length(unique(DT2$PROJECT_ID)))){
  print(i)
  # select the species cover data from a given project
  species.long <- DT2[DT2$PROJECT_ID==i,]
  # turn the long format into a RELEVE_NR by TaxonName matrix
  species <- reshape2::acast(species.long, RELEVE_NR~TaxonName, value.var ="Cover_Perc", fill=0)
  # select the header data from a given project
  env <- header[header$PROJECT_ID==i,]
  RELEVE_NR=env$RELEVE_NR
  year <- as.matrix(env$YEAR)
  dimnames(year)[[1]] <- env$RELEVE_NR
  dimnames(year)[[2]] <- "year"
  # RS_PLOT holds a unique (within the site) code of the resurveyed plot; 
  # it is used to pair observations from different times recorded in 
  # the same plot; gives a unique identifier for the resurveyed plot or 
  # set of plots in time if combined with RS_PROJECT (see RS_PROJECT_PLOT). 
  # Several plots in the same year might have the same RS_PLOT code if they have 
  # to be summarised for temporal comparisons.
  RS_PLOT <- env$RS_PLOT
  # plot.list holds the list of all RS_PLOT
  plot.list <- sort(names(table(RS_PLOT))[table(RS_PLOT)>1])
  # produce a presence/absence version of the RELEVE_NR by TaxonName matrix
  species_pa <- species
  species_pa[species_pa > 0] <- 1
  
  # link vegetation data to header data for this project
  index2 <- match(RELEVE_NR,env$RELEVE_NR)
  # loop for all RS_PLOT
  for (j in 1:length(plot.list)){
    target.year <- year[RS_PLOT %in% plot.list[j],]
    if (length(unique(target.year))>1){
      # Check whether the plot series actually has more than one data
      # this might happen when a previous plot was recorded several times
      # but no new record was made
      
      # produce a subset of the plot by species matrix
      species.target.plot <- species[RS_PLOT %in% plot.list[j],]
      # Remove empty cols, i.e. species that do not occur in the subset
      species.target.plot <- species.target.plot[,colSums(species.target.plot)!=0, drop=F]
      # produce a presence/absence version of the species.target.plot matrix
      species.target.plot.pa <- species.target.plot 
      species.target.plot.pa[species.target.plot.pa>0] <- 1
      # calculate plot species richness
      richness.target.plot <- rowSums(species.target.plot.pa)
      names(richness.target.plot) <- target.year
      # calculate plot Shannon diversity
      if (dim(species.target.plot)[[2]]>1){
        shannon.target.plot <- vegan::diversity(species.target.plot, index = "shannon")
      } else {
        shannon.target.plot <- c(0,0)
      }
      names(shannon.target.plot) <- target.year
      # calculate plot Pilou evenness
      evenness.target.plot <- shannon.target.plot/log(richness.target.plot )
      
      # take the mean of replicated plots
      species.target.plot <- aggregate(species.target.plot, by=list(target.year),FUN=mean)
      rownames(species.target.plot) <- species.target.plot$Group.1
      species.target.plot <- species.target.plot[,-1, drop=F] 
      # remove the first column from the aggregation function
      species.target.plot <- species.target.plot[order(as.numeric(rownames(species.target.plot)),decreasing=T),,drop=F]
      # reorder the plot records from the latest year being the first row
      
      # take the mean of richness
      richness.target.plot <- aggregate(richness.target.plot, by=list(target.year),FUN=mean)
      rownames(richness.target.plot) <- richness.target.plot$Group.1
      richness.target.plot <- richness.target.plot[,-1, drop=F]
      richness.target.plot <- richness.target.plot[order(as.numeric(rownames(richness.target.plot)),decreasing=T),,drop=F]
      
      # take the mean of shannon
      shannon.target.plot <- aggregate(shannon.target.plot, by=list(target.year),FUN=mean)
      rownames(shannon.target.plot) <- shannon.target.plot$Group.1
      shannon.target.plot <- shannon.target.plot[,-1, drop=F]
      shannon.target.plot <- shannon.target.plot[order(as.numeric(rownames(shannon.target.plot)),decreasing=T),,drop=F]
      # take the mean of eveness
      evenness.target.plot <- aggregate(evenness.target.plot, by=list(target.year),FUN=mean)
      rownames(evenness.target.plot) <- evenness.target.plot$Group.1
      evenness.target.plot <- evenness.target.plot[,-1, drop=F]
      evenness.target.plot <- evenness.target.plot[order(as.numeric(rownames(evenness.target.plot)),decreasing=T),,drop=F]
      
      # calculate differences in rank abundance curves
      # according to Avolio, M.L., Carroll, I.T., Collins, S.L., Houseman, G.R., 
      # Hallett, L.M., Isbell, F., Koerner, S.E., Komatsu, K.J., Smith, M.D., 
      # Wilcox, K.R., 2019. A comprehensive approach to analyzing community 
      # dynamics using rank abundance curves 10: e02881. 10.1002/ecs2.2881
      species.target.plot2 <- sweep(species.target.plot, 1, apply(species.target.plot,1,FUN=sum), FUN="/")
      species.target.plot2 <- species.target.plot2[order(as.numeric(rownames(species.target.plot2)),decreasing=T),,drop=F]
      # species.target.plot2 holds relative abundance values for each plot
      # the species are ranked per row using the matrixStats package
      species.target.plot.ranks2 <- as.matrix(rowRanks(as.matrix(-species.target.plot2, ties.method="average")))
      dimnames(species.target.plot.ranks2)[[1]] <- row.names(species.target.plot2)
      dimnames(species.target.plot.ranks2)[[2]] <- colnames(species.target.plot2)
      # turn ranks into relative ranks
      species.target.plot.ranks2 <- sweep(species.target.plot.ranks2, 1, apply(species.target.plot.ranks2,1,FUN=max), FUN="/")
      
      # calculate delta curve (curve.diff) according to Avolio et al. (2019)
      unique.relative.ranks <- c(0,unique(sort(species.target.plot.ranks2)))
      species.target.plot.ranks3 <- apply(species.target.plot.ranks2,2,FUN=cut,breaks=unique.relative.ranks, labels=F)
      dimnames(species.target.plot.ranks3)[[1]] <- rownames(species.target.plot2)
      curve.diff <- 0
      y1sum <- 0
      y2sum <- 0
      for (k in 2: length(unique.relative.ranks)){
        y1 <- species.target.plot2[1,species.target.plot.ranks3[1,]==k-1]
        y2 <- species.target.plot2[2,species.target.plot.ranks3[2,]==k-1]
        y1sum <- y1sum + ifelse(length(y1)>0,sum(y1),0)
        y2sum <- y2sum + ifelse(length(y2)>0,sum(y2),0)
        curve.diff <- curve.diff+y1sum-y2sum
      }
      curve.diff
      
      # calculate rank change (rank.change) according to Avolio et al. (2019)
      species.target.plot.ranks4 <- species.target.plot.ranks2[1,]-species.target.plot.ranks2[2,]
      rank.change <- sum(abs(species.target.plot.ranks4))/ncol(species.target.plot.ranks2)
      
      # calculate rank change separately for pos. and neg. change (rank.change.sign)
      rank.change.neg <- mean(species.target.plot.ranks4[species.target.plot.ranks4<0])
      rank.change.pos <- mean(species.target.plot.ranks4[species.target.plot.ranks4>0])
      
      rank.change.neg.sum <- sum(species.target.plot.ranks4[species.target.plot.ranks4<0])
      rank.change.pos.sum <- sum(species.target.plot.ranks4[species.target.plot.ranks4>0])
      
      # turn cover values in relative change values in this plot 
      # relative to max cover of the whole time series for that plot
      # currently not used in the paper, but was used in a previous version
      # and gave similar results
      species.target.plot3 <- sweep(species.target.plot, 2, apply(species.target.plot[,,drop=F],2,FUN=max), FUN="/")
      
      # calculate mean and median relative cover
      species.target.plot4 <- species.target.plot
      species.target.plot4[species.target.plot4==0] <- NA
      mean.cover <- apply(species.target.plot4,1, FUN=mean, na.rm=T)
      median.cover <- apply(species.target.plot4,1, FUN=median, na.rm=T)
      
      # Number of years in the time series
      n <- table(target.year)[order(as.numeric(names(table(target.year))),decreasing=T)]
      # Compare subsequent changes in a time series
      for (k in 1:(length(n)-1)){
        # diff.absolute.change is the difference in absolute cover values
        diff.absolute.change <- species.target.plot[k,,drop=F] - species.target.plot[k+1,,drop=F]
        diff.absolute.change[species.target.plot[k,]==0 & species.target.plot[k+1,]==0] <- NA
        # diff.relative.change is the difference in absolute cover values
        diff.relative.change <- species.target.plot3[k,,drop=F] - species.target.plot3[k+1,,drop=F]
        diff.relative.change[species.target.plot3[k,]==0 & species.target.plot3[k+1,]==0] <- NA
        # diff.relative.rank.change is the difference in relative ranks
        diff.relative.rank.change <- species.target.plot.ranks2[k,] - species.target.plot.ranks2[k+1,]
        diff.relative.rank.change[species.target.plot.ranks2[k,]==0 & species.target.plot.ranks2[k+1,]==0] <- NA
        # diff.absolute.change.colonizer is diff.absolute.change only for 
        # new colonizers in the interval
        diff.absolute.change.colonizer <- diff.absolute.change
        # only keep records that are new and did not occur in Year 1
        diff.absolute.change.colonizer[species.target.plot[k+1,]!=0] <- NA
        # diff.absolute.change.extinct is diff.absolute.change only for 
        # species that went extinct in the interval
        diff.absolute.change.extinct <- diff.absolute.change
        # only keep records that went extinct and did not occur in Year 2
        diff.absolute.change.extinct[species.target.plot[k,]!=0] <- NA
        # calculate number of species that increased or decreased in an interval
        diff.rel.cover <- species.target.plot[k,] - species.target.plot[k+1,]
        n.gains <- length(diff.rel.cover[diff.rel.cover>0])
        n.losses <- length(diff.rel.cover[diff.rel.cover<0])
        
        # collect all metrics for resurvey ID x species x time interval combinations
        
        ## FMS: Switching to LIST increases the speed of this loop substantially (since it avoids R overwriting a vector with 10^5 rows at every epoch)
        species.change.list <-  rbind(species.change.list, data.frame(
          RS_PROJECT=Projects$RS_PROJECT[i],
          RS_PLOT=plot.list[j],
          from.n=as.numeric(n[k+1]), to.n=as.numeric(n[k]), from=as.numeric(names(n)[k+1]),
          to=as.numeric(names(n)[k]),
          species=names(diff.relative.change),
          absolute.change=as.numeric(diff.absolute.change),
          relative.change=as.numeric(diff.relative.change),
          relative.rank.change=as.numeric(diff.relative.rank.change),
          absolute.change.colonizer=as.numeric(diff.absolute.change.colonizer),
          absolute.change.extinct=as.numeric(diff.absolute.change.extinct)))
        
        # collect all metrics for resurvey ID x time interval combinations
        # and calculate log response ratios for all metrics
        richness.change.list <- rbind(richness.change.list, data.frame(
          RS_PROJECT=Projects$RS_PROJECT[i],RS_PLOT=plot.list[j],
          from.n=as.numeric(n[k+1]), to.n=as.numeric(n[k]), from=as.numeric(names(n)[k+1]),
          to=as.numeric(names(n)[k]),
          log.richness.change=log(richness.target.plot[k,]/richness.target.plot[k+1,]),
          log.shannon.change=log(shannon.target.plot[k,]/shannon.target.plot[k+1,]),
          log.evenness.change=log(evenness.target.plot[k,]/evenness.target.plot[k+1,]),
          curve.diff=curve.diff,
          log.rank.change=log(rank.change.pos/-rank.change.neg),
          log.rank.change.sum=log(rank.change.pos.sum/-rank.change.neg.sum),
          mean.cover.change=log(mean.cover[k]/mean.cover[k+1]),
          median.cover.change=log(median.cover[k]/median.cover[k+1]),
          diff.gains.losses=n.gains-n.losses))
      }
    }
  }
}


richness.change.list$RS_PROJECT_PLOT <- paste(richness.change.list$RS_PROJECT,richness.change.list$RS_PLOT, sep="_")
# Make a unique ID from RS_PROJECT and RS_PLOT 

# remove empty rows from the resurvey ID x species x time interval combinations
str(species.change.list) #579449 obs. of  13 variables:
species.change.list2 <- species.change.list[!is.na(species.change.list$absolute.change),]
# add some variables for analysis
species.change.list2$mid.year <- ceiling((species.change.list2$to-species.change.list2$from)/2+species.change.list2$from)
species.change.list2$decade <- floor(species.change.list2$mid.year/10)*10
species.change.list2$end.decade <- floor(species.change.list2$to/10)*10
species.change.list2$slope <- species.change.list2$relative.change/(species.change.list2$to-species.change.list2$from)
species.change.list2$RS_PROJECT_PLOT <- paste(species.change.list2$RS_PROJECT,
                                              species.change.list2$RS_PLOT, sep="_")
str(species.change.list2) #458311 obs. of  17 variables, in_MS Results, Methods


species.change.list2 <- data.table(species.change.list2)
change <- species.change.list2[,{n_pos <- length(absolute.change[absolute.change>0]);
n_neg <- length(absolute.change[absolute.change<0]);
n_all <- ifelse(n_pos+n_neg>0,n_pos+n_neg,1);
est <- binom.test(n_pos,n_all);
list(n=length(absolute.change),
     pos=length(absolute.change[absolute.change>0]),
     equal=length(absolute.change[absolute.change==0]), neg=length(absolute.change[absolute.change<0]),
     est.binom=est$estimate,
     conf.binom.minus=est$conf.int[1],
     conf.binom.plus=est$conf.int[2],
     p.values.binom=est$p.value,
     mean.absolute.change=mean(absolute.change, na.rm=T))}, by=species]
change # 1794 species
writexl::write_xlsx(change,path = "resurvey_species_change.xlsx")


# all spcies even the ones with no sign change
species.change.list2$Species <- as.factor(species.change.list2$species)

change_list <- species.change.list2 %>%
  group_by(Species) %>%
  summarise(n_records       = n(),    
    relative.change = mean(relative.change, na.rm = TRUE),
    slope = mean(slope, na.rm = TRUE),
    .groups = 'drop')

writexl::write_xlsx(change_list,path = "resurvey_species_change_list.xlsx")


