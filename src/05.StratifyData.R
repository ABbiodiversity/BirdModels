# ---
# title: ABMI models - stratify data
# author: Elly Knight
# created: May 10, 2024
# ---

#NOTES################################

#PURPOSE: This code picks one survey per location per block of region and time interval for each bootstrap, checks that the list of species to model have sufficient data, and filters counts to the 99.9% quantile to remove outliers

#This code uses several lookups to compile the list of birds that are actually modelled.

#FUTURE VERSION: count outlier removal may make more sense in a previous step.
#FUTURE VERSION: Bird list lookups are a relic of previous modeling versions and could probably be combined for parsimony.
#FUTURE VERSION: The list of species modelled could probably be expanded based on exploration of sample sizes and model assumptions.

#PREAMBLE############################

#1. Load packages----
library(tidyverse) #basic data wrangling

#2. Set root path for data on google drive----
root <- "G:/Shared drives/ABMI_ECKnight/Projects/BirdModels"

#3. Restrict scientific notation----
options(scipen = 99999)

#4. Load wrangled data set----
load(file.path(root, "Data", "Wrangled.Rdata"))

#BLOCKING UNITS#####

#1. Make bins----

#number of spatial bins is same as other taxa, has been updated to UTM
#temporal bins are determined by quantiles
tbins <- quantile(covs$year, probs=seq(0, 1, by=0.2))

bins <- covs |>
  mutate(blockx = cut(Easting, c(180000, 407000, 634000, 861000)),
         blocky = cut(Northing, c(5420000, 5730000, 6040000, 6350000, 6660000)),
         blockt = cut(year, tbins),
         blockxy = interaction(blockx, blocky, sep="::"),
         block = interaction(blockxy, blockt, sep="::"))

ftable(bins$blockt, bins$blockx, bins$blocky)

#2. Add location id----
#round UTM to nearest 20 m
tolerance <- 20

loc <- bins |>
  mutate(Xr = round(Easting/tolerance)*tolerance,
         Yr = round(Northing/tolerance)*tolerance) |>
  group_by(Xr, Yr) |>
  mutate(locationid = cur_group_id()) |>
  ungroup()

#3. Get unique locations----
covsu <- loc |>
  dplyr::select(locationid, Xr, Yr) |>
  unique()

#4. Add quantiles for witholding data----
set.seed(1234)
covsu$q <- sample.int(100, nrow(covsu), replace=TRUE)

#5. Join back to all surveys----
covsbins <- loc |>
  left_join(covsu) |>
  arrange(surveyid)

#MAKE BOOTSTRAPS##########

#1. Set parameters----
#Minimum sample size
NMIN <- 20

#Number of bootstraps
B <- 100

#2. Set up loop----
boot <- covsbins |>
  dplyr::filter(q > 10) |>
  dplyr::select(block, locationid) |>
  unique() |>
  dplyr::select(-block, -locationid) |>
  data.frame()

for(i in 1:B){
  
  #3. Set seed----
  set.seed(i)
  
  #4. Get random sample of one location per block----
  boot.i <- covsbins |>
    dplyr::filter(q > 10) |>
    group_by(block, locationid) %>%
    mutate(rowid = row_number(),
           use = sample(1:max(rowid), 1)) %>%
    ungroup() |>
    dplyr::filter(rowid==use)
  
  #5. Bind to boot object----
  boot[,i] <- boot.i$surveyid
  
  cat("Finished", i, "of", B, "\n")
  
}

#6. Name columns---
colnames(boot) <- c(1:B)

#SPECIES LIST####

#1. Make to do list----
loop <- expand.grid(boot = 1:B,
                    region = c("north", "south"))

#2. Set up loop----
count.list <- list()
for(i in 1:nrow(loop)){
  
  #3. Filter survey data to region and bootstrap----
  if(loop$region[i]=="north"){
    survey.i <- covs[covs$surveyid %in% boot[,loop$boot[i]],] |>
      dplyr::filter(useNorth==TRUE)
  }
  
  if(loop$region[i]=="south"){
    survey.i <- covs[covs$surveyid %in% boot[,loop$boot[i]],] |>
      dplyr::filter(useSouth==TRUE)
  }
  
  #4. Filter bird data and count detections----
  count.list[[i]] <- bird[bird$surveyid %in% survey.i$surveyid,] |>
    dplyr::select(-surveyid, -gisid) |>
    pivot_longer(ALFL:YRWA, names_to="species", values_to="count") |>
    dplyr::filter(count > 0) |>
    group_by(species) |>
    summarize(detections = n()) |>
    ungroup() |>
    mutate(region = loop$region[i],
           boot = loop$boot[i])
  
}

#5. Summarize across bootstraps----
#Remove birds for regions with a mean # of detections less than 20
nmin <- 20

#Get the bird wish list (list of species to run, same as previous years)
birdlist <- read.csv(file.path(root, "Data", "lookups", "birds-v2024.csv")) |>
  left_join(read.csv(file.path(root, "Data", "lookups", "birdlist.csv")) |> 
              rename(common = species,
                     species = code)) |> 
  mutate(north = ifelse(show %in% c("c", "n"), TRUE, FALSE),
         south = ifelse(show %in% c("c", "s"), TRUE, FALSE)) |> 
  dplyr::select(species, north, south) |> 
  pivot_longer(north:south, names_to="region", values_to="model") |> 
  dplyr::filter(model==TRUE) |> 
  dplyr::select(-model) |> 
  left_join(do.call(rbind, count.list) |>
              group_by(region, species) |>
              summarize(det.mean = round(mean(detections)),
                        det.min = min(detections)) |>
              ungroup()) |> 
  dplyr::filter(det.mean >= nmin)

#6. Filter counts to 99.9% quantile----
bird.q <- bird
for(i in 1:nrow(birdlist)){
  
  #Get just the presence points
  bird.i <- bird[,birdlist$species[i]][bird[,birdlist$species[i]]>0]
  
  #Get the quantile
  q.i <- ceiling(quantile(bird.i, 0.995))
  
  #Replace
  bird.q[,birdlist$species[i]] <- ifelse(bird[,birdlist$species[i]] > q.i,
                                         q.i, bird[,birdlist$species[i]])
  
}

#Rename and tidy
bird <- bird.q |> 
  dplyr::select(c("surveyid", all_of(unique(birdlist$species))))
rm(bird.q)

#PACKAGE######

#1. Add train/test to covs object & repackage----
covs <- covsbins |>
  mutate(use = ifelse(q <= 10, "test", "train")) |>
  dplyr::select(-blockx, -blocky, -blockt, -blockxy, -Xr, -Yr, -q)

#2. Save----
save(covs, bird, off, boot, birdlist, file=file.path(root, "Data", "Stratified.Rdata"))
