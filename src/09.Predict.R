# ---
# title: ABMI models - make model predictions
# author: Elly Knight
# created: Sept 25, 2024
# ---

#NOTES################################

#Adapted from Brandon Allen's prediction script.

#PREAMBLE############################

#1. Load packages----
library(tidyverse) #basic data wrangling
library(Matrix) #sparse matrices
library(abind) #binding arrays
library(sf) #shapefiles
library(gridExtra) #plotting

#2. Set root path for data on google drive----
root <- "G:/Shared drives/ABMI_ECKnight/Projects/BirdModels"

#3. Restrict scientific notation----
options(scipen = 99999)

#4. Load the coefficients----
load(file.path(root, "Results", "Birds2024.RData"))

#5. Define link functions----
inv.link  <- function (eta) {pmin(pmax(exp(eta), .Machine$double.eps), .Machine$double.xmax)}
link <- poisson()$linkfun

#6. Load the kgrid & backfill----
load(file.path(root, "Data", "gis", "kgrid_2.2.Rdata"))
load(file.path(root, "Data", "gis", "backfillV7_w2w_2021HFI.Rdata"))

#7. Load BA's functions----
source("src/00.Functions.R")

#8. Load lookups----
veg.lookup <- read.csv(file.path(root, "Data", "lookups", "lookup-veg-hf-age-v2020.csv"))
soil.lookup <- read.csv(file.path(root, "Data", "lookups", "lookup-soil-hf-v2020.csv"))

#WRANGLING##########

#1. Add intercept to climate kgrid----
kgrid$Intercept <- 1

#2. Add easting and northing derived columns
kgrid$Easting2 <- kgrid$Easting^2
kgrid$Northing2 <- kgrid$Northing^2
kgrid$EastingNorthing <- kgrid$Easting*kgrid$Northing

#2. Standardize veg & soil kgrid, convert to proportions, and add climate vars----

#Current veg
veg.cur.raw <- landscape_hf_summary(data.in = as.data.frame(as.matrix(d.wide$veg.current)),
                                    landscape.lookup = veg.lookup,
                                    class.in = "ID", class.out = "UseInAnalysis_Simplified")

veg.cur <- data.frame(veg.cur.raw/rowSums(veg.cur.raw)) |>
  mutate(LinkID = row.names(veg.cur.raw)) |>
  left_join(kgrid, by= "LinkID")
row.names(veg.cur) <- veg.cur$LinkID

#Current soil
soil.cur.raw <- landscape_hf_summary(data.in = as.data.frame(as.matrix(d.wide$soil.current)),
                                     landscape.lookup = soil.lookup,
                                     class.in = "ID", class.out = "UseInAnalysis_Simplified")

soil.cur <- data.frame(soil.cur.raw/rowSums(soil.cur.raw)) |>
  mutate(LinkID = row.names(soil.cur.raw)) |>
  left_join(kgrid, by="LinkID")
row.names(soil.cur) <- soil.cur$LinkID

#3. Simplify bird coefficients for kgrid predictions----
#Take mean of the treed fen coeffs in north
#Take mean of urban & industrial for both
#Take mean of mine and minev for both

#Get the offending coeffs
treedfen.all <- birds$north$joint[,c("TreedFenR", "TreedFen1", "TreedFen2", "TreedFen3", "TreedFen4", "TreedFen5", "TreedFen6", "TreedFen7", "TreedFen8"),]
urbind.n.all <- birds$north$joint[,c("Urban", "Industrial"),]
urbind.s.all <- birds$south$joint[,c("Urban", "Industrial"),]

#Take the means for the north
means.n <- array(0, c(dim(birds$north$joint)[1], 2, dim(birds$north$joint)[3]))
for(i in 1:dim(treedfen.all)[3]){
  means.n[,1,i] <- rowMeans(treedfen.all[,,i])
  means.n[,2,i] <- rowMeans(urbind.n.all[,,i])
}

#Take the means for the south
means.s <- array(0, c(dim(birds$south$joint)[1], 1, dim(birds$south$joint)[3]))
for(i in 1:dim(urbind.s.all)[3]){
  means.s[,1,i] <- rowMeans(urbind.s.all[,,i])
}

#Rename
dimnames(means.n) <- list(dimnames(birds$north$joint)[[1]], c("TreedFen", "UrbInd"), dimnames(birds$north$joint)[[3]])

dimnames(means.s) <- list(dimnames(birds$south$joint)[[1]], c("UrbInd"), dimnames(birds$south$joint)[[3]])

#Remove the unwanted
removed.n <- birds$north$joint[,!dimnames(birds$north$joint)[[2]] %in% c("TreedFenR", "TreedFen1", "TreedFen2", "TreedFen3", "TreedFen4", "TreedFen5", "TreedFen6", "TreedFen7", "TreedFen8", "Urban", "Industrial", "MineV"),]

removed.s <- birds$south$joint[,!dimnames(birds$south$joint)[[2]] %in% c("Urban", "Industrial", "MineV"),]

#Put back together
northjoint <- abind(removed.n, means.n, along = 2)
southjoint <- abind(removed.s, means.s, along = 2)

#4. Some renaming----
dimnames(northjoint)[[2]] <- gsub("GrassHerb", "Grass", dimnames(northjoint)[[2]])
dimnames(southjoint)[[2]] <- gsub("GrassHerb", "Grass", dimnames(southjoint)[[2]])

#MAKE PREDICTIONS##########

#1. List of spp----
#use the median Bootstrap from the birdlist object
spp <- birds$species |>
  dplyr::filter(!(ModelNorth==FALSE & ModelSouth==FALSE)) |>
  mutate(label=paste0(Comments, "_", Bootstrap))

#2. Set up loop----
predictions <- list()
for(i in 1:nrow(spp)){
  
  #3. Determine if the spp is north, south, or both----
  mod.i <- case_when(spp$ModelNorth[i]==TRUE & spp$ModelSouth[i]==TRUE ~ "both",
                     spp$ModelNorth[i]==TRUE & spp$ModelSouth[i]==FALSE ~ "north",
                     spp$ModelNorth[i]==FALSE & spp$ModelSouth[i]==TRUE ~ "south")
  
  #4. Get the coefficients----
  if(mod.i %in% c("both", "north")){
    climate.i <- birds$north$marginal[spp$SpeciesID[i],,spp$Bootstrap[i]]
  }
  
  if(mod.i=="south"){
    climate.i <- birds$south$marginal[spp$SpeciesID[i],,spp$Bootstrap[i]]
  }
  
  #5. Make predictions----
  
  #Get the climate variables
  use.climate <- kgrid |>
    dplyr::select(names(climate.i)) |>
    as.matrix()
  rownames(use.climate) <- kgrid$LinkID
  
  #make the climate predictions
  raw.climate <- matrix(inv.link(drop(use.climate %*% climate.i)),
                        ncol=1,
                        dimnames = list(rownames(use.climate), "Climate"))
  
  q99 <- quantile(raw.climate, 0.99)
  pred.climate <- ifelse(raw.climate > q99, q99, raw.climate)
  
  #get the veg coefficients and make veg predictions
  if(mod.i %in% c("both", "north")){
    
    #Get the coefficients
    vegclim.i <- northjoint[spp$SpeciesID[i],,spp$Bootstrap[i]]
    veg.i <- vegclim.i[names(vegclim.i)!="Climate"]
    
    #Get the kgrid veg data, but take out climate
    use.veg <- veg.cur |>
      dplyr::select(names(veg.i)) |>
      as.matrix()
    
    #predict the joint climate contribution
    joint.climate.veg <- matrix(pred.climate * vegclim.i["Climate"],
                                ncol = ncol(use.veg), nrow = nrow(use.veg))
    
    #Create climate adjusted veg covariates
    veg.coef <- t(t(joint.climate.veg) + veg.i)
    colnames(veg.coef) <- names(veg.i)
    
    #prediction
    veg.cur$Veg <- rowSums(use.veg * inv.link(veg.coef))
    
    if(mod.i=="north"){soil.cur$Soil <- 0}
    
  }
  
  #get the soil coefficients and make soil predictions
  if(mod.i %in% c("both", "south")){
    
    # #Get the coefficients
    soilclim.i <- southjoint[spp$SpeciesID[i],,spp$Bootstrap[i]]
    soil.i <- soilclim.i[!names(soilclim.i) %in% c("Climate", "pAspen")]
    
    #Get the kgrid soil data
    use.soil <- soil.cur |>
      dplyr::select(names(soil.i)) |>
      as.matrix()
    
    #get the paspen
    paspen.pred <- kgrid$pAspen
    
    #predict the joint climate contribution
    joint.climate.soil <- matrix(pred.climate * soilclim.i["Climate"],
                                 ncol = ncol(use.soil), nrow = nrow(use.soil))
    
    #predict the joint paspen contribution
    joint.paspen.soil <- matrix(paspen.pred * soilclim.i["pAspen"],
                                ncol = ncol(use.soil), nrow = nrow(use.soil))
    
    #put together
    joint.soil <- joint.climate.soil + joint.paspen.soil
    
    #Create climate adjusted soil covariates
    soil.coef <- t(t(joint.soil) + soil.i)
    
    #prediction
    soil.cur$Soil <- rowSums(use.soil * inv.link(soil.coef))
    
    if(mod.i=="south"){veg.cur$Veg <- 0}
    
  }
  
  #6. Combine & store----
  q <- spp$PlotQuantile[i]
  predictions[[i]] <- try(data.frame(Climate = pred.climate,
                                     Veg = veg.cur$Veg,
                                     Soil = soil.cur$Soil,
                                     wN = veg.cur$wN,
                                     wS = veg.cur$wS,
                                     LinkID = veg.cur$LinkID) |>
                            mutate(Veg = ifelse(is.infinite(Veg), median(Veg), Veg),
                                   Soil = ifelse(is.infinite(Soil), median(Soil), Soil)) |>
                            mutate(Provincial = Veg*wN + Soil*wS,
                                   Soil = ifelse(wS==0, 0, Soil),
                                   Veg = ifelse(wN==0, 0, Veg),
                                   species = spp$Comments[i],
                                   speciesname = spp$SpeciesID[i],
                                   Bootstrap = spp$Bootstrap[i]) |>
                            mutate(Provincial = ifelse(Provincial > quantile(Provincial, q), quantile(Provincial, q), Provincial)))
  
  cat("Finished", spp$Comments[i], "predictions :", i, "of", nrow(spp), "\n")
  
}

names(predictions) <- spp$Comments

#7. Save----
save(predictions, file=file.path(root, "Results", "Predictions2024.RData"))
