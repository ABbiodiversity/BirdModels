# ---
# title: ABMI models - run landcover models
# author: Elly Knight
# created: June 6, 2024
# ---

#PREAMBLE############################

#1. Load packages----
print("* Loading packages on master *")
library(tidyverse) #basic data wrangling
library(parallel) #parallel computing

#2. Determine if testing and on local or cluster----
test <- FALSE
cc <- FALSE

#3. Set nodes for local vs cluster----
if(cc){ nodes <- 32}
if(!cc | test){ nodes <- 4}

#4. Create and register clusters----
print("* Creating clusters *")
cl <- makePSOCKcluster(nodes, type="PSOCK")

#5. Set root path----
print("* Setting root file path *")
if(cc){root <- "/scratch/ecknight/ABModels"}
if(!cc){root <- "G:/Shared drives/ABMI_ECKnight/Projects/BirdModels"}

tmpcl <- clusterExport(cl, c("root"))

#6. Load packages on clusters----
print("* Loading packages on workers *")
tmpcl <- clusterEvalQ(cl, library(tidyverse))
tmpcl <- clusterEvalQ(cl, library(AICcmodavg))
tmpcl <- clusterEvalQ(cl, library(MuMIn))

#7. Load data package----
print("* Loading data on master *")

load(file.path(root, "Data", "Stratified.Rdata"))

#8. Load model script----
if(cc){source("00.NorthModels.R")
  source("00.SouthModels.R")}
if(!cc){source("modelling2.0/00.NorthModels.R")
  source("modelling2.0/00.SouthModels.R")}

#9. Load data objects----
print("* Loading data on workers *")

tmpcl <- clusterExport(cl, c("bird", "off", "covs", "boot", "birdlist", "modelsnorth", "modelssouth"))

#WRITE FUNCTION##########

model_landcover <- function(i){
  
  #2. Loop settings----
  boot.i <- loop$bootstrap[i]
  species.i <- as.character(loop$species[i])
  region.i <- loop$region[i]
  
  #3. Get the data----
  covs.i <- covs[covs$surveyid %in% boot[,boot.i],]
  bird.i <- bird[bird$surveyid %in% boot[,boot.i], species.i]
  off.i <- off[off$surveyid %in% boot[,boot.i], species.i]
  
  #4. Load the climate model predictions----
  climate.i <- read.csv(loop$path[i])
  
  #5. Put together----
  dat.i <- data.frame(count = bird.i, offset = off.i, climate = climate.i$x) |>
    cbind(covs.i)
  
  #6. Subset the data and rename the weights----
  if(region.i=="north"){use.i <- dplyr::filter(dat.i, useNorth==TRUE) |>
    rename(weight = vegw)}
  if(region.i=="south") {use.i <- dplyr::filter(dat.i, useSouth==TRUE) |>
    rename(weight = soilw)}
  
  #7. Rename the model script---
  if(region.i=="north"){models <- modelsnorth}
  if(region.i=="south"){models <- modelssouth}
  
  #8. Make model output list---
  lc.list <- list()
  
  #9. Run null model----
  lc.list[[1]] <- glm(count ~ climate + offset(offset),
                      data = use.i,
                      family = "poisson",
                      weights = use.i$weight)
  
  #10. Set up model stage loop----
  for(j in 1:length(models)){
    
    #11. Get models----
    mods.j <- models[j][[1]]
    
    #12. Set up model loop----
    full.list <- list()
    for(k in 1:length(mods.j)){
      
      full.list[[k]] <- try(update(lc.list[[j]], formula=mods.j[[k]]))
      
    }
    
    #13. Remove try errors----
    full.list <- full.list[!sapply(full.list, inherits, "try-error")]
    
    if(length(full.list) > 0){
      
      #13. Get the BIC table----
      bictable <- AICcmodavg::bictab(cand.set = full.list, sort = F)
      
      #14. Pick best model----
      bestmodeltable <- bictable |>
        dplyr::filter(K > 1, Delta_BIC <=2) |>
        dplyr::filter(K == min(K))
      
      #15. Save that model for the next stage----
      lc.list[[j+1]] <- full.list[[as.numeric(str_sub(bestmodeltable$Modnames[1], 4, 5))]]
      
    }
    
  }
  
  #16. Get the final model----
  bestmodel <- lc.list[[length(lc.list)]]
  
  #17. Save it----
  if(region.i=="north"){ save(bestmodel, file = file.path(root, "Results", "LandcoverModels", "Models", "north", paste0("NorthModel_", species.i, "_", boot.i, ".Rdata"))) }
  
  if(region.i=="south"){ save(bestmodel, file = file.path(root, "Results", "LandcoverModels", "Models", "south", paste0("SouthModel_", species.i, "_", boot.i, ".Rdata"))) }
  
  #18. Save the coefficients----
  coef <- data.frame(name = names(bestmodel$coefficients),
                     value = as.numeric(bestmodel$coefficients)) |>
    mutate(name = ifelse(name=="(Intercept)", "Intercept", name))
  
  if(region.i=="north"){ write.csv(coef, file = file.path(root, "Results", "LandcoverModels", "Coefficients", "north", paste0("NorthModel_", species.i, "_", boot.i, ".csv")),  row.names = FALSE) }
  
  if(region.i=="south"){ write.csv(coef, file = file.path(root, "Results", "LandcoverModels", "Coefficients", "south", paste0("SouthModel_", species.i, "_", boot.i, ".csv")), row.names = FALSE) }
  
  #19. Tidy up----
  rm(lc.list, bestmodel, bictable, bestmodeltable)
  
}

#RUN MODELS###############

#1. Get list of climate models----
climate <- data.frame(path = list.files(file.path(root, "Results", "ClimateModels", "Predictions"),  full.names = TRUE, pattern="*.csv"),
                      file = list.files(file.path(root, "Results", "ClimateModels", "Predictions"), pattern="*.csv")) |>
  separate(file, into=c("model", "species", "bootstrap", "filetype")) |>
  dplyr::select(-model, -filetype)

#3. Make list----
todo <- birdlist |>
  dplyr::select(species, region) |>
  inner_join(climate, multiple="all")

#4. Check against models already run----
done <- data.frame(file = list.files(file.path(root, "Results", "LandCoverModels", "Models"), pattern="*.Rdata", recursive = TRUE)) |>
  separate(file, into=c("region", "f1", "species", "bootstrap", "f2")) |>
  dplyr::select(region, species, bootstrap) |>
  inner_join(data.frame(file = list.files(file.path(root, "Results", "LandCoverModels", "Coefficients"), pattern="*.csv", recursive = TRUE)) |>
               separate(file, into=c("region", "f3", "species", "bootstrap", "f4")) |>
               dplyr::select(region, species, bootstrap))

loop <- anti_join(todo, done)

if(nrow(loop) > 0){
  
  #For testing
  if(test) {loop <- loop[1:nodes,]}
  
  print("* Loading model loop on workers *")
  tmpcl <- clusterExport(cl, c("loop"))
  
  #10. Run BRT function in parallel----
  print("* Fitting models *")
  mods <- parLapply(cl,
                    X=1:nrow(loop),
                    fun=model_landcover)
  
}

#CONCLUDE####

#1. Close clusters----
print("* Shutting down clusters *")
stopCluster(cl)

if(cc){ q() }
