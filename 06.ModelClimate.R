# ---
# title: ABMI models - run climate models
# author: Elly Knight
# created: June 5, 2024
# ---

#NOTES################################

#TO DO: MAYBE CHANGE BIRDLIST?

#TO DO: MAKE THE MODEL SAVE LIGHTER########

#FUTURE VERSION: Try using climate variables that are smoothed via a large focal moving window (e.g., 20 km) to smooth the response and prediction surfaces

#PREAMBLE############################

#1. Load packages----
print("* Loading packages on master *")
library(tidyverse) #basic data wrangling
library(AICcmodavg) # For model selection
library(MuMIn) #Model averaging
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
if(cc){source("00.ClimateModels.R")}
if(!cc){source("modelling2.0/00.ClimateModels.R")}

#9. Number of climate models----
n <- length(modelsclimate)

#10. Bird list from previous version----
birdlist <- read.csv(file.path(root, "Data", "BirdListStatic.csv"))

#11. Link functions----
inv.link  <- function (eta) {pmin(pmax(exp(eta), .Machine$double.eps), .Machine$double.xmax)}
link <- poisson()$linkfun

#1. Load data objects----
print("* Loading data on workers *")

tmpcl <- clusterExport(cl, c("bird", "off", "covs", "boot", "birdlist", "modelsclimate", "n", "inv.link", "link"))

#WRITE FUNCTION##########

model_climate <- function(i){
  
  #1. Loop settings----
  boot.i <- loop$bootstrap[i]
  species.i <- as.character(loop$species[i])
  
  #2. Get the data----
  covs.i <- covs[covs$surveyid %in% boot[,boot.i],]
  bird.i <- bird[bird$surveyid %in% boot[,boot.i], species.i]
  off.i <- off[off$surveyid %in% boot[,boot.i], species.i]
  
  dat.i <- data.frame(count = bird.i, offset = off.i) |>
    cbind(covs.i)
  
  #3. Make model list---
  climate.list <- list()
  
  #4. Run null model----
  climate.list[[1]] <- try(glm(count ~ 1  + offset(offset),
                               data = dat.i,
                               family = "poisson",
                               y=FALSE,
                               model=FALSE))
  
  if(!inherits(climate.list[[1]], "try-error")){
    
    #5. Run the other climate models----
    for(j in 1:n){climate.list[[j + 1]] <- try(update(climate.list[[1]], formula=modelsclimate[[j]]))}
    
    #6. Add XY----
    for(j in 1:n){climate.list[[j + n + 1]] <- try(update(climate.list[[j + 1]], . ~ . + Easting + Northing + Easting:Northing))}
    
    for(j in 1:n){climate.list[[j + 2*n + 1]] <- try(update(climate.list[[j + n + 1]], . ~ . + I(Easting^2) + I(Northing^2)))}
    
    #Take out any try-errors
    climate.list <- climate.list[sapply(climate.list, function(x) !inherits(x, "try-error"))]
    
    #7. Model averaging----
    averagemodel <- model.avg(climate.list , rank = "AICc")
    
    #8. Get predictions----
    averageprediction <- inv.link(predict(averagemodel, type="link", full=TRUE))
    
    #9. Get coefficients----
    averagecoefficients <- data.frame(coef = averagemodel$coefficients["full",])
    
    #10. Save some things----
    save(averagemodel, file = file.path(root, "Results", "ClimateModels", "Models", paste0("ClimateModel_", species.i, "_", boot.i, ".Rdata")))
    
    write.csv(averageprediction, file = file.path(root, "Results", "ClimateModels", "Predictions", paste0("ClimateModelPrediction_", species.i, "_", boot.i, ".csv")), row.names = FALSE)
    
    write.csv(averagecoefficients, file = file.path(root, "Results", "ClimateModels", "Coefficients", paste0("ClimateModelCoefficients_", species.i, "_", boot.i, ".csv")), row.names = TRUE)
  }
  
  if(inherits(climate.list[[1]], "try-error")){
    
    error <- data.frame(class = class(climate.list[[1]]),
                        error = climate.list[[1]][[1]],
                        species = species.i,
                        boot = boot.i)
    
    write.csv(error, file = file.path(root, "Results", "ClimateModels", "Try-errors", paste0("ClimateModelErrors_", species.i, "_", boot.i, ".csv")), row.names = FALSE)
    
  }
  
  #11. Clean up----
  rm(climate.list, averagemodel, averageprediction, averagecoefficients, dat.i)
  
}

#RUN MODELS###############

#1. Set species list----
spp <- unique(birdlist$species)

#2. Set bootstrap list----
b <- c(1:50)

#3. Make todo list----
todo <- expand.grid(species = spp, bootstrap = b) |>
  arrange(species) |>
  #    dplyr::filter(species %in% c("ALFL", "BOCH", "CAWA", "PIWO", "BOBO", "VEER", "WIWR")) |>
  dplyr::filter(species %in% c("ALFL"))

#4. Check against models already run----
done <- data.frame(file = list.files(file.path(root, "Results", "ClimateModels", "Predictions"), pattern="*.csv")) |>
  separate(file, into=c("f1", "species", "bootstrap", "f2")) |>
  dplyr::select(species, bootstrap) |>
  mutate(bootstrap = as.numeric(bootstrap)) |>
  inner_join(data.frame(file = list.files(file.path(root, "Results", "ClimateModels", "Coefficients"), pattern="*.csv")) |>
               separate(file, into=c("f3", "species", "bootstrap", "f4")) |>
               dplyr::select(species, bootstrap) |>
               mutate(bootstrap = as.numeric(bootstrap)))

#5. Check against try-error models----
error <- data.frame(file = list.files(file.path(root, "Results", "ClimateModels", "Try-errors"), pattern="*.csv")) |>
  separate(file, into=c("f1", "species", "bootstrap", "f2")) |>
  dplyr::select(species, bootstrap) |>
  mutate(bootstrap = as.numeric(bootstrap))

#6. Make modelling list----
loop <- anti_join(todo, done) |>
  anti_join(error)

if(nrow(loop) > 0){
  
  #For testing
  if(test) {loop <- loop[1:nodes,]}
  
  print("* Loading model loop on workers *")
  tmpcl <- clusterExport(cl, c("loop"))
  
  #7. Run BRT function in parallel----
  print("* Fitting models *")
  mods <- parLapply(cl,
                    X=1:nrow(loop),
                    fun=model_climate)
  
}

#CONCLUDE####

#1. Close clusters----
print("* Shutting down clusters *")
stopCluster(cl)

if(cc){ q() }
