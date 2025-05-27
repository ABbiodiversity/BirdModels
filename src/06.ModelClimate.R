# ---
# title: ABMI models - run climate models
# author: Elly Knight
# created: June 5, 2024
# ---

#NOTES################################

#PURPOSE: This script runs the first stage of the modelling process, which is the climate models that are run with all data for the entire province.

#This script is written to be run on the alliance canada cluster in parallel. Messages are therefore printed in the output for each step in the script to help with debugging. Use the Climate.sh object to run the script on the cluster. The cedar cluster is recommended.

#There are two objects at the beginning of the script that can be set to TRUE if running tests. One controls whether you are running on a subset of model iterations (i.e., testing), and one controls whether you are running on your local machine or alliance canada. Ideally you would:
#1. set test to true and cc to false (test on local)
#2. set cc to true (test on alliance canada)
#3. set test to false (run full model set on alliance canada)

#The script is set to inventory the models already run and remove them from the to-do list ("loop" object) every time the script is run. This means you can submit it as a [relatively] small job on alliance canada instead of requesting enough resources for the entire list of models (which is a lot). Just make sure you keep the "results" folder with all the output in it until you are finished running everything, because this is what is being used to inventory the models already run.

#The steps for running on alliance canada (for this and other scripts are):
#1. Transfer this script and your data object between local computer and alliance canada's servers using Globus Connect.
#2. Create any folders you need for saving output into (as per your script). You'll want to work in the scratch space, given the volume of output we will produce.
#3. Log in to alliance canada in the terminal with the ssh function.
#4. Figure out which modules you need to run your R packages (if making any changes) and load them.
#5. Load an instance of R in the test node and install the packages you need. Close it.
#6. Use the nano function to create a shell script that tells the slurm the resources you need, the modules you need, and the script to run. NOTE: You have to do this via the terminal or with a Linux machine. The slurm can't read the encoding if you create it on a windows machine and transfer it over via Globus.
#7. Use the cd function to navigate to the write working directory within copute canada's file servers (i.e., where you put your files with globus connect)
#8. Run your script with sbatch and the name of your sh file!
#9. Use squeue -u and your username to check on the status of your job.
#10. Rerun the script as many times as you need to work through all the required models (see note above).

#see medium.com/the-nature-of-food/how-to-run-your-r-script-with-compute-canada-c325c0ab2973 for the most straightforward tutorial I found for alliance canada

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
if(cc){ nodes <- 48}
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
source("00.ClimateModels.R")

#9. Number of climate models----
n <- length(modelsclimate)

#10. Link functions----
inv.link  <- function (eta) {pmin(pmax(exp(eta), .Machine$double.eps), .Machine$double.xmax)}
link <- poisson()$linkfun

#11. Load data objects----
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
    
    #Make a species folder in predictions
    if(!(file.exists(file.path(root, "Results", "ClimateModels", "Predictions", species.i)))){
      dir.create(file.path(root, "Results", "ClimateModels", "Predictions", species.i))
    }
    
    #Save the prediction
    write.csv(averageprediction, file = file.path(root, "Results", "ClimateModels", "Predictions", species.i, paste0("ClimateModelPrediction_", species.i, "_", boot.i, ".csv")), row.names = FALSE)
    
    #Make a species folder in coefficients
    if(!(file.exists(file.path(root, "Results", "ClimateModels", "Coefficients", species.i)))){
      dir.create(file.path(root, "Results", "ClimateModels", "Coefficients", species.i))
    }
    
    #Save the coefficients
    write.csv(averagecoefficients, file = file.path(root, "Results", "ClimateModels", "Coefficients", species.i, paste0("ClimateModelCoefficients_", species.i, "_", boot.i, ".csv")), row.names = TRUE)
  }
  
  #Save any try-errors
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
b <- c(1:100)

#3. Make todo list----
todo <- expand.grid(species = spp, bootstrap = b) |>
  arrange(species)

#4. Check against models already run----
done <- data.frame(file = list.files(file.path(root, "Results", "ClimateModels", "Predictions"), pattern="*.csv", recursive = TRUE)) |>
  separate(file, into=c("f1", "f2", "species", "bootstrap", "f3")) |> 
  dplyr::select(species, bootstrap) |>
  mutate(bootstrap = as.numeric(bootstrap)) |>
  inner_join(data.frame(file = list.files(file.path(root, "Results", "ClimateModels", "Coefficients"), pattern="*.csv", recursive=TRUE)) |>
               separate(file, into=c("f4", "f5", "species", "bootstrap", "f6")) |>
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
  
  #7. Run model function in parallel----
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
