# ---
# title: ABMI models - package model coefficients
# author: Elly Knight
# created: June 17, 2024
# ---

#NOTES################################

#The translation of the model coefficients to the coefficients that are standardized across taxa relies on a lookup table. The proportions of the age coefficients in the raw models that contribute to the age stages in the standardized coefficients are area under the curve.

#PREAMBLE############################

#1. Load packages----
library(tidyverse) #basic data wrangling
library(purrr) #functional programming
library(vegan) #for picking median bootstrap
library(pROC) #AUC

#2. Set root path for data on google drive----
root <- "G:/Shared drives/ABMI_ECKnight/Projects/BirdModels"

#3. Restrict scientific notation----
options(scipen = 99999)

#4. Get the list of models run----
mods <- data.frame(path = list.files(file.path(root, "Results", "LandcoverModels", "Coefficients"), full.names = TRUE, recursive = TRUE, pattern="*.csv"),
                   file = list.files(file.path(root, "Results", "LandcoverModels", "Coefficients"), recursive = TRUE, pattern="*.csv")) |>
  separate(file, into=c("region", "f1", "f2", "species", "bootstrap")) |>
  mutate(bootstrap = as.numeric(bootstrap)) |>
  dplyr::select(-f1, -f2)

#5. Load model scripts----
source("00.NorthModels.R")
source("00.SouthModels.R")

#6. Load data----
load(file.path(root, "Data", "Stratified.Rdata"))

#7. Function to fix interaction term names----
#make formula term names sorted and predictable, i.e. always A:B instead of B:A
fix_names <- function(x, sep=":") {
  unlist(lapply(x, function(z) {
    paste(sort(strsplit(z, sep)[[1]]), collapse=sep)
  }))
}

#8. Function to get model terms-----
#Formula from Peter's approach
get_terms <- function(mods, type=c("formula", "list"), intercept=TRUE) {
  type <- match.arg(type)
  x <- unlist(lapply(unlist(mods), function(z) as.character(z)[3]))
  #    x <- unname(substr(x, 5, nchar(x)))
  x <- gsub(". + ", "", x, fixed=TRUE)
  x <- unlist(strsplit(x, "+", fixed=TRUE))
  x <- unlist(strsplit(x, "*", fixed=TRUE))
  if (type == "list")
    x <- unlist(strsplit(x, ":", fixed=TRUE))
  x <- sub("^[[:space:]]*(.*?)[[:space:]]*$", "\\1", x, perl=TRUE)
  x <- unique(x)
  if (type == "formula") {
    x <- paste("~", paste(x, collapse=" + ", sep=""))
    if (!intercept)
      x <- paste(x, "- 1")
    x <- as.formula(x)
  }
  x
}

#9. Load veg age lookup----
load(file.path(root, "Data", "lookups", "Xn-veg-v2024.Rdata"))
colnames(age) <- fix_names(colnames(age))

#CLIMATE MODEL COEFFICIENTS##########

#1. Get list of climate models----
climate <- data.frame(path = list.files(file.path(root, "Results",  "ClimateModels", "Coefficients"),  full.names = TRUE, pattern="*.csv", recursive = TRUE),
                      file = list.files(file.path(root, "Results",  "ClimateModels", "Coefficients"), pattern="*.csv", recursive = TRUE)) |>
  separate(file, into=c("f1", "model", "species", "bootstrap", "filetype")) |>
  dplyr::select(-model, -filetype, -f1) |>
  mutate(bootstrap = as.numeric(bootstrap)) |>
  inner_join(mods |>
               dplyr::select(species, bootstrap) |>
               unique())

#2. Get list of species----
spp <- unique(climate$species)

#3. Get max # of bootstraps----
bmax <- max(climate$bootstrap)

#4. Set the column names----
coefnames <- c("(Intercept)", "MAP", "FFP", "TD", "CMD", "EMT", "Easting", "Northing", "I(Easting^2)", "I(Northing^2)", "Easting:Northing")
coefrename <- c("Intercept", "MAP", "FFP", "TD", "CMD", "EMT", "Easting", "Northing", "Easting2", "Northing2", "EastingNorthing")

#5. Set up loop-----
for(i in 1:length(spp)){
  
  #6. List of models for that species----
  climate.i <- dplyr::filter(climate, species==spp[i],
                             bootstrap %in% c(1:bmax))
  
  #7. Read them in----
  #sort and fix names
  coef.list <- list()
  for(j in 1:nrow(climate.i)){
    coef.list[[j]] <- read.csv(climate.i$path[j]) |>
      pivot_wider(names_from="X", values_from = "coef") |>
      dplyr::select(all_of(coefnames))
    colnames(coef.list[[j]]) <- coefrename
  }
  
  #8. Make the array if first species----
  if(i==1){
    marginal <- array(0, c(length(spp), length(coefrename), bmax))
    dimnames(marginal) <- list(spp, coefrename, paste0("b", c(1:bmax)))
  }
  
  #9. Add species to array----
  marginal[spp[i],,] <- as.matrix(t(do.call(rbind, coef.list)))
  
  
  cat("Finished species", i, "of", length(spp), "\n")
  
}

#LANDCOVER MODEL COEFFICIENTS - NORTH#############

#1. Get the list of north models----
north <- mods |>
  dplyr::filter(region=="north")

#2. Get list of species----
spp <- unique(north$species)

#3. Get the full list of potential continuous covariates----
#don't need the categorical ones because they're in every model
names.north <- do.call(rbind, unlist(modelsnorth)) |>
  data.frame() |>
  separate(X3, into=c("X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10", "X11", "X12", "X13", "X14", "X15", "X16", "X17", "X18", "X19", "X20", "X21", "X22", "X23"), sep = " +") |>
  dplyr::select(-X1, -X2, -X3) |>
  pivot_longer(cols=X4:X23, names_to="position", values_to="name") |>
  dplyr::filter(!is.na(name), name!="+") |>
  dplyr::select(-position) |>
  unique() |>
  mutate(name = fix_names(name)) |>
  dplyr::filter(!name %in% c("method", "vegc"))

#4. Get the full north data----
covs.n <- dplyr::filter(covs, useNorth==TRUE)

#5. Get the model matrix----
Xn <- model.matrix(get_terms(modelsnorth), covs.n)
colnames(Xn) <- fix_names(c("Intercept", colnames(Xn)[-1]))

#6. Set up loop----
for(i in 1:length(spp)){
  
  #7. List of models for that species----
  north.i <- dplyr::filter(north, species==spp[i])
  
  coef.list <- list()
  for(j in 1:nrow(north.i)){
    
    #8. Get the raw coefficients-----
    #Add zeros for any covariates that weren't in the best model
    coef.j <- read.csv(north.i$path[j]) |>
      mutate(name = fix_names(name)) |>
      full_join(data.frame(name = names.north$name)) |>
      mutate(value = ifelse(is.na(value), 0, value)) |>
      suppressMessages()
    raw.j <- coef.j$value
    names(raw.j) <- coef.j$name
    
    #9. Translate to standardized & transformed coefficients for the veg coefficients----
    mu <- drop(age %*% raw.j[colnames(age)])
    lam.j <- exp(drop(age %*% raw.j[colnames(age)]))
    
    #10. Adjust the linear feature coefficients for competing models----
    #replace others with msoft if it is non-zero
    #because mSoft is competed with the others in the model set
    hf.j <- coef.j |>
      dplyr::filter(name %in% c("mWell", "mSoft", "mEnSft", "mTrSft", "mSeism")) |>
      mutate(value = exp(value)) |>
      pivot_wider(names_from=name, values_from=value) |>
      mutate(mEnSft = ifelse(mSoft!=0 & mEnSft==0, mSoft, mEnSft),
             mTrSft = ifelse(mSoft!=0 & mTrSft==0, mSoft, mTrSft),
             mSeism = ifelse(mSoft!=0 & mSeism==0, mSoft, mSeism)) |>
      pivot_longer(mWell:mSeism, names_to="name", values_to="value") |>
      data.frame()
    
    #11. Adjust the linear feature and well coefficients----
    
    #Human-modified landcover types that we don't want interfering with adjustment
    hfc <- c("Crop", "Industrial", "Mine", "RoughP", "Rural", "TameP", "Urban")
    
    #Dataframe of variables to adjust
    vars <- data.frame(var = c("mWell", "mEnSft", "mTrSft", "mSeism"))
    
    for(k in 1:nrow(vars)){
      
      #Make a mock dataframe for the variable
      zero.k <- data.frame(var = 0)
      colnames(zero.k) <- vars$var[k]
      
      #Get the rows of interest from the model matrix
      #> 0 proportion of variable, no human-modified landover type, no harvest
      rows.k <- covs.n |>
        mutate(rowid = row_number()) |>
        anti_join(zero.k) |>
        dplyr::filter(!vegc %in% hfc &
                        fcc2==0) |>
        suppressMessages()
      
      Xn.k <- Xn[rows.k$rowid, colnames(age)]
      
      #Make predictions from those rows using the raw coefficient
      lam.k <- exp(Xn.k %*% raw.j[colnames(Xn.k)])
      
      #Multiply transformed coefficients by those predictions and take the mean
      vars$est[k] <- mean(lam.k * hf.j[hf.j$name==vars$var[k],]$value)
      
      #Add some extra tracking information
      vars$n[k] <- nrow(rows.k)
      vars$original[k] <- exp(raw.j[vars$var[k]])
      
    }
    
    #12. Cap open habitat values----
    #Get maximum lambda for open habitat types
    lam.open <- max(lam.j[c(names(lam.j)[endsWith(names(lam.j), "R")],
                            "GrassHerb", "Shrub", "GraminoidFen", "Marsh")])
    
    #Get maximum lambda for open human footprint types
    lam.hf <- max(lam.j[c("Industrial", "Rural", "Urban")])
    
    #Overall cap value
    lam.max <- max(lam.open, lam.hf)
    
    #Cap the linear features----
    linear <- vars |>
      rowwise() |>
      mutate(capped = min(est, lam.max)) |>
      ungroup() |>
      mutate(cap_open = lam.open,
             cap_hf = lam.hf,
             cap_total = lam.max,
             species = spp[i],
             boot = j,
             name = c("Wellsites", "EnSoftLin", "TrSoftLin", "EnSeismic"))
    
    linear.j <- linear$capped
    names(linear.j) <- linear$name
    
    #13. Put together----
    lam.out <- c(Climate = exp(raw.j["climate"]),
                 lam.j[names(lam.j)!="Mine"],
                 linear.j,
                 HardLin = 0,
                 Water = 0,
                 Bare = 0,
                 SnowIce = 0,
                 Mine = 0,
                 MineV = unname(lam.j["Mine"]))
    names(lam.out) <- gsub("Spruce", "WhiteSpruce", names(lam.out))
    names(lam.out) <- gsub("Decid", "Deciduous", names(lam.out))
    names(lam.out) <- gsub("Climate.climate", "Climate", names(lam.out))
    names(lam.out) <- gsub("TreedBog", "BlackSpruce", names(lam.out))
    
    #14. Transform back and cap values----
    lam.final <- log(lam.out)
    lam.final[lam.final > 10^4] <- 10^4
    lam.final[lam.final < -10^4] <- -10^4
    
    coef.list[[j]] <- lam.final
    
  }
  
  #15. Make the array if first species----
  if(i==1){
    joint.n <- array(0, c(length(spp), length(coef.list[[j]]), nrow(north.i)))
    dimnames(joint.n) <- list(spp, names(coef.list[[j]]), paste0("b", c(1:nrow(north.i))))
  }
  
  #16. Add species to array----
  joint.n[spp[i],,] <- as.matrix(do.call(cbind, coef.list))
  
  cat("Finished species", i, "of", length(spp), "\n")
  
}

#17. Truncate NESP climate coef because it blows up the relative abundance prediction----
joint.n["NESP","Climate",] <- joint.n["NESP","Climate",]*0.3

#LANDCOVER MODEL COEFFICIENTS - SOUTH#############

#1. Get the list of north models----
south <- mods |>
  dplyr::filter(region=="south")

#2. Get list of species----
spp <- unique(south$species)

#3. Get the full list of potential continuous covariates----
#don't need the categorical ones because they're in every model
names.south <- do.call(rbind, unlist(modelssouth)) |>
  data.frame() |>
  separate(X3, into=c("X3", "X4", "X5", "X6", "X7", "X8", "X9"), sep = " +") |>
  dplyr::select(-X1, -X2, -X3) |>
  pivot_longer(cols=X4:X9, names_to="position", values_to="name") |>
  dplyr::filter(!is.na(name), name!="+") |>
  dplyr::select(-position) |>
  unique() |>
  mutate(name = fix_names(name)) |>
  dplyr::filter(!name %in% c("method", "soilc"))

#4. Get the full south data----
covs.s <- dplyr::filter(covs, useSouth==TRUE, !is.na(soilc))

#5. Get the model matrix----
#use formula from north section above
Xs <- model.matrix(get_terms(modelssouth), covs.s)
colnames(Xs) <- fix_names(c("Intercept", colnames(Xs)[-1]))

#6. Set up loop----
for(i in 1:length(spp)){
  
  #7. List of models for that species----
  south.i <- dplyr::filter(south, species==spp[i])
  
  coef.list <- list()
  for(j in 1:nrow(south.i)){
    
    #9. Get the raw coefficients-----
    #Add zeros for any covariates that weren't in the best model
    coef.j <- read.csv(south.i$path[j]) |>
      mutate(name = fix_names(name)) |>
      full_join(data.frame(name = names.south$name)) |>
      mutate(value = ifelse(is.na(value), 0, value)) |>
      suppressMessages()
    raw.j <- coef.j$value
    names(raw.j) <- coef.j$name
    
    #10. Transform the soil estimates and handle intercept----
    names.soil <- c("soilcBlowout", "soilcClaySub", "soilcCrop",
                    "soilcIndustrial", "soilcMine", "soilcOther", "soilcRapidDrain",
                    "soilcRoughP", "soilcRural", "soilcSandyLoam",  "soilcSoilUnknown",
                    "soilcTameP", "soilcThinBreak", "soilcUrban", "soilcWater", "soilcWellsites")
    
    lam.soil <- exp(c(raw.j[1], raw.j[1] + raw.j[names.soil]))
    names(lam.soil) <- levels(covs.s$soilc)
    
    #11. Adjust the linear feature and well coefficients----
    
    #Human-modified landcover types that we don't want interfering with adjustment
    hfc <- c("Crop", "Industrial", "Mine", "RoughP", "Rural", "TameP", "Urban")
    
    #Dataframe of variables to adjust
    vars <- data.frame(var = c("mWell", "mSoft"))
    
    for(k in 1:nrow(vars)){
      
      #Make a mock dataframe for the variable
      zero.k <- data.frame(var = 0)
      colnames(zero.k) <- vars$var[k]
      
      #Get the rows of interest from the model matrix
      #> 0 proportion of variable, no human-modified landover type, no harvest
      rows.k <- covs.s |>
        mutate(rowid = row_number()) |>
        anti_join(zero.k) |>
        dplyr::filter(!vegc %in% hfc) |>
        suppressMessages()
      
      Xs.k <- Xs[rows.k$rowid,]
      
      #Make predictions from those rows using the raw coefficient
      lam.k <- exp(Xs.k %*% raw.j[colnames(Xs)])
      
      #Multiply transformed coefficient by those predictions and take the mean
      vars$est[k] <- mean(lam.k * exp(raw.j[vars$var[k]]))
      
      #Add some extra tracking information
      vars$n[k] <- nrow(rows.k)
      vars$original[k] <- exp(raw.j[vars$var[k]])
      
    }
    
    #12. Cap open habitat values----
    #Get maximum lambda for open habitat types
    lam.open <- max(lam.soil[c("Loamy", "Blowout", "ClaySub", "RapidDrain", "SandyLoam", "ThinBreak", "Other")])
    
    #Get maximum lambda for open human footprint types
    lam.hf <- max(lam.soil[c("Industrial", "Rural", "Urban")])
    
    #Overall cap value
    lam.max <- max(lam.open, lam.hf)
    
    #Cap the linear features----
    linear <- vars |>
      rowwise() |>
      mutate(capped = min(est, lam.max)) |>
      ungroup() |>
      mutate(cap_open = lam.open,
             cap_hf = lam.hf,
             cap_total = lam.max,
             species = spp[i],
             boot = j,
             name = c("Wellsites", "EnSoft"))
    
    linear.j <- linear$capped
    names(linear.j) <- linear$name
    
    #13. Put together----
    lam.out <- c(Climate = exp(raw.j["climate"]),
                 lam.soil[!names(lam.soil) %in% c("Mine", "Water", "Wellsites")],
                 pAspen = exp(unname(raw.j["paspen"])),
                 Wellsites = linear$capped[1],
                 EnSeismic = linear$capped[2],
                 EnSoftLin = linear$capped[2],
                 TrSoftLin = linear$capped[2],
                 HardLin = 0,
                 Water = 0,
                 Mine = 0,
                 MineV = unname(lam.soil["Mine"]))
    names(lam.out) <- gsub("Climate.climate", "Climate", names(lam.out))
    
    #14. Transform back and cap values----
    lam.final <- log(lam.out)
    lam.final[lam.final > 10^4] <- 10^4
    lam.final[lam.final < -10^4] <- -10^4
    
    coef.list[[j]] <- t(lam.final)
    
  }
  
  #15. Make the array if first species----
  if(i==1){
    joint.s <- array(0, c(length(spp), length(coef.list[[j]]), nrow(south.i)))
    dimnames(joint.s) <- list(spp, names(lam.final), paste0("b", c(1:nrow(south.i))))
  }
  
  #16. Add species to array----
  joint.s[spp[i],,] <- as.matrix(do.call(cbind, coef.list))
  
  cat("Finished species", i, "of", length(spp), "\n")
  
}

#17. Truncate NESP climate coef because it blows up the relative abundance prediction----
joint.n["NESP","Climate",] <- joint.n["NESP","Climate",]*0.3

#MEDIAN BOOTSTRAP SELECTION########

#1. List of species*region----
regions <- mods |>
  dplyr::select(species, region) |>
  unique()

todo <- regions |>
  dplyr::select(species) |>
  unique()

#2. Set up loop----
boots <- data.frame()
for(i in 1:nrow(todo)){
  
  #3. Settings----
  spp.i <- todo$species[i]
  regions.i <- regions |>
    dplyr::filter(species==spp.i)
  region.i <- ifelse(nrow(regions.i)==2, "both", unique(regions.i$region))
  
  #4. Get the right matrix----
  if("north" %in% regions.i$region){
    north.i <- t(as.matrix(joint.n[spp.i,,]))
  }
  if("south" %in% regions.i$region){
    south.i <- t(as.matrix(joint.s[spp.i,,]))
  }
  
  if(region.i=="both"){coefs.i <- cbind(north.i, south.i)}
  if(region.i=="north"){coefs.i <- north.i}
  if(region.i=="south"){coefs.i <- south.i}
  
  #5. Run an NMDS----
  nmds.i <- metaMDS(coefs.i, k=2, distance="euclidean", trace=0)
  
  #6. Get the site scores----
  scores.i <- scores(nmds.i, display="sites")
  
  #7. Calculate centroid----
  centroid.i <- colMeans(scores.i)
  
  #8. Get distances to centroid----
  dist.i <- data.frame(scores.i) |>
    mutate(distance = sqrt((NMDS1-centroid.i[1])^2 + (NMDS2-centroid.i[2])^2)) |>
    rownames_to_column() |>
    rename(bootstrap = rowname) |>
    mutate(bootstrap = as.numeric(str_sub(bootstrap, 2, 3)))
  
  #9. Get closest bootstrap----
  min.i <- dist.i |>
    dplyr::filter(distance==min(distance))
  
  boots <- rbind(boots,
                 min.i |>
                   mutate(species = spp.i,
                          region = region.i))
  
  cat("Finished nmds", i, "of", nrow(todo), "\n")
  
  
}

#MODEL EVALUATION#########

#1. Peter's AUC functions----
simple_roc <- function(labels, scores){
  Labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(
    TPR=cumsum(Labels)/sum(Labels),
    FPR=cumsum(!Labels)/sum(!Labels),
    Labels=Labels)
}

simple_auc <- function(ROC) {
  ROC$inv_spec <- 1-ROC$FPR
  dx <- diff(ROC$inv_spec)
  sum(dx * ROC$TPR[-1]) / sum(dx)
}

#2. Get list of median models----
loop <- data.frame(path = list.files(file.path(root, "Results", "LandcoverModels", "Models"),  full.names = TRUE, pattern="*.Rdata", recursive = TRUE),
                   file = list.files(file.path(root, "Results", "LandcoverModels", "Models"), pattern="*.Rdata", recursive = TRUE)) |> 
  separate(file, into=c("region", "sp1", "region2", "species", "bootstrap", "filetype")) |>
  mutate(bootstrap = as.numeric(bootstrap)) |> 
  dplyr::select(-filetype, -region2, -sp1) |>
  inner_join(rbind(boots |> 
                dplyr::filter(region %in% c("north", "south")),
              boots |> 
                dplyr::filter(region=="both") |> 
                mutate(region="north"),
              boots |> 
                dplyr::filter(region=="both") |> 
                mutate(region="south")),
             by=c("species", "region", "bootstrap"))

#3. Set up loop----
loop$AUC_binary <- NA
loop$AUC_poisson <- NA
for(i in 1:nrow(loop)){
  
  #4. Load the model----
  load.i <-try(load(loop$path[i]))
  if(inherits(load.i, "try-error")){ next }
  
  #5. Make predictions on training data----
  pred.i <- predict(bestmodel, type="response")
  
  #6. AUC for binary response----
  bestmodel$binary <- ifelse(bestmodel$y > 1, 1, bestmodel$y)
  roc.i <- try(roc(response=bestmodel$binary, predictor=pred.i, quiet=TRUE))
  if(inherits(roc.i, "try-error")){loop$AUC_binary[i] <- NA} else {loop$AUC_binary[i] <- auc(roc.i)}
  
  #7. AUC for count response----
  roc.p <- simple_roc(labels=bestmodel$y, scores=pred.i)
  loop$AUC_poisson[i] <- simple_auc(roc.p)
  
  cat("Finished model ", i, "of", nrow(loop), "\n")
  
  
}

#8. Tidy----
auc.out <- loop
rm(loop)

#PACKAGING##########

#1. Change array names from code to common name----

#Get lookup table
birdnames <- read.csv(file.path(root, "Data", "lookups", "birdlist.csv"))

#Match
birdnames.m <- birdnames |>
  arrange(code) |>
  dplyr::filter(code %in% dimnames(marginal)[[1]])

birdnames.n <- birdnames |>
  arrange(code) |>
  dplyr::filter(code %in% dimnames(joint.n)[[1]])

birdnames.s <- birdnames |>
  arrange(code) |>
  dplyr::filter(code %in% dimnames(joint.s)[[1]])

#Replace
dimnames(marginal) <- list(birdnames.m$sppid, dimnames(marginal)[[2]], dimnames(marginal)[[3]])

dimnames(joint.n) <- list(birdnames.n$sppid, dimnames(joint.n)[[2]], dimnames(marginal)[[3]])

dimnames(joint.s) <- list(birdnames.s$sppid, dimnames(joint.s)[[2]], dimnames(marginal)[[3]])

#2. Bird look up table-----

#Get the show details
birds2024 <- read.csv(file.path(root, "Data", "lookups", "birds-v2024.csv")) |>
  rename(species = common) |>
  dplyr::select(species, show)

#Get the occurrence totals
occurrence <- bird |>
  pivot_longer(ALFL:YRWA, names_to="code", values_to="count") |>
  dplyr::filter(count > 0) |>
  group_by(code) |>
  summarize(Occurrences = n()) |>
  ungroup()

#Get the AUC values
auc <- auc.out |>
  rename(code = species) |>
  group_by(code, region) |>
  summarize(AUC = mean(AUC_poisson)) |>
  ungroup() |>
  pivot_wider(names_from="region", values_from=c("AUC"))

#Wrangle the median bootstrap
bootuse <- boots |>
  dplyr::select(species, bootstrap) |>
  rename(code = species)

#set species to remove
spp.remove <- c("VEER", "LEYE")

#set different plot truncation quantile species
q990 <- c("NOFL")
q999 <- c("BLBW", "MOCH", "WCSP", "CAWA")

#Make the table
birdtable <- birdnames |>
  dplyr::filter(sppid %in% dimnames(joint.n)[[1]] | sppid %in% dimnames(joint.s)[[1]]) |>
  left_join(birds2024) |>
  left_join(occurrence) |>
  left_join(auc) |>
  left_join(bootuse) |>
  mutate(ModelNorth = ifelse(show %in% c("c", "n"), TRUE, FALSE),
         ModelSouth = ifelse(show %in% c("c", "s"), TRUE, FALSE),
         LinkHabitat = "log",
         LinkSpclim = "log",
         Group = "birds",
         Nonnative = FALSE,
         Use = ifelse(code %in% spp.remove, FALSE, TRUE),
         PlotQuantile = case_when(code %in% q990 ~ 0.990,
                                  code %in% q999 ~ 0.999,
                                  !is.na(code) & Use==TRUE ~ 0.995)) |>
  rename(SpeciesID = sppid,
         ScientificName = scinam,
         CommonName = species,
         Comments = code) |> 
  rename(AUCNorth = north, AUCSouth = south, Bootstrap = bootstrap) |>
  dplyr::select(SpeciesID, ScientificName, CommonName, ModelNorth, ModelSouth, Occurrences, Nonnative, LinkHabitat, LinkSpclim, AUCNorth, AUCSouth, Bootstrap, Comments, Group, Use, PlotQuantile)

#3. Put together----
birds <- list(
  north = list(marginal = marginal, joint = joint.n),
  south = list(marginal = marginal, joint = joint.s),
  species = birdtable
)

#4. Save----
save(birds, file=file.path(root, "Results", "Birds2024.RData"))
