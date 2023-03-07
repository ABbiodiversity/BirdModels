# ---
# title: ABMI models - package results for allinone
# author: Peter Solymos, updated by Elly Knight
# created: Sep 2020, updated Feb 22 2023
# ---

#NOTES################################


#PREAMBLE############################

#1. Load packages----
library(tidyverse) #basic data wrangling
library(mefa4) #Solymos data wrangling
library(intrval) #For QPAD offsets

#2. Set root path for data on google drive----
root <- "G:/My Drive/ABMI/Projects/BirdModels/"

#3. Load functions----
source("script/00.Functions.R")

#4. Restrict scientific notation----
options(scipen=99999)

#A. LOAD THINGS#####

#1. North dataset----
en <- new.env()
load(file.path(root, "Data", "3Packaged-North.Rdata"), envir=en)

#2. South dataset----
es <- new.env()
load(file.path(root, "Data", "3Packaged-South.Rdata"), envir=es)

#3. Model matrices----
Xn <- get_model_matrix(en$DAT, en$mods)
Xs <- get_model_matrix(es$DAT, es$mods)

#4. Veglookup----
Xage <- as.matrix(read.csv(file.path(root, "Data", "lookups", "Xn-veg-v2020.csv")))
colnames(Xage) <- colnames(Xn)[match(colnames(Xage), make.names(colnames(Xn)))]

#5. Get number of bootstraps----
#From 04A.north.R or 04.B.south.R scripts
B <- 8*32

#6. Full taxa list----
ee <- new.env()
load(file.path(root, "Data", "2Wrangled.Rdata"), envir=ee)
bt <- ee$tax
row.names(bt) <- bt$code
rm(ee)

#7. Model species lists----
SPPn <- substr(list.files(file.path(root, "Results", "CCOutput", "out2023-QPADV3", "north")), 1, 4)
SPPs <- substr(list.files(file.path(root, "Results", "CCOutput", "out2023-QPADV3", "south")), 1, 4)

#Add full common names
names(SPPn) <- bt[SPPn, "sppid"]
names(SPPs) <- bt[SPPs, "sppid"]

#8. Filter to species for modelling----
# c=N+S combo
# n=N model only
# s=S model only
# u=use avail (no model)
# o=exclude (passing through, extinct, bogus)
blist <- read.csv(file.path(root, "Data", "lookups", "birds-v2023.csv"))
compare_sets(blist$common, bt$species)
blist$id <- bt$sppid[match(blist$common, bt$species)]
rownames(blist) <- blist$id

SPPn <- SPPn[names(SPPn) %in% rownames(blist)[blist$show %in% c("c", "n")]]
SPPs <- SPPs[names(SPPs) %in% rownames(blist)[blist$show %in% c("c", "s")]]

#9. Calculate AUC----
#North
AUCNorth <- list()
for (spp in SPPn) {
  cat(spp, "N\n");flush.console()
  resn <- load_species(file.path(root, "Results", "CCOutput", "out2023-QPADV3", "north", paste0(spp, ".RData")))
  yn <- as.numeric(en$YY[,spp])
  off <- if (spp %in% colnames(en$OFF))
    en$OFF[,spp] else en$OFFmean
  lamn <- exp(predict_with_SSH(resn, Xn, en$SSH, stage="Space") + off)
  rocn <- simple_roc(ifelse(yn > 0, 1, 0), rowMeans(lamn))
  aucn <- simple_auc(rocn)
  AUCNorth[[spp]] <- aucn
}
AUCNorth <- unlist(AUCNorth)
save(AUCNorth, file=file.path(root, "Results", "BIRDS-North-QPADV3-AUC.RData"))
load(file.path(root, "Results", "BIRDS-North-QPADV3-AUC.RData"))

#South
AUCSouth <- list()
for (spp in SPPs) {
    cat(spp, "N\n");flush.console()
    ress <- load_species(file.path(root, "Results", "CCOutput", "out2023-QPADV3", "south", paste0(spp, ".RData")))
    ys <- as.numeric(es$YY[,spp])
    off <- if (spp %in% colnames(es$OFF))
        es$OFF[,spp] else es$OFFmean
    lams <- exp(predict_with_SSH(ress, Xs, es$SSH, stage="Space") + off)
    rocs <- simple_roc(ifelse(ys > 0, 1, 0), rowMeans(lams))
    aucs <- simple_auc(rocs)
    AUCSouth[[spp]] <- aucs
}
AUCSouth <- unlist(AUCSouth)
save(AUCSouth, file=file.path(root, "Results", "BIRDS-South-QPADV3-AUC.RData"))
load(file.path(root, "Results", "BIRDS-South-QPADV3-AUC.RData"))

#10. Calculate coefficients---
#North
cfn <- list()
for (spp in SPPn) {
  cat(spp, "\n")
  flush.console()
  
  res <- load_species(file.path(root, "Results", "CCOutput", "out2023-QPADV3", "north", paste0(spp, ".RData")))
  
  est1 <- suppressWarnings(get_coef(res, Xn, stage="ARU", na.out=FALSE))
  est2 <- suppressWarnings(get_coef(res, Xn, stage="Space", na.out=FALSE))
  
  BB <- min(B, nrow(est1))
  cf1 <- sapply(1:BB, function(i) get_coef_north(est1, subset=i))
  cf2 <- sapply(1:BB, function(i) get_coef_north(est2, subset=i))
  
  cfn[[spp]] <- list(estARU=est1, estSpace=est2, coefARU=cf1, coefSpace=cf2)
}
save(cfn, file=file.path(root, "Results", "BIRDS-North-QPADV3-coefs.RData"))
load(file=file.path(root, "Results", "BIRDS-North-QPADV3-coefs.RData"))

#South
cfs <- list()
for (spp in SPPs) {
    cat(spp, "\n")
    flush.console()

    res <- load_species(file.path(root, "Results", "CCOutput", "out2023-QPADV3", "south", paste0(spp, ".RData")))
    
    est1 <- suppressWarnings(get_coef(res, Xs, stage="ARU", na.out=FALSE))
    est2 <- suppressWarnings(get_coef(res, Xs, stage="Space", na.out=FALSE))

    BB <- min(B, nrow(est1))
    cf1 <- sapply(1:BB, function(i) get_coef_south(est1, subset=i))
    cf2 <- sapply(1:BB, function(i) get_coef_south(est2, subset=i))

    cfs[[spp]] <- list(estARU=est1, estSpace=est2, coefARU=cf1, coefSpace=cf2)
}
save(cfs, file=file.path(root, "Results", "BIRDS-South-QPADV3-coefs.RData"))
load(file.path(root, "Results", "BIRDS-South-QPADV3-coefs.RData"))

#11. Align species lists----
TAX <- read.csv(file.path(root, "Data", "lookups", "birdlist.csv"))
row.names(TAX) <- TAX$code

compare_sets(blist$common, TAX$species)
TAX$type <- blist$show[match(TAX$species, blist$common)]

#New species added in 2023 for experimental purposes
#North:YERA RTHA NSWO NOPI HOGR CLNU COLO BADO
#South: CONI RTHA NOPI

SPPn <- names(ALLBIRDSPP$north[match(names(cfn), ALLBIRDSPP$north)])
SPPs <- names(ALLBIRDSPP$south[match(names(cfs), ALLBIRDSPP$south)])
rownames(TAX) <- TAX$sppid
TAX <- TAX[sort(union(SPPn,SPPs)),]
names(cfs) <- SPPs
names(cfn) <- SPPn

#12. Get coefficient names----
#Load previous version of coefficients
library(allinone)
dir <- getOption("allinone")$dir
fn <- file.path(dir, "COEFS.RData")
load(fn)

#Get coeff names
cns <- rownames(COEFS$mites$south[1,,])
cns0 <- cns[1:(which(cns=="Intercept")-2)]
cns1 <- cns[(which(cns=="Intercept")-1):length(cns)]
cnn <- rownames(COEFS$mites$north[1,,])
cnn0 <- cnn[1:(which(cnn=="Intercept")-1)]
cnn1 <- cnn[which(cnn=="Intercept"):length(cnn)]

cnb <- c("pWater_KM", "pWater2_KM", "xPET",
         "xMAT", "xAHM", "xFFP", "xMAP", "xMWMT", "xMCMT", "xY", "xX",
         "xY2", "xX2", "xFFP:xMAP", "xMAP:xPET", "xAHM:xMAT", "xX:xY")

#13. Set bootstrap maximum----
BMAX <- 250

Bs <- min(BMAX, nrow(cfs[[spp]]$estSpace))
Bn <- min(BMAX, nrow(cfn[[spp]]$estSpace))

#14. Organize north coefficients----
set.seed(1234)
for (spp in names(cfn)) {
  
  #Get marginal and joint estimates
  if(nrow(cfn[[spp]]$estARU) >= BMAX){
    cfnm <- rbind(log(cfn[[spp]]$coefARU)[,1:Bn])
    cfnj <- rbind(log(cfn[[spp]]$coefSpace)[,1:Bn], t(cfn[[spp]]$estSpace[1:Bn,cnb]))
  }
  if(nrow(cfn[[spp]]$estARU) < BMAX){
    cfnm <- rbind(log(cfn[[spp]]$coefARU))
    cfnj <- rbind(log(cfn[[spp]]$coefSpace), t(cfn[[spp]]$estSpace[,cnb]))
  }

  #Set reasonable values on estimates
  cfnm[cfnm > 10^4] <- 10^4
  cfnm[cfnm < -10^4] <- -10^4
  cfnj[cfnj > 10^4] <- 10^4
  cfnj[cfnj < -10^4] <- -10^4
  
  #Subsample bootstraps if needed
  if (ncol(cfnm) < BMAX) {
    b <- sample.int(ncol(cfnm), BMAX-ncol(cfnm), replace=TRUE)
    cfnm <- cbind(cfnm, cfnm[,b])
    cfnj <- cbind(cfnj, cfnj[,b])
  }
  
  #Make arrays if is first species
  if(spp==names(cfn)[1]){
    CFnm <- array(0, c(length(cfn), nrow(cfnm), BMAX))
    dimnames(CFnm) <- list(names(cfn), rownames(cfnm), NULL)
    CFnj <- array(0, c(length(cfn), nrow(cfnj), BMAX))
    dimnames(CFnj) <- list(names(cfn), rownames(cfnj), NULL)
  }
  
  #Put together
  CFnm[spp,,] <- cfnm
  CFnj[spp,,] <- cfnj
  
  print(spp)
}

#15. Organize south coefficients----
set.seed(1234)
for (spp in names(cfs)) {
  
  #Get marginal and joint estimates
  cfsm <- rbind(log(cfs[[spp]]$coefARU)[,1:Bs], pAspen=cfs[[spp]]$estARU[1:Bs,"pAspen"])
  cfsj <- rbind(log(cfs[[spp]]$coefSpace)[,1:Bs], t(cfs[[spp]]$estSpace[1:Bs,c("pAspen", cnb)]))
  
  #Set reasonable values on estimates
  cfsm[cfsm > 10^4] <- 10^4
  cfsm[cfsm < -10^4] <- -10^4
  cfsj[cfsj > 10^4] <- 10^4
  cfsj[cfsj < -10^4] <- -10^4
  
  #Subsample bootstraps if needed
  if (ncol(cfsm) < BMAX) {
    b <- sample.int(ncol(cfsm), BMAX-ncol(cfsm), replace=TRUE)
    cfsm <- cbind(cfsm, cfsm[,b])
    cfsj <- cbind(cfsj, cfsj[,b])
  }
  
  #Make arrays if is first species
  if(spp==names(cfs)[1]){
    CFsm <- array(0, c(length(cfs), nrow(cfsm), BMAX))
    dimnames(CFsm) <- list(names(cfs), rownames(cfsm), NULL)
    CFsj <- array(0, c(length(cfs), nrow(cfsj), BMAX))
    dimnames(CFsj) <- list(names(cfs), rownames(cfsj), NULL)
  }
  
  #Put together
  CFsm[spp,,] <- cfsm
  CFsj[spp,,] <- cfsj
  
  print(spp)
}

#16. Compare coeff names----
compare_sets(rownames(cfnm),dimnames(COEFS$lichens$north)[[2]])
compare_sets(rownames(cfsm),dimnames(COEFS$lichens$south)[[2]])

#17. Put together bird look up table----
edat <- new.env()
load(file=file.path(root, "Data", "2Wrangled.Rdata"), envir=edat)
uu <- edat$tax %>% 
  dplyr::filter(!code %in% c("GRAJ", "MCLO"))
yy <- as.matrix(edat$yy)

uu$keep <- as.character(uu$species) %in%
  as.character(blist[blist$show != "o","common"])
uu <- uu[uu$keep,]
uu$show <- blist$show[match(uu$species, blist$common)]
uu <- uu[!is.na(uu$show),]
table(uu$show)

uu$AUCn <- AUCNorth[match(uu$code, names(AUCNorth))]
uu$AUCs <- AUCSouth[match(uu$code, names(AUCSouth))]
uu$pocc <- colSums(yy[,rownames(uu)] > 0) 

#put it together
spb <- data.frame(
  SpeciesID=uu$sppid,
  ScientificName=uu$scinam,
  TSNID=NA,
  CommonName=uu$species,
  ModelNorth=uu$show %in% c("c", "n"),
  ModelSouth=uu$show %in% c("c", "s"),
  UseavailNorth=NA,
  UseavailSouth=NA,
  Occurrences=uu$pocc,
  nSites=NA,
  SizeNorth=NA,
  SizeSouth=NA,
  Nonnative=FALSE,
  LinkHabitat="log",
  LinkSpclim="log",
  AUCNorth=uu$AUCn,
  AUCSouth=uu$AUCs,
  R2North=NA,
  R2South=NA,
  Comments=uu$code,
  Group="birds")

#Change everything to character
for (j in 1:ncol(spb))
  if (is.factor(spb[,j]))
    spb[,j] <- as.character(spb[,j])

#Make the rownames
rownames(spb) <- spb$SpeciesID

#18. Make a bird object----
birds <- list(
  north=list(marginal=CFnm, joint=CFnj),
  south=list(marginal=CFsm, joint=CFsj),
  species=spb)

#19. Some checks----

#Check bird list
spb20 <- COEFS$birds$species
spb23 <- dplyr::filter(spb, ModelNorth==FALSE, ModelSouth==FALSE)
head(COEFS$birds$species, 2)
head(birds$species, 2)

#Check Number of species in species list
nrow(COEFS$birds$species)
nrow(birds$species)
#number of north birds
dim(COEFS$birds$north$joint)[1]
dim(birds$north$joint)[1]
#number of south birds
dim(COEFS$birds$south$joint)[1]
dim(birds$south$joint)[1]

#Check number of coefficients
length(COEFS$birds$north$marginal[1,,1])
length(birds$north$marginal[1,,1])
length(COEFS$birds$north$joint[1,,1])
length(birds$north$joint[1,,1])

#Check coefficients are there
head(COEFS$birds$north$marginal[1,,1])
head(birds$north$marginal[1,,1])
head(COEFS$birds$north$joint[1,,1])
head(birds$north$joint[1,,1])

#Check linear coefficients
COEFS$birds$north$marginal[1,c(90:93),1:10]
birds$north$marginal[1,c(90:93),1:10]
COEFS$birds$north$joint[1,c(90:93),1:10]
birds$north$joint[1,c(90:93),1:10]

COEFS$birds$south$marginal[1,,1:10]
birds$south$marginal[1,,1:10]
COEFS$birds$south$joint[1,,1:10]
birds$south$joint[1,,1:10]

#Check number of bootstraps
length(COEFS$birds$north$marginal[1,1,])
length(birds$north$marginal[1,1,])
length(COEFS$birds$south$marginal[1,1,])
length(birds$south$marginal[1,1,])

# #Fix number of bootstraps
# birds <- list(
#   north=list(marginal=CFnm[,,1:100], joint=CFnj[,,1:100]),
#   south=list(marginal=CFsm[,,1:100], joint=CFsj[,,1:100]),
#   species=spb)
# 
# #Check number of bootstraps again
# length(COEFS$birds$north$marginal[1,1,])
# length(birds$north$marginal[1,1,])
# length(COEFS$birds$south$marginal[1,1,])
# length(birds$south$marginal[1,1,])

#20. Write out----
save(birds, file=file.path(root, "Results", "Birds2023-QPADV3.RData"))
