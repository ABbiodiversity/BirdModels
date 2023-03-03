# ---
# title: ABMI models - package data for running on compute canada
# author: Peter Solymos, updated by Elly Knight
# created: Sep 23 2020, updated Feb 8 2023
# ---

#NOTES################################
#In this script we filter, transform, mutate, and get some offsets calculated.

#PREAMBLE############################

#1. Load packages----
library(tidyverse) #basic data wrangling
library(mefa4) #Solymos data wrangling
library(QPAD) #For QPAD offsets
library(raster) #For QPAD offsets
library(intrval) #For QPAD offsets
library(maptools) #For QPAD offsets
library(opticut) #For lorenz curve in model

#2. Set root path for data on google drive----
root <- "G:/My Drive/ABMI/Projects/BirdModels/"

#3. Load harmonized set----
load(file.path(root, "Data", "2Wrangled.Rdata"))

#4. Load functions----
source("script/00.Functions.R")

#A. REMOVE MISSING DATA#########
#Logic summary
# - `ROAD`: NAs are expected because we need hard linear amount for non-BBS surveys
# - `X` and `Y`: should be dropped
# - climate: out-of-AB-bound issue, to be dropped

#1. Check NAs first----
#summary of NAs in veg+soil+hf summaries
dd$vshf <- ifelse(rowSums(is.na(vc1)) > 0, NA, 0)
(aa <- data.frame(n_of_NAs=colSums(is.na(dd))))

#2. Drop missing lat/lons----
dd <- droplevels(dd[!is.na(dd$X) & !is.na(dd$vshf),])

#3. Drop missing climate vars----
dd <- dd[!is.na(dd$pAspen) & !is.na(dd$NRNAME),]

#4. Remove data out of year range----
dd <- dd[dd$YEAR %[]% c(1993, 2022),]

#5. Remove missing dataes----
dd <- dd[!is.na(dd$DATE),]

#6. Remove missing duration or distance method----
dd <- dd[!is.na(dd$MAXDUR),]

#7. Constrain date and time to make QPAD expectations----
dd$JULIAN <- as.POSIXlt(dd$DATE)$yday
dd$start <- as.POSIXlt(dd$DATI)$hour + as.POSIXlt(dd$DATI)$min / 60
keep <- is.na(dd$DATI) # these will be constant phi
keep[dd$JULIAN %[]% c(125, 200)] <- TRUE
keep[dd$start %[]% c(3, 12)] <- TRUE
dd <- droplevels(dd[keep,])

#8. Check NAs again----
(aa <- data.frame(n_of_NAs=colSums(is.na(dd))))

#B. ADD SOME COLUMNS FOR QPAD####

#1. Normalized ordinal day----
dd$JDAY <- dd$JULIAN / 365

#2. Normalized time since local sunrise----
Coor <- as.matrix(dd[,c("X", "Y")])
JL <- as.POSIXct(dd$DATI, tz="America/Edmonton")
subset <- rowSums(is.na(Coor))==0 & !is.na(JL)
sr <- sunriset(Coor[subset,], JL[subset], direction="sunrise", POSIXct.out=FALSE) * 24
dd$srise <- NA
dd$srise[subset] <- sr
dd$TSSR <- (dd$start - dd$srise) / 24

#3. Quadratic terms----
dd$JDAY2 <- dd$JDAY^2
dd$TSSR2 <- dd$TSSR^2

#4. Maximum counting distance in 1m units----
table(dd$MAXDIS)

#C. SUBSET SPECIES DATA####
# Keep species with at least 20 detections

#1. Set threshold----
spp.n <- 20

#2. Remove species below threshold----
yy <- yy[rownames(dd), colSums(yy > 0) >= spp.n]

#D. DERIVE PREDICTORS####

#1. Year relative to start year 1993----
dd$YR <- dd$YEAR - min(dd$YEAR)

#2. Point level intersection (veg+soil+HF)----
dd <- data.frame(dd, vs0[rownames(dd),])

#3. Read lookup tables----
#3a. Vegetation
tv0 <- read.csv(file.path(root, "Data", "lookups", "lookup-veg-hf-age-v61.csv"))
rownames(tv0) <- tv0[,1]
tv <- read.csv(file.path(root, "Data", "lookups", "lookup-veg-hf-age-v2020.csv"))
rownames(tv) <- tv[,1]
tv0 <- tv0[rownames(tv),]
tv <- data.frame(tv, tv0)
tv$UseInAnalysisNoAge <- gsub("CC", "", gsub("[[:digit:]]", "", as.character(tv$UseInAnalysis)))
tv$UseInAnalysisNoAge[endsWith(tv$UseInAnalysisNoAge, "R")] <-
  substr(tv$UseInAnalysisNoAge[endsWith(tv$UseInAnalysisNoAge, "R")], 1,
         nchar(tv$UseInAnalysisNoAge[endsWith(tv$UseInAnalysisNoAge, "R")])-1)

#3b. Soil
ts0 <- read.csv(file.path(root, "Data", "lookups", "lookup-soil-hf-v61.csv"))
rownames(ts0) <- ts0[,1]
ts <- read.csv(file.path(root, "Data", "lookups", "lookup-soil-hf-v2020.csv"))
rownames(ts) <- ts[,1]
ts0 <- ts0[rownames(ts),]
ts <- data.frame(ts, ts0)

#4. Calculate 7 ha (150 m radius) level proportions----
dd$pWater <- rowSums(row_std(vc1[rownames(dd),])[,tv[colnames(vc1), "is_water"]])
dd$pRoad <- rowSums(row_std(vc1[rownames(dd),])[,tv[colnames(vc1), "is_road"]])
dd$pRoadVeg <- rowSums(row_std(vc1[rownames(dd),])[,tv[colnames(vc1), "is_road_veg"]])
dd$pClosed <- rowSums(row_std(vc1[rownames(dd),])[,tv[colnames(vc1), "is_closed"]])
dd$pHarest <- rowSums(row_std(vc1[rownames(dd),])[,tv[colnames(vc1), "is_harvest"]])
dd$pHF <- rowSums(row_std(vc1[rownames(dd),])[,tv[colnames(vc1), "is_HF"]])

#5. Summarize 4 and 2 level land cover classification for QPAD offsets----
lcc4 <- row_std(groupSums(vc1[rownames(dd),], 2, tv[colnames(vc1),"LCC4"]))
tmp <- find_max(lcc4)
dd$LCC4 <- tmp$index
dd$LCC2 <- dd$LCC4
levels(dd$LCC2) <- c("OpenWet", "Forest", "Forest", "OpenWet")

#6. Calculate 1 km2 (564 m radius) level proportions----
dd$pWater_KM <- rowSums(row_std(vc2[rownames(dd),])[,tv[colnames(vc2), "is_water"]])
dd$pWater2_KM <- dd$pWater_KM^2
dd$pWet_KM <- rowSums(row_std(vc2[rownames(dd),])[,tv[colnames(vc2), "is_wet"]])
dd$pWetWater_KM <-dd$pWater_KM + dd$pWet_KM

#7. Placeholder for Surrounding Suitable Habitat (SSH) in 1km2----
dd$SSH_KM <- 0
dd$SSH05_KM <- sqrt(dd$SSH_KM)

#8. Surrounding footprint in 1km2----
#no abandoned or rough pasture here
dd$THF_KM <- rowSums(row_std(vc1[rownames(dd),])[,tv[colnames(vc2), "is_HF"]])
dd$Lin_KM <- rowSums(row_std(vc1[rownames(dd),])[,tv[colnames(vc2), "is_HF"] & tv[colnames(vc2), "is_linear"]])
dd$Cult_KM <- rowSums(row_std(vc1[rownames(dd),])[,c("CultivationCrop","CultivationTamePasture", "HighDensityLivestockOperation")])
dd$Alien_KM <- rowSums(row_std(vc1[rownames(dd),])[,tv[colnames(vc2), "is_HF"] & tv[colnames(vc2), "is_alien"]])
dd$Nonlin_KM <- dd$THF_KM - dd$Lin_KM
dd$Noncult_KM <- dd$THF_KM - dd$Cult_KM
dd$Succ_KM <- dd$THF_KM - dd$Alien_KM
dd$THF2_KM <- dd$THF_KM^2
dd$Succ2_KM <- dd$Succ_KM^2
dd$Alien2_KM <- dd$Alien_KM^2
dd$Noncult2_KM <- dd$Noncult_KM^2
dd$Nonlin2_KM <- dd$Nonlin_KM^2

#9. Transform climate variable and lat/lon----
dd <- data.frame(dd, transform_clim(dd))

#10. Define road----
#based on hard linear threshold where it is missing (non-BBS)
dd$ROAD[is.na(dd$ROAD) & dd$pRoad > 0.04] <- 1
dd$ROAD[is.na(dd$ROAD)] <- 0
table(dd$ROAD, dd$pRoad > 0.04, useNA="a")

#11. Define ARU----
#this version replaces the use of the "CMETHOD" (riverforks vs ARU vs human) parameter in model set #4 with an "SM2" parameter that accounts for the smaller EDR of the SM2 model relative to human observers and other recorder types (see Yip et al. 2017 in ACE-ECO). The differences in availability for detection between humans and ARUs are now directly incorporated into QPAD V4.
dd$SM2 <- ifelse(dd$EQUIP=="SM2", "SM2", ifelse(dd$EQUIP=="unknown", "unknown", "human"))

#E. RECLASSIFY LANDCOVER####

#1. Vegetation----
vc1r <- row_std(groupSums(vc1[rownames(dd),], 2, tv[colnames(vc1),"UseInAnalysisNoAge"]))
vr1r <- row_std(groupSums(vr1[rownames(dd),], 2, tv[colnames(vr1),"UseInAnalysisNoAge"]))
vc2r <- row_std(groupSums(vc2[rownames(dd),], 2, tv[colnames(vc2),"UseInAnalysisNoAge"]))
vr2r <- row_std(groupSums(vr2[rownames(dd),], 2, tv[colnames(vr2),"UseInAnalysisNoAge"]))

#2. Soil----
sc1r <- row_std(groupSums(sc1[rownames(dd),], 2, ts[colnames(sc1),"UseInAnalysis"]))
sr1r <- row_std(groupSums(sr1[rownames(dd),], 2, ts[colnames(sr1),"UseInAnalysis"]))
sc2r <- row_std(groupSums(sc2[rownames(dd),], 2, ts[colnames(sc2),"UseInAnalysis"]))
sr2r <- row_std(groupSums(sr2[rownames(dd),], 2, ts[colnames(sr2),"UseInAnalysis"]))

#3. Reclassify categorical point value----
#(w/o forest age)
dd$vegpt <- as.factor(tv$UseInAnalysisNoAge[match(dd$VEGHFAGEclass, rownames(tv))])
dd$soilpt <- as.factor(ts$UseInAnalysis[match(dd$SOILHFclass, rownames(ts))])

#4. Find dominant land cover type (veg+hf)

#4a. Remove classes we don't use
#Calculate proportion without classes we do not use, because:
# - it is not a stratum we are interested in (water, snow/ice is 0 density),
# - its is too small of a feature to make up a full 7-ha buffer.
tmp <- find_max(vc1r[,colnames(vc1r) %nin% c("Water","HWater", "SnowIce", "Bare",
                                             "EnSeismic", "HardLin", "TrSoftLin", "EnSoftLin", "Well")])
#4b. Keep levels consistent
dd$vegc <- factor(as.character(tmp$index), levels(dd$vegpt))
dd$vegv <- tmp$value

#5. Identify harvest areas----

#5a. Classify harvest area
cc <- paste0(ifelse(tv[colnames(vc1), "is_harvest"], "CC", ""),
             as.character(tv[colnames(vc1), "UseInAnalysisNoAge"]))
tmp <- row_std(groupSums(vc1[rownames(dd),], 2, cc))
tmp <- find_max(tmp[,colnames(tmp) %ni% c("Water","HWater", "SnowIce", "Bare",
                                          "EnSeismic", "HardLin", "TrSoftLin", "EnSoftLin", "Well")])
dd$vegccc <- factor(as.character(tmp$index), c(levels(dd$vegc),
                                               paste0("CC", c(c("Spruce","Decid","Mixedwood","Pine")))))
dd$vegvcc <- tmp$value
table(cc=dd$vegccc, not=dd$vegc)

#5b. Add indicator variable for harvest area
# Weight is a function of proportion:
# - 0 below 0.25 (too small to be considered dominant)
# - 1 above 0.75 (it is large enough to consider it dominant)
# - 0-1 in between
dd$isCC <- startsWith(as.character(dd$vegccc), "CC")
tmp <- dd$vegccc
levels(tmp) <- gsub("CC", "", levels(dd$vegccc))
dd$vegc[dd$isCC] <- tmp[dd$isCC]
dd$vegw[dd$isCC] <- dd$vegvcc[dd$isCC]
dd$vegw <- pmax(0, pmin(1, 2*dd$vegv-0.5))

#5c. Check how often we get the dominant class at the center
a <- table(pt=dd$vegpt,bf=droplevels(dd$vegc))
a <- a[colnames(a),]
sum(diag(a))/sum(a) # pretty good

#5d. Finally: drop unused levels and relevel to have Decid as reference
dd$vegc <- droplevels(dd$vegc)
dd$vegc <- relevel(dd$vegc, "Decid")
dd$vegccc <- droplevels(dd$vegccc)
dd$vegccc <- relevel(dd$vegccc, "Decid")
data.frame(table(dd$vegc))

#6. Indicator variables for stand types----
#don't use age for TreedSwamp
dd$isMix <- ifelse(dd$vegc == "Mixedwood", 1L, 0L)
dd$isWSpruce <- ifelse(dd$vegc == "Spruce", 1L, 0L)
dd$isPine <- ifelse(dd$vegc == "Pine", 1L, 0L)
dd$isBog <- ifelse(dd$vegc == "TreedBog", 1L, 0L)
dd$isFen <- ifelse(dd$vegc == "TreedFen", 1L, 0L)
dd$isBogFen <- ifelse(dd$vegc %in% c("TreedBog", "TreedFen"), 1L, 0L)
dd$isUpCon <- ifelse(dd$vegc %in% c("Spruce", "Pine"), 1L, 0L)
dd$isCon <- ifelse(dd$vegc %in% c("TreedBog", "TreedFen",
                                  "Spruce", "Pine"), 1L, 0L)
#F. WEIGHTED AGE CALCULATION####

#1. Get age classes of each column----
tv <- tv[colnames(vc1),]
ac <- as.character(tv[, "AGE"])
ac[is.na(ac)] <- ""
ac[tv$UseInAnalysisNoAge == "TreedSwamp"] <- ""

#2. Summarize?----
vc1age <- row_std(groupSums(vc1[rownames(dd),], 2, ac))
AgePtCr <- t(vc1age[,c("R", "1", "2", "3", "4", "5", "6", "7", "8", "9")])

#3. Assign year values----
## exclude unknown (0) and non-forest (blank)
AgeMin <- structure(c(0,10,20,40,60,80,100,120,140,160)/200,
                    names=c("R", "1", "2", "3", "4", "5", "6", "7", "8", "9"))

#4. Weight----
dd$wtAge <- colSums(AgePtCr * AgeMin) / colSums(AgePtCr)
dd$wtAge[is.na(dd$wtAge)] <- 0
dd$isFor <- dd$vegc %in% c("Spruce","Decid","Mixedwood","Pine","TreedBog", "TreedFen")
dd$wtAge[!dd$isFor] <- 0

#5. Square and sqrt----
dd$wtAge2 <- dd$wtAge^2
dd$wtAge05 <- sqrt(dd$wtAge)

#G. FORESTRY CONVERGENCE####

#1. fCC1: linear----
MAXFOR <- 50/200 #Is this age when forestry converges with fire???
dd$fCC1 <- 0
dd$fCC1[dd$isCC==1] <- pmax(0, 1 - (dd$isCC * dd$wtAge/MAXFOR)[dd$isCC==1])
plot(fCC1 ~ wtAge, dd[dd$isCC==1,])

#2. fCC2: Dave Huggard's recovery trajectories----
age <- c(0, 1:20*4)/200
conif <- 1-c(0, 1.3, 4.7, 10, 17.3, 26, 35.5, 45.3, 54.6, 63.1, 70.7, 77.3,
             82.7, 87, 90.1, 92.3, 94, 95.3, 96.7, 98.2, 100)/100
decid <- 1-c(0, 6.5, 15.1, 25.2, 36.1, 47.2, 57.6, 66.7, 74.3, 80.4, 85,
             88.3, 90.5, 92, 93, 94, 95.1, 96.4, 97.6, 98.8, 100)/100

dd$fCC2 <- 0
tmp1 <- approxfun(age, decid)(dd$wtAge)
tmp1[is.na(tmp1)] <- 0 # this happens when wtAge > 0.4 (out of range of approx)
ii <- dd$isFor & !dd$isCon & dd$isCC==1
dd$fCC2[ii] <- tmp1[ii]
tmp2 <- approxfun(age, conif)(dd$wtAge)
tmp2[is.na(tmp2)] <- 0 # this happens when wtAge > 0.4 (out of range of approx)
ii <- dd$isCon & dd$isCC==1
dd$fCC2[ii] <- tmp2[ii]

#3. Visualize----
plot(dd$wtAge, dd$fCC1, col=2, pch=".")
points(dd$wtAge, dd$fCC2, col=4, pch=".")
sum(is.na(dd$fCC2))
by(dd$wtAge*200, list(veg=interaction(dd$isCC,dd$vegc,drop=TRUE)), fstat, level=1)

#H. MODIFIER VARIABLES####

#These variables will modify the effects given some adjacent habitat around them.
#Adjacent habitat is the dominant land cover, e.g. deciduous, or crop around a wellpad.

#1. Modifiers used only in the north----
dd$mEnSft <- vc1r[,"EnSoftLin"]
dd$mTrSft <- vc1r[,"TrSoftLin"]
dd$mSeism <- vc1r[,"EnSeismic"]

#2. Modifiers used in the north and the south----
dd$mWell <- vc1r[,"Well"]
dd$mHard <- vc1r[,"HardLin"] # optional, use ROAD instead
dd$mSoft <- dd$mSeism + dd$mEnSft + dd$mTrSft

#I. IDENTIFY NORTH REGION####

#1. Data subset to be used in the north----
dd$useNorth <- TRUE

#2. do not consider the Grassland natural region----
dd$useNorth[dd$NRNAME == "Grassland"] <- FALSE

#3. do not consider sites that are mostly open water----
#makes dominant land cover designation questionable
dd$useNorth[dd$useNorth & (dd$pWater > 0.5 | dd$vegpt == "Water")] <- FALSE

#4. drop observations that were assigned 0 weight---- 
#will not contribute to likelihood
dd$useNorth[dd$vegw == 0] <- FALSE

#5. Surrounding habitat classification----
dd$vegca <- dd$vegc
levels(dd$vegca) <- c(levels(dd$vegca), "SpruceO","DecidO","MixedwoodO","PineO","TreedBogO", "TreedFenO")
ii <- dd$vegc %in% c("Spruce","Pine","TreedBog", "TreedFen") & dd$wtAge*200 >= 80
dd$vegca[ii] <- paste0(as.character(dd$vegc[ii]), "O")
ii <- dd$vegc %in% c("Decid","Mixedwood") & dd$wtAge*200 >= 50
dd$vegca[ii] <- paste0(as.character(dd$vegc[ii]), "O")
table(dd$vegca, dd$isCC)

ao <- paste0(as.character(tv[colnames(vc1), "UseInAnalysisNoAge"]),
             ifelse(tv[colnames(vc1), "MatureOld"], "O", ""))
SSH_veg <- row_std(groupSums(vc2[rownames(dd),], 2, ao))
SSH_veg <- SSH_veg[,colnames(SSH_veg) %ni% c("Water","HWater", "SnowIce", "Bare",
                                             "EnSeismic", "HardLin", "TrSoftLin", "EnSoftLin", "Well")]
compare_sets(levels(dd$vegca), colnames(SSH_veg))

#J. SOIL SOUTH RECLASS####

#1. Dominant land cover type: soil+HF---- 
#similar to the veg counterpart
tmp <- find_max(sc1r[,colnames(sc1r) %ni% c("SoilWater", "SoilUnknown", "HWater",
                                            "EnSeismic", "EnSoftLin", "TrSoftLin", "HardLin", "Well")])
dd$soilc <- factor(as.character(tmp$index), levels(dd$soilpt))
dd$soilv <- tmp$value

#2. Use backfilled soil type where HFor is dominant----
ii <- dd$soilc == "HFor"
tmp2 <- find_max(sr1r[,colnames(sr1r) %ni% c("SoilWater", "SoilUnknown", "HWater")])
dd$soilc[ii] <- tmp2$index[ii]
dd$soilv[ii] <- tmp2$value[ii]
dd$soilw <- pmax(0, pmin(1, 2*dd$soilv-0.5))
a <- table(pt=dd$soilpt,bf=droplevels(dd$soilc))
a <- a[colnames(a),]
sum(diag(a))/sum(a) # pretty good
dd$soilc <- droplevels(dd$soilc)
dd$soilc <- relevel(dd$soilc, "Loamy")

#K. IDENTIFY SOUTH REGION####

#1. Data subset to be used in the south----
dd$useSouth <- FALSE

#2. use the grassland and Parkland natural regions----
dd$useSouth[dd$NRNAME %in% c("Grassland", "Parkland")] <- TRUE

#3. plus the Dry Mixedwood subregion within the Boreal----
dd$useSouth[dd$NSRNAME %in% c("Dry Mixedwood")] <- TRUE

#4. but only below the magical latitude limit of 56.7 degrees----
dd$useSouth[dd$useSouth & dd$Y > 56.7] <- FALSE

#5. do not consider sites that are mostly open water----
#makes dominant land cover designation questionable
dd$useSouth[dd$useSouth & (dd$pWater > 0.5 | dd$vegpt == "Water")] <- FALSE

#6. do not consider sites where we have no soil info----
dd$useSouth[dd$useSouth & sr1r[,"SoilUnknown"] > 0] <- FALSE

#7. drop observations that were assigned 0 weight----
#will not contribute to likelihood
dd$useSouth[dd$soilw == 0] <- FALSE

#8. Surrounding land cover in the south----
SSH_soil <- sc2r[rownames(dd),]

#L. QPAD OFFSETS####

#1. Load version 4 of estimates----
load_BAM_QPAD(version=3)

#2. Set WD to qpad-offsets package----
root.qpad <- "C:/Users/Elly Knight/Documents/BAM/Projects/QPAD/qpad-offsets"

#3. Read raster data----
rlcc <- raster(file.path(root.qpad, "data", "lcc.tif"))
rtree <- raster(file.path(root.qpad, "data", "tree.tif"))
rtz <- raster(file.path(root.qpad, "data", "utcoffset.tif"))
rd1 <- raster(file.path(root.qpad, "data", "seedgrow.tif"))
crs <- proj4string(rtree)

#4. Source functions----
source(file.path(root.qpad, "functions.R"))

#5. Make prediction object---
x <- dd %>% 
  mutate(time=str_sub(as.character(DATI), 12, 16)) %>% 
  rename(date=DATE,
         lon=X,
         lat=Y,
         dur=MAXDUR,
         dis=MAXDIS,
         tagmethod = TAGMETHOD) %>% 
  dplyr::select(time, date, lon, lat, dur, dis, tagmethod) %>% 
  make_x(tz="local")

#6. Replace LCC with backfill-derived values----
x$LCC2 <- dd$LCC2
x$LCC4 <- dd$LCC4
x$TREE <- dd$pClosed

#7. Get species list----
sppp <- intersect(colnames(yy), getBAMspecieslist())

#8. Set up output----
off <- matrix(0, nrow(x), length(sppp))
colnames(off) <- sppp
rownames(off) <- rownames(dd)

#9. Make offsets----
for (spp in sppp) {
  cat(spp, "\n")
  flush.console()
  o <- make_off(spp, x, useMethod="y")
  off[,spp] <- round(o$offset, 4)
}

#10. sanity checks----
(Ra <- apply(off, 2, range))
summary(t(Ra))
which(!is.finite(Ra[1,]))
which(!is.finite(Ra[2,]))
off_mean <- log(rowMeans(exp(off)))

#M. BOOTSTRAP BLOCKING UNITS####
#Spatial blocking units are based on unique locations (SS), we aim for a well balanced data set

#1. Get visit variables----
cn=c("PCODE", "SS", "SSYR", "PKEY", "YEAR", "DATE", "DATI", "MAXDUR",
     "MAXDIS", "CMETHOD", "ROAD", "X", "Y", "NRNAME", "NSRNAME", "LUF_NAME", "useNorth", "useSouth")
ddd=dd[,cn]

#2. Plot----
with(ddd, plot(X, Y, col=factor(NRNAME), pch="."))

#3. Identify bins----
tmp <- nonDuplicated(ddd, SS, TRUE)
cx <- cut(tmp$X, c(-121, -116, -112,-109))
cy <- cut(tmp$Y, c(48, 51, 54, 57, 61))
ct <- cut(tmp$YEAR, c(1992, 2001, 2009, 2013, 2016, 2019, 2022))
table(cy, cx)
table(ct)
ftable(ct, cy, cx)

#4. Classify into bins----
dd$BLOCK_X <- cut(dd$X, c(-121, -116, -112,-109))
dd$BLOCK_Y <- cut(dd$Y, c(48, 51, 54, 57, 61))
dd$BLOCK_T <- cut(dd$YEAR, c(1992, 2001, 2009, 2013, 2016, 2019, 2022))
dd$BLOCK_XY <- interaction(droplevels(dd$BLOCK_X), droplevels(dd$BLOCK_Y), sep="::", drop=TRUE)
dd$BLOCK_XYT <- interaction(dd$BLOCK_XY, dd$BLOCK_T, sep="::", drop=TRUE)
ftable(dd$BLOCK_T, dd$BLOCK_Y, dd$BLOCK_X)

#5. Random quantiles----
#these are also based on SS
set.seed(1)
tmp$RND <- sample.int(100, nrow(tmp), replace=TRUE)
dd$RND <- tmp$RND[match(dd$SS, tmp$SS)]

#N. CREATE MODEL SUBSETS - SOUTH####

#1. Set parameters----

#1a. Minimum sample size
NMIN <- 20

#1b. Number of bootstraps
B <- 256

#2. Read in the models----
source("script/00.models-soil.R")
setdiff(get_terms(mods_soil, "list"), colnames(dd))

#3. Subset covariate object----
#relevant terms only
#use south only
cn2 <- c(cn, get_terms(mods_soil, "list"), "soilw")
DAT <- droplevels(dd[dd$useSouth & dd$RND > 10, cn2])

#4. Subset detection object----
YY <- yy[rownames(DAT),]
YY <- YY[,colSums(YY>0) >= NMIN]

#5. Subset offset object----
OFF <- off[rownames(DAT), intersect(colnames(off), colnames(YY))]
OFFmean <- off_mean[rownames(DAT)]

#6. Subset surrounding habitat object----
mods <- mods_soil
SSH <- SSH_soil[rownames(DAT),]

#7. Select data rows for each bootstrap----
set.seed(1234)
BB <- pbapply::pbsapply(1:B, bfun, DAT$SS, DAT$BLOCK_XYT)

#8. Some checks----
nrow(DAT)
max(BB)
(lu <- length(unique(as.numeric(BB))))
stopifnot(all(BB <= nrow(DAT)))
stopifnot(lu <= nrow(DAT))
nrow(BB)

#9. Test model-----
z <- run_path1(1, "AMRO", mods, CAICalpha=1, wcol="soilw", ssh_class="soilc", ssh_fit="Space")
z$timer
cat("Estimate for", ncol(YY), "species and", B, "runs is", ceiling(unname(ncol(YY)*B*z$timer[3])/(60*60)), "hrs\n")

#10. Save out----
#10a. Save to google drive for archive
save(DAT, YY, OFF, OFFmean, SSH, BB, mods, file=file.path(root, "Data", "3Packaged-South.Rdata"))

#10b. Save to local for compute canada
save(DAT, YY, OFF, OFFmean, SSH, BB, mods, file=file.path("script", "04.ComputeCanada", "data", "3Packaged-South.Rdata"))

#O. CREATE MODEL SUBSETS - NORTH####

rm(DAT, YY, OFF, OFFmean, SSH, BB, mods)

#1. Set parameters----

#1a. Minimum sample size
NMIN <- 20

#1b. Number of bootstraps
B <- 256

#2. Read in the models----
source("script/00.models-veg.R")
setdiff(get_terms(mods_veg, "list"), colnames(dd))

#3. Subset covariate object----
#relevant terms only
#use north only
cn2 <- c(cn, get_terms(mods_veg, "list"), "vegw", "vegca")
DAT <- dd[dd$useNorth & dd$RND > 10, cn2]

#4. Subset detection object----
YY <- yy[rownames(DAT),]
YY <- YY[,colSums(YY>0) >= NMIN]

#5. Subset offset object----
OFF <- off[rownames(DAT), intersect(colnames(off), colnames(YY))]
OFFmean <- off_mean[rownames(DAT)]

#6. Subset surrounding habitat object----
mods <- mods_veg
SSH <- SSH_veg[rownames(DAT),]

set.seed(1234)
BB <- pbapply::pbsapply(1:B, bfun, DAT$SS, DAT$BLOCK_XYT)

#8. Some checks----
nrow(DAT)
max(BB)
(lu <- length(unique(as.numeric(BB))))
stopifnot(all(BB <= nrow(DAT)))
stopifnot(lu <= nrow(DAT))
nrow(BB)

#9. Test model----
z <- run_path1(1, "AMRO", mods, CAICalpha=1, wcol="vegw", ssh_class="vegca", ssh_fit="Space")
z$timer
cat("Estimate for", ncol(YY), "species and", B, "runs is", ceiling(unname(ncol(YY)*B*z$timer[3])/(60*60)), "hrs\n")

#10. Save out-----

#10a. Save to google drive for archive
save(DAT, YY, OFF, OFFmean, SSH, BB, mods, file=file.path(root, "Data", "3Packaged-North.Rdata"))

#10b. Save to local for compute canada
save(DAT, YY, OFF, OFFmean, SSH, BB, mods, file=file.path("script", "04.ComputeCanada", "data", "3Packaged-North.Rdata"))
