# ---
# title: ABMI models - wrangle data
# author: Elly Knight
# created: May 10, 2024
# ---

#NOTES################################

#PURPOSE: This script wrangles the gis extraction output into the covariates used in the modelling workflow and calculates qpad offsets.

#PREAMBLE############################

#1. Load packages----
library(tidyverse) #basic data wrangling
library(mefa4) #Solymos data wrangling
library(lubridate) #date wrangling
library(wildrtrax) #offset calculation
library(sf) #coordinate reprojection
library(QPAD) #for offsets

#2. Set root path for data on google drive----
root <- "G:/Shared drives/ABMI_ECKnight/Projects/BirdModels"

#3. Restrict scientific notation----
options(scipen = 99999)

#4. Load harmonized set----
load(file.path(root, "Data", "Harmonized.Rdata"))

#5. Load gis data extractions----

#climate data
load(file.path(root, "Data", "gis", "bird-point-climate_2024.Rdata"))

#point data - load this one before the next one so you get the correct version of d.wide.150 and d.wide.564
load(file.path(root, "Data", "gis", "birds-point-sites_1993-2023.Rdata"))

#buffer data
load(file.path(root, "Data", "gis", "birds-point-sites-simplified_1993-2023.Rdata"))

#tidy up things we don't need
rm(d.long.150, d.long.564)

#6. Load missing uid lookup table----
missingid <- read.csv(file.path(root, "Data", "lookups", "MissingID_lookup.csv")) |>
  dplyr::select(UID, gisid)

#7. Fix the ids----
climate.id <- climate.data |>
  left_join(missingid,
            multiple= "all") |>
  mutate(gisid = ifelse(is.na(gisid), UID, gisid))

d.long.pts.id <- d.long.pts |>
  left_join(missingid,
            multiple= "all") |>
  mutate(gisid = ifelse(is.na(gisid), UID, gisid))

backfill.id <- data.frame(d.wide.150$veg.current) |>
  rownames_to_column(var = "UID") |>
  left_join(missingid,
            multiple= "all") |>
  mutate(gisid = ifelse(is.na(gisid), UID, gisid))

#REMOVE SURVEYS WITH MISSING DATA#########

#1. Get unique IDs that are present in all datasets----
ids <- location |>
  dplyr::select(gisid) |>
  dplyr::filter(gisid %in% climate.id$gisid,
                gisid %in% d.long.pts.id$gisid,
                gisid %in% backfill.id$gisid)

#2. Filter survey object----
#standardize time zone 
#remove surveys with "NONE" transcription method
#remove missing lat/lons
#remove missing climate variables
#remove missing date/time
#remove missing duration or distance method
#remove data out of year range
#remove data out of assumptions for QPAD range
#remove surveys without covariate extraction

complete <- use |>
  dplyr::select(-locationfix, -missingcomments) |> #TAKE THIS OUT ONCE NAMES FIXED
  dplyr::filter(gisid %in% ids$gisid,
                task_method!="None",
                !is.na(latitude),
                !is.na(longitude),
                !is.na(date_time),
                !is.na(duration),
                !is.na(distance),
                year(date_time) >= 1993,
                year(date_time) <= 2024,
                yday(date_time) >= 125,
                yday(date_time) <= 200,
                hour(date_time) >= 3,
                hour(date_time) <= 12) |>
  arrange(gisid, date_time) |>
  mutate(surveyid = row_number())

#CALCULATE OFFSETS####

#1. Load QPAD estimates----
load_BAM_QPAD(3)

#2. Make species list----
#Remove species without enough data
spmin <- 20

spp.use <- complete |>
  dplyr::select(ALFL:YRWA) |>
  mutate_all(~ifelse(.>0, 1, 0)) |>
  pivot_longer(ALFL:YRWA, names_to="species", values_to="detection") |>
  dplyr::filter(detection > 0) |>
  group_by(species) |>
  summarize(detections = n()) |>
  ungroup() |>
  dplyr::filter(detections >= spmin)

bird <- complete |>
  dplyr::select(gisid, all_of(spp.use$species))

#Add GRAJ manually
spp <- sort(c(intersect(getBAMspecieslist(), colnames(bird)), "GRAJ"))

#3. WildTRax login----
source("00.WTlogin.R")
wt_auth()

#4. Wrangle data----
complete.x <- complete |> 
  dplyr::select(surveyid, organization, project_id, location, longitude, latitude, date_time, duration, distance, task_method) |>
  rename(recording_date_time = date_time,
         task_duration = duration,
         task_distance = distance)

#5. Make offsets----
off <- wt_qpad_offsets(complete.x, species=spp, version=3, together=FALSE)

#6. Put together----
offset.calc <- data.frame(surveyid = complete.x$surveyid) |>
  cbind(data.frame(off)) |>
  arrange(surveyid) |>
  rename(CAJA = GRAJ)

#7. Get remaining species----
spp.mean <- colnames(bird |> dplyr::select(-gisid)) |>
  setdiff(c(spp, "CAJA"))

#8. Calculate mean offset for them----
offset.mean <- rowMeans(offset.calc[,-1])
offset.df <- data.frame(replicate(offset.mean, n=length(spp.mean)))
colnames(offset.df) <- spp.mean

#9. Put together----
offset.spp <- cbind(offset.calc |>
                      dplyr::select(-surveyid),
                    offset.df)
offsets <- cbind(offset.calc |>
                   dplyr::select(surveyid),
                 offset.spp[, order(colnames(offset.spp))])

# DERIVE PREDICTORS####

#1. Get lookup tables----
#These tables dictate how the different landcover types from the gis extraction should be handled in the script. They should be reviewed and updated each time the model is run.

#Veg
tv <- read.csv(file.path(root, "Data", "lookups", "lookup-veg-hf-age-birds-v2024.csv"))
row.names(tv) <- tv$ID

#Soil
ts <- read.csv(file.path(root, "Data", "lookups", "lookup-soil-hf-birds-v2024.csv"))
row.names(ts) <- ts$ID

#2. Initial wrangling----
#We manipulate a couple variables and add some column values from the point-level gis extraction and climate extraction
covs1 <- complete |>
  dplyr::select(surveyid, gisid, source:distance) |>
  mutate(year = year(date_time) - min(year(date_time))+1) |>
  left_join(d.long.pts.id |>
              dplyr::select(gisid, NRNAME, NSRNAME),
            multiple = "all") |>
  left_join(climate.id |>
              dplyr::select(gisid, paspen, MAP, TD, CMD, FFP, EMT)) |>
  unique()

#3. Scale and format GIS data----
#We divide by the area of each buffer size to scale from 0 to 1
#We join to survey data frame to link covariates (extracted at the year * location level) to surveys, which can be multiple in one year
vc1 <- data.frame(d.wide.150$veg.current) |>
  mutate_all(~./70680) |>
  rownames_to_column(var = "UID") |>
  left_join(missingid,
            multiple= "all") |>
  mutate(gisid = ifelse(is.na(gisid), UID, gisid)) |>
  inner_join(covs1 |> dplyr::select(gisid, surveyid), multiple="all") |>
  remove_rownames() |>
  column_to_rownames(var = "surveyid") |>
  dplyr::select(-gisid, -UID)

vc2 <- data.frame(d.wide.564$veg.current) |>
  mutate_all(~./999291) |>
  rownames_to_column(var = "UID") |>
  left_join(missingid,
            multiple= "all") |>
  mutate(gisid = ifelse(is.na(gisid), UID, gisid)) |>
  inner_join(covs1 |> dplyr::select(gisid, surveyid), multiple="all") |>
  remove_rownames() |>
  column_to_rownames(var = "surveyid") |>
  dplyr::select(-gisid, -UID)

sc1 <- data.frame(d.wide.150$soil.current) |>
  mutate_all(~./70680) |>
  rownames_to_column(var = "UID") |>
  left_join(missingid,
            multiple= "all") |>
  mutate(gisid = ifelse(is.na(gisid), UID, gisid)) |>
  inner_join(covs1 |> dplyr::select(gisid, surveyid), multiple="all") |>
  remove_rownames() |>
  column_to_rownames(var = "surveyid") |>
  dplyr::select(-gisid, -UID)

sr1 <- data.frame(d.wide.150$soil.reference) |>
  mutate_all(~./70680) |>
  rownames_to_column(var = "UID") |>
  left_join(missingid,
            multiple= "all") |>
  mutate(gisid = ifelse(is.na(gisid), UID, gisid)) |>
  inner_join(covs1 |> dplyr::select(gisid, surveyid), multiple="all") |>
  remove_rownames() |>
  column_to_rownames(var = "surveyid") |>
  dplyr::select(-gisid, -UID)

#4. Calculate water proportions----
covs1$pWater <- rowSums(vc1[,tv[colnames(vc1), "is_water"]])
covs1$pWater_KM <- rowSums(vc2[,tv[colnames(vc2), "is_water"]])
covs1$pWater2_KM <- covs1$pWater_KM^2

#5. Dominant landcover type, weight, & indicator variable----
#This section identifies the landcover type that comprises the largest proportion of the 150 m radius buffer and the proportion of that cover type, and then calculates a weighting factor for use in the glm

#We calculate proportion without some classes because:
# - it is not a stratum we are interested in (water, snow/ice is 0 density),
# - its is too small of a feature to make up a full 7-ha buffer.
#We therefore lose a few surveys here that are entirely water or bare or mostly linear feature

covs2 <- vc1 |>
  rownames_to_column("surveyid") |>
  pivot_longer(DeciduousR:CCPine4, names_to="vegage", values_to="proportion") |>
  left_join(tv |>
              rename(vegage=ID, vegc = UseInAnalysisNoAge) |>
              dplyr::select(vegage, vegc)) |>
  dplyr::filter(!vegage %in% c("Water","HWater", "SnowIce", "Bare", "EnSeismic", "HardLin", "TrSoftLin", "EnSoftLin", "Well"),
                proportion > 0) |>
  group_by(surveyid, vegc) |>
  summarize(vegv = sum(proportion)) |>
  group_by(surveyid) |>
  dplyr::filter(vegv==max(vegv)) |>
  ungroup() |>
  mutate(vegw = pmax(0, pmin(1, 2*vegv-0.5)),
         surveyid = as.numeric(surveyid),
         isMix = ifelse(vegc=="Mixedwood", 1L, 0L),
         isWSpruce = ifelse(vegc=="Spruce", 1L, 0L),
         isPine = ifelse(vegc=="Pine", 1L, 0L),
         isBogFen = ifelse(vegc %in% c("TreedBog", "TreedFen"), 1L, 0L),
         isUpCon = ifelse(vegc %in% c("Spruce", "Pine"), 1L, 0L),
         isCon = ifelse(vegc %in% c("TreedBog", "TreedFen", "Spruce", "Pine"), 1L, 0L)) |>
  left_join(covs1)

#6. Weighted Age----
#This section identifies the proportion of each polygon of the 150m radius of each age class and then calculates a weighted age factor by multipling those proportions by a weighting factor, summing them, and dividing by the total proportion.
#We do this only forested landcover types

# weighting values: exclude unknown (0) and non-forest (blank)
AgeMin <- data.frame(wt = c(0,10,20,40,60,80,100,120,140,160)/200,
                     age = c("R", "1", "2", "3", "4", "5", "6", "7", "8", "9"))

covs3 <- vc1 |>
  rownames_to_column("surveyid") |>
  pivot_longer(DeciduousR:CCPine4, names_to="vegage", values_to="proportion") |>
  left_join(tv |>
              rename(vegage=ID, vegc = UseInAnalysisNoAge, age=AGE) |>
              dplyr::select(vegage, vegc, age)) |>
  dplyr::filter(vegc %in% c("Spruce","Decid","Mixedwood","Pine","TreedBog", "TreedFen"),
                proportion > 0) |>
  group_by(surveyid, age) |>
  summarize(proportion = sum(proportion)) |>
  left_join(AgeMin) |>
  group_by(surveyid) |>
  summarize(wtAge = sum(proportion*wt)/sum(proportion)) |>
  ungroup() |>
  mutate(surveyid = as.numeric(surveyid)) |>
  right_join(covs2) |>
  mutate(wtAge = ifelse(is.na(wtAge), 0, wtAge),
         wtAge2 = wtAge^2,
         wtAge05 = sqrt(wtAge))

#7. Forestry convergence----
#In this section we determine which surveys are predemoninantly harvest within the 150m buffer, and create a variable that accounts for differences in forestry convergence of those areas relative to a natural burn. These trajectories are expert-defined and differ between deciduous and coniferous stands.

#Dave Huggard's recovery trajectories----
age <- c(0, 1:20*4)/200
conif <- 1-c(0, 1.3, 4.7, 10, 17.3, 26, 35.5, 45.3, 54.6, 63.1, 70.7, 77.3,
             82.7, 87, 90.1, 92.3, 94, 95.3, 96.7, 98.2, 100)/100
decid <- 1-c(0, 6.5, 15.1, 25.2, 36.1, 47.2, 57.6, 66.7, 74.3, 80.4, 85,
             88.3, 90.5, 92, 93, 94, 95.1, 96.4, 97.6, 98.8, 100)/100

covs4 <- vc1 |>
  rownames_to_column("surveyid") |>
  pivot_longer(DeciduousR:CCPine4, names_to="vegage", values_to="proportion") |>
  left_join(tv |>
              rename(vegage=ID) |>
              mutate(vegccc = ifelse(is_harvest, paste0("CC", UseInAnalysisNoAge), UseInAnalysisNoAge)) |>
              dplyr::select(vegage, vegccc)) |>
  dplyr::filter(!vegage %in% c("Water","HWater", "SnowIce", "Bare", "EnSeismic", "HardLin", "TrSoftLin", "EnSoftLin", "Well"),
                proportion > 0) |>
  group_by(surveyid, vegccc) |>
  summarize(vegccv = sum(proportion)) |>
  group_by(surveyid) |>
  dplyr::filter(vegccv==max(vegccv)) |>
  ungroup() |>
  mutate(isCC = ifelse(str_sub(vegccc, 1, 2)=="CC", 1L, 0L),
         surveyid = as.numeric(surveyid)) |>
  left_join(covs3) |>
  mutate(fcc2 = case_when(isCC & isCon ~ approxfun(age, conif)(wtAge),
                          isCC & !isCon ~ approxfun(age, decid)(wtAge)),
         fcc2 = ifelse(is.na(fcc2), 0, fcc2))

#Visualize
ggplot(covs4 |> dplyr::filter(isCC==TRUE)) +
  geom_line(aes(x=wtAge, y=fcc2, colour=factor(isCon)))

#8. Dominant soil type----
#In this section, we find the dominant soil type within the 150 m buffer, as per veg type, except that we use the dominant soil type from the backfill layer (reference) when HFor is the dominant class.
#We also calculate weights, as per the veg categories

#Calculate reference values
srsum <- sr1 |>
  rownames_to_column("surveyid") |>
  pivot_longer(ClaySub:Water, names_to="soilrc", values_to="proportion") |>
  dplyr::filter(!soilrc %in% c("SoilWater", "SoilUnknown", "HWater"),
                proportion > 0) |>
  group_by(surveyid, soilrc) |>
  summarize(soilrv = sum(proportion)) |>
  group_by(surveyid) |>
  dplyr::filter(soilrv == max(soilrv)) |>
  ungroup()

covs5 <- sc1 |>
  rownames_to_column("surveyid") |>
  pivot_longer(ClaySub:Wellsites, names_to="soilc", values_to="proportion") |>
  dplyr::filter(!soilc %in% c("Water","HWater", "SnowIce", "Bare", "EnSeismic", "HardLin", "TrSoftLin", "EnSoftLin", "Well"),
                proportion > 0) |>
  group_by(surveyid, soilc) |>
  summarize(soilv = sum(proportion)) |>
  group_by(surveyid) |>
  dplyr::filter(soilv == max(soilv)) |>
  ungroup() |>
  left_join(srsum) |>
  mutate(soilc = ifelse(soilc=="HFor", soilrc, soilc),
         soilv = ifelse(soilc=="HFor", soilrv, soilv),
         soilw = pmax(0, pmin(1, 2*soilv-0.5)),
         surveyid = as.numeric(surveyid)) |>
  inner_join(covs4)


#9. Add nuisance variables, modifier variables, and tidy----
#Surveys that are BBS or surveys with > 0.4 proportion "is_road"
#Fix factor levels for vegc and soilc and method
#Add proportions for the linear features & well sites
#Make hybrid surveys 1SPM

covs6 <- covs5 |>
  mutate(road = case_when(project_id %in% c(1384, 1385) ~ 1,
                          pWater > 0.4 ~ 1,
                          !is.na(pWater) ~ 0),
         method = case_when(source=="eBird" ~ "eBird",
                            task_method=="1SPM Audio/Visual hybrid" ~ "1SPM",
                            !is.na(task_method) ~ task_method)) |>
  left_join(vc1 |>
              rownames_to_column("surveyid") |>
              mutate(surveyid = as.numeric(surveyid)) |>
              dplyr::select(surveyid, EnSoftLin, HardLin, TrSoftLin, EnSeismic, Wellsites)) |>
  mutate(mSoft = EnSeismic + EnSoftLin + TrSoftLin,
         vegc = relevel(as.factor(vegc), "Decid"),
         soilc = relevel(as.factor(soilc), "Loamy"),
         method = relevel(as.factor(method), "PC")) |>
  rename(mWell = Wellsites,
         mEnSft = EnSoftLin,
         mTrSft = TrSoftLin,
         mSeism = EnSeismic)

#10. Identify north and south regions-----
#We convert lat lon to EPSG 3400 for 10TM UTM
covs7 <- covs6 |>
  st_as_sf(coords=c("longitude", "latitude"), crs=4326, remove=FALSE) |>
  st_transform(crs=3400) |>
  st_coordinates() |>
  cbind(covs6) |>
  mutate(useNorth = case_when(pWater > 0.5 ~ FALSE,
                              vegw==0 ~ FALSE,
                              NRNAME=="Grassland" ~ FALSE,
                              !is.na(surveyid) ~ TRUE),
         useSouth = case_when(NRNAME %in% c("Grassland", "Parkland") ~ TRUE,
                              NSRNAME=="Dry Mixedwood" ~ TRUE,
                              latitude > 56.7 ~ FALSE,
                              soilw==0 ~ FALSE,
                              pWater > 0.5 ~ FALSE,
                              soilc=="SoilUnknown" ~ FALSE,
                              !is.na(surveyid) ~ FALSE)) |>
  dplyr::filter(!(useNorth==FALSE & useSouth==FALSE))

#Visualize
ggplot(covs7) +
  geom_point(aes(x=X, y=Y, colour=useNorth))

ggplot(covs7) +
  geom_point(aes(x=X, y=Y, colour=useSouth))

#PACKAGE####

#1. Get final covariates----
covs <- covs7 |>
  dplyr::select(organization, project_id, surveyid, gisid, useNorth, useSouth, X, Y, year, date_time, vegc, soilc, vegw, soilw, wtAge, wtAge2, wtAge05, isCon, isUpCon, isBogFen, isMix, isPine, isWSpruce, fcc2, road, mWell, mSoft, mEnSft, mTrSft, mSeism, method, pWater_KM, pWater2_KM, paspen, MAP, TD, CMD, FFP, EMT) |>
  rename(Easting = X, Northing = Y)

bird <- complete |>
  dplyr::select(surveyid, gisid, all_of(colnames(offsets[,-1]))) |>
  dplyr::filter(surveyid %in% covs$surveyid)

off <- offsets |>
  dplyr::filter(surveyid %in% covs$surveyid)

save(covs, bird, off, file=file.path(root, "Data", "Wrangled.Rdata"))
