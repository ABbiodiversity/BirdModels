# ---
# title: ABMI models - get data
# author: Elly Knight
# created: November 29, 2022
# ---

#NOTES################################

#The "projectInstructions.csv" file is a list of all projects currently in WildTrax should not be used in ABMI models (instructions=="DO NOT USE"). This file should be updated for each iteration of in collaboration with Erin Bayne. Future versions of this spreadsheet can hopefully be derived by a combination of organization and a google form poll for consent from other organizations.

#The "projectInstructions.csv" file also contains information on which ARU projects are processed for a single species or taxa (instructions=="DO NOT USE") and therefore those visits should only be used for models for the appropriate taxa. This file should be updated for each iteration of the national models in collaboration with Erin Bayne. These projects are currently not included in the models.

#There are a handful of projects that are not downloading properly via wildRtrax. An issue is open on this. These projects are listed in the error.log object. These files should be downloaded manually.

#Riverforks data is stored in the ABMI Oracle database and should be downloaded  using dBeaver. Use the UNRESTRICTED_ACCESS acount to ensure full retrieval of species at risk records. Contact Joan (qfang@ualberta.ca) for access information. This dataset will hopefully be incorporated into WildTrax in the future to avoid this step.

#raw eBird data is downloaded from the eBird interface at https://ebird.org/data/download/ebd prior to wrangling with the "auk" package and will require a request for access. Use the custom download tool to download only the datasets for Canada and the US instead of the global dataset. Note you will also need the global sampling file to use the auk package for zero filling.

#raw eBird data omits Great Grey Owl & Northern Hawk Owl as sensitive species (https://support.ebird.org/en/support/solutions/articles/48000803210?b_id=1928&_gl=1*xq054u*_ga*ODczMTUyMjcuMTY2OTE0MDI4Ng..*_ga_QR4NVXZ8BM*MTY2OTE0MDI4NS4xLjEuMTY2OTE0MDM3OC4zNS4wLjA.&_ga=2.147122167.150058226.1669140286-87315227.1669140286) and should not be used for modelling these two species.

#wrangling eBird data with the auk package requires installation of AWK on windows computers. Please see #https://cornelllabofornithology.github.io/auk/articles/auk.html.

#eBird data has not been zerofilled because there was no species filtering done and we are assuming that all stationary counts have at least 1 bird observed.

#The replace TMTTs script will be replaced by a wildRtrax function in the near future.

#The location csv output of this script (line 582) should be sent to Eric for GIS covariate extraction. The output Eric provides will be used as input in the next script.

#PREAMBLE############################

#1. Load packages----

library(tidyverse) #basic data wrangling
library(wildRtrax) #to download data from wildtrax
library(data.table) #for binding lists into dataframes
library(lubridate) #date wrangling
library(auk) #eBird wrangling
library(downloader) #download zipped provincial boundary shp
library(sf) #shapefile wrangling
library(terra) #raster wrangling

#2. Set root path for data on google drive----

root <- "G:/My Drive/ABMI/Projects/BirdModels/"

#3. Login to WildTrax----
config <- "script/00.WTlogin.R"
source(config)

#4. Authenticate----
wt_auth()

#A. DOWNLOAD DATA FROM WILDTRAX#######################

#1. Get list of projects from WildTrax----
#sensor = PC gives all ARU and point count projects
projects <- wt_get_download_summary(sensor_id = 'PC')

#2. Filter out projects that shouldn't be used----
#nothing in BU training & all "DO NOT USE" projects in projectInventory file
#filter out 'NONE' method ARU projects later after this field is parsed out

instructions <- read.csv(file.path(root, "Data", "projectInventory", "projectInstructions.csv"))

projects.use <- projects %>% 
  left_join(instructions) %>%
  mutate(instruction = ifelse(is.na(instruction), "NO RESTRICTION", instruction)) %>% 
  dplyr::filter(organization!="BU-TRAINING",
                !instruction %in% c("DO NOT USE", "SINGLE SPECIES"))

#3. Loop through projects to download data----
aru.list <- list()
pc.list <- list()
error.log <- data.frame()
for(i in 1:nrow(projects.use)){
  
  #authenticate each time because this loop takes forever
  wt_auth()
  
  #Do each sensor type separately because the reports have different columns and we need different things for each sensor type
  if(projects$sensor[i]=="ARU"){
    
    dat.try <- try(wt_download_report(project_id = projects.use$project_id[i], sensor_id = projects.use$sensor[i], weather_cols = F, report = "main"))
    
    if(class(dat.try)=="data.frame"){
      aru.list[[i]] <- dat.try
    }
    
  }
  
  if(projects$sensor[i]=="PC"){
    
    dat.try <- try(wt_download_report(project_id = projects.use$project_id[i], sensor_id = projects.use$sensor[i], weather_cols = F, report="main"))
    
    if(class(dat.try)=="data.frame"){
      pc.list[[i]] <- dat.try
    }
    
  }
  
  #Log projects that error
  if(class(dat.try)!="data.frame"){
    error.log <- rbind(error.log, 
                       projects.use[i,])
    
  }
  
  print(paste0("Finished dataset ", projects.use$project[i], " : ", i, " of ", nrow(projects.use), " projects"))
  
}

#4. Go download error projects from wildtrax.ca----
error.log %>% 
  inner_join(projects.use %>% 
              mutate(i = row_number())) %>% 
  View()

#5. Read in error projects----
error.files.aru <- list.files(file.path(root, "Data", "WildTrax", "errorFiles", "ARU"), full.names = TRUE)

aru.error <- data.frame()
for(i in 1:length(error.files.aru)){
  aru.error <- read.csv(error.files.aru[i]) %>% 
    rbind(aru.error)
}

error.files.pc <- list.files(file.path(root, "Data", "WildTrax", "errorFiles", "PC"), full.names = TRUE)

pc.error <- data.frame()
for(i in 1:length(error.files.pc)){
  pc.error <- read.csv(error.files.pc[i]) %>% 
    rbind(pc.error)
}

#6. Collapse lists----
#Take out ARU projects with "None" method
aru.wt <- rbindlist(aru.list[], fill=TRUE)  %>% 
  rbind(aru.error, fill=TRUE) %>% 
  left_join(projects %>% 
              dplyr::rename(project_status = status)) %>% 
  dplyr::filter(task_method!="None")

pc.wt <- rbindlist(pc.list[], fill=TRUE)  %>% 
  rbind(pc.error, fill=TRUE) %>% 
  left_join(projects %>% 
              dplyr::rename(project_status = status))

#7. Save date stamped data & project list----
save(aru.wt, pc.wt, projects.use, error.log, file=paste0(root, "/Data/WildTrax/wildtrax_raw_", Sys.Date(), ".Rdata"))

#B. GET EBIRD DATA##########################

#1. Set ebd path----
auk_set_ebd_path(file.path(root, "Data/ebd/ebd_CA-AB_relOct-2022"), overwrite=TRUE)

#2. Define filters----
filters <- auk_ebd(file="ebd_CA-AB_relOct-2022.txt") %>% 
  auk_protocol("Stationary") %>% 
  auk_duration(c(0, 10)) %>% 
  auk_complete()

#3. Filter data----
#select columns to keep
filtered <- auk_filter(filters, file=file.path(root, "Data/ebd_data_filtered.txt"), overwrite=TRUE,
                       keep = c("group identifier", "sampling_event_identifier", "scientific name", "common_name", "observation_count", "latitude", "longitude", "locality_type", "observation_date", "time_observations_started", "observer_id", "duration_minutes"))

#C. HARMONIZE###############################

#1. Get list of bird species----
library(QPAD)
load_BAM_QPAD(3)
spp <- QPAD::getBAMspecieslist()

#2. Set desired columns----
colnms <- c("source", "organization", "project", "sensor", "tagMethod", "equipment", "location", "buffer", "lat", "lon", "year", "date", "observer", "duration", "distance", "species", "abundance", "isSeen", "isHeard")

#3. Load WildTrax data-----
load(file.path(root, "Data", "WildTrax", "wildtrax_raw_2023-11-21.Rdata"))

#4. Wrangle WT ARU data----
#To do filter to QPAD list
#Fix replace_tmtt error

aru.use <- aru.wt %>% 
  wt_tidy_species() %>% 
  wt_replace_tmtt() %>% 
  wt_make_wide() %>% 
  mutate(source="WildTrax") %>% 
  mutate(buffer = as.numeric(ifelse(!buffer %in% c("50", "10000"), str_sub(buffer, -100, -2), buffer)),
         buffer = ifelse(is.na(buffer), 0, buffer),
         location = case_when(str_sub(location, 1, 4)=="1577" ~ paste0("1577B", str_sub(location, 5, 100)),
                              str_sub(location, 1, 7)=="OG-1600" ~ paste0("OG-ABMI-1600", str_sub(location, 8, 100)),
                              !is.na(location) ~ location))

#5. Wrangle WT PC data----
#wrangle distance and duration maximums
#remove counts with unknown duration and distance
pc.use <- pc.wt %>% 
  dplyr::filter(durationMethod!="UNKNOWN",
                distanceMethod!="UNKNOWN") %>% 
  wt_tidy_species() %>% 
  wt_make_wide() %>% 
  mutate(source="WildTrax") %>% 

  
  
  
  
  
  
  
  
  
  




#2b. Get the point count data----

pc.wt.meth <- dat.wt %>% 
  dplyr::filter(sensor=="PC") %>% 
  dplyr::select(durationMethod, distanceMethod) %>% 
  unique() %>% 
  rowwise() %>% 
  mutate(durationMethod = ifelse(str_sub(durationMethod, -1, -1)=="+", str_sub(durationMethod, -100, -2), durationMethod),
         chardur = str_locate_all(durationMethod, "-"),
         chardurmax = max(chardur),
         duration = as.numeric(str_sub(durationMethod, chardurmax+1, -4)),
         chardis = str_locate_all(distanceMethod, "-"),
         chardismax = max(chardis),
         distance1 = str_sub(distanceMethod, chardismax+1, -2),
         distance = ifelse(distance1 %in% c("AR", "IN"), Inf, as.numeric(distance1))) %>% 
  dplyr::select(distanceMethod, durationMethod, distance, duration)

pc.wt <- dat.wt %>% 
  dplyr::filter(sensor=="PC") %>% 
  mutate(tagMethod = ifelse(distanceMethod=="0m-INF-ARU", "1SPT", "PC"),
         equipment = "Human") %>% 
  left_join(pc.wt.meth) %>% 
  dplyr::select(all_of(colnms)) %>% 
  data.frame()

#2c. Get the aru data----
#filter ARU data to first detection of each individual
#wrangle duration
#remove 'NONE' method
#wrangle equipment type
aru.wt.equip <- dat.wt %>% 
  dplyr::filter(sensor=="ARU",
                !is.na(equipment)) %>% 
  dplyr::select(equipment) %>% 
  unique() %>% 
  rowwise() %>% 
  mutate(equipmentend = str_locate(equipment, "\\("),
         equipmentlong = str_sub(equipment, 1, equipmentend[,1]-1),
         equipmentshort = case_when(equipmentlong=="?+" ~ "unknown",
                                    str_detect(equipmentlong, "SM2")==TRUE ~ "SM2",
                                    str_detect(equipmentlong, "SM3")==TRUE ~ "SM3",
                                    str_detect(equipmentlong, "SM4")==TRUE ~ "SM4",
                                    str_detect(equipmentlong, "mini")==TRUE ~ "mini",
                                    str_detect(equipmentlong, "Mini")==TRUE ~ "mini",
                                    is.na(equipmentlong) ~ "unknown"))

aru.wt <- dat.wt %>% 
  dplyr::filter(sensor=="ARU") %>% 
  separate(method, into=c("duration", "tagMethod"), remove=TRUE) %>% 
  dplyr::filter(tagMethod %in% c("1SPM", "1SPT")) %>% 
  mutate(duration = as.numeric(str_sub(duration, -100, -2))/60,
         distance = Inf) %>% 
  group_by(source, organization, project, sensor, tagMethod, equipment, location, buffer, lat, lon, year, date, observer, duration, distance, species, abundance, isSeen, isHeard, individual_appearance_order) %>%
  mutate(first_tag = min(tag_start_s)) %>%
  ungroup() %>%
  dplyr::filter(tag_start_s == first_tag) %>% 
  left_join(aru.wt.equip) %>% 
  dplyr::mutate(equipment = ifelse(is.na(equipmentshort), "unknown", equipmentshort),
                isSeen = "f",
                isHeard = "t") %>% 
  dplyr::select(all_of(colnms))

#2d. Replace TMTTs with predicted abundance----
tmtt <- read.csv("C:/Users/elly/Documents/ABMI/WildTrax/TMTT/data/tmtt_predictions_mean.csv") %>% 
  rename(species = species_code)

tmtt.wt <- aru.wt %>% 
  dplyr::filter(abundance=="TMTT") %>%
  mutate(species = ifelse(species %in% tmtt$species, species, "species"),
         observer_id = as.integer(ifelse(observer %in% tmtt$observer_id, observer, 0))) %>% 
  data.frame() %>% 
  left_join(tmtt) %>% 
  mutate(abundance = round(pred)) %>% 
  dplyr::select(colnames(aru.wt))

#2e. Put back together----
#summarize abundance
use.wt <- aru.wt %>% 
  dplyr::filter(abundance!="TMTT") %>% 
  rbind(tmtt.wt) %>% 
  rbind(pc.wt)
  
#3. Wrangle Riverforks data----

#3a. Read it in-----
det.rf <- read.csv(file.path(root, "Data", "Riverforks", "M_RT_BIRD_COUNT_202301161601.csv"), header=TRUE)
loc.rf <- read.csv(file.path(root, "Data", "Riverforks", "A_RT_SITE_PHYCHAR_202301161554.csv"), header=TRUE)

#3b. Put together----
#concatenate site name with point count # (there are 9 per site) before joining
raw.rf <- det.rf %>%
  mutate(SITE = paste0(SITE, "-", TBB_POINT_COUNT)) %>% 
  left_join(loc.rf %>% 
              mutate(SITE = paste0(SITE, "-", TSFG_POINT_COUNT)) %>% 
              dplyr::select(ROTATION, SITE, YEAR, SITE_LATITUDE, SITE_LONGITUDE) %>% 
              unique())

#3c. Get the taxonomy lookup----
tax.wt <- read.csv(file.path(root, "Data", "lookups", "lu_species.csv")) %>% 
  mutate(SCIENTIFIC_NAME = paste(species_genus, species_name)) %>% 
  rename(species = species_code, COMMON_NAME = species_common_name) %>% 
  dplyr::select(COMMON_NAME, SCIENTIFIC_NAME, species) %>% 
  unique() %>% 
  dplyr::filter(nchar(species)==4)

#3d. Harmonize----
use.rf <- raw.rf %>% 
  rename(location=SITE, lat = SITE_LATITUDE, lon = SITE_LONGITUDE, year = YEAR, observer = ANALYST) %>% 
  mutate(source = "riverforks",
         organization="ABMI",
         project="ABMI-Riverforks",
         sensor="ARU",
         tagMethod = "1SPM",
         equipment="riverforks",
         buffer = 5500,
         duration = 3,
         distance = Inf,
         abundance = 1,
         date = ymd_hm(paste0(ADATE, " ", TBB_START_TIME)),
         isSeen = "f",
         isHeard = "t") %>% 
  left_join(tax.wt) %>% 
  dplyr::select(all_of(colnms))

#4. Wrangle ebird data----
#Note this assumes observations with "X" individuals are 1s
#Filter out hotspots
#Replace common name with alpha code

raw.ebd <- read_ebd(file.path(root, "Data", "ebd", "ebd_data_filtered.txt"))

tax.wt <- read.csv(file.path(root, "Data", "lookups", "lu_species.csv")) %>% 
  mutate(scientific_name = paste(species_genus, species_name)) %>% 
  rename(species = species_code, common_name = species_common_name) %>% 
  dplyr::select(scientific_name, common_name, species) %>% 
  unique() %>% 
  dplyr::filter(nchar(species)==4)

use.ebd <- raw.ebd %>% 
  dplyr::filter(locality_type!="H") %>% 
  mutate(source = "eBird",
         organization = "eBird",
         project="eBird",
         sensor="PC",
         tagMethod="PC",
         equipment="human",
         singlesp="n",
         buffer=0,
         date = ymd_hms(paste0(observation_date, time_observations_started)),
         year = year(date),
         distance = Inf,
         abundance = as.numeric(ifelse(observation_count=="X", 1, observation_count)),
         isSeen = NA,
         isHeard = NA) %>% 
  rename(lat = latitude,
         lon = longitude,
         observer = observer_id,
         duration = duration_minutes,
         location = checklist_id) %>% 
  left_join(tax.wt) %>% 
  dplyr::select(all_of(colnms))

#E. PUT TOGETHER############################

#1. Put everything together----
use <- rbind(use.wt, use.rf, use.ebd)

#2. Clip by provincial boundaries & filter to AB----
use.visit <- use %>% 
  dplyr::select(-species, -abundance, -isSeen, -isHeard) %>% 
  unique() %>% 
  dplyr::filter(!is.na(lon),
                !is.na(lat))

#Download provincial boundaries shapefile
temp <- tempfile()
download("https://www12.statcan.gc.ca/census-recensement/2021/geo/sip-pis/boundary-limites/files-fichiers/lpr_000b21a_e.zip", temp)
unzip(zipfile=temp, exdir=file.path(root, "Data", "gis"))
unlink(temp)

#Filter to Alberta
shp <- read_sf(file.path(root, "Data", "gis", "lpr_000b21a_e.shp")) %>% 
  dplyr::filter(PRNAME=="Alberta") %>% 
  vect()

#Create rasters (much faster than from polygon)
r <- rast(ext(shp), resolution=1000, crs=crs(shp))
ab <- rasterize(x=shp, y=r, field="PRNAME")

#Extract raster value
visit.ab <- use.visit %>% 
  st_as_sf(coords=c("lon", "lat"), crs=4326) %>% 
  st_transform(crs=crs(ab)) %>% 
  vect() %>% 
  extract(x=ab) %>% 
  cbind(use.visit) %>% 
  dplyr::filter(PRNAME=="Alberta") %>% 
  dplyr::select(colnames(use.visit))

#apply to bird data and clip again by shp for precision
use.ab <- use %>% 
  inner_join(visit.ab) %>% 
  st_intersection(shp)

#3. Remove duplicate surveys----

#3a. Investigate----
#Take out buffered locations
visit.dup <- visit.ab %>% 
  mutate(latr = round(lat, 4),
         lonr = round(lon, 4)) %>% 
  group_by(date, latr, lonr) %>% 
  summarize(n = n()) %>% 
  ungroup() %>% 
  dplyr::filter(n > 1) %>% 
  left_join(visit.ab %>% 
              mutate(latr = round(lat, 4),
                     lonr = round(lon, 4)))
#no eBird duplicates
#some point count data entry errors, but really only 2 sources of duplicates:
#1. JOSM points that have human & ARU data & eBird data
#2. Duplicate processing of recordings

#3b. Identify JOSM points with human & ARU data----
#Keep the ARU ones because frankly they're probably better data
visit.josm <- visit.dup %>% 
  dplyr::filter(project %in% c("JOSM ECCC Cause and Effect Monitoring for Landbirds 2013",
                               "Oil Sands Monitoring CWS Prairie Region 2012",
                               "Oil Sands Monitoring CWS Prairie Region 2013",
                               "Oilsands Monitoring CWS Prairie Region 2014"))

visit.josm.remove <- visit.josm %>% 
  dplyr::filter(sensor=="PC")

#3c. Identify recordings that have been processed twice----
#Use recording with longer duration processing, random if processing is equal
visit.duprec <- visit.dup %>% 
  dplyr::filter(sensor=="ARU") %>% 
  left_join(visit.ab %>% 
              mutate(latr = round(lat, 4),
                     lonr = round(lon, 4))) %>% 
  group_by(date, lat, lon, location) %>% 
  summarize(n=n(),
            maxdur = max(duration),
            mindur = min(duration)) %>% 
  left_join(visit.ab)

visit.mindur.remove <- visit.duprec %>% 
  dplyr::filter(mindur!=maxdur,
                duration==mindur)
  
set.seed(1234)
visit.eqdur.remove <- visit.duprec %>% 
  dplyr::filter(mindur==maxdur) %>% 
  group_by(date, lat, lon, location) %>% 
  sample_n(1) %>% 
  ungroup() %>% 
  left_join(visit.ab)

#3d. Look at remaining duplicates----
visit.dup2 <- visit.ab %>% 
  anti_join(visit.josm.remove) %>% 
  anti_join(visit.mindur.remove) %>% 
  anti_join(visit.eqdur.remove) %>% 
  mutate(latr = round(lat, 4),
         lonr = round(lon, 4)) %>% 
  group_by(date, latr, lonr) %>% 
  summarize(n = n()) %>% 
  ungroup() %>% 
  dplyr::filter(n > 1) %>% 
  left_join(visit.ab %>% 
              mutate(latr = round(lat, 4),
                     lonr = round(lon, 4)))

#3e. Take out the ebird duplicates----
visit.ebird.remove <- visit.dup2 %>% 
  dplyr::filter(source=="eBird")

#3f. Look at duplicates again----
visit.dup3 <- visit.ab %>% 
  anti_join(visit.josm.remove) %>% 
  anti_join(visit.mindur.remove) %>% 
  anti_join(visit.eqdur.remove) %>% 
  anti_join(visit.ebird.remove) %>% 
  mutate(latr = round(lat, 4),
         lonr = round(lon, 4)) %>% 
  group_by(date, latr, lonr) %>% 
  summarize(n = n()) %>% 
  ungroup() %>% 
  dplyr::filter(n > 1) %>% 
  left_join(visit.ab %>% 
              mutate(latr = round(lat, 4),
                     lonr = round(lon, 4)))

set.seed(1234)
visit.dup3.keep <- visit.dup3 %>% 
  group_by(latr, lonr, date, n) %>% 
  sample_n(1) %>% 
  ungroup()

visit.dup3.remove <- visit.dup3 %>% 
  anti_join(visit.dup3.keep)

#3g. Check again----
visit.dup4 <- visit.ab %>% 
  anti_join(visit.josm.remove) %>% 
  anti_join(visit.mindur.remove) %>% 
  anti_join(visit.eqdur.remove) %>% 
  anti_join(visit.ebird.remove) %>% 
  anti_join(visit.dup3.remove) %>% 
  mutate(latr = round(lat, 4),
         lonr = round(lon, 4)) %>% 
  group_by(date, latr, lonr) %>% 
  summarize(n = n()) %>% 
  ungroup() %>% 
  dplyr::filter(n > 1) %>% 
  left_join(visit.ab %>% 
              mutate(latr = round(lat, 4),
                     lonr = round(lon, 4)))

#3h. Filter out duplicates from dataset----
#remove surveys before 1993 (first year with substantial data)
#create unique ID
visit <- visit.ab %>% 
  anti_join(visit.josm.remove) %>% 
  anti_join(visit.mindur.remove) %>% 
  anti_join(visit.eqdur.remove) %>% 
  anti_join(visit.ebird.remove) %>% 
  anti_join(visit.dup3.remove) %>% 
  dplyr::filter(year >= 1993) %>% 
  mutate(gisid = paste0(location, "_", year))

bird <- use.ab %>% 
  anti_join(visit.josm.remove) %>% 
  anti_join(visit.mindur.remove) %>% 
  anti_join(visit.eqdur.remove) %>% 
  anti_join(visit.ebird.remove) %>% 
  anti_join(visit.dup3.remove) %>% 
  dplyr::filter(year >= 1993) %>% 
  mutate(gisid = paste0(location, "_", year))

#F. IDENTIFY LOCATIONS FOR COVARIATE EXTRACTION####

#1. Identify projects with buffered locations----
#Any project with buffer of 55000
secret1 <- visit %>% 
  dplyr::filter(organization=="ABMI", buffer==5500) %>% 
  select(organization, project) %>% 
  unique()

#projects that have the same lat lon for more than one location
secret2 <- visit %>% 
  dplyr::filter(organization=="ABMI") %>% 
  dplyr::select(organization, project, location, lat, lon) %>% 
  unique() %>% 
  group_by(organization, project, lat, lon) %>% 
  summarize(n=n()) %>% 
  ungroup() %>% 
  dplyr::filter(n > 1) %>% 
  anti_join(secret1) %>% 
  dplyr::select(organization, project) %>% 
  unique()

secret <- rbind(secret1, secret2) %>% 
  mutate(topsecret = 1)

#2. Filter to just unique combinations of year & location----
location <- visit %>% 
  dplyr::select(gisid, source, organization, sensor, project, buffer, location, lat, lon, year) %>% 
  unique() %>% 
  left_join(secret) %>% 
  mutate(topsecret = ifelse(is.na(topsecret), 0, topsecret))

#3. Check Riverforks for inconsistency with GIS data----
gis <- read.csv(file.path(root, "Data", "gis", "topsecret_inventory.csv"))

check <- location %>% 
  dplyr::filter(topsecret==1,
                project=="ABMI-Riverforks") %>% 
  full_join(gis %>% 
              rename(gisid = match_,
                     year = year_) %>% 
              dplyr::filter(gisid!=""))
summary(check$NameFixed) #No NAs - good

#4. Check that location name fixes worked----
location.fix <- gis %>% 
  dplyr::filter(NameFixed=="yes") %>% 
  left_join(location)
summary(location.fix$NameFixed) #No NAs - good

write.csv(location, file.path(root, "Data", "gis", "birds_ab_locations.csv"), row.names = FALSE)

#G. SAVE!#############################
save(location, visit, bird, file=file.path(root, "Data", "1Harmonized.Rdata"))
  
#H. COMPARE############################
load(file.path(root, "data/ab-birds-all-2020-09-23.Rdata"))
load(file.path(root, "data/Harmonized.Rdata"))

nrow(dd)
nrow(visit)
table(visit$source)
dd.n <- data.frame(table(dd$PCODE)) %>% 
  arrange(-Freq)
visit.n <- data.frame(table(visit$project)) %>% 
  arrange(-Freq)

missing <- dd %>% 
  rename(lat = Y, lon = X) %>% 
  mutate(date = as.character(DATI),
         latr = round(lat, 4),
         lonr = round(lon, 4)) %>% 
  dplyr::select(date, latr, lonr) %>% 
  anti_join(visit %>% 
               mutate(date = as.character(date),
                      latr = round(lat, 4),
                      lonr= round(lon, 4))) %>% 
  dplyr::filter(!is.na(latr), 
                !is.na(date)) %>% 
  left_join(dd %>% 
              rename(lat = Y, lon = X) %>% 
              mutate(date = as.character(DATI),
                     latr = round(lat, 4),
                     lonr = round(lon, 4)))

year <- data.frame(table(dd$YEAR)) %>% 
  rename(year = Var1, dd.n = Freq) %>% 
  full_join(data.frame(table(visit$year)) %>% 
              rename(year = Var1, visit.n = Freq)) %>% 
  mutate(dd.n=ifelse(is.na(dd.n), 0, dd.n),
         visit.n=ifelse(is.na(visit.n), 0, visit.n)) %>% 
  dplyr::filter(as.numeric(as.character(year)) >= 1993,
                as.numeric(as.character(year)) <= 2022)

ggplot(year) +
  geom_point(aes(x=dd.n, y=visit.n, colour=year)) +
  geom_abline(aes(intercept=0, slope=1)) +
  xlab("# visits in previous dataset") +
  ylab("# visits in current dataset")

ggsave(filename=file.path(root, "Figures", "datasetversionN.jpeg"), width=6, height=4)

proj.dd <- dd %>% 
  dplyr::filter(YEAR %in% c(2012:2017)) %>% 
  group_by(PCODE, YEAR) %>% 
  summarize(n=n()) %>% 
  ungroup() %>% 
  pivot_wider(values_from=n, names_from="YEAR", names_prefix="Year", names_sort=TRUE) %>% 
  arrange(PCODE)

write.csv(proj.dd, "data/PreviousVersionProjects2012-2017.csv", row.names=FALSE)

proj.visit <- visit %>% 
  dplyr::filter(year %in% c(2012:2017)) %>% 
  group_by(organization, project, year) %>% 
  summarize(n=n()) %>% 
  ungroup() %>% 
  pivot_wider(values_from=n, names_from="year", names_prefix="Year", names_sort=TRUE) %>% 
  arrange(organization, project)

write.csv(proj.visit, "data/CurrentVersionProjects2012-2017.csv", row.names=FALSE)

#I. EXTRAS####
#1. Peak at multiyear data----
multiyear <- visit %>% 
  dplyr::select(location, lat, lon, year) %>% 
  unique() %>% 
  group_by(location, lat, lon) %>% 
  summarize(n=n()) %>% 
  ungroup() %>% 
  dplyr::filter(n > 1) %>% 
  left_join(visit %>% 
              dplyr::select(location, lat, lon, year) %>% 
              unique())

table(multiyear$n)

ggplot(multiyear) +
  geom_point(aes(x=lon, y=lat, colour=n)) +
  scale_colour_viridis_c()

#2. Plot all the visits----
ggplot(location) +
  geom_point(aes(x=lon, y=lat, colour=source)) +
  facet_wrap(~sensor)

ggsave(filename=file.path(root, "Figures", "Locations.jpeg"), width=6, height=4)
