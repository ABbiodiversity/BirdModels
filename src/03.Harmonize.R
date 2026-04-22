# ---
# title: ABMI models - harmonize data
# author: Elly Knight
# created: May 6, 2024
# ---

#NOTES################################

#PURPOSE: This script harmonizes data from wildtrax and ebird and formats for subsequent modelling steps.

#Riverforks data is stored in the ABMI Oracle database and should be downloaded  using dBeaver. Use the UNRESTRICTED_ACCESS acount to ensure full retrieval of species at risk records. Contact Joan (qfang@ualberta.ca) for access information.

#The location csv output of this script should be sent to Eric and/or Brandon for GIS covariate extraction. The output Eric provides will be used as input in the next script.

#FUTURE VERSIONS: Riverforks data will be available in wildtrax and will not require extraction from the ABMI Oracle database.

#PREAMBLE############################

#1. Load packages----

library(tidyverse) #basic data wrangling
library(wildrtrax) #to download data from wildtrax
library(auk) #to handle eBird data
library(data.table) #for binding lists into dataframes
library(lubridate) #date wrangling
library(downloader) #download zipped provincial boundary shp
library(sf) #shapefile wrangling
library(terra) #raster wrangling
library(hms)

#Source temp fixes for wildrtrax.
source('temp_functions/wt_replace_tmtt_fixed.R')
source('temp_functions/wt_tidy_species_fixed.R')
source('temp_functions/wt_make_wide_fixed.R')

#2. Set root path for data on google drive----

root <- "G:/Shared drives/ABMI_RHedley/Projects/BirdModels"

#3. Login to WildTrax----
config <- "src/00.WTlogin.R"
source(config)

#4. Authenticate----
wt_auth()

#5. Pointer to most recent RData file:
most_recent_file <- file.path(root, "Data", "WildTrax", "wildtrax_raw_2026-03-10.Rdata")

#HARMONIZE###############################

#1. Get list of bird species----
birdcodes <- wt_get_species() |> 
  rename(common = species_common_name) |> 
  dplyr::select(common, species_code) |> 
  mutate(common = case_when(common == "Bald eagle" ~ "Bald Eagle",
                            common == "Northern Pygmy-owl" ~ "Northern Pygmy-Owl",
                            common == "Black-crowned Night Heron" ~ "Black-crowned Night-Heron",
                            !is.na(common) ~ common))

###################################################
# Manual step. Update the list of birds observed in the projects.
# Specifically, there are often species splits, creating duplicates of 
# the same species with multiple codes. Then replace the birds-vYYYY.csv file.
# I went through and manually vetted the list, and added a preferred_synonym
# column for species with name changes. I also removed all species
# in birds-v2024 that had show=='o'. Those could be revisited in future and
# added back in if their situation changes.
###################################################
spp <- read.csv(file.path(root, "Data", "lookups", "birds-v2026.csv")) 
spp <- spp[order(spp$common),]
spp <- left_join(spp, birdcodes)

#Subset species with aliases due to name changes.
aliases <- spp[nchar(spp$preferred_synonym) == 4, ]

#Filter spp to remove species with aliases.
spp <- spp[nchar(spp$preferred_synonym) != 4, ]

rm(birdcodes) #birdcodes not used hereafter. In fact, now not needed at all.

#2. Set desired columns----
colnms <- c("source", "organization", "project_id", "sensor", "task_method", "location", "buffer", "latitude", "longitude", "date_time", "duration", "distance")

#3. Load WildTrax data and fix a few species codes-----
load(most_recent_file)

#Where a species_code exists in the data, replace with preferred_synonym.
for(i in 1:nrow(aliases)) {
  aru.wt$species_code[aru.wt$species_code == aliases$species_code[i]] <- aliases$preferred_synonym[i]
  pc.wt$species_code[pc.wt$species_code == aliases$species_code[i]] <- aliases$preferred_synonym[i]
}

#Set options for lubridate, to enable checking of date parsing.
options(lubridate.verbose = TRUE)

#4. Wrangle WT ARU data----
#fix a few names to match the GIS
#remove detections > 180s for hybrid method and fix duration
use.aru <- aru.wt |> 
  dplyr::filter(!(task_method=="1SPM Audio/Visual hybrid" & as.numeric(detection_time) > 180)) |> 
  wt_tidy_species_fixed_aru(remove=c("mammal", "amphibian", "abiotic", "insect", "human", "unknown")) |> 
  mutate(abundance = ifelse(species_code == 'NONE', 0, abundance)) |>
  mutate(abundance = ifelse(abundance == 'CI 1', 1, abundance)) |>
  mutate(abundance = ifelse(abundance == 'CI 2', 1, abundance)) |>
  mutate(abundance = ifelse(abundance == 'CI 3', 1, abundance)) |>
  mutate(abundance = ifelse(abundance == 'N/A', 1, abundance)) |>
  wt_replace_tmtt_fixed() |>
  wt_make_wide_fixed() |> 
  mutate(source="WildTrax",
         sensor="ARU",
         distance=Inf,
         date_time = ymd_hms(recording_date_time, tz="America/Edmonton", truncated=3),
         #duration = round(as.numeric(str_sub(task_duration, -100, -2))),
         duration = (round(task_duration / 60) * 60),
         task_method = ifelse(is.na(task_method), 'Unknown', task_method),
         duration = ifelse(task_method=="1SPM Audio/Visual hybrid", 180, duration),
         location_buffer_m = ifelse(is.na(location_buffer_m), 0, location_buffer_m)) |> 
  mutate(tod = as_hms(date_time)) |>
  rename(buffer = location_buffer_m) |>
  filter(duration > 0, !(duration %in% c(240, 360, 420, 480))) #Remove a few weird durations.

#Richard: further rounding to the nearest minute, to remove some weird task_durations.
table(use.aru$task_duration)
table(use.aru$duration, useNA = 'always')


#5. Wrangle WT PC data----
#remove counts with unknown duration and distance
use.pc <- pc.wt |> 
  dplyr::filter(survey_duration_method!="UNKNOWN",
                survey_distance_method!="UNKNOWN") |> 
  wt_tidy_species_fixed_pc(remove=c("mammal", "amphibian", "abiotic", "insect", "human", "unknown")) |> 
  wt_make_wide_fixed() |> 
  mutate(source="WildTrax",
         sensor="PC",
         task_method="PC",
         date_time = ymd_hms(survey_date, tz="America/Edmonton", truncated=3)) |> 
  rowwise() |>
  mutate(durationMethod = ifelse(str_sub(survey_duration_method, -1, -1)=="+",
                                 str_sub(survey_duration_method, -100, -2),
                                 survey_duration_method),
         chardur = str_locate_all(durationMethod, "-"),
         chardurmax = max(chardur),
         duration = as.numeric(str_sub(durationMethod, chardurmax+1, -4))*60,
         chardis = str_locate_all(survey_distance_method, "-"),
         chardismax = max(chardis),
         distance1 = str_sub(survey_distance_method, chardismax+1, -2),
         distance = ifelse(distance1 %in% c("AR", "IN"), Inf, as.numeric(distance1))) |>
  ungroup() |> 
  rename(buffer = location_buffer_m)

table(use.pc$survey_duration_method, useNA = 'always')
table(use.pc$duration, useNA = 'always')

#6. Wrangle Riverforks data----

#6a. Read it in-----
det.rf <- read.csv(file.path(root, "Data", "Riverforks", "RT_BIRD_COUNT_BL.csv"), header=TRUE)
loc.rf <- read.csv(file.path(root, "Data", "Riverforks", "A_RT_SITE_PHYCHAR_202301161554.csv"), header=TRUE)

#6b. Put together----
#concatenate site name with point count # (there are 9 per site) before joining
raw.rf <- det.rf |>
  mutate(SITE = paste0(SITE, "-", TBB_POINT_COUNT),
         YEAR = as.integer(YEAR)) |> 
  left_join(loc.rf |> 
              mutate(SITE = paste0(SITE, "-", TSFG_POINT_COUNT)) |> 
              dplyr::select(ROTATION, SITE, YEAR, SITE_LATITUDE, SITE_LONGITUDE) |> 
              unique())

#6c. Wrangle----
#remove visits with no time data
#fix some names to match the gis
#fix some dates to match the gis
use.rf <- raw.rf |> 
  rename(location=SITE, latitude = SITE_LATITUDE, longitude = SITE_LONGITUDE, observer = ANALYST, common = COMMON_NAME) |> 
  mutate(source = "riverforks",
         organization="ABMI",
         project_id=9998,
         sensor="ARU",
         task_method = "1SPM",
         buffer = 5500,
         duration = 3*60,
         distance = Inf,
         abundance = 1,
         date_time = dmy_hm(paste0(ADATE, " ", TBB_START_TIME), tz="America/Edmonton")) |> 
  dplyr::filter(!is.na(date_time)) |> 
  mutate(date_time = case_when(str_sub(location, 1, 4)=="1525" ~ dmy_hm(paste0(str_sub(ADATE, -100, -3), "07", " ", TBB_START_TIME)),
                               !is.na(date_time) ~ date_time)) |> 
  left_join(spp |> 
              dplyr::select(species_code, common)) |> 
  dplyr::select(all_of(colnms), species_code, abundance) |> 
  pivot_wider(names_from=species_code, values_from=abundance, values_fn=sum, values_fill=0)

#7. Wrangle ebird data----
#Note this assumes observations with "X" individuals are 1s
#Filter out hotspots
#Replace common name with alpha code

#7a. Read it in----
raw.ebd <- read_ebd(file.path(root, "Data", "ebd", "ebd_data_filtered.txt"))

#7b. Wrangle----
use.ebd <- raw.ebd |> 
  dplyr::filter(locality_type!="H") |> 
  mutate(source = "eBird",
         organization = "eBird",
         project_id=99999,
         sensor="PC",
         task_method="PC",
         buffer=0,
         date_time = ymd_hms(paste0(observation_date, time_observations_started), tz="America/Edmonton"),
         distance = Inf,
         abundance = as.numeric(ifelse(observation_count=="X", 1, observation_count)),
         duration = duration_minutes*60) |> 
  rename(observer = observer_id,
         location = checklist_id,
         common = common_name) |> 
  left_join(spp |> 
              dplyr::select(species_code, common)) |> 
  dplyr::select(all_of(colnms), species_code, abundance) |> 
  pivot_wider(names_from=species_code, values_from=abundance, values_fn=sum, values_fill=0)

#PUT TOGETHER############################

#1. Put everything together----
use <- data.table::rbindlist(list(use.aru, use.pc, use.rf, use.ebd), fill=TRUE) |> 
  dplyr::select(all_of(c(colnms, spp$species_code))) |> 
  dplyr::filter(!is.na(latitude)) |> 
  mutate(across(c(ALFL:YRWA), \(x) replace_na(0)))

#2. Clip by provincial boundaries & filter to AB----

#Download provincial boundaries shapefile - only need to do this one
# temp <- tempfile()
# download("https://www12.statcan.gc.ca/census-recensement/2021/geo/sip-pis/boundary-limites/files-fichiers/lpr_000b21a_e.zip", temp)
# unzip(zipfile=temp, exdir=file.path(root, "Data", "gis"))
# unlink(temp)

#2a. Get boundaries----
shp <- read_sf(file.path(root, "Data", "gis", "lpr_000b21a_e.shp")) |> 
  dplyr::filter(PRNAME=="Alberta")

#2b. Filter points in Alberta----
use.ab <- use |> 
  dplyr::filter(!is.na(longitude)) |> 
  st_as_sf(coords=c("longitude", "latitude"), crs=4326, remove = FALSE) |> 
  st_transform(crs=st_crs(shp)) |>
  st_filter(st_geometry(shp)) |> 
  st_drop_geometry() |> 
  dplyr::select(colnames(use)) |> 
  rbind(use |> 
          dplyr::filter(is.na(longitude)))

#3. Remove duplicate surveys----

#Flag column for surveys to remove.
use.ab$remove_dup <- FALSE
#Round lat and lon to 4 digits (11 meters at the equator, about 6 meters wide in AB).
use.ab$latr <- round(use.ab$latitude,4)
use.ab$lonr <- round(use.ab$longitude,4)
#Column for unique lat/lon/date/time.
use.ab$point_in_time <- paste(use.ab$latr, use.ab$lonr, use.ab$date_time, sep = '_')
#Flag column for secret locations.
use.ab$secret <- FALSE

#Identify projects with at least one buffer > 0, and use that to flag secret locations.
buffered_projects <- aggregate(buffer ~ project_id, data = use.ab, FUN = max)
buffered_projects <- buffered_projects[buffered_projects$buffer > 0,]

#use buffered projects to set secret == TRUE
use.ab$secret[use.ab$project_id %in% buffered_projects$project_id] <- TRUE

#Identify duplicates surveys (same location and date/time).
duplicates <- table(use.ab$point_in_time)
duplicates <- duplicates[duplicates > 1]

set.seed(1234)
for(i in 1:length(duplicates)) {
  sub <- use.ab[use.ab$point_in_time == names(duplicates)[i],]
  
  # #If different surveys have different location names, flag all as secret locations.
  # if(length(unique(sub$location)) == nrow(sub)) {
  #   sub$secret <- TRUE
  # } 
  
  #If two rows, one with ARU and PC, select ARU.
  if(nrow(sub) == 2 & sum(sub$sensor == 'ARU') == 1 &
     sum(sub$sensor == 'PC') == 1) {
    sub$remove_dup[sub$sensor == 'PC'] <- TRUE
  }
  
  #If the duplicates are from multiple sources, remove eBird.
  if(length(unique(sub$source)) > 1 & 'eBird' %in% sub$source) {
    sub$remove_dup[sub$source == 'eBird'] <- TRUE
  }
  
  #If the duplicates have the same location name, select randomly among the 
  #longest surveys.
  if(length(unique(sub$location)) == 1) {
    longest <- which(sub$duration == max(sub$duration))
    selected <- sample(longest, 1)
    sub$remove_dup <- TRUE
    sub$remove_dup[selected] <- FALSE
  }
  
  #If the above procedure has not accounted for any duplicates, a message will be printed.
  if(sum(!sub$remove_dup) != 1 & !all(sub$secret)) {
    #if(sum(sub$remove_dup) == 0) {
    message(i,' Projects: ', paste(unique(sub$project_id), collapse = ' '), ' nrow: ', nrow(sub),
            ' locations: ', paste(unique(sub$location), collapse = ' '))
  }
  
  # #Otherwise select longest survey.
  # sub <- sub[which.max(sub$duration),]
  # 
  # if(nrow())
  
}

use.dup <- use.ab |> 
  mutate(latr = round(latitude, 4),
         lonr = round(longitude, 4)) |> 
  group_by(date_time, latr, lonr) |> #Looking for same date_time, latr, lonr.
  summarize(n = n()) |> 
  ungroup()

#3a. Investigate----
#Take out buffered locations
use.dup <- use.ab |> 
  mutate(latr = round(latitude, 4),
         lonr = round(longitude, 4)) |> 
  group_by(date_time, latr, lonr) |> #Looking for same date_time, latr, lonr.
  summarize(n = n()) |> 
  ungroup() |> 
  dplyr::filter(n > 1) |> 
  left_join(use.ab |> 
              dplyr::select(all_of(colnms)) |> 
              mutate(latr = round(latitude, 4),
                     lonr = round(longitude, 4)),
            multiple="all")
#no eBird duplicates
#some point count data entry errors, but really only 2 sources of duplicates:
#1. JOSM points that have human & ARU data & eBird data
#2. Duplicate processing of recordings

#3b. Identify JOSM points with human & ARU data----
#Keep the ARU ones because frankly they're probably better data
use.josm <- use.dup |> 
  left_join(projects.use) |> 
  dplyr::filter(project %in% c("JOSM ECCC Cause and Effect Monitoring for Landbirds 2013",
                               "Oil Sands Monitoring CWS Prairie Region 2012",
                               "Oil Sands Monitoring CWS Prairie Region 2013",
                               "Oilsands Monitoring CWS Prairie Region 2014"))

use.josm.remove <- use.josm |> 
  dplyr::filter(sensor=="PC")

#3c. Identify recordings that have been processed twice----
#Use recording with longer duration processing, random if processing is equal
use.duprec <- use.dup |> 
  dplyr::filter(sensor=="ARU") |> 
  left_join(use.ab |> 
              mutate(latr = round(latitude, 4),
                     lonr = round(longitude, 4)) |> 
              dplyr::select(all_of(colnms)),
            multiple="all") |> 
  group_by(date_time, latitude, longitude, location) |> 
  summarize(n=n(),
            maxdur = max(duration),
            mindur = min(duration)) |> 
  left_join(use.ab,
            multiple="all")

use.mindur.remove <- use.duprec |> 
  dplyr::filter(mindur!=maxdur,
                duration==mindur)

set.seed(1234)
use.eqdur.remove <- use.duprec |> 
  dplyr::filter(mindur==maxdur) |> 
  group_by(date_time, latitude, longitude, location) |> 
  sample_n(1) |> 
  ungroup() |> 
  left_join(use.ab,
            multiple="all")

#3d. Look at remaining duplicates----
use.dup2 <- use.ab |> 
  dplyr::filter(!is.na(longitude)) |> 
  anti_join(use.josm.remove) |> 
  anti_join(use.mindur.remove) |> 
  anti_join(use.eqdur.remove) |> 
  mutate(latr = round(latitude, 4),
         lonr = round(longitude, 4)) |> 
  group_by(date_time, latr, lonr) |> 
  summarize(n = n()) |> 
  ungroup() |> 
  dplyr::filter(n > 1) |> 
  left_join(use.ab  |> 
              mutate(latr = round(latitude, 4),
                     lonr = round(longitude, 4)),
            multiple="all")

#3e. Take out the ebird duplicates----
use.ebird.remove <- use.dup2 |> 
  dplyr::filter(source=="eBird")

#3f. Look at duplicates again----
use.dup3 <- use.ab |> 
  dplyr::filter(!is.na(longitude)) |> 
  anti_join(use.josm.remove) |> 
  anti_join(use.mindur.remove) |> 
  anti_join(use.eqdur.remove) |> 
  anti_join(use.ebird.remove) |> 
  mutate(latr = round(latitude, 4),
         lonr = round(longitude, 4)) |> 
  group_by(date_time, latr, lonr) |> 
  summarize(n = n()) |> 
  ungroup() |> 
  dplyr::filter(n > 1) |> 
  left_join(use.ab |> 
              mutate(latr = round(latitude, 4),
                     lonr = round(longitude, 4)),
            multiple="all") |> 
  arrange(latr, lonr, date_time) |> 
  left_join(projects.use)

set.seed(1234)
use.dup3.keep <- use.dup3 |> 
  group_by(latr, lonr, date_time, n) |> 
  sample_n(1) |> 
  ungroup()

use.dup3.remove <- use.dup3 |> 
  anti_join(use.dup3.keep)

#3g. Check again----
use.dup4 <- use.ab |> 
  dplyr::filter(!is.na(longitude)) |> 
  anti_join(use.josm.remove) |> 
  anti_join(use.mindur.remove) |> 
  anti_join(use.eqdur.remove) |> 
  anti_join(use.ebird.remove) |> 
  anti_join(use.dup3.remove) |> 
  mutate(latr = round(latitude, 4),
         lonr = round(longitude, 4)) |> 
  group_by(date_time, latr, lonr) |> 
  summarize(n = n()) |> 
  ungroup() |> 
  dplyr::filter(n > 1) |> 
  left_join(use.ab |> 
              dplyr::select(all_of(colnms)) |> 
              mutate(latr = round(latitude, 4),
                     lonr = round(longitude, 4)),
            multiple="all")

#3h. Filter out duplicates from dataset----
#remove surveys before 1993 (first year with substantial data)
#remove ABMI surveys that can't resolve with GIS data
#create unique ID
#Fix location naming errors

fix <- read.csv(file.path(root, "Data", "gis", "birds_mismatches_lookup.csv")) |> 
  dplyr::filter(locationfix!="")

use <- use.ab |> 
  anti_join(use.josm.remove) |> 
  anti_join(use.mindur.remove) |> 
  anti_join(use.eqdur.remove) |> 
  anti_join(use.ebird.remove) |> 
  anti_join(use.dup3.remove) |> 
  mutate(year = year(date_time)) |> 
  dplyr::filter(year >= 1993) |>  
  dplyr::filter(!(location=="701-1" & year==2004)) |> 
  mutate(gisid = paste0(location, "_", project_id, "_", year)) |> 
  left_join(fix) |> 
  mutate(location = ifelse(!is.na(locationfix), locationfix, location))

#IDENTIFY LOCATIONS FOR COVARIATE EXTRACTION####

#1. Identify projects with buffered locations----
#Any project with buffer of 55000
secret1 <- use |> 
  dplyr::filter(organization=="ABMI", buffer==5500) |> 
  select(organization, project_id) |> 
  unique() |> 
  left_join(projects.use |> 
              dplyr::select(project, project_id))

#projects that have the same lat lon for more than one location
secret2 <- use |> 
  dplyr::filter(organization=="ABMI") |> 
  dplyr::select(organization, project_id, location, latitude, longitude) |> 
  unique() |> 
  group_by(organization, project_id, latitude, longitude) |> 
  summarize(n=n()) |> 
  ungroup() |> 
  dplyr::filter(n > 1) |> 
  anti_join(secret1) |> 
  dplyr::select(organization, project_id) |> 
  unique() |> 
  left_join(projects.use |> 
              dplyr::select(project, project_id))

secret <- rbind(secret1, secret2) |> 
  mutate(topsecret = 1,
         project = ifelse(is.na(project), "ABMI-Riverforks", project))

#2. Filter to just unique combinations of year & location----
location <- use |> 
  dplyr::select(gisid, source, organization, sensor, project_id, buffer, location, latitude, longitude, year) |> 
  left_join(projects.use |> 
              dplyr::select(project, project_id)) |> 
  mutate(project = case_when(project_id==9998 ~ "ABMI-Riverforks",
                             project_id==99999 ~ "eBird", 
                             !is.na(project) ~ project)) |> 
  unique() |> 
  left_join(secret) |> 
  mutate(topsecret = ifelse(is.na(topsecret), 0, topsecret)) 

#3. Check against GIS inventory----
gis <- read.csv(file.path(root, "Data", "gis", "OffGridCamARU_Analysis_camaru_20231208.csv"))

missing <- read.csv(file.path(root, "Data", "gis", "birds_mismatches_lookup.csv")) |> 
  dplyr::filter(locationfix=="")

gis.location <- location |> 
  dplyr::filter(topsecret==1) |> 
  left_join(gis |> 
              rename(location = Site_ID,
                     year = survey_year) |> 
              mutate(gis=1,
                     location = case_when(str_sub(location, 1, 5)=="OG-EI" ~ Site,
                                          !is.na(location) ~ location)),
            multiple="all") |> 
  mutate(gis = ifelse(is.na(gis), 0, 1))

gis.na <- dplyr::filter(gis.location, gis==0) |> 
  arrange(project, location) |> 
  dplyr::select(project, location, year) |> 
  left_join(missing)

write.csv(gis.na, (file.path(root, "Data", "gis", "birds_gis_missing.csv")), row.names = FALSE)

#4. Take out the mismatches for now----
location <- location |> 
  anti_join(gis.na |> 
              dplyr::select(project, location, year))

use <- use |> 
  anti_join(gis.na |> 
              dplyr::select(project, location, year))

#5. Save----
write.csv(location, file.path(root, "Data", "gis", "birds_ab_locations.csv"), row.names = FALSE)

#SAVE!#############################
save(location, use, file=file.path(root, "Data", "Harmonized.Rdata"))