# ---
# title: ABMI models - get data
# author: Elly Knight
# created: November 29, 2022
# ---

#NOTES################################

#The "BAMProjects_WildTrax.csv" file is a list of all projects currently in WildTrax that can be used in ABMI models. This file should be updated for each iteration of in collaboration with Erin Bayne. Future versions of this spreadsheet can hopefully be derived by a combination of organization and a google form poll for consent from other organizations.

#The "BAMProjects_WildTrax.csv" file also contains information on which ARU projects are processed for a single species or taxa and therefore those visits should only be used for models for the appropriate taxa. This file should be updated for each iteration of the national models in collaboration with Erin Bayne.

#There are a handful of projects that are not downloading properly via wildRtrax. An issue is open on this. These projects are listed in the error.log object.

#BAM patch was compiled by Melina Houle in BAM to include datasets that are not yet loaded into WildTrax. Future iterations of the national models should not need this patch, with the exception of the BBS data which is likely too large to ever be uploaded in WT.

#raw eBird data is downloaded from the eBird interface at https://ebird.org/data/download/ebd prior to wrangling with the "auk" package and will require a request for access. Use the custom download tool to download only the datasets for Canada and the US instead of the global dataset. Note you will also need the global sampling file to use the auk package for zero filling.

#raw eBird data omits Great Grey Owl & Northern Hawk Owl as sensitive species (https://support.ebird.org/en/support/solutions/articles/48000803210?b_id=1928&_gl=1*xq054u*_ga*ODczMTUyMjcuMTY2OTE0MDI4Ng..*_ga_QR4NVXZ8BM*MTY2OTE0MDI4NS4xLjEuMTY2OTE0MDM3OC4zNS4wLjA.&_ga=2.147122167.150058226.1669140286-87315227.1669140286) and should not be used for modelling these two species.

#wrangling eBird data with the auk package requires installation of AWK on windows computers. Please see #https://cornelllabofornithology.github.io/auk/articles/auk.html.

#eBird data has not been zerofilled because there was no species filtering done and we are assuming that all stationary counts have at least 1 bird observed.

#The column "sensor" currently only differentiates between ARU & human point count data types. Future versions should consider differentiating between SM2 ARU data types and other ARU data types due to differences in the perceptibility of these two approaches, either via QPAD or a correction factor.

#The replace TMTTs script will be replaced by a wildRtrax function in the near future.

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

rootabmi <- "G:/My Drive/ABMI/Projects/BirdModels/"

#For access to the bam patch
rootbam <- "G:/.shortcut-targets-by-id/0B1zm_qsix-gPbkpkNGxvaXV0RmM/BAM.SharedDrive/RshProjs/PopnStatus/NationalModelsV4.1/PointCount/"

#A. DOWNLOAD DATA FROM WILDTRAX#######################

#1. Login to WildTrax----
config <- "script/login.R"
source(config)

#1. Get list of projects from WildTrax----
wt_auth()

#sensor = PC gives all ARU and point count projects
project.list <- wt_get_download_summary(sensor_id = 'PC')

#2. Convert to a plain dataframe----
projects <- data.frame(project = as.character(project.list$project),
                       project_id = as.numeric(project.list$project_id),
                       sensorId = as.character(project.list$sensorId),
                       tasks = as.numeric(project.list$tasks),
                       status = as.character(project.list$status))

#3. Filter by list of projects we can use----
#Fix names to match report
project.names <- data.frame(project_id = c(325, 432, 856, 977, 978),
                            newname = c("Tłı̨chǫ Winter Road CWS Northern Region 2019",
                                        "Ts’udé Nilįné Tuyeta PA CWS Northern Region 2020",
                                        "Ts’udé Nilįné Tuyeta PA CWS Northern Region 2021",
                                        "Nááts'įhch'oh NPR 2018 CWS Northern Region",
                                        "Nááts'įhch'oh NPR 2019 CWS Northern Region"))

use <- read.csv(file.path("data/ABMIProjects_WildTrax.csv"))

projects.wt <- projects %>% 
  dplyr::filter(project_id %in% use$project_id) %>% 
  left_join(project.names) %>% 
  mutate(project = ifelse(is.na(newname), project, newname)) %>% 
  dplyr::select(-newname)

#4. Loop through projects to download data----
dat.list <- list()
error.log <- data.frame()
for(i in 1:nrow(projects.wt)){
  
  #Do each sensor type separately because the reports have different columns and we need different things for each sensor type
  if(projects.wt$sensorId[i]=="ARU"){
    
    #Get summary report
    report.try <- try(wt_download_report(project_id = projects$project_id[i], sensor_id = projects$sensorId[i], weather_cols = F, report = "summary"))
    
    #Get task report for ARU model
    task.try <- try(wt_download_report(project_id = projects$project_id[i], sensor_id = projects$sensorId[i], weather_cols = F, report = "task"))
    
    if(class(report.try)=="data.frame"){
      dat.try <- report.try %>% 
        left_join(task.try)
    }
  }
  
  if(projects.wt$sensorId[i]=="PC"){
    
    dat.try <- try(wt_download_report(project_id = projects.wt$project_id[i], sensor_id = projects.wt$sensorId[i], weather_cols = F, report="report"))
    
  }
  
  if(class(dat.try)=="data.frame"){
    dat.list[[i]] <- dat.try
  }
  
  #Log projects that error
  if(class(dat.try)!="data.frame"){
    error.log <- rbind(error.log, 
                       projects.wt[i,])
    
  }
  
  print(paste0("Finished dataset ", projects.wt$project[i], " : ", i, " of ", nrow(projects.wt), " projects"))
  
}

#5. Collapse list----
#standardize column names between sensor types
#fix special character project names
report.names <- data.frame(project=c("EdÃ©hzhÃ­e National Wildlife Area CWS Northern Region 2016",
                                     "EdÃ©hzhÃ­e National Wildlife Area CWS Northern Region 2019",
                                     "TÅ‚Ä±Ì¨chÇ« Winter Road CWS Northern Region 2019",
                                     "Tsâ€™udÃ© NilÄ¯nÃ© Tuyeta PA CWS Northern Region 2020",
                                     "ForÃªt Montmorency long-term bird survey 2000",
                                     "RÃ©gularisation du Lac KÃ©nogami EIA 2001",
                                     "Tsâ€™udÃ© NilÄ¯nÃ© Tuyeta PA CWS Northern Region 2021",
                                     "NÃ¡Ã¡ts'Ä¯hch'oh NPR 2018 CWS Northern Region",
                                     "NÃ¡Ã¡ts'Ä¯hch'oh NPR 2019 CWS Northern Region",
                                     "EdÃ©hzhÃ­e National Wildlife Area CWS Northern Region 2021",
                                     "Thaidene NÃ«nÃ© PA CWS Northern Region 2022"),
                           newname=c("Edéhzhíe National Wildlife Area CWS Northern Region 2016",
                                     "Edéhzhíe National Wildlife Area CWS Northern Region 2019",
                                     "Tchicho Winter Road CWS Northern Region 2019",
                                     "Ts’udé Niliné Tuyeta PA CWS Northern Region 2020",
                                     "Forêt Montmorency long-term bird survey 2000",
                                     "Régularisation du Lac Kénogami EIA 2001",
                                     "Ts’udé Niliné Tuyeta PA CWS Northern Region 2021",
                                     "Nááts'ihch'oh NPR 2018 CWS Northern Region",
                                     "Nááts'ihch'oh NPR 2019 CWS Northern Region",
                                     "Edéhzhíe National Wildlife Area CWS Northern Region 2021",
                                     "Thaidene Nëné PA CWS Northern Region 2022"))

raw.wt <- rbindlist(dat.list, fill=TRUE) %>%
  mutate(project = ifelse(is.na(project), project_name, project),
         speciesCode = ifelse(is.na(speciesCode), species_code, speciesCode),
         date = ifelse(is.na(date), recording_date, date)) %>%
  dplyr::select(-project_name, -recording_date, -species_code) %>% 
  left_join(report.names) %>% 
  mutate(project = ifelse(is.na(newname), project, newname)) %>% 
  dplyr::select(-newname)

#6. Save date stamped data & project list----
save(raw.wt, projects.wt, error.log, file=paste0(rootabmi, "Data/wildtrax_raw_", Sys.Date(), ".Rdata"))

#B. GET PATCH DATA###############################

#1. BAM patch----
raw.bam <- readRDS(file.path(rootbam, "pc_patch.rds"))

#C. GET EBIRD DATA##########################

#1. Set ebd path----
auk_set_ebd_path(file.path(rootabmi, "Data/ebd_CA-AB_relOct-2022"), overwrite=TRUE)

#2. Define filters----
filters <- auk_ebd(file="ebd_CA-AB_relOct-2022.txt") %>% 
  auk_protocol("Stationary") %>% 
  auk_duration(c(0, 10)) %>% 
  auk_complete()

#3. Filter data----
#select columns to keep
filtered <- auk_filter(filters, file=file.path(rootabmi, "Data/ebd_data_filtered.txt"), overwrite=TRUE,
                       keep = c("group identifier", "sampling_event_identifier", "scientific name", "common_name", "observation_count", "latitude", "longitude", "locality_type", "observation_date", "time_observations_started", "observer_id", "duration_minutes"))

#D. HARMONIZE###############################

#1. Set desired columns----
colnms <- c("source", "organization", "project", "sensor", "equipment", "singlesp", "location", "buffer", "lat", "lon", "year", "date", "observer", "duration", "distance", "species", "abundance", "isSeen", "isHeard")

#2. Wrangle wildtrax data-----

load(file.path(rootabmi, "data/wildtrax_raw_2022-12-05.Rdata"))

#2a. A bit of prep----
dat.wt <- raw.wt %>% 
  rename(lat = latitude, lon = longitude, species = speciesCode, equipment = equipment_used) %>% 
  full_join(projects.wt %>% 
              rename(sensor = sensorId) %>% 
              dplyr::select(project_id, project, sensor)) %>% 
  mutate(source = "WildTrax", 
         date = ymd_hms(date),
         year = year(date),
         buffer=ifelse(is.na(bufferRadius.m.), buffer, bufferRadius.m.)) %>% 
  separate(buffer, into=c("buffer"), sep=" ", remove=TRUE) %>% 
  mutate(buffer = as.numeric(ifelse(!buffer %in% c("50", "10000"), str_sub(buffer, -100, -2), buffer)),
         buffer = ifelse(is.na(buffer), 0, buffer))

#2b. Get the point count data----
#wrangle distance and duration maximums
#remove counts with unknown duration and distance
pc.wt.meth <- dat.wt %>% 
  dplyr::filter(sensor=="PC",
                durationMethod!="UNKNOWN",
                distanceMethod!="UNKNOWN") %>% 
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
  mutate(singlesp = "n") %>% 
  left_join(pc.wt.meth) %>% 
  dplyr::select(all_of(colnms)) %>% 
  data.frame()

#2c. Get the aru data----
#filter ARU data to first detection of each individual
#wrangle duration
#identify datasets that should only be used for single species

ssp <- read.csv("data/ABMIProjects_WildTrax.csv") %>% 
  dplyr::filter(ssp=="y")

aru.wt <- dat.wt %>% 
  dplyr::filter(sensor=="ARU") %>% 
  mutate(singlesp = ifelse(project_id %in% ssp$project_id, "y", "n")) %>% 
  separate(method, into=c("duration", "method"), remove=TRUE) %>% 
  mutate(duration = as.numeric(str_sub(duration, -100, -2))/60,
         distance = Inf) %>% 
  group_by(source, project, sensor, singlesp, location, buffer, lat, lon, year, date, observer, duration, distance, species, abundance, individual_appearance_order) %>%
  mutate(first_tag = min(tag_start_s)) %>%
  ungroup() %>%
  dplyr::filter(tag_start_s == first_tag) %>% 
  dplyr::select(all_of(colnms))

#2d. Replace TMTTs with predicted abundance----
tmtt <- read.csv("C:/Users/Elly Knight/Documents/ABMI/Projects/TMTT/data/tmtt_predictions_mean.csv") %>% 
  rename(species = species_code)

user <- read.csv("C:/Users/Elly Knight/Documents/ABMI/Projects/TMTT/data/app_user.csv") %>% 
  rename(observer = user_name) %>% 
  dplyr::select(observer, user_id)

tmtt.wt <- aru.wt %>% 
  dplyr::filter(abundance=="TMTT") %>% 
  left_join(user) %>% 
  mutate(user_id = ifelse(is.na(user_id), observer, user_id))%>% 
  mutate(species = ifelse(species %in% tmtt$species, species, "species"),
         user_id = as.integer(ifelse(user_id %in% tmtt$user_id, user_id, 0))) %>% 
  data.frame() %>% 
  left_join(tmtt) %>% 
  mutate(abundance = round(pred)) %>% 
  dplyr::select(colnames(aru.wt))

#2e. Put back together----
use.wt <- aru.wt %>% 
  dplyr::filter(abundance!="TMTT") %>% 
  rbind(tmtt.wt) %>% 
  rbind(pc.wt)

#3. Wrangle BAM patch data----
#wrangle distance and duration maximums
#remove counts with odd duration method entries
#replace all unknown dates (MN-BBATLAS, NEFBMP2012-19) with June 15 of 2012
bam.meth <- raw.bam %>% 
  dplyr::select(distanceMethod, durationMethod) %>% 
  unique() %>% 
  dplyr::filter(!durationMethod %in% c("", " during the 10 minutes.", " seemed to bring food then brooded; assume nestlings stil", " timeperiod C")) %>% 
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

use.bam <- raw.bam %>% 
  mutate(source = "BAM",
         sensor = "PC",
         equipment = NA,
         singlesp = "n",
         date = ymd_hms(date),
         year = year(date)) %>% 
  rename(buffer = 'bufferRadius(m)', lat = latitude, lon = longitude, species=speciesCode) %>%
  mutate(buffer = ifelse(is.na(buffer), 0, buffer)) %>% 
  dplyr::filter(!is.na(date)) %>% 
  left_join(bam.meth) %>% 
  dplyr::select(all_of(colnms))

#5. Wrangle ebird data----
#Note this assumes observations with "X" individuals are 1s
#Filter out hotspots
#Replace common name with alpha code

raw.ebd <- read_ebd(file.path(rootabmi, "Data/ebd_data_filtered.txt"))

tax.wt <- read.csv("data/lu_species.csv") %>% 
  mutate(scientific_name = paste(species_genus, species_name)) %>% 
  rename(species = species_code) %>% 
  dplyr::select(scientific_name, species)

use.ebd <- raw.ebd %>% 
  dplyr::filter(locality_type=="H") %>% 
  mutate(source = "eBird",
         organization = "eBird",
         project=NA,
         sensor="PC",
         equipment=NA,
         singlesp="n",
         location=NA,
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
         duration = duration_minutes) %>% 
  left_join(tax.wt) %>% 
  dplyr::select(all_of(colnms))

#E. PUT TOGETHER############################

#1. Put everything together----
use <- rbind(use.wt, use.bam, use.ebd)

#2. Clip by provincial boundaries & filter to AB----
use.visit <- use %>% 
  dplyr::select(-species, -abundance, -isSeen, -isHeard) %>% 
  unique() %>% 
  dplyr::filter(!is.na(lon),
                !is.na(lat))

#Download provincial boundaries shapefile
temp <- tempfile()
download("https://www12.statcan.gc.ca/census-recensement/2021/geo/sip-pis/boundary-limites/files-fichiers/lpr_000b21a_e.zip", temp)
unzip(zipfile=temp, exdir=file.path(rootabmi, "Data"))
unlink(temp)

#Filter to Alberta
shp <- read_sf(file.path(rootabmi, "Data/lpr_000b21a_e.shp")) %>% 
  dplyr::filter(PRNAME=="Alberta")

#Create rasters (much faster than from polygon)
r <- rast(ext(shp), resolution=1000, crs=crs(shp))
ab <- rasterize(x=shp, y=r, field="PRNAME")

#Extract raster value
#Identify topsecret locations
visit.ab <- use.visit %>% 
  st_as_sf(coords=c("lon", "lat"), crs=4326) %>% 
  st_transform(crs=crs(ab)) %>% 
  vect() %>% 
  extract(x=ab) %>% 
  cbind(use.visit) %>% 
  dplyr::filter(PRNAME=="Alberta") %>% 
  dplyr::select(colnames(use.visit)) %>% 
  mutate(topsecret = ifelse((buffer > 0 & organization=="ABMI"), 1, 0))

#apply to bird data
use.ab <- use %>% 
  inner_join(visit.ab)

#4. Remove duplicate surveys----

#4a. Investigate----
#Take out buffered locations
visit.dup <- visit.ab %>% 
  dplyr::filter(buffer==0) %>% 
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
#1. JOSM points that have human & ARU data
#2. Duplicate processing of recordings

#4b. Identify JOSM points with human & ARU data----
#Keep the ARU ones because frankly they're probably better data
visit.josm <- visit.dup %>% 
  dplyr::filter(project %in% c("JOSM ECCC Cause and Effect Monitoring for Landbirds 2013",
                               "Oil Sands Monitoring CWS Prairie Region 2012",
                               "Oil Sands Monitoring CWS Prairie Region 2013",
                               "Oilsands Monitoring CWS Prairie Region 2014"))

visit.josm.remove <- visit.josm %>% 
  dplyr::filter(sensor=="PC")

#4c. Identify recordings that have been processed twice----
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

#4d. Filter out duplicates from dataset----
#remove surveys before 1993 (first year with substantial data)
visit <- visit.ab %>% 
  anti_join(visit.josm.remove) %>% 
  anti_join(visit.mindur.remove) %>% 
  anti_join(visit.eqdur.remove) %>% 
  dplyr::filter(year >= 1993)

bird <- use.ab %>% 
  anti_join(visit.josm.remove) %>% 
  anti_join(visit.mindur.remove) %>% 
  anti_join(visit.eqdur.remove)

#F. IDENTIFY LOCATIONS FOR COVARIATE EXTRACTION####
location <- visit %>% 
  dplyr::select(source, sensor, location, buffer, lat, lon, year, topsecret) %>% 
  unique()

write.csv(location, "data/birds_ab_locations_V2.csv", row.names = FALSE)

#G. SAVE!#############################
save(location, visit, bird, file=file.path(rootabmi, "data/Harmonized.Rdata"))
  
#H. COMPARE############################
load(file.path(rootabmi, "data/ab-birds-all-2020-09-23.Rdata"))

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

ggsave(filename="figs/datasetversionN.jpeg", width=6, height=4)

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
