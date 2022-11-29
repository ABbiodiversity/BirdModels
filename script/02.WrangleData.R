library(tidyverse) #basic data wrangling
library(lubridate) #date wrangling
library(downloader) #download files
library(sf) #for spatial wrangling
library(terra) #for spatial wrangling

#Load previous dataset----
load("data/ab-birds-all-2022-03-09.RData")

#1. Load data downloaded from WildTrax----
load("data/wildtrax_data_2022-11-15.RData")

#2. Wrangle----
#tidy column names between sensor types
#tidy method data
#add sensor type
#wrangle date
dat <- raw  %>% 
  data.frame() %>% 
  # mutate(project = ifelse(is.na(project), project_name, project),
  #        speciesCode = ifelse(is.na(speciesCode), species_code, speciesCode), 
  #        date = ifelse(is.na(date), recording_date, date)) %>% 
  # dplyr::select(-project_name, -species_code, -recording_date) %>% 
  separate(method, into=c("duration", "method"), remove=TRUE) %>% 
  mutate(duration = ifelse(!method %in% c("1SPM", "1SPT", "None"), method, duration),
         method = ifelse(!method %in% c("1SPM", "1SPT", "None"), "None", method),
         date = ymd_hms(date),
         year = year(date)) %>% 
  left_join(projects %>% 
              rename(sensor = sensorId) %>% 
              dplyr::select(project, sensor))

#3. Filter data----
#remove unstructured
#remove year = 1900
dat.method <- dat %>% 
  mutate(durationMethod = ifelse(is.na(durationMethod), "ARU", durationMethod)) %>% 
  dplyr::filter(durationMethod!="UNKNOWN",
                !(method=="None" & sensor=="ARU"),
                year!=1900)

#4. Identify unique survey locations----
visit <- dat.method %>% 
  dplyr::select(organization, project, location, latitude, longitude, buffer, year) %>% 
  unique() %>% 
  dplyr::filter(!is.na(latitude))

#5. Remove data outside of alberta----

#5a. Get provincial shapefile----
# temp <- tempfile()
# download("https://extranet.gov.ab.ca/srd/geodiscover/srd_pub/boundaries/census/Alberta_Census_Boundaries_SHP.zip", temp)
# unzip(zipfile=temp, exdir="data/ABshp")
# unlink(temp)

#5b.Read in & wrangle shapefile----
ab <- read_sf("data/ABshp/Data/AB_CD_2021.shp") %>% 
  st_union() %>% 
  st_transform(crs=4326) %>% 
  st_as_sf() %>% 
  mutate(prov = "AB") 

#5c. Make survey data spatial----
visit.sf <- visit %>% 
  st_as_sf(coords=c("longitude", "latitude"), crs=4326)

#5d. Intersect with shapefile----
visit.ab <- visit.sf %>% 
  st_intersection(ab)

#6. Identify buffered points & tidy to save----
#remove buffered points that aren't ABMI
visit.use <- visit.ab %>% 
  st_coordinates() %>% 
  data.frame() %>% 
  data.frame(visit.ab) %>% 
  dplyr::select(-geometry) %>% 
  mutate(buffer = ifelse(buffer=="0m (Buffered location)", NA, buffer)) %>% 
  mutate(topsecret = ifelse((!is.na(buffer) & organization=="ABMI"), 1, 0)) %>% 
  dplyr::filter(!(!is.na(buffer) & topsecret!=1))

#7. Save for GIS extraction----
write.csv(visit.use, "data/birds_ab_locations.csv", row.names = FALSE)
