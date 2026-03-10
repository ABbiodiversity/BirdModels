# ---
# title: ABMI models - get WildTRax data
# author: Elly Knight
# created: November 29, 2022
# modified: March 9, 2026
# modified by: Richard Hedley
# ---

#NOTES################################

#PURPOSE: This script downloads WildTrax data for modelling. This script should be run once per modelling version and the date-stamped raw object output from this script should be archived for reproducibility.

#The "projectInstructions.csv" file is a list of all projects currently in WildTrax should not be used in ABMI models (instructions=="DO NOT USE"). This file should be updated for each iteration of in collaboration with Erin Bayne. Future versions of this spreadsheet can hopefully be derived by a combination of organization and a google form poll for consent from other organizations.

#The "projectInstructions.csv" file also contains information on which ARU projects are processed for a single species or taxa (instructions=="DO NOT USE") and therefore those visits should only be used for models for the appropriate taxa. This file should be updated for each iteration of the national models in collaboration with Erin Bayne. These projects are currently not included in the models.

#There are a handful of projects that are not downloading properly via wildRtrax. An issue is open on this. These projects are listed in the error.log object. These files should be downloaded manually.

#FUTURE VERSION: Some components of this workflow may require revision to comply with the current version of wildRtrax and associated data cleaning functions.

#PREAMBLE############################

#1. Load packages----

library(tidyverse) #basic data wrangling
library(wildrtrax) #to download data from wildtrax
library(data.table) #for binding lists into dataframes

#2. Set root path for data on google drive----

root <- "G:/Shared drives/ABMI_RHedley/Projects/BirdModels"

#3. Login to WildTrax----
config <- "src/00.WTlogin.R"
source(config)
rm(config)

#4. Authenticate----
wt_auth()

#A. DOWNLOAD DATA FROM WILDTRAX#######################

#1. Get list of projects from WildTrax----
#sensor = PC gives all ARU and point count projects
projects <- rbind(wt_get_projects(sensor = 'PC'), wt_get_projects(sensor = 'ARU'))

#2. Filter out projects that shouldn't be used----
#nothing in BU training & all "DO NOT USE" projects in projectInventory file
#filter out 'NONE' method ARU projects later after this field is parsed out

instructions <- read.csv(file.path(root, "Data", "projectInventory", "Projects_WildTrax_Use_CheckedByEMB_2026-02-15.csv")) %>%  
  dplyr::filter(Richard.To.Use=="Yes")

projects.use <- projects %>% 
  dplyr::filter(project_id %in% instructions$project_id)

rm(instructions)

#3. Loop through projects to download data----
aru.list <- list()
pc.list <- list()
error.log <- data.frame()
for(i in 1:nrow(projects.use)){
  
  #authenticate each time because this loop takes forever
  wt_auth()
  
  #Do each sensor type separately because the reports have different columns and we need different things for each sensor type
  dat.try <- try(wt_download_report(project_id = projects.use$project_id[i], sensor_id = projects.use$project_sensor[i], reports = "main"))
  
  #If the function call returned a data frame, store the results, otherwise log the error.
  if("data.frame" %in% class(dat.try)){
    if(projects.use$project_sensor[i]=="ARU"){
      
      #If the sensor is an ARU, store in aru.list.
      aru.list[[i]] <- as.data.frame(dat.try)
      
    }
    if(projects.use$project_sensor[i]=="PC"){
      
      #If the sensor is a PC, store in pc.list.
      pc.list[[i]] <- as.data.frame(dat.try)
      
    }
  } else {
    
    #Construct error log.
    error.log <- rbind(error.log, 
                       projects.use[i,])
    
  }
  
  print(paste0("Finished dataset ", projects.use$project[i], " : ", i, " of ", nrow(projects.use), " projects"))
  
  rm(dat.try)
}

#DEAL WITH ERRORED PROJECT############

#1. Go download error projects from wildtrax.ca----
error.log %>% 
  inner_join(projects.use %>% 
               mutate(i = row_number())) %>% 
  View()

#2. Read in error projects----
error.files.aru <- list.files(file.path(root, "Data", "WildTrax", "errorFiles", "ARU"), full.names = TRUE, pattern = '*.csv')

aru.error <- data.frame()
for(i in 1:length(error.files.aru)){
  aru.error <- read.csv(error.files.aru[i]) %>% 
    rbind(aru.error)
}

error.files.pc <- list.files(file.path(root, "Data", "WildTrax", "errorFiles", "PC"), full.names = TRUE, pattern = '*.csv')

pc.error <- data.frame()
for(i in 1:length(error.files.pc)){
  pc.error <- read.csv(error.files.pc[i]) %>% 
    rbind(pc.error)
}

#Fix issue that arises due to combining character datetimes with
#posixct type in wt_download_report.
aru.error$recording_date_time <- as.POSIXct(aru.error$recording_date_time)
pc.error$recording_date_time <- as.POSIXct(pc.error$recording_date_time)

#PUT TOGETHER##############

#1. Collapse lists----
#Take out ARU projects with "None" method
aru.wt <- rbindlist(aru.list[], fill=TRUE)  %>% 
  rbind(aru.error, fill=TRUE) %>%
  left_join(projects.use) %>% 
  dplyr::filter(task_method!="None")

pc.wt <- rbindlist(pc.list[], fill=TRUE)  %>% 
 rbind(pc.error, fill=TRUE) %>%
  left_join(projects.use)

#2. Save date stamped data & project list----
save(aru.wt, pc.wt, projects.use, error.log, file=paste0(root, "/Data/WildTrax/wildtrax_raw_", Sys.Date(), ".Rdata"))

#2. Save date stamped data & project list----
save(aru.wt, pc.wt, projects.use, error.log, file=file.path(root, "Data", "WildTrax", "wildtrax_raw_2023-11-21.Rdata"))
