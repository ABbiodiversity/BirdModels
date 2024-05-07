# ---
# title: ABMI models - get WildTRax data
# author: Elly Knight
# created: November 29, 2022
# ---

#NOTES################################

#The "projectInstructions.csv" file is a list of all projects currently in WildTrax should not be used in ABMI models (instructions=="DO NOT USE"). This file should be updated for each iteration of in collaboration with Erin Bayne. Future versions of this spreadsheet can hopefully be derived by a combination of organization and a google form poll for consent from other organizations.

#The "projectInstructions.csv" file also contains information on which ARU projects are processed for a single species or taxa (instructions=="DO NOT USE") and therefore those visits should only be used for models for the appropriate taxa. This file should be updated for each iteration of the national models in collaboration with Erin Bayne. These projects are currently not included in the models.

#There are a handful of projects that are not downloading properly via wildRtrax. An issue is open on this. These projects are listed in the error.log object. These files should be downloaded manually.

#PREAMBLE############################

#1. Load packages----

library(tidyverse) #basic data wrangling
library(wildRtrax) #to download data from wildtrax
library(data.table) #for binding lists into dataframes

#2. Set root path for data on google drive----

root <- "G:/My Drive/ABMI/Projects/BirdModels"

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

instructions <- read.csv(file.path(root, "Data", "projectInventory", "projectInstructions.csv")) %>%  
  dplyr::filter(instruction!="ABMI ONLY")

projects.use <- projects %>% 
  dplyr::filter(organization!="BU-TRAINING",
                !project_id %in% instructions$project_id)

#3. Loop through projects to download data----
aru.list <- list()
pc.list <- list()
error.log <- data.frame()
for(i in 1:nrow(projects.use)){
  
  #authenticate each time because this loop takes forever
  wt_auth()
  
  #Do each sensor type separately because the reports have different columns and we need different things for each sensor type
  if(projects.use$sensor[i]=="ARU"){
    
    dat.try <- try(wt_download_report(project_id = projects.use$project_id[i], sensor_id = projects.use$sensor[i], weather_cols = F, report = "main"))
    
    if(class(dat.try)=="data.frame"){
      aru.list[[i]] <- dat.try
    }
    
  }
  
  if(projects.use$sensor[i]=="PC"){
    
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

#DEAL WITH ERRORED PROJECT############

#1. Go download error projects from wildtrax.ca----
error.log %>% 
  inner_join(projects.use %>% 
              mutate(i = row_number())) %>% 
  View()

#2. Read in error projects----
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

#PUT TOGETHER##############

#1. Collapse lists----
#Take out ARU projects with "None" method
aru.wt <- rbindlist(aru.list[], fill=TRUE)  %>% 
  rbind(aru.error, fill=TRUE) %>%
  left_join(projects.use %>% 
              dplyr::rename(project_status = status)) %>% 
  dplyr::filter(task_method!="None")

pc.wt <- rbindlist(pc.list[], fill=TRUE)  %>% 
 rbind(pc.error, fill=TRUE) %>%
  left_join(projects %>% 
              dplyr::rename(project_status = status))

#2. Save date stamped data & project list----
save(aru.wt, pc.wt, projects.use, error.log, file=paste0(root, "/Data/WildTrax/wildtrax_raw_", Sys.Date(), ".Rdata"))