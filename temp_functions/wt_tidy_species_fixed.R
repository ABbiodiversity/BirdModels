
#' Filter species from a report
#'
#' @description This function filters the species provided in WildTrax reports to only the groups of interest. The groups available for filtering are mammal, bird, amphibian, abiotic, insect, and unknown. Zero-filling functionality is available to ensure all surveys are retained in the dataset if no observations of the group of interest are available.
#'
#' @param data WildTrax main report or tag report from the `wt_download_report()` function.
#' @param remove Character; groups to filter from the report ("mammal", "bird", "amphibian", "abiotic", "insect", "human", "unknown"). Defaults to retaining bird group only.
#' @param zerofill Logical; indicates if zerofilling should be completed. If TRUE, unique surveys with no observations after filtering are added to the dataset with "NONE" as the value for species_code and/or species_common_name. If FALSE, only surveys with observations of the retained groups are returned. Default is TRUE.
#'
#' @import dplyr
#' @export
#'
#' @examples
#' \dontrun{
#' dat.tidy <- wt_tidy_species(dat, remove=c("mammal", "unknown"), zerofill = T)
#' }
#' @return A dataframe identical to input with observations of the specified groups removed.

wt_tidy_species_fixed_aru <- function(data,
                                  remove = "",
                                  zerofill = TRUE) {
  
  if (is.null(remove)) {
    stop('Not removing any species')
  }
  
  if (any(!(remove %in% c("mammal", "bird", "amphibian", "abiotic", "insect", "human", "unknown")))) {
    stop("Select remove options from bird, mammal, amphibian, abiotic, insect, human or unknown.")
  }
  
  if('bird' %in% remove){
    message('Note: By removing birds, you will not be able to use wt_qpad_offsets since QPAD offsets are only available for birds.')
  }
  
  #Convert to the sql database labels for species class
  remove <- case_when(remove=="mammal" ~ "MAMMALIA",
                      remove=="amphibian" ~ "AMPHIBIA",
                      remove=="abiotic" ~ "ABIOTIC",
                      remove=="insect" ~ "INSECTA",
                      remove=="bird" ~ "AVES",
                      remove=="human" ~ "HUMAN ACTIVITY",
                      remove=="unknown" ~ "unknown",
                      remove=="" ~ remove)
  
  .species <- wt_get_species()
  
  #Get the species codes for what you want to filter out
  species.remove <- .species |>
    filter(species_class %in% remove)
  
  #Add the unknowns if requested
  if("unknown" %in% remove){
    species.remove <- .species %>%
      filter(substr(species_common_name, 1, 12) == "Unidentified") %>%
      rbind(species.remove)
  }
  
  #Remove those codes from the data
  filtered <- data |>
    filter(!(species_code %in% species.remove$species_code))
  
  #if you don't need nones, remove other NONEs & return the filtered object
  if(zerofill==FALSE){
    
    filtered.sp <- filter(filtered, species_code!="NONE")
    
    return(filtered.sp)
  }
  
  #if you do need nones, add them
  if(zerofill==TRUE){
    
    #first identify the unique visits (task_id) ensure locations are included for proper join
    visit <- data |>
      select(organization, project_id, location, latitude, longitude, location_id, recording_date_time, task_id, task_duration) |>
      distinct()
    
    #see if there are any that have been removed
    none <- suppressMessages(anti_join(visit, filtered)) |>
      mutate(species_code = "NONE",
             species_common_name = "NONE",
             species_scientific_name = NA) #Replaced with NA to align with other NONE tags.
    
    #add to the filtered data
    filtered.none <- suppressMessages(full_join(filtered, none)) |>
      arrange(organization, project_id, location, recording_date_time)
    
    #return the filtered object with nones added
    return(filtered.none)
    
  }
  
}


#' Filter species from a report
#'
#' @description This function filters the species provided in WildTrax reports to only the groups of interest. The groups available for filtering are mammal, bird, amphibian, abiotic, insect, and unknown. Zero-filling functionality is available to ensure all surveys are retained in the dataset if no observations of the group of interest are available.
#'
#' @param data WildTrax main report or tag report from the `wt_download_report()` function.
#' @param remove Character; groups to filter from the report ("mammal", "bird", "amphibian", "abiotic", "insect", "human", "unknown"). Defaults to retaining bird group only.
#' @param zerofill Logical; indicates if zerofilling should be completed. If TRUE, unique surveys with no observations after filtering are added to the dataset with "NONE" as the value for species_code and/or species_common_name. If FALSE, only surveys with observations of the retained groups are returned. Default is TRUE.
#'
#' @import dplyr
#' @export
#'
#' @examples
#' \dontrun{
#' dat.tidy <- wt_tidy_species(dat, remove=c("mammal", "unknown"), zerofill = T)
#' }
#' @return A dataframe identical to input with observations of the specified groups removed.

wt_tidy_species_fixed_pc <- function(data,
                                  remove = "",
                                  zerofill = TRUE) {
  
  if (is.null(remove)) {
    stop('Not removing any species')
  }
  
  if (any(!(remove %in% c("mammal", "bird", "amphibian", "abiotic", "insect", "human", "unknown")))) {
    stop("Select remove options from bird, mammal, amphibian, abiotic, insect, human or unknown.")
  }
  
  #Rename survey_id to task_id, survey_date to recording_date_time,
  #and add empty task_duration
  data <- data |>
    select(-recording_date_time) |>
    rename(task_id=survey_id,
           recording_date_time=survey_date)
  #If task_duration does not exist, add placeholder with NA.
  #This is needed because task_duration was added to list of columns below.
  if(!'task_duration' %in% colnames(data)) {data$task_duration <- NA}
  
  if('bird' %in% remove){
    message('Note: By removing birds, you will not be able to use wt_qpad_offsets since QPAD offsets are only available for birds.')
  }
  
  #Convert to the sql database labels for species class
  remove <- case_when(remove=="mammal" ~ "MAMMALIA",
                      remove=="amphibian" ~ "AMPHIBIA",
                      remove=="abiotic" ~ "ABIOTIC",
                      remove=="insect" ~ "INSECTA",
                      remove=="bird" ~ "AVES",
                      remove=="human" ~ "HUMAN ACTIVITY",
                      remove=="unknown" ~ "unknown",
                      remove=="" ~ remove)
  
  .species <- wt_get_species()
  
  #Get the species codes for what you want to filter out
  species.remove <- .species |>
    filter(species_class %in% remove)
  
  #Add the unknowns if requested
  if("unknown" %in% remove){
    species.remove <- .species %>%
      filter(substr(species_common_name, 1, 12) == "Unidentified") %>%
      rbind(species.remove)
  }
  
  #Remove those codes from the data
  filtered <- data |>
    filter(!(species_code %in% species.remove$species_code))
  
  #if you don't need nones, remove other NONEs & return the filtered object
  if(zerofill==FALSE){
    
    filtered.sp <- filter(filtered, species_code!="NONE")
    
    #Translate point count field names back
    filtered.sp <- filtered.sp |>
      rename(survey_id=task_id,
             survey_date = recording_date_time)
    
    return(filtered.sp)
  }
  
  #if you do need nones, add them
  if(zerofill==TRUE){
    
    #first identify the unique visits (task_id) ensure locations are included for proper join
    visit <- data |>
      select(organization, project_id, location, latitude, longitude, location_id, recording_date_time, task_id, survey_duration_method) |>
      distinct()
    
    #see if there are any that have been removed
    none <- suppressMessages(anti_join(visit, filtered)) |>
      mutate(species_code = "NONE",
             species_common_name = "NONE",
             species_scientific_name = NA) #Replaced with NA to align with other NONE tags.
    
    #add to the filtered data
    filtered.none <- suppressMessages(full_join(filtered, none)) |>
      arrange(organization, project_id, location, recording_date_time)
    
    #Translate point count field names back
    filtered.none <- filtered.none |>
      rename(survey_id=task_id,
             survey_date = recording_date_time)
    
    #return the filtered object with nones added
    return(filtered.none)
    
  }
}




