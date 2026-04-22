
#Version taken from Elly's Feb 2 PR, which was subsequently over-written.


#' Convert to a wide survey by species dataframe
#'
#' @description This function converts a long-formatted report into a wide survey by species dataframe of abundance values.
#'
#' @param data WildTrax main report or tag report from the `wt_download_report()` function.
#' @param sound Character; vocalization type(s) to retain ("all", "Song", "Call", "Non-vocal"). Can be used to remove certain types of detections. Defaults to "all" (i.e., no filtering).
#'
#' @import dplyr
#' @importFrom tidyr pivot_wider
#' @export
#'
#' @examples
#' \dontrun{
#' dat.tidy <- wt_tidy_species(dat)
#' dat.tmtt <- wt_replace_tmtt(dat.tidy)
#' dat.wide <- wt_make_wide(dat.tmtt, sound="all")
#' }
#' @return A dataframe identical to input with observations of the specified groups removed.

wt_make_wide_fixed <- function(data, sound="all"){
  
  #Steps for ARU data
  if(!"survey_url" %in% colnames(data)){
    
    #Filter to first detection per individual
    summed <- data |>
      group_by(organization, project_id, location, recording_date_time, task_method, task_is_complete, observer_id, species_code, species_common_name, individual_order) |>
      mutate(first = min(detection_time)) |>
      ungroup() |>
      filter((species_code != "NONE" & detection_time == first) | species_code == "NONE")
    
    #Remove undesired sound types
    if(!"all" %in% sound){
      sound <- gsub("\\b(\\w)", "\\U\\1", tolower(sound), perl = TRUE)
      summed <- filter(summed, vocalization %in% sound)
    }
    
    #Make it wide
    wide <- summed |>
      mutate(abundance = case_when(is.na(abundance) & species_code == "NONE" ~ "0", grepl("^C",  abundance) ~ NA_character_, TRUE ~ as.character(abundance)) |> as.numeric()) |>
      pivot_wider(id_cols = organization:task_method,
                  names_from = "species_code",
                  values_from = "abundance",
                  values_fn = sum,
                  values_fill = 0,
                  names_sort = TRUE)
    
  }
  
  #Steps for point count data
  if("survey_url" %in% colnames(data)){
    
    #Make it wide and return field names to point count format
    wide <- data |>
      rename(abundance = individual_count) |>
      mutate(abundance = case_when(is.na(abundance) & species_code == "NONE" ~ "0", grepl("^C",  abundance) ~ NA_character_, TRUE ~ as.character(abundance)) |> as.numeric()) |>
      pivot_wider(id_cols = organization:survey_duration_method,
                  names_from = "species_code",
                  values_from = "abundance",
                  values_fn = sum,
                  values_fill = 0,
                  names_sort = TRUE)
    
  }
  
  return(wide)
  
}