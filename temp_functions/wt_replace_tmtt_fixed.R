

#' Replace 'TMTT' abundance with model-predicted values
#'
#' @description This function uses a lookup table of model-predicted values to replace 'TMTT' entries in listener-processed ARU data from WildTrax. The model-predicted values were produced using estimated abundances for 'TMTT' entries in mixed effects model with a Poisson distribution and random effects for species and observer.
#'
#' @param data Dataframe of WildTrax observations, for example the summary report.
#' @param calc Character; method to convert model predictions to integer ("round", "ceiling", or "floor"). See `?round()` for details.
#'
#' @import dplyr
#' @export
#'
#' @examples
#' \dontrun{
#' dat.tmtt <- wt_replace_tmtt(dat, calc="round")
#' }
#' @return A dataframe identical to input with 'TMTT' entries in the abundance column replaced by integer values.

wt_replace_tmtt_fixed <- function(data, calc="round"){
  
  if(!"recording_date_time" %in% colnames(data)){
    stop("The `wt_replace_tmtt` function only works on data from the ARU sensor")
  }
  
  check_none <- data |>
    select(species_code) |>
    distinct() |>
    pull()
  
  if (length(check_none) == 1 && check_none == 'NONE') {
    stop('There are no species in this project')
  }
  
  .tmtt <- readRDS(system.file("extdata", "tmtt_predictions.rds", package="wildrtrax"))
  
  dat.tmtt <- mutate(data, id = row_number())
  
  # only TMTT rows for replacement
  dat.tmt <- dat.tmtt |> filter(abundance %in% c("TMTT", "TNPE"))
  
  if(nrow(dat.tmt) > 0){
    dat.tmt <- dat.tmt |>
      mutate(species_code_temp = species_code) |>
      mutate(
        species_code = ifelse(species_code %in% .tmtt$species_code, species_code, "species"),
        observer_id = as.integer(ifelse(observer_id %in% .tmtt$observer_id, observer_id, 0))
      ) |>
      inner_join(.tmtt |> select(species_code, observer_id, pred),
                 by = c("species_code", "observer_id")) |>
      mutate(
        abundance = case_when(
          calc == "round"   ~ round(pred),
          calc == "ceiling" ~ ceiling(pred),
          calc == "floor"   ~ floor(pred),
          TRUE ~ NA_real_
        )
      ) |>
      mutate(species_code = species_code_temp) |>
      select(-pred, -species_code_temp)
  }
  
  # replace TMTT rows with predictions
  
  dat.tmtt <- dat.tmtt |>
    mutate(abundance = case_when(abundance %in% c("TMTT", "TNPE") ~ NA, TRUE ~ abundance)) |>
    mutate(abundance = as.numeric(abundance)) |>
    rows_update(dat.tmt, by = c("id")) |>
    select(-id)
  
  return(dat.tmtt)
}


