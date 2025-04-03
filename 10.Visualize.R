# ---
# title: ABMI models - plot predictions
# author: Elly Knight
# created: Sept 25, 2024
# ---

#NOTES################################

#PREAMBLE############################

#1. Load packages----
library(tidyverse) #basic data wrangling
library(Matrix) #sparse matrices
library(terra) #rasters
library(sf) #polygons
library(gridExtra) #plotting
library(ebirdst) #ebird predictions
library(auk) #ebird codes

#2. Set root path for data on google drive----
root <- "G:/Shared drives/ABMI_ECKnight/Projects/BirdModels"

#3. Restrict scientific notation----
options(scipen = 99999)

#4. Load the new coefficients----
new <- new.env()
load(file.path(root, "Results", "Birds2024.RData"), envir=new)

#5. Load the old coefficients----
old <- new.env()
load(file.path(root, "Results", "Archive", "2023", "Birds2023.RData"), envir=old)

#6. Load the data----
load(file.path(root, "Data", "Stratified.Rdata"))

#7. Load the bird kgrid----
load(file.path(root, "Data", "gis", "kgrid_2.2_birds.Rdata"))

#KGRID PREDICTIONS#######

#1. Get the eBird codes----
codes <- get_ebird_taxonomy()

#2. Species list----
spp <- inner_join(new$birds$species |>
                    dplyr::filter(ModelNorth==TRUE | ModelSouth==TRUE) |>
                    dplyr::select(Comments),
                  old$birds$species |>
                    dplyr::filter(ModelNorth==TRUE | ModelSouth==TRUE) |>
                    dplyr::select(Comments, ScientificName)) |>
  rename(species = Comments) |>
  mutate(scientific_name = ScientificName,
         ScientificName = str_replace_all(ScientificName, " ", "-"),
         scientific_name = case_when(scientific_name=="Ammodramus bairdii" ~ "Centronyx bairdii",
                                     scientific_name=="Picoides pubescens" ~ "Dryobates pubescens",
                                     scientific_name=="Picoides villosus" ~ "Dryobates villosus",
                                     scientific_name=="Ammodramus nelsoni" ~ "Ammospiza nelsoni",
                                     scientific_name=="Oreothlypis celata" ~ "Leiothlypis celata",
                                     scientific_name=="Carpodacus purpureus" ~ "Haemorhous purpureus",
                                     scientific_name=="Regulus calendula" ~ "Corthylio calendula",
                                     scientific_name=="Oreothlypis peregrina" ~ "Leiothlypis peregrina",
                                     !is.na(scientific_name) ~ scientific_name)) |>
  left_join(codes)

#3. Get the eBird predictions----
#only need to do this once
# for(i in 1:nrow(spp)){
#
#     ebirdst_download_status(species = spp$species_code[i],
#                             download_abundance = TRUE,
#                             download_ranges = FALSE,
#                             force=TRUE,
#                             pattern = "abundance_seasonal_mean_3km")
# }

#4. Get the new predictions & kgrid----
load(file.path(root, "Results", "Predictions2024.RData"))
load(file.path(root, "Data", "gis", "kgrid_2.2.Rdata"))

#5. Get an alberta border for cropping----
ab <- read_sf(file.path(root, "Data", "gis", "lpr_000b21a_e.shp")) |>
  dplyr::filter(PRNAME=="Alberta") |>
  st_transform(crs(load_raster(spp$species_code[1], product="abundance", period="seasonal", resolution="3km"))) |>
  vect()

#6. Set up species loop----
out <- data.frame()
for(i in 1:nrow(spp)){
  
  #7. Load the old raster----
  pred.old <- try(rast(file.path(root, "PreviousModelRasters", paste0(spp$ScientificName[i], ".tif")))$Current)
  
  #8. Get the new rasters----
  pred.i <- predictions[[spp[i,]$species]]
  
  if(class(pred.i)=="try-error"){next}
  
  #9. Get the landcover prediction----
  pred.new <- pred.i |>
    left_join(kgrid, by="LinkID") |>
    dplyr::select(X, Y, Provincial) |>
    raster::rasterFromXYZ(crs=3400) |>
    rast()
  
  #10. Get the climate prediction----
  clim <- pred.i |>
    left_join(kgrid, by="LinkID") |>
    dplyr::select(X, Y, Climate) |>
    raster::rasterFromXYZ(crs=3400) |>
    rast()

  #11. Load eBird and crop----
  if(class(pred.old)!="try-error"){
    pred.ebd.raw <-load_raster(spp$species_code[i], period="seasonal", resolution="3km", product="abundance")
    if("breeding" %in% names(pred.ebd.raw)){
      pred.ebd <- pred.ebd.raw[["breeding"]] |>
        crop(ab, mask=TRUE) |>
        project(crs(pred.old)) |>
        resample(pred.old, method="bilinear")
    }
    if("resident" %in% names(pred.ebd.raw)){
      pred.ebd <- pred.ebd.raw[["resident"]] |>
        crop(ab, mask=TRUE) |>
        project(crs(pred.old)) |>
        resample(pred.old, method="bilinear")
    }
  }
  
  #12. Stack them----
  if(class(pred.old)!="try-error"){
    pred <- c(pred.old, pred.new, clim, pred.ebd)
    names(pred) <- c("old", "new", "climate", "ebird")
  } else {
    pred <- c(pred.new, clim)
    names(pred) <- c("new", "climate")
  }

  #14. Pearson correlation----
  cor <- layerCor(pred, "cor")
  
  #15. Difference in population size----
  pop <- global(pred, fun="sum", na.rm=TRUE)
  
  #16. Wrangle----
  if(class(pred.old)!="try-error"){
    out <- data.frame(species = spp$species[i],
                      cor_oldnew = cor$correlation["old", "new"],
                      cor_oldebd = cor$correlation["old", "ebird"],
                      cor_newebd = cor$correlation["new", "ebird"],
                      pop.old = pop["old",],
                      pop.new = pop["new",]) |>
      rbind(out)
  } else {
    out <- data.frame(species = spp$species[i],
                      cor_oldnew = NA,
                      cor_oldebd = NA,
                      cor_newebd = NA,
                      pop.old = NA,
                      pop.new = pop["new",]) |>
      rbind(out)
  }

  #17. Get a bootstrap of data points to plot----
  pts <- bird |>
    dplyr::filter(surveyid %in% boot$'1') |>
    dplyr::select(surveyid, spp$species[i]) |>
    left_join(covs |>
                dplyr::select(surveyid, Easting, Northing))
  colnames(pts) <- c("surveyid", "count", "Easting", "Northing")
  
  #17. Plot----
  pred.plotdf <- as.data.frame(pred, xy=TRUE) |>
    pivot_longer(-c(x, y), names_to="layer", values_to="value")
  
  if(class(pred.old)!="try-error"){
    old.plot <- ggplot(pred.plotdf |> dplyr::filter(layer=="old")) +
      geom_raster(aes(x=x, y=y, fill=value)) +
      scale_fill_viridis_c() +
      theme(legend.position = "bottom",
            axis.title = element_blank(),
            axis.text = element_blank()) +
      ggtitle("Old prediction")
    
    ebd.plot <- ggplot(pred.plotdf |> dplyr::filter(layer=="ebird")) +
      geom_raster(aes(x=x, y=y, fill=value)) +
      scale_fill_viridis_c() +
      theme(legend.position = "bottom",
            axis.title = element_blank(),
            axis.text = element_blank()) +
      ggtitle("eBird prediction")
    
  }

  new.plot <- ggplot(pred.plotdf |> dplyr::filter(layer=="new")) +
    geom_raster(aes(x=x, y=y, fill=value)) +
    scale_fill_viridis_c() +
    theme(legend.position = "bottom",
          axis.title = element_blank(),
          axis.text = element_blank()) +
    ggtitle("New prediction - packaged")
  
  clim.plot <- ggplot(pred.plotdf |> dplyr::filter(layer=="climate")) +
    geom_raster(aes(x=x, y=y, fill=value)) +
    scale_fill_viridis_c() +
    theme(legend.position = "bottom",
          axis.title = element_blank(),
          axis.text = element_blank()) +
    ggtitle("Climate prediction - packaged")
  
  pts.plot <- ggplot(pred.plotdf |> dplyr::filter(layer=="ebird")) +
    geom_raster(aes(x=x, y=y), fill="white") +
    geom_point(data=dplyr::filter(pts, count > 0),
               aes(x=Easting, y=Northing, colour=count)) +
    scale_colour_viridis_c() +
    theme(legend.position = "bottom",
          axis.title = element_blank(),
          axis.text = element_blank()) +
    ggtitle("Detections")
  
  if(class(pred.old)!="try-error"){
    ggsave(grid.arrange(clim.plot, new.plot, old.plot, ebd.plot, pts.plot,
                        ncol=5, nrow=1, top=spp$species[i]),
           filename = file.path(root, "Results", "Plots", "Comparisons", paste0(spp$species[i], ".jpeg")),
           width = 14, height = 6)
  } else {
    ggsave(grid.arrange(clim.plot, new.plot, pts.plot,
                        ncol=3, nrow=1, top=spp$species[i]),
           filename = file.path(root, "Results", "Plots", "Comparisons", paste0(spp$species[i], ".jpeg")),
           width = 9, height = 6)
  }

  cat("Finished", spp$species[i], "predictions :", i, "of", nrow(spp), "\n")
  
}

write.csv(out, file.path(root, "Results", "ModelComparison.csv"), row.names = FALSE)
