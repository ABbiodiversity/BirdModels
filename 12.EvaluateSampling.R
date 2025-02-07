# ---
# title: ABMI models - extrapolation analysis of covariates
# author: Elly Knight
# created: Dec 5, 2024
# ---

#NOTES################################

#Adapted from Anna Drake's script from the BAM modelling workflow
# This code uses dsextra to quantify extrapolation
# This package calls the "ExDet" function which needs it's tolerances adjusted to avoid singular matrices (rounding low values down)

#PREAMBLE############################

#1. Load packages----
library(tidyverse) #basic data wrangling
library(dsmextra) #extrapolation analysis

#2. Set root path for data on google drive----
root <- "G:/My Drive/ABMI/Projects/BirdModels"

#3. Restrict scientific notation----
options(scipen = 99999)

#4. Load data----
load(file.path(root, "Data", "Stratified.Rdata"))

#5. Get the kgrid & backfill----
load(file.path(root, "Data", "gis", "kgrid_2.2.Rdata"))
load(file.path(root, "Data", "gis", "backfillV7_w2w_2021HFI.Rdata"))
load(file.path(root, "Data", "gis", "kgrid_2.2_birds.Rdata"))

#6. Define projection system----
crs.all <- sp::CRS("EPSG:3400")

#CLIMATE MODELS###########

#1. Get list of covariates----
source("modelling2.0/00.ClimateModels.R")
vars.c <- modelsclimate |>
  unlist() |>
  lapply(function(f) all.vars(f)[-1]) |>
  unlist() |>
  unique() |>
  data.frame() |>
  purrr::set_names("var") |>
  mutate(mic = row_number())

#2. Get the samples----
samp.c <- covs |>
  dplyr::select("Easting", "Northing", all_of(vars.c$var)) |>
  rename(x=Easting, y=Northing)

#3. Wrangle kgrid-----
grid.c <- kgrid |>
  dplyr::select("Easting", "Northing", all_of(vars.c$var)) |>
  rename(x=Easting, y=Northing)

#4. Calculate extrapolation----
extrap.c <- compute_extrapolation(samples = samp.c,
                                  covariate.names = vars.c$var,
                                  prediction.grid = grid.c,
                                  coordinate.system = crs.all,
                                  resolution = 1000)

#5. Wrangle----
out.c <- extrap.c$data$all |>
  left_join(vars.c) |>
  dplyr::select(x, y, ExDet, mic_combinatorial, var) |>
  rename(Type2_climate = mic_combinatorial,
         Type1_climate = var)

#6. Plot univariate extrapolation----
ggplot(out.c) +
  geom_raster(aes(x=x, y=y, fill=Type1))

#7. Plot multivariate extrapolation----
ggplot(out.c) +
  geom_raster(aes(x=x, y=y, fill=Type2))

#NORTH MODELS#############

#1. Get list of covariates----
source("modelling2.0/00.NorthModels.R")
vars.n <- modelsnorth |>
  unlist() |>
  lapply(function(f) all.vars(f)[-1]) |>
  unlist() |>
  unique() |>
  data.frame() |>
  purrr::set_names("var") |>
  mutate(mic = row_number()) |>
  dplyr::filter(!var %in% c("method", "vegc"))

#2. Get the samples-----
samp.n <- covs |>
  dplyr::select("Easting", "Northing", all_of(vars.n$var)) |>
  rename(x=Easting, y=Northing)

#3. Wrangle kgrid-----
grid.n <- kgrid_birds |>
  dplyr::select("Easting", "Northing", all_of(vars.n$var)) |>
  rename(x=Easting, y=Northing)

#4. Calculate extrapolation----
extrap.n <- compute_extrapolation(samples = samp.n,
                                  covariate.names = vars.n$var,
                                  prediction.grid = grid.n,
                                  coordinate.system = crs.all,
                                  resolution = 1000)

#5. Wrangle----
out.n <- extrap.n$data$all |>
  left_join(vars.n) |>
  dplyr::select(x, y, ExDet, mic_combinatorial, var) |>
  rename(Type2_north = mic_combinatorial,
         Type1_north = var)

#6. Plot univariate extrapolation----
ggplot(out.n) +
  geom_raster(aes(x=x, y=y, fill=Type1))

#7. Plot multivariate extrapolation----
ggplot(out.n) +
  geom_raster(aes(x=x, y=y, fill=Type2))

#PACKAGE#########

#1. Put together----
out <- full_join(out.n, out.c)

#2. Save----
write.csv(out, "Extrapolation_North&Climate.csv", row.names = FALSE)
