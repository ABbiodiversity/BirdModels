# ---
# title: ABMI models - filter eBird data
# author: Elly Knight
# created: May 6, 2024
# ---

#NOTES################################

#PURPOSE: This script wrangles raw eBird datasets for combination with wildtrax data.

#raw eBird data is downloaded from the eBird interface at https://ebird.org/data/download/ebd prior to wrangling with the "auk" package and will require a request for access. Use the custom download tool to download only the datasets for Canada and the US instead of the global dataset. Note you will also need the global sampling file to use the auk package for zero filling.

#raw eBird data omits Great Grey Owl & Northern Hawk Owl as sensitive species (https://support.ebird.org/en/support/solutions/articles/48000803210?b_id=1928&_gl=1*xq054u*_ga*ODczMTUyMjcuMTY2OTE0MDI4Ng..*_ga_QR4NVXZ8BM*MTY2OTE0MDI4NS4xLjEuMTY2OTE0MDM3OC4zNS4wLjA.&_ga=2.147122167.150058226.1669140286-87315227.1669140286) and should not be used for modelling these two species.

#wrangling eBird data with the auk package requires installation of AWK on windows computers. Please see #https://cornelllabofornithology.github.io/auk/articles/auk.html.

#eBird data has not been zerofilled because there was no species filtering done and we are assuming that all stationary counts have at least 1 bird observed.

#PREAMBLE############################

#1. Load packages----

library(tidyverse) #basic data wrangling
library(auk) #eBird wrangling

#2. Set root path for data on google drive----

root <- "G:/Shared drives/ABMI_RHedley/Projects/BirdModels"

#FILTER DATA###############

#1. Set ebd path----
auk_set_ebd_path(file.path(root, "Data/ebd/ebd_CA-AB_relOct-2022"), overwrite=TRUE)

#2. Define filters----
filters <- auk_ebd(file="ebd_CA-AB_relOct-2022.txt") %>%
    auk_protocol("Stationary") %>%
    auk_duration(c(0, 10)) %>%
    auk_complete()

#3. Filter data----
#select columns to keep
filtered <- auk_filter(filters, file=file.path(root, "Data/ebd_data_filtered.txt"), overwrite=TRUE,
                       keep = c("group identifier", "sampling_event_identifier", "scientific name", "common_name", "observation_count", "latitude", "longitude", "locality_type", "observation_date", "time_observations_started", "observer_id", "duration_minutes"))

#Temporary: read from Elly's BAM version while I await eBird data request.
library(sf)
filtered <- data.table::fread(file = file.path(root, "Data/ebd/03_ebd_filtered_CA_Jan-2026.txt"))
filtered <- st_as_sf(filtered, coords=c("LONGITUDE", "LATITUDE"), crs=4326)
ab <- st_read(file.path(root, "Data/gis/lpr_000b21a_e.shp"))
ab <- ab[ab$PRNAME == 'Alberta',]
ab <- st_transform(ab, st_crs(filtered))
filtered <- st_filter(filtered, ab)
filtered <- st_drop_geometry(filtered)
write.table(filtered, file.path(root, "Data/ebd/ebd_data_filtered_temp.txt"), row.names = F)
