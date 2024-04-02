# ---
# title: ABMI models - wrangle data
# author: Elly Knight
# created: February 3, 2023
# ---

#NOTES################################
#The veghfsoil package will need to be downloaded from github to convert raw backfill GIS output to usable covariate data. See https://github.com/ABbiodiversity/veghfsoil.

#This script formats data as per Solymos inputs and is separate from the package data script to retain Peter's workflow.

#PREAMBLE############################

#1. Load packages----
library(tidyverse) #basic data wrangling
library(data.table) #for binding lists into dataframes
library(lubridate) #date wrangling
library(RSQLite) #read in GIS SQL data
library(veghfsoil) #wrangle GIS data
library(mefa4) #wrangle into previous data format

#2. Set root path for data on google drive----
root <- "G:/My Drive/ABMI/Projects/BirdModels/"

#3. Load harmonized set----
load(file.path(root, "Data", "1Harmonized.Rdata"))
#load(file.path(root, "Data", "ab-birds-all-2020-09-23.Rdata"))

#A. LOAD GIS EXTRACTION####

#1. Climate data----
load(file.path(root, "Data", "gis", "bird-point-climate_2024.Rdata"))

#2. Point data----
load(file.path(root, "Data", "gis", "birds-point-sites_1993-2023.Rdata"))

#3. Buffer data----
load(file.path(root, "Data", "gis", "birds-point-sites-simplified_1993-2023.Rdata"))


#B. WRANGLE GIS DATA####

#1. Make smaller buffer size object ----
#Smaller buffer size (150 M) excludes rows with buffer bands that are larger than 100 m
df1 <- df %>% 
  dplyr::filter(BUFF_DIST < 100)

#2. Make long summary----
d_long1 <- make_veghf_long(d=df1,
                          col.label = "UID_old",
                          col.veg = "Combined_ChgByCWCS",
                          col.baseyear = 1993,
                          col.hfyear = "YEAR",
                          col.soil = "Soil_Type_1",
                          unround = FALSE,
                          hf_fine = TRUE)

d_long2 <- make_veghf_long(d=df,
                           col.label = "UID_old",
                           col.veg = "Combined_ChgByCWCS",
                           col.baseyear = 1993,
                           col.hfyear = "YEAR",
                           col.soil = "Soil_Type_1",
                           unround = FALSE,
                           hf_fine = TRUE)

#3. Make wide summary----
d_wide1 <- make_veghf_wide(d=df1,
                          long_output = d_long1,
                          col.label = "UID_old",
                          col.area = "Shape_Area",
                          hf_fine = TRUE,
                          tol = 0,
                          sparse = TRUE)

d_wide2 <- make_veghf_wide(d=df,
                           long_output = d_long2,
                           col.label = "UID_old",
                           col.area = "Shape_Area",
                           hf_fine = TRUE,
                           tol = 0,
                           sparse = TRUE)


#4. Make point summary----
d_long <- make_veghf_long(d=df.pt,
                          col.label = "UID",
                          col.veg = "Combined_ChgByCWCS",
                          col.baseyear = 1993,
                          col.hfyear = "YEAR",
                          col.soil = "Soil_Type_1",
                          unround = FALSE,
                          hf_fine = TRUE)

#remove large objects
#rm(df, df1)

#C. PACKAGE LIKE PREVIOUS DATASET####

#1. Load previous dataset----
e1 <- new.env()
load(file.path(root, "Data", "Archive", "ab-birds-all-2020-09-23.RData"), envir = e1)
#`yy`: occurence data (matrix)
#`dd`: survey metadata (data frame)
#`vs0`: veg+soil current+reference at point level (data frame)
#`vc1`: veg current 150m buffer (matrix)
#`vc2`: veg current 564m buffer (matrix)
#`vr1`: veg reference 150m buffer (matrix)
#`vr2`: veg reference 564m buffer (matrix)
#`sc1`: soil current 150m buffer (matrix)
#`sc2`: soil current 564m buffer (matrix)
#`sr1`: soil reference 150m buffer (matrix)
#`sr2`: soil reference 564m buffer (matrix)

#2. Remove sites not in the GIS dataset----
visit <- location %>% 
  dplyr::filter(gisid %in% row.names(d.wide.pts$veg.current))

#check for duplicates
nrow(visit %>% 
       group_by(visitid) %>% 
       summarize(n=n()) %>% 
       dplyr::filter(n > 1))

#4. Occurrence data (yy)----
#get species list for filtering
spp <- e1$tax$code

#summarize abundance per species per survey
use$visitid <- paste0(use$gisid, "_", as.character(use$date_time))

bird.yy <- use %>% 
  dplyr::select(ALFL:YTVI)

#make wide sparse matrix
yy <- Xtab(abundance ~ visitid + species, bird.yy)

#remove none column
yy <- yy[,colnames(yy)!="NONE"]

#5. Survey data (dd)----
#add in needed columns from gis extraction
#create road column
#create cmethod column
dd <- visit %>% 
  left_join(dcli %>% 
              dplyr::select(gisid, AHM, Eref, FFP, MAP, MAT, MCMT, MWMT, Populus_tremuloides_brtpred_nofp) %>% 
              unique()) %>% 
  left_join(df.pt %>%
              dplyr::select(gisid, NRNAME, NSRNAME, LUF_NAME) %>% 
              unique()) %>% 
  rename(PCODE = project,
         SS = location,
         SSYR = gisid,
         PKEY = visitid,
         YEAR = year,
         DATI = date,
         MAXDUR = duration,
         MAXDIS = distance,
         X = lon,
         Y= lat,
         TAGMETHOD = tagMethod,
         EQUIP = equipment,
         pAspen = Populus_tremuloides_brtpred_nofp,
         PET = Eref)  %>% 
  mutate(DATE = as.Date(str_sub(as.character(DATI), 1, 10)),
         ROAD = ifelse(PCODE %in% c("BAM-BBS-1987-2006", "BAM-BBS-2007-2022"), 1, NA),
         CMETHOD = case_when(TAGMETHOD=="PC" ~ "HS",
                             source=="riverforks" ~ "RF",
                             !is.na(TAGMETHOD) ~ "SM")) %>% 
  dplyr::select(colnames(e1$dd), TAGMETHOD, EQUIP) %>% 
  arrange(PKEY)
rownames(dd) <- dd$PKEY

#compare to previous - should be 2 new columns (tagmeth, equip)
compare_sets(colnames(e1$dd), colnames(dd))

#check rows against yy
compare_sets(rownames(yy), rownames(dd))
  
#6. Point level covariates (vs0)----
vs0 <- visit %>% 
  dplyr::select(visitid, gisid) %>% 
  left_join(d_long) %>% 
  dplyr::select(colnames(e1$vs0), visitid) %>% 
  unique() %>% 
  arrange(visitid)
rownames(vs0) <- vs0$visitid
vs0 <- vs0[,colnames(vs0)!="visitid"]

#compare to previous
compare_sets(colnames(e1$vs0), colnames(vs0))

#check rows against dd
compare_sets(rownames(dd), rownames(vs0))

#7. Veg & soil objects-----

#7a. Veg current
vc1 <- as.data.frame(as.matrix(d_wide1$veg_current)) %>% 
  mutate()

vc1 <- visit %>% 
  dplyr::select(visitid, gisid) %>% 
  left_join(dsite) %>% 
  left_join(as.data.frame(as.matrix(d_wide1$veg_current)) %>% 
              mutate(UID = rownames(d_wide1$veg_current))) %>% 
  dplyr::select(colnames(e1$vc1), visitid) %>% 
  unique() %>% 
  arrange(visitid)
rownames(vc1) <- vc1$visitid
vc1 <- as(as.matrix(vc1[,colnames(vc1)!="visitid"]), "dgCMatrix")

vc2 <- as.data.frame(as.matrix(d_wide2$veg_current)) %>% 
  mutate()

vc2 <- visit %>% 
  dplyr::select(visitid, gisid) %>% 
  left_join(dsite) %>% 
  left_join(as.data.frame(as.matrix(d_wide2$veg_current)) %>% 
              mutate(UID = rownames(d_wide2$veg_current))) %>% 
  dplyr::select(colnames(e1$vc2), visitid) %>% 
  unique() %>% 
  arrange(visitid)
rownames(vc2) <- vc2$visitid
vc2 <- as(as.matrix(vc2[,colnames(vc2)!="visitid"]), "dgCMatrix")

#compare to previous
compare_sets(colnames(e1$vc1), colnames(vc1))
compare_sets(colnames(e1$vc2), colnames(vc2))

#check rows against dd
compare_sets(rownames(dd), rownames(vc1))
compare_sets(rownames(dd), rownames(vc2))

#7b. Veg reference
vr1 <- as.data.frame(as.matrix(d_wide1$veg_reference)) %>% 
  mutate()

vr1 <- visit %>% 
  dplyr::select(visitid, gisid) %>% 
  left_join(dsite) %>% 
  left_join(as.data.frame(as.matrix(d_wide1$veg_reference)) %>% 
              mutate(UID = rownames(d_wide1$veg_reference))) %>% 
  dplyr::select(colnames(e1$vr1), visitid) %>% 
  unique() %>% 
  arrange(visitid)
rownames(vr1) <- vr1$visitid
vr1 <- as(as.matrix(vr1[,colnames(vr1)!="visitid"]), "dgCMatrix")

vr2 <- as.data.frame(as.matrix(d_wide2$veg_reference)) %>% 
  mutate()

vr2 <- visit %>% 
  dplyr::select(visitid, gisid) %>% 
  left_join(dsite) %>% 
  left_join(as.data.frame(as.matrix(d_wide2$veg_reference)) %>% 
              mutate(UID = rownames(d_wide2$veg_reference))) %>% 
  dplyr::select(colnames(e1$vr2), visitid) %>% 
  unique() %>% 
  arrange(visitid)
rownames(vr2) <- vr2$visitid
vr2 <- as(as.matrix(vr2[,colnames(vr2)!="visitid"]), "dgCMatrix")

#compare to previous
compare_sets(colnames(e1$vr1), colnames(vr1))
compare_sets(colnames(e1$vr2), colnames(vr2))

#check rows against dd
compare_sets(rownames(dd), rownames(vr1))
compare_sets(rownames(dd), rownames(vr2))

#7c. Soil current
sc1 <- as.data.frame(as.matrix(d_wide1$soil_current)) %>% 
  mutate()

sc1 <- visit %>% 
  dplyr::select(visitid, gisid) %>% 
  left_join(dsite) %>% 
  left_join(as.data.frame(as.matrix(d_wide1$soil_current)) %>% 
              mutate(UID = rownames(d_wide1$soil_current))) %>% 
  dplyr::select(colnames(e1$sc1), visitid) %>% 
  unique() %>% 
  arrange(visitid)
rownames(sc1) <- sc1$visitid
sc1 <- as(as.matrix(sc1[,colnames(sc1)!="visitid"]), "dgCMatrix")

sc2 <- as.data.frame(as.matrix(d_wide2$soil_current)) %>% 
  mutate()

sc2 <- visit %>% 
  dplyr::select(visitid, gisid) %>% 
  left_join(dsite) %>% 
  left_join(as.data.frame(as.matrix(d_wide2$soil_current)) %>% 
              mutate(UID = rownames(d_wide2$soil_current))) %>% 
  dplyr::select(colnames(e1$sc2), visitid) %>% 
  unique() %>% 
  arrange(visitid)
rownames(sc2) <- sc2$visitid
sc2 <- as(as.matrix(sc2[,colnames(sc2)!="visitid"]), "dgCMatrix")

#compare to previous
compare_sets(colnames(e1$sc1), colnames(sc1))
compare_sets(colnames(e1$sc2), colnames(sc2))

#check rows against dd
compare_sets(rownames(dd), rownames(sc1))
compare_sets(rownames(dd), rownames(sc2))

#7d. Soil reference
sr1 <- as.data.frame(as.matrix(d_wide1$soil_reference)) %>% 
  mutate()

sr1 <- visit %>% 
  dplyr::select(visitid, gisid) %>% 
  left_join(dsite) %>% 
  left_join(as.data.frame(as.matrix(d_wide1$soil_reference)) %>% 
              mutate(UID = rownames(d_wide1$soil_reference))) %>% 
  dplyr::select(colnames(e1$sr1), visitid) %>% 
  unique() %>% 
  arrange(visitid)
rownames(sr1) <- sr1$visitid
sr1 <- as(as.matrix(sr1[,colnames(sr1)!="visitid"]), "dgCMatrix")

sr2 <- as.data.frame(as.matrix(d_wide2$soil_reference)) %>% 
  mutate()

sr2 <- visit %>% 
  dplyr::select(visitid, gisid) %>% 
  left_join(dsite) %>% 
  left_join(as.data.frame(as.matrix(d_wide2$soil_reference)) %>% 
              mutate(UID = rownames(d_wide2$soil_reference))) %>% 
  dplyr::select(colnames(e1$sr2), visitid) %>% 
  unique() %>% 
  arrange(visitid)
rownames(sr2) <- sr2$visitid
sr2 <- as(as.matrix(sr2[,colnames(sr2)!="visitid"]), "dgCMatrix")

#compare to previous
compare_sets(colnames(e1$sr1), colnames(sr1))
compare_sets(colnames(e1$sr2), colnames(sr2))

#check rows against dd
compare_sets(rownames(dd), rownames(sr1))
compare_sets(rownames(dd), rownames(sr2))

#8. Package and save----
save(tax, yy, dd, vs0, vc1, vr1, sc1, sr1, vc2, vr2, sc2, sr2,
     file=file.path(root, "Data", "2Wrangled.Rdata"))
