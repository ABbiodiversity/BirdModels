# ---
# title: ABMI models - compare versions
# author: Elly Knight
# created: Sept 25, 2024
# ---

#NOTES################################

#PREAMBLE############################

#1. Load packages----
library(tidyverse) #basic data wrangling
library(Matrix) #sparse matrices

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

#6. Output from visualization script----
out <- read.csv(file.path(root, "Results", "ModelComparison.csv"))

#MODEL EVALUATION####

#1. Wrangle----
auc <- old$birds$species |>
  dplyr::select(Comments, AUCNorth, AUCSouth) |>
  mutate(version = "old") |>
  rbind(new$birds$species |>
          dplyr::select(Comments, AUCNorth, AUCSouth) |>
          mutate(version = "new")) |> 
  dplyr::filter(!(is.na(AUCNorth) & is.na(AUCSouth))) |> 
  mutate(model = case_when(is.na(AUCNorth) ~ "south",
                           is.na(AUCSouth) ~ "north",
                           !is.na(AUCNorth) & !is.na(AUCSouth) ~ "both")) |> 
  pivot_longer(cols = AUCNorth:AUCSouth, names_to="Region", values_to="AUC", names_prefix = "AUC") |>
  pivot_wider(names_from=version, values_from=AUC) |> 
  dplyr::filter(!is.na(old) & !is.na(new))

#2. Plot against each other----
ggplot(auc) +
  geom_text(aes(x=old, y=new, label=Comments, colour = model)) +
  geom_abline(aes(intercept=0, slope=1), linetype="dashed") +
  facet_wrap(~Region)

#3. Plot totals----
ggplot(pivot_longer(auc, old:new, names_to="version", values_to="AUC")) +
  geom_violin(aes(x=version, y=AUC))

ggplot(pivot_longer(auc, old:new, names_to="version", values_to="AUC")) +
  geom_violin(aes(x=version, y=AUC, fill=model))

#POPULATION SIZE###########

#1. Define species with unreasonable population sizes in old version ----
remove <- c("MAGO", "LASP", "SAVS", "CCLO", "BANS", "LBCU", "BRBL")

#2. Correlation----
ggplot(out |> dplyr::filter(!species %in% remove)) +
  geom_text(aes(x=pop.old, y=pop.new, label=species, colour=cor_oldnew)) +
  geom_abline(aes(intercept=0, slope=1), linetype="dashed") +
  scale_colour_viridis_c()

#3. Correlation log scale----
ggplot(out |> dplyr::filter(!species %in% remove)) +
  geom_text(aes(x=log(pop.old), y=log(pop.new), label=species, colour=cor_oldnew)) +
  geom_abline(aes(intercept=0, slope=1), linetype="dashed") +
  scale_colour_viridis_c()

#4. Histograms----
ggplot(out |> dplyr::filter(!species %in% remove)) +
  geom_histogram(aes(x=pop.old))

ggplot(out |> dplyr::filter(!species %in% remove)) +
  geom_histogram(aes(x=pop.new))

#5. Boxplot----
ggplot(out |> pivot_longer(pop.old:pop.new, names_to="version", values_to="males") |>
         dplyr::filter(!species %in% remove)) +
  geom_boxplot(aes(x=version, y=log(males)))

#SPATIAL CORRELATION###########

#1. Put together with AUC----
out_auc <- auc |> 
  rename(species = Comments) |> 
  inner_join(out) |> 
  mutate(auc_diff = new - old)

#2. Old vs new vs AUC ----
ggplot(out_auc) + 
  geom_text(aes(x=cor_oldnew, y=auc_diff, label=species))

#3. EBird correlation ----
ggplot(out_auc) +
  geom_text(aes(x=cor_oldebd, y=cor_newebd, label=species)) +
  geom_abline(aes(intercept=0, slope=1), linetype="dashed")

ggplot(out_auc |> pivot_longer(cor_oldebd:cor_newebd, names_to="version", values_to="corr")) +
  geom_boxplot(aes(x=version, y=corr))
