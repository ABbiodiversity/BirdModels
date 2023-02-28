# ---
# title: ABMI models - evaluate results
# author: Elly Knight
# created: Feb 26 2023
# ---

#NOTES################################


#PREAMBLE############################

#1. Load packages----
library(tidyverse) #basic data wrangling
library(allinone) #previous results

#2. Set root path for data on google drive----
root <- "G:/My Drive/ABMI/Projects/BirdModels/"

#3. Load 2023 coefficients----
load(file.path(root, "results", "Birds2023.RData"))

#4. Load 2022 coefficients----
load(file.path(getOption("allinone")$dir, "COEFS.RData"))

#COMPARE VERSIONS##########

#1. AUC----
auc22 <- COEFS$birds$species %>% 
  dplyr::select(Comments, AUCNorth, AUCSouth, Occurrences) %>% 
  rename(aucn22 = AUCNorth, aucs22 = AUCSouth, n22 = Occurrences)
auc23 <- birds$species %>% 
  dplyr::select(Comments, AUCNorth, AUCSouth, Occurrences) %>% 
  rename(aucn23 = AUCNorth, aucs23 = AUCSouth, n23 = Occurrences)
auc <- full_join(auc22, auc23)

auc.n <- ggplot(auc) +
  geom_abline(intercept = 0, slope = 1)+
  geom_point(aes(x=aucn22, y=aucn23)) +
  xlab("North AUC 2022") +
  ylab("North AUC 2023") +
  theme_bw()
auc.n

auc.s <- ggplot(auc) +
  geom_abline(intercept = 0, slope = 1)+
  geom_point(aes(x=aucs22, y=aucs23)) +
  xlab("South AUC 2022") +
  ylab("South AUC 2023") +
  theme_bw()
auc.s

ggsave(gridExtra::grid.arrange(auc.n, auc.s, ncol=2, nrow=1), filename = file.path(root, "Figures", "AUC.jpeg"), width = 10, height = 4)

ggplot(auc) +
  geom_abline(intercept = 0, slope = 1)+
  geom_point(aes(x=n22, y=n23))

#2. Coefficients----
plot(exp(COEFS$birds$north$marginal[,1,1]), exp(birds$north$marginal[,1,1]))
