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

options(scipen=99999)

#COMPARE VERSIONS##########

#1. AUC----
auc22 <- COEFS$birds$species %>% 
  dplyr::select(Comments, SpeciesID, AUCNorth, AUCSouth, Occurrences) %>% 
  rename(aucn22 = AUCNorth, aucs22 = AUCSouth, n22 = Occurrences)
auc23 <- birds$species %>% 
  dplyr::select(Comments, SpeciesID, AUCNorth, AUCSouth, Occurrences) %>% 
  rename(aucn23 = AUCNorth, aucs23 = AUCSouth, n23 = Occurrences)
auc <- full_join(auc22, auc23) %>% 
  mutate(aucn = aucn23 - aucn22,
         aucs = aucs23 - aucs22)

auc.n <- ggplot(auc) +
  geom_abline(intercept = 0, slope = 1)+
  geom_text(aes(x=aucn22, y=aucn23, label=Comments)) +
  xlab("North AUC 2022") +
  ylab("North AUC 2023") +
  theme_bw()
auc.n

auc.s <- ggplot(auc) +
  geom_abline(intercept = 0, slope = 1)+
  geom_text(aes(x=aucs22, y=aucs23, label=Comments)) +
  xlab("South AUC 2022") +
  ylab("South AUC 2023") +
  theme_bw()
auc.s

ggsave(gridExtra::grid.arrange(auc.n, auc.s, ncol=2, nrow=1), filename = file.path(root, "Figures", "AUC.jpeg"), width = 10, height = 4)

#2. Coefficients----

#South####
coef22.s <- rowMeans(COEFS$birds$south$joint, dims=2) %>%
  data.frame() %>%
  rownames_to_column("species") %>%
  rename(xMAP.xPET = xMAP:xPET, xAHM.xMAT = xAHM:xMAT, xX.xY=xX:xY, xFFP.xMAP=xFFP:xMAP) %>%
  pivot_longer(cols = Loamy:xX.xY, names_to="covariate", values_to="coef22")

coef23.s <- rowMeans(birds$south$joint, dims=2) %>%
  data.frame() %>%
  rownames_to_column("species") %>%
  rename(xMAP.xPET = xMAP:xPET, xAHM.xMAT = xAHM:xMAT, xX.xY=xX:xY, xFFP.xMAP=xFFP:xMAP) %>%
  pivot_longer(cols = Loamy:xX.xY, names_to="covariate", values_to="coef23")

coef.s <- full_join(coef22.s, coef23.s) %>% 
  dplyr::filter(!is.na(coef22), !is.na(coef23)) %>% 
  mutate(same = ifelse(coef22==coef23, 1, 0),
         diff = abs(coef22 - coef23))

table(coef.s$covariate, coef.s$same)

#Plot the coefficients against each other-----
ggplot(coef.s %>% 
         dplyr::filter(!covariate %in% c("EnSeismic", "EnSoftLin", "TrSoftLin", "Wellsites"),
                       coef22 > -10000,
                       coef23 > -10000,
                       !species %in% c("BrewersSparrow", "LarkSparrow", "LarkBunting", "Bobolink"))) +
  geom_point(aes(x=coef22, y=coef23, colour=species)) +
  geom_abline(aes(intercept=0, slope = 1)) +
  geom_hline(aes(yintercept=0)) +
  geom_vline(aes(xintercept=0))

#Correlation by species----
spp <- unique(coef.s$species)
resid.s <- data.frame()
cor.s <- data.frame()
for(i in 1:length(spp)){
  
  coef.i <- coef.s %>% 
    dplyr::filter(species==spp[i],
                  !covariate %in% c("EnSeismic", "EnSoftLin", "TrSoftLin", "Wellsites"))
  
  lm.i <- lm(coef22 ~ coef23, coef.i)
  
  coef.i$resid <- residuals(lm.i)
  
  resid.s <- rbind(resid.s, coef.i)
  
  cor.s <- rbind(cor.s, 
                 data.frame(species = spp[i],
                      cor = cor(coef.i$coef22, coef.i$coef23)))
  
}

#Plot the residuals
ggplot(resid.s) +
  geom_boxplot(aes(x=covariate, y=resid)) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5))

ggplot(resid.s) +
  geom_boxplot(aes(x=species, y=resid)) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5))

#Correlation by covariate----
cov.s <- unique(coef.s$covariate)
cor.cov.s <- data.frame()
for(i in 1:length(cov.s)){
  
  coef.i <- coef.s %>% 
    dplyr::filter(covariate==cov.s[i])
  
  cor.cov.s <- rbind(cor.cov.s,
                     data.frame(covariate=cov.s[i],
                                cor = cor(coef.i$coef22, coef.i$coef23)))
}

cor.cov.s %>% 
  arrange(desc(cor))

#Correlation by species and covariate-----
spp <- as.character(unlist(unique(dimnames(birds$south$joint))[[1]]))
cov <- as.character(unlist(unique(dimnames(birds$south$joint))[[2]]))
loop <- expand.grid(spp=spp, cov=cov, stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE) %>% 
  dplyr::filter(spp %in% as.character(unlist(unique(dimnames(COEFS$birds$south$joint))[[1]])))

t.s <- data.frame()
for(i in 1:nrow(loop)){
  
  coef22.i <- COEFS$birds$south$joint[loop$spp[i],loop$cov[i],]
  coef23.i <- birds$south$joint[loop$spp[i],loop$cov[i],1:100]
  
  t.i <- try(t.test(coef22.i, coef23.i), silent = TRUE)
  
  if(class(t.i)!="try-error"){
    t.s <- rbind(t.s, data.frame(spp = loop$spp[i],
                                 cov = loop$cov[i],
                                 p = t.i$p.value,
                                 t = t.i$statistic,
                                 mean22 = mean(coef22.i),
                                 mean23 = mean(coef23.i)) %>% 
                   mutate(meand = mean22-mean23))
  }
  
}

t.s %>% 
  dplyr::filter(p < 0.05) %>% 
  group_by(cov) %>% 
  summarize(n=n()) %>% 
  ungroup() %>% 
  arrange(-n) %>% 
  View()

ggplot(t.s %>% 
         dplyr::filter(p < 0.05)) +
  geom_boxplot(aes(x=cov, y=p)) +
  theme(axis.text.x = element_text(angle=90))

ggplot(t.s) +
  geom_boxplot(aes(x=cov, y=log(meand))) +
  theme(axis.text.x = element_text(angle=90)) +
  ylab("log of mean 2023 - mean 2022 coefficient")

ggsave(filename = file.path(root, "Figures", "Coefficient_diff_south.jpeg"), width = 8, height = 5)

ggplot(t.s) +
  geom_point(aes(x=mean22, y=mean23)) +
  geom_abline(intercept=0, slope=1) +
  facet_wrap(~cov, scales="free")

ggsave(filename = file.path(root, "Figures", "Coefficients_south.jpeg"), width = 20, height = 15)

#North####
coef22.n <- rowMeans(COEFS$birds$north$joint, dims=2) %>%
  data.frame() %>%
  rownames_to_column("species") %>%
  rename(xMAP.xPET = xMAP:xPET, xAHM.xMAT = xAHM:xMAT, xX.xY=xX:xY, xFFP.xMAP=xFFP:xMAP) %>%
  pivot_longer(cols = WhiteSpruceR:xX.xY, names_to="covariate", values_to="coef22")

coef23.n <- rowMeans(birds$north$joint, dims=2) %>%
  data.frame() %>%
  rownames_to_column("species") %>%
  rename(xMAP.xPET = xMAP:xPET, xAHM.xMAT = xAHM:xMAT, xX.xY=xX:xY, xFFP.xMAP=xFFP:xMAP) %>%
  pivot_longer(cols = WhiteSpruceR:xX.xY, names_to="covariate", values_to="coef23")

coef.n <- full_join(coef22.n, coef23.n) %>% 
  dplyr::filter(!is.na(coef22), !is.na(coef23)) %>% 
  mutate(same = ifelse(coef22==coef23, 1, 0))

ggplot(coef.n %>% 
         dplyr::filter(!covariate %in% c("EnSeismic", "EnSoftLin", "TrSoftLin", "Wellsites"),
                       coef22 > -10000,
                       coef23 > -10000)) +
  geom_point(aes(x=coef22, y=coef23, colour=species)) +
  geom_abline(aes(intercept=0, slope = 1)) +
  geom_hline(aes(yintercept=0)) +
  geom_vline(aes(xintercept=0)) +
  facet_wrap(~covariate, scales="free") +
  theme(legend.position = "none")

table(coef.n$covariate, coef.n$same)
#Just the linear ones, all good

spp <- unique(coef.n$species)
resid.n <- data.frame()
cor.n <- data.frame()
for(i in 1:length(spp)){
  
  coef.i <- coef.n %>% 
    dplyr::filter(species==spp[i])
  
  lm.i <- lm(coef22 ~ coef23, coef.i)
  
  coef.i$resid <- residuals(lm.i)
  
  resid.n <- rbind(resid.n, coef.i)
  
  cor.n <- rbind(cor.n, 
                 data.frame(species = spp[i],
                            cor = cor(coef.i$coef22, coef.i$coef23)))
  
}

resid.n.nl <- dplyr::filter(resid.n, !covariate %in% c("EnSeismic", "EnSoftLin", "TrSoftLin", "Wellsites"))

ggplot(dplyr::filter(resid.n, !covariate %in% c("EnSeismic", "EnSoftLin", "TrSoftLin", "Wellsites"))) +
  geom_boxplot(aes(x=covariate, y=resid)) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5))

ggplot(resid.n) +
  geom_boxplot(aes(x=covariate, y=resid)) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5))

#Correlation by covariate
cov.n <- unique(coef.n$covariate)
cor.cov.n <- data.frame()
for(i in 1:length(cov.n)){
  coef.i <- coef.n %>% 
    dplyr::filter(covariate==cov.n[i])
  cor.cov.n <- rbind(cor.cov.n,
                     data.frame(covariate=cov.n[i],
                                cor = cor(coef.i$coef22, coef.i$coef23)))
}

#Correlation by species and covariate-----
spp <- as.character(unlist(unique(dimnames(birds$north$joint))[[1]]))
cov <- as.character(unlist(unique(dimnames(birds$north$joint))[[2]]))
loop <- expand.grid(spp=spp, cov=cov, stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE) %>% 
  dplyr::filter(spp %in% as.character(unlist(unique(dimnames(COEFS$birds$north$joint))[[1]])))

t.n <- data.frame()
for(i in 1:nrow(loop)){
  
  coef22.i <- COEFS$birds$north$joint[loop$spp[i],loop$cov[i],]
  coef23.i <- birds$north$joint[loop$spp[i],loop$cov[i],1:100]
  
  t.i <- try(t.test(coef22.i, coef23.i), silent = TRUE)
  
  if(class(t.i)!="try-error"){
    t.n <- rbind(t.n, data.frame(spp = loop$spp[i],
                                 cov = loop$cov[i],
                                 p = t.i$p.value,
                                 t = t.i$statistic,
                                 mean22 = mean(coef22.i),
                                 mean23 = mean(coef23.i)) %>% 
                   mutate(meand = mean22-mean23))
  }

}

t.n %>% 
  dplyr::filter(p < 0.05) %>% 
  group_by(cov) %>% 
  summarize(n=n()) %>% 
  ungroup() %>% 
  arrange(-n) %>% 
  View()

ggplot(t.n %>% 
         dplyr::filter(p < 0.05)) +
  geom_boxplot(aes(x=cov, y=p)) +
  theme(axis.text.x = element_text(angle=90))

ggplot(t.n) +
  geom_boxplot(aes(x=cov, y=log(meand))) +
  theme(axis.text.x = element_text(angle=90))

ggsave(filename = file.path(root, "Figures", "Coefficient_diff_north.jpeg"), width = 8, height = 5)

ggplot(t.n) +
  geom_point(aes(x=mean22, y=mean23)) +
  geom_abline(intercept=0, slope=1) +
  facet_wrap(~cov, scales="free") +
  theme_bw() +
  xlab("Mean coef 2022") +
  ylab("Mean coef 2023")

ggsave(filename = file.path(root, "Figures", "Coefficients_north.jpeg"), width = 20, height = 15)

#3. Correlation vs delta AUC----
cor.auc <- cor.n %>% 
  rename(corn = cor) %>% 
  full_join(cor.s %>% 
              rename(cors=cor)) %>% 
  full_join(auc %>% 
              rename(species = SpeciesID)) %>% 
  mutate(n=n23-n22)

ggplot(cor.auc) +
  geom_point(aes(x=corn, y=aucn))

ggplot(cor.auc) +
  geom_point(aes(x=cors, y=aucs))

ggplot(cor.auc) +
  geom_point(aes(x=n, y=corn))

ggplot(cor.auc) +
  geom_point(aes(x=n, y=cors))

ggplot(cor.auc) +
  geom_point(aes(x=n, y=aucn))

ggplot(cor.auc) +
  geom_point(aes(x=n, y=aucs))
