
library(ggplot2)  
library(phytotools)
library(plyr)
library(tidyr)
library(dplyr)
library(reshape2)
library(rstatix)
library(cowplot)
library(ggpubr)
rm(list=ls())

#calculating RLC parameters using phytotools (pratt 1980)
setwd("/Volumes/CBP_Students/Jessica Bellworthy/Jessica Bellworthy/NSF-BSF Genomic/R/Github")
d <- read.csv('rETR_insitu_2022.csv')
View(d)


rlc.data <- d[c(2,3,1)] # Subset data and ensure that column order is: PAR - ETR - ID/COLONY
rlc.data<- setNames(rlc.data, c("par","etr", "id")) # Set column names
rlc.data$id <- as.factor(rlc.data$id)
head(rlc.data)

rlc.data$etr <- na_if(rlc.data$etr, 0) # Change any zeros in etr to NA
View(rlc.data)

ncurves <- length(unique(rlc.data$id)) # number of unique ids in the data 
ids <- unique(rlc.data$id) # store the unique ids 


rlc.parameters <- data.frame(
  id = ids, 
  alpha = 0, 
  beta = 0, 
  ETRmax = 0, 
  Ek = 0, 
  ps = 0
)

# Run fitPGH.R script

for (i in 1:ncurves){
  
  temp.id = ids[i] # extract the id of the curve to be fitted
  
  print(paste("Now fitting curve ", as.character(temp.id))) # to keep track what's happening if the data has many curves
  
  temp.rlc.data <- rlc.data[rlc.data$id==temp.id,] # extract the the data of a single curve into a temporary variable
  PAR = temp.rlc.data$par 
  ETR = temp.rlc.data$etr
  
  fit = fitPGH(PAR, ETR, fitmethod = "Port") # for more options and explanation see package phytotools manual
  
  # store the fitted RLC values into temporary variables
  alpha.rlc = fit$alpha[1]
  beta.rlc = fit$beta[1]
  ps.rlc = fit$ps[1]
  
  # store the parameters
  rlc.parameters$id[i] <- temp.id
  rlc.parameters$alpha[i] <- alpha.rlc
  rlc.parameters$beta[i] <- beta.rlc
  rlc.parameters$ps[i] <- ps.rlc
  
  # calculate ETRmax and Ek for the PGH model (see e.g.Ralph & Gademann 2005 Aquatic Botany 82 (3): 222 - 237). 
  # Note that the equation depends on the model fitted, the code below applies only to the PGH model! 
  # Model equations are documented in the phytotools package code examples (and in the original papers): https://cran.r-project.org/web/packages/phytotools/phytotools.pdf
  
  ETRmax = ps.rlc*(alpha.rlc/(alpha.rlc + beta.rlc))*(beta.rlc/(alpha.rlc+beta.rlc))^(beta.rlc/alpha.rlc)
  Ek = ETRmax/alpha.rlc 
  
  # store the variables
  rlc.parameters$ETRmax[i] <- ETRmax
  rlc.parameters$Ek[i] <- Ek
  
  #plotting the curve and fitted model into a tiff file. By default the file name is the id of the curve. 
  tiff(file=paste0(temp.id, ".tiff"), compression="lzw")
  
  # plot the data, 
  plot(x=PAR, y=ETR, main=temp.id) 
  
  # plot the model fit
  with(fit, {
    P <- ps.rlc*(1-exp(-1*alpha.rlc*PAR/ps.rlc))*exp(-1*beta.rlc*PAR/ps.rlc) # the PGH model equation
    lines(PAR,P)
  }
  ) # end of with
  dev.off() #close the plotting divide. if this is not done, the next run of the loop will override the plot. 
  
}

# now the data frame rlc.parameters contains the fitted values for each curve. Tiff plots should be in current working directory. 
rlc.parameters
View(rlc.parameters)

warnings()

write.csv(file="rlc.parameters.insitu_60days_StyloSpat2022.csv",rlc.parameters)


#### check statistical differences
library(lme4)
library(nlme)
library(lmerTest)
library(MASS)
library(car)
library(predictmeans)
library(ggpubr)

### manually copy the FVFM to the rlc parameter data file and and the metadata e.g. treatment

in.data <- read_csv('rlc.parameters_metadata_TranslocationStyloSpat_2021.2022.csv')

head(in.data)
str(in.data)
in.data$chamber <- as.factor(in.data$chamber)
in.data$treatment <- as.factor(in.data$treatment)
in.data$year <- as.factor(in.data$year)
in.data$age <- as.factor(in.data$age)


## summary of results 
Sum_all <- in.data %>% 
  group_by(treatment) %>% 
  summarise_each(funs(mean(., na.rm=TRUE), n = sum(!is.na(.)), max(., na.rm = TRUE), min(., na.rm = TRUE),
                      se = sd(., na.rm=TRUE)/sqrt(sum(!is.na(.)))), alpha:fvfm)

View(Sum_all)
write.csv(file="rlc.summarytable.InSitu_60days.csv", Sum_all)


### Assumption testing for in situ data by day
#Subset data by treatment group
SD_8 = in.data %>% subset(age == '8' & treatment == "SD")
SS_8 = in.data %>% subset(age == '8' & treatment == "SS")
SD_60 = in.data %>% subset(age == '60' & treatment == "SD")
SS_60 = in.data %>% subset(age == '60' & treatment == "SS")

in.data.Ek = subset(in.data, Ek < 90) # remove two high anomalies from Ek data
# subset new Ek data table
SD_8eK = in.data.Ek %>% subset(age == '8' & treatment == "SD")
SS_8eK = in.data.Ek %>% subset(age == '8' & treatment == "SS")
SD_60eK = in.data.Ek %>% subset(age == '60' & treatment == "SD")
SS_60eK = in.data.Ek %>% subset(age == '60' & treatment == "SS")

# 8 days - assumption OK, unless written otherwise
shapiro.test(SD_8$ETRmax) 
shapiro.test(SS_8$ETRmax) 
leveneTest(ETRmax~treatment,d=(subset(in.data, age == '8')))

shapiro.test(SD_8$alpha) 
shapiro.test(SS_8$alpha) 
leveneTest(alpha~treatment,d=(subset(in.data, age == '8')))

shapiro.test(SD_8eK$Ek) 
shapiro.test(SS_8eK$Ek) 
leveneTest(Ek~treatment,d=(subset(in.data.Ek, age == '8')))

shapiro.test(SD_8$fvfm) 
shapiro.test(SS_8$fvfm) 
leveneTest(fvfm~treatment,d=(subset(in.data, age == '8')))


# 60 days - assumption OK, unless written otherwise
shapiro.test(SD_60$ETRmax) # not normal, both only 2021 and both years combined, log/ sqrt doesnt help, sq good   === SQ
shapiro.test(SS_60$ETRmax) # all years data are not normal, log transform ok, sqrt transform ok, sq not good    ==== SQRT
leveneTest(ETRmax.log~treatment,d=(subset(in.data, age == '60')))

shapiro.test(SD_60$alpha) 
shapiro.test(SS_60$alpha) # not normal, both only 2021 and both years combined, log/ srt/sq donest help
leveneTest(alpha~treatment,d=(subset(in.data, age == '60')))# all years data are not homogenous, log/ sqrt doesnt help, sq good   ==== SQ

shapiro.test(SD_60$Ek) # all years data are not normal, log doesnt help, ====== SQ
shapiro.test(SS_60$Ek) # log doesnt help, sq doesnt help
leveneTest(Ek~treatment,d=(subset(in.data, age == '60')))

shapiro.test(SD_60$fvfm) #  not normal, both only 2021 and both years combined, log/sqrt/sq doesnt help
shapiro.test(SS_60$fvfm) 
leveneTest(fvfm~treatment,d=(subset(in.data, age == '60')))

# # # Multiple breaches of assumptions in 60 days data - will need to use non parametric statistics

# STATS for graphics
stat.etr <- in.data%>%
  group_by(age) %>%
  wilcox_test(ETRmax ~ treatment) %>%
  add_significance() %>%
  p_round(digits = 3) %>%
  add_xy_position(x = "treatment")
stat.etr

stat.alpha <- in.data %>%
  group_by(age) %>%
  wilcox_test(alpha ~ treatment) %>%
  add_significance() %>%
  p_round(digits = 3) %>%
  add_xy_position(x = "treatment")
stat.alpha

stat.eK <- in.data.Ek %>%
  group_by(age) %>%
  wilcox_test(Ek ~ treatment) %>%
  add_significance() %>%
  p_round(digits = 3) %>%
  add_xy_position(x = "treatment")
stat.eK

stat.fvfm <- in.data %>%
  group_by(age) %>%
  wilcox_test(fvfm ~ treatment) %>%
  add_significance() %>%
  p_round(digits = 3) %>%
  add_xy_position(x = "treatment")
stat.fvfm



## Plot all in situ data from 2021 and 2022 together

mytheme = theme_classic()+
  theme(axis.title.x = element_blank(), axis.text = element_text(colour = "black", size = 8), axis.title = element_text(size = 10))

### GRAPHS ###

in.etr = ggplot(in.data, aes(y = ETRmax, x = treatment)) +
  geom_boxplot(outlier.shape = NA, aes(fill= treatment, alpha = 0.8), fatten = 0.5, alpha = 0.7, lwd = 0.2)+
  # geom_jitter(aes(shape= year), position = position_jitter(width = .25), fill = "black", size = 0.4)+
  geom_jitter(position = position_jitter(width = .25), fill = "black", size = 0.25)+
  facet_wrap(facets = "age", nrow = 1)+
  mytheme+
  scale_fill_manual(values = c("#f0f921", "#FCA510"))+
  scale_shape_manual(values = c(19, 17))+
  scale_y_continuous(limits = c(0, 50, 10))+
  theme(axis.title.x = element_blank(), axis.text = element_text(colour = "black"), strip.background = element_rect(size = 0.5))+
  stat_pvalue_manual(stat.etr, label = "{p.signif} {p}", tip.length = 0.005, hide.ns = FALSE, size = 2.5)+
  guides(fill = "none", alpha = "none")
in.etr

in.alpha= ggplot(in.data, aes(y = alpha, x = treatment)) +
  geom_boxplot(outlier.shape = NA, aes(fill= treatment, alpha = 0.8), fatten = 0.5, alpha = 0.7, lwd = 0.2)+
  #  geom_jitter(aes(shape= year), position = position_jitter(width = .25), size = 0.4,fill = "black")+
  geom_jitter(position = position_jitter(width = .25), fill = "black", size = 0.25)+
   facet_wrap(facets = "age", nrow = 1)+
  mytheme+
  scale_fill_manual(values = c("#f0f921", "#FCA510"))+
  scale_shape_manual(values = c(17, 19))+
  scale_y_continuous(limits = c(0.2, 0.9, 0.2))+
  theme(axis.title.x = element_blank(), axis.text = element_text(colour = "black"), strip.background = element_rect(size = 0.5))+
stat_pvalue_manual(stat.alpha, label = "{p.signif} {p}", tip.length = 0.005, hide.ns = FALSE, size = 2.5)+
  guides(fill = "none", alpha = "none")
in.alpha


in.eK= ggplot(in.data.Ek, aes(y = Ek, x = treatment)) +
 geom_boxplot(outlier.shape = NA, aes(fill= treatment, alpha = 0.8), fatten = 0.5, alpha = 0.7, lwd = 0.2)+
  #  geom_jitter(aes(shape= year), position = position_jitter(width = .25), size = 0.4, fill = "black")+
  geom_jitter(position = position_jitter(width = .25), fill = "black", size = 0.25)+
  facet_wrap(facets = "age", nrow = 1)+
mytheme+
  scale_fill_manual(values = c("#f0f921", "#FCA510"))+
  scale_shape_manual(values = c(17, 19))+
  scale_y_continuous(limits = c(0, 100, 20))+
  theme(axis.title.x = element_blank(), axis.text = element_text(colour = "black"), strip.background = element_rect(size = 0.5))+
 stat_pvalue_manual(stat.eK, label = "{p.signif} {p}", tip.length = 0.005, hide.ns = FALSE, size = 2.5)+
  guides(fill = "none", alpha = "none")
in.eK

in.fvfm= ggplot(in.data, aes(y = fvfm, x = treatment)) +
  geom_boxplot(outlier.shape = NA, aes(fill= treatment, alpha = 0.8), fatten = 0.5, alpha = 0.7, lwd = 0.2)+
  # geom_jitter(aes(shape= year), position = position_jitter(width = .25), size = 0.7, fill = "black")+
  geom_jitter(position = position_jitter(width = .25), fill = "black", size = 0.25)+
  facet_wrap(facets = "age", nrow = 1)+
mytheme+
  labs(y= ~F[V] ~F[M])+
  scale_fill_manual(values = c("#f0f921", "#FCA510"))+
  scale_shape_manual(values = c(17, 19))+
  scale_y_continuous(limits = c(0.3, 0.7, 0.2))+
  theme(axis.title.x = element_blank(), axis.text = element_text(colour = "black"), strip.background = element_rect(size = 0.5))+
stat_pvalue_manual(stat.fvfm, label = "{p.signif} {p}", tip.length = 0.005, hide.ns = FALSE, size = 2.5)+
  guides(fill = "none", alpha = "none")
in.fvfm

# combined plots to grid
in.plots <- plot_grid(in.etr, in.alpha,in.eK, in.fvfm, labels = c('A', 'B', 'C', 'D'),label_x = 0.08,
                       label_y = 0.985, label_size = 14,ncol = 1, align = "v", byrow = F, hjust =2)

in.plots
ggsave("PAM_insitu_byday_stat_StyloSpat_reclour.jpeg", plot = in.plots, width = 12, height = 18,dpi=300, 
       units = "cm")

ggsave("PAM_insitu_byday_stat_StyloSpat_reclour.pdf", plot = in.plots, width = 12, height = 18,dpi=300, 
       units = "cm")




# PCA Photophysiology

library(stats)
library(tidyverse)
library(factoextra)

a = read.csv("PAM_PCA.csv")
head(a)
str(a)
a$id = factor(a$id, levels = c("SS8","SS60","SD8","SD60"))

pca2 = prcomp(a[,c(1:4)], scale. = TRUE)
pca2

fviz_eig(pca2)

fviz_pca_var(pca2,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

fviz_pca_biplot(pca2, repel = TRUE,
                mean.point = FALSE, # remove mean point
                pointshape = 19, pointsize = 2,
                label = "var", # labels to display
                col.var = "black", # Variables color
                col.ind = a$id,# Individuals color
)+
  scale_color_manual(values = c("#f0f921", "#FCA510", "lightblue", "blue"))+
  theme(legend.title=element_blank())

