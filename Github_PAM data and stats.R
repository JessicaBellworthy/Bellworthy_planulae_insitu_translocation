
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

write.csv(file="rlc.parameters.insitu_60days_StyloSpat2021.csv",rlc.parameters)


#### check statistical differences
library(lme4)
library(nlme)
library(lmerTest)
library(MASS)
library(car)
library(predictmeans)
library(ggpubr)

### manually copy the FVFM to the rlc parameter data file and and the metadata e.g. treatment

d <- read_csv('rlc.parameters.insitu_60days_StyloSpat2021.csv')
com = read_csv('rlc.parameters_ExSitu2020&2021combined.csv')
com$treatment <- as.factor(com$treatment)
com$age <- as.factor(com$age)


head(d)
str(d)
d$id <- as.factor(d$id)
d$chamber <- as.factor(d$chamber)
d$treatment <- as.factor(d$treatment)

sh = subset(d, treatment == "SS")
deep = subset(d, treatment == "SD")


# # # # # # # # #  CAUTION - these tests group together the two time points # # # # # # # # # # 
# rETR
gghistogram(d, x = "ETRmax", rug = TRUE, fill = "treatment", bins = 10)
etrMAX <- lmer(ETRmax~treatment + (1|chamber), data=d)
summary(etrMAX)
anova(etrMAX) 
# No sig diff p = 0.6922

#normality assumptions - OK
shapiro.test(sh$ETRmax)
shapiro.test(deep$ETRmax)
leveneTest(ETRmax~treatment,d=d)


###alpha
# need to remove one low outlier from the SD group, row 4
a = subset(d, alpha > 0.45)
Alpha <- lmer(alpha~treatment + (1|chamber), data=a)
anova(Alpha)
# No sig diff p = 0.05188

#normality assumptions - OK
shapiro.test(sh$alpha)
shapiro.test(deep$alpha)
leveneTest(alpha~treatment,d=a)

###eK
eK = lmer(Ek~treatment + (1|chamber), data=d)
summary(eK)
anova(eK)
# No sig diff p = 0.6588

#normality assumptions -ok
shapiro.test(sh$Ek)
shapiro.test(deep$Ek)
leveneTest(Ek~treatment,d=d)

###FV/FM
fvfm = lmer(fvfm~treatment + (1|chamber), data=d)
anova(fvfm)
# significant  F = 77.962, p = 0.000000001977

#normality assumptions
shapiro.test(sh$fvfm)
shapiro.test(deep$fvfm)
leveneTest(fvfm~treatment,d=d)


###beta
b = lmer(beta~treatment + (1|chamber), data=d)
anova(b)
# No sig diff p = 0.677


#normality assumptions
shapiro.test(sh$beta) # not normal
shapiro.test(deep$beta) #not normal
leveneTest(beta~treatment,d=d)

## summary of results 
library(plotrix)

Sum_all <- d %>% 
  group_by(treatment) %>% 
  summarise_each(funs(mean(., na.rm=TRUE), n = sum(!is.na(.)), max(., na.rm = TRUE), min(., na.rm = TRUE),
                      se = sd(., na.rm=TRUE)/sqrt(sum(!is.na(.)))), alpha:fvfm)

View(Sum_all)
write.csv(file="rlc.summarytable.InSitu_60days.csv", Sum_all)


#### Graphics in Box plots #########
#a = read.csv("rlc.parameters.insitu_8days_StyloSpat2021.csv")

com$treatment = factor(com$treatment, levels = c("DD", "SD","DS", "SS"))
com = subset(com, age != 7)
View(com)
#d$treatment = factor(d$treatment, levels = c("DD", "SD","DS", "SS"))
#d$treatment = factor(d$treatment, levels = c("SD", "SS"))

p.etr = ggplot(com, aes(y = ETRmax, x = treatment)) +
  geom_boxplot(outlier.shape = NA, aes(fill= treatment, alpha = 0.8), fatten = 0.5, alpha = 0.7, lwd = 0.2)+
  geom_jitter(position = position_jitter(width = .25), size = 0.25)+
facet_wrap(facets = "age", nrow = 1)+
  theme_classic()+
 scale_fill_viridis_d()+
  #scale_fill_manual(values = c('#31688EFF', '#FDE725FF'))+
  theme(axis.title.x = element_blank(), axis.text = element_text(colour = "black"), strip.background = element_rect(size = 0.5))+
 # stat_pvalue_manual(stat.etr, label = "{p.signif} {p}", tip.length = 0.005, hide.ns = TRUE, size = 1.7)+
  guides(fill = FALSE, alpha = FALSE)
p.etr

p.alpha= ggplot(com, aes(y = alpha, x = treatment)) +
  geom_boxplot(outlier.shape = NA, aes(fill= treatment, alpha = 0.8), fatten = 0.5, alpha = 0.7, lwd = 0.2)+
  geom_jitter(position = position_jitter(width = .25), size = 0.25)+
 facet_wrap(facets = "age", nrow = 1)+
  theme_classic()+
  scale_fill_viridis_d()+
  #scale_fill_manual(values = c('#31688EFF', '#FDE725FF'))+
  theme(axis.title.x = element_blank(), axis.text = element_text(colour = "black"), strip.background = element_rect(size = 0.5))+
  guides(fill = FALSE, alpha = FALSE)
p.alpha

p.beta= ggplot(com, aes(y = beta, x = treatment)) +
  geom_boxplot(outlier.shape = NA, aes(fill= treatment, alpha = 0.8), fatten = 0.5, alpha = 0.7, lwd = 0.2)+
  geom_jitter(position = position_jitter(width = .25), size = 0.25)+
  facet_wrap(facets = "age", nrow = 1)+
  theme_classic()+
  scale_fill_viridis_d()+
  #scale_fill_manual(values = c('#31688EFF', '#FDE725FF'))+
  theme(axis.title.x = element_blank(), axis.text = element_text(colour = "black"), strip.background = element_rect(size = 0.5))+
  guides(fill = FALSE, alpha = FALSE)
p.beta

p.eK= ggplot(com, aes(y = Ek, x = treatment)) +
  geom_boxplot(outlier.shape = NA, aes(fill= treatment, alpha = 0.8), fatten = 0.5, alpha = 0.7, lwd = 0.2)+
  geom_jitter(position = position_jitter(width = .25), size = 0.25)+
  facet_wrap(facets = "age", nrow = 1)+
  theme_classic()+
  scale_fill_viridis_d()+
  scale_y_continuous(limits = c(0, 400))+
  labs(y= ~Ek ~(?mol ~m^-2 ~s^-1))+
# scale_fill_manual(values = c('#31688EFF', '#FDE725FF'))+
  theme(axis.title.x = element_blank(), axis.text = element_text(colour = "black"), strip.background = element_rect(size = 0.5))+
  guides(fill = FALSE, alpha = FALSE)
p.eK

p.fvfm= ggplot(com, aes(y = fvfm, x = treatment)) +
  geom_boxplot(outlier.shape = NA, aes(fill= treatment, alpha = 0.8), fatten = 0.5, alpha = 0.7, lwd = 0.2)+
  geom_jitter(position = position_jitter(width = .25), size = 0.25)+
 facet_wrap(facets = "age", nrow = 1)+
  theme_classic()+
  labs(y= ~F[V] ~F[M])+
scale_fill_viridis_d()+
  #scale_fill_manual(values = c('#31688EFF', '#FDE725FF'))+
  theme(axis.title.x = element_blank(), axis.text = element_text(colour = "black"), strip.background = element_rect(size = 0.5))+
  guides(fill = FALSE, alpha = FALSE)
p.fvfm


library(cowplot)

rlc.plots <- plot_grid(p.etr, p.alpha, p.eK, p.fvfm, labels = c('A', 'B', 'C', 'D'),label_x = 0.07,
                       label_y = 0.985, label_size = 14,ncol = 1, align = "v", byrow = F, hjust =2)

rlc.plots
ggsave("PAM_exsitu_publication_StyloSpat.jpeg", plot = rlc.plots, width = 15, height = 20,dpi=300, 
       units = "cm")
ggsave("PAM_exsitu_publication_StyloSpat.pdf", plot = rlc.plots, width = 15, height = 20,dpi=300, 
       units = "cm")


ggsave("PAM_insitu_60days_StyloSpat.jpeg", plot = rlc.plots, width = 4, height = 15,dpi=300, 
       units = "cm")
ggsave("PAM_insitu_60days_StyloSpat.pdf", plot = rlc.plots, width = 4, height = 15,dpi=300, 
       units = "cm")

# reduced graph for report
rlc.plots2 <- plot_grid(p.alpha, p.fvfm, labels = c('A', 'B'),label_x = 0.18,
                       label_y = 0.985, label_size = 10,ncol = 1, align = "v", byrow = F, hjust =2)

ggsave("PAM_insitu_8days_StyloSpat_report.jpeg", plot = rlc.plots2, width = 4, height = 9,dpi=300, 
       units = "cm")

# 8 days, 2021, in situ
###alpha
kruskal.test(data = a, alpha~ treatment)
t.test(data = a, alpha~ treatment)
###fvfm
kruskal.test(data = a, fvfm~ treatment)
t.test(data = a, fvfm~ treatment)
###b
kruskal.test(data = a, beta~ treatment)
t.test(data = a, beta~ treatment)
###etr
kruskal.test(data = a, ETRmax~ treatment)
t.test(data = a, ETRmax~ treatment)
###fvfm
kruskal.test(data = a, Ek~ treatment)
t.test(data = a, Ek~ treatment)



## Plot all in situ data from 2019 and 2021 and 2022 together (7, 60, and 120 days)

mytheme = theme_classic()+
  theme(axis.title.x = element_blank(), axis.text = element_text(colour = "black", size = 8), axis.title = element_text(size = 10))

d <- read.csv('rlc.parameters_TranslocationStyloSpat_alldata_edit.csv')
View(d)

in.data = subset(d, in_ex == "in")
in.data = subset(in.data, year != "2019")
in.data$year = as.factor(in.data$year)
in.data$age = as.factor(in.data$age)
in.data$treatment = factor(in.data$treatment, levels=c("SS", "SD"))

View(in.data)


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

in.beta= ggplot(in.data, aes(y = beta, x = treatment)) +
  geom_boxplot(outlier.shape = NA, aes(fill= treatment, alpha = 0.8), fatten = 0.5, alpha = 0.7, lwd = 0.2)+
  #  geom_jitter(aes(shape= year), position = position_jitter(width = .25), size = 0.4, fill = "black")+
  geom_jitter(position = position_jitter(width = .25), fill = "black", size = 0.25)+
    facet_wrap(facets = "age", nrow = 1)+
mytheme+
  scale_fill_manual(values = c("#f0f921", "#FCA510"))+
  scale_shape_manual(values = c(17, 19))+
  scale_y_continuous(limits = c(-0.2, 10, 2.5))+
  theme(axis.title.x = element_blank(), axis.text = element_text(colour = "black"), strip.background = element_rect(size = 0.5))+
  stat_pvalue_manual(stat.beta, label = "{p.signif} {p}", tip.length = 0.005, hide.ns = FALSE, size = 2.5)+
  guides(fill = "none", alpha = "none")
in.beta

in.data.Ek = subset(in.data, Ek < 90)

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


### Assumption testing for in situ all data by day
#Subset data by group
SD_8 = in.data %>% subset(age == '8' & treatment == "SD")
SS_8 = in.data %>% subset(age == '8' & treatment == "SS")
SD_60 = in.data %>% subset(age == '60' & treatment == "SD")
SS_60 = in.data %>% subset(age == '60' & treatment == "SS")
SD_120 = in.data %>% subset(age == '120' & treatment == "SD")
SS_120 = in.data %>% subset(age == '120' & treatment == "SS")

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

shapiro.test(SD_8$beta) # not normal, both only 2021 and both years combined
shapiro.test(SS_8$beta) # not normal,  both only 2021 and both years combined
leveneTest(beta~treatment,d=(subset(in.data, age == '8')))

shapiro.test(SD_8eK$Ek) 
shapiro.test(SS_8eK$Ek) 
leveneTest(Ek~treatment,d=(subset(in.data.Ek, age == '8')))

shapiro.test(SD_8$fvfm) 
shapiro.test(SS_8$fvfm) 
leveneTest(fvfm~treatment,d=(subset(in.data, age == '8')))


# 60 days - assumption OK, unless written otherwise
shapiro.test(SD_60$ETRmax.log) # not normal, both only 2021 and both years combined, log/ sqrt doesnt help, sq good   === SQ
shapiro.test(SS_60$ETRmax.log) # all years data are not normal, log transform ok, sqrt transform ok, sq not good    ==== SQRT
leveneTest(ETRmax.log~treatment,d=(subset(in.data, age == '60')))

shapiro.test(SD_60$alpha.log) 
shapiro.test(SS_60$alpha.log) # not normal, both only 2021 and both years combined, log/ srt/sq donest help
leveneTest(alpha.log~treatment,d=(subset(in.data, age == '60')))# all years data are not homogenous, log/ sqrt doesnt help, sq good   ==== SQ

shapiro.test(SD_60$beta.log) # not normal, both only 2021 and both years combined, log/ sqrt/sq doesnt help
shapiro.test(SS_60$beta.log) # not normal, both only 2021 and both years combined, log/ sqrt/sq doesnt help
leveneTest(beta.log~treatment,d=(subset(in.data, age == '60'))) # not homogeneous for 2021 alone, log doesnt help, sq good   ==== SQ

shapiro.test(SD_60$Ek.log) # all years data are not normal, log doesnt help, ====== SQ
shapiro.test(SS_60$Ek.log) # log doesnt help, sq doesnt help
leveneTest(Ek.log~treatment,d=(subset(in.data, age == '60')))

shapiro.test(SD_60$fvfm.log) #  not normal, both only 2021 and both years combined, log/sqrt/sq doesnt help
shapiro.test(SS_60$fvfm.log) 
leveneTest(fvfm.log~treatment,d=(subset(in.data, age == '60')))



# 120 days - assumption OK, unless written otherwise
shapiro.test(SD_120$ETRmax) 
shapiro.test(SS_120$ETRmax) 
leveneTest(ETRmax~treatment,d=(subset(in.data, age == '120')))

shapiro.test(SD_120$alpha) 
shapiro.test(SS_120$alpha) 
leveneTest(alpha~treatment,d=(subset(in.data, age == '120')))

shapiro.test(SD_120$beta) # not normal, both only 2021 and both years combined
shapiro.test(SS_120$beta) # not normal, both only 2021 and both years combined
leveneTest(beta~treatment,d=(subset(in.data, age == '120')))

shapiro.test(SD_120$Ek) 
shapiro.test(SS_120$Ek) 
leveneTest(Ek~treatment,d=(subset(in.data, age == '120')))

shapiro.test(SD_120$fvfm) 
shapiro.test(SS_120$fvfm) 
leveneTest(fvfm~treatment,d=(subset(in.data, age == '120'))) # not homogeneous, both only 2021 and both years combined

in.data.22 = subset(in.data, year == "2022")
in.data.21 = subset(in.data, year == "2021")

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

stat.beta <- in.data %>%
  group_by(age) %>%
  wilcox_test(beta ~ treatment) %>%
  add_significance() %>%
  p_round(digits = 3) %>%
  add_xy_position(x = "treatment")
stat.beta

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






#### Printing ex situ graphics in Box plots in larger size and quality for presentations #########
setwd("/Volumes/CBP_Students/Jessica Bellworthy/Jessica Bellworthy/R/RapidLightCurves/TranslocationStyloSpat")
com = read_csv('rlc.parameters_ExSitu2020&2021combined.csv')
com$treatment <- as.factor(com$treatment)
com$age <- as.factor(com$age)

com$treatment = factor(com$treatment, levels = c("DD", "SD","DS", "SS"))
com = subset(com, age != 7)
View(com)


p.etr = ggplot(com, aes(y = ETRmax, x = treatment)) +
  geom_boxplot(outlier.shape = NA, aes(fill= treatment, alpha = 0.8), fatten = 0.5, alpha = 0.8, lwd = 0.6)+
  geom_jitter(position = position_jitter(width = .25), size = 0.6)+
  facet_wrap(facets = "age", nrow = 1)+
  theme_classic()+
  scale_fill_viridis_d()+
  theme(axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.title.x = element_blank(), axis.text = element_text(colour = "black", size = 20), strip.background = element_rect(size = 2), axis.title = element_text(size = 24), strip.text = element_text(size = 20))+
  guides(fill = "none", alpha = "none")
p.etr

p.alpha= ggplot(com, aes(y = alpha, x = treatment)) +
  geom_boxplot(outlier.shape = NA, aes(fill= treatment, alpha = 0.8), fatten = 0.5, alpha = 0.8, lwd = 0.6)+
  geom_jitter(position = position_jitter(width = .25), size = 0.6)+
  facet_wrap(facets = "age", nrow = 1)+
  theme_classic()+
  scale_fill_viridis_d()+
  theme(axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.title.x = element_blank(), axis.text = element_text(colour = "black", size = 20), strip.background = element_rect(size = 2), axis.title = element_text(size = 24), strip.text = element_text(size = 20))+
  guides(fill = "none", alpha = "none")
p.alpha

p.beta= ggplot(com, aes(y = beta, x = treatment)) +
  geom_boxplot(outlier.shape = NA, aes(fill= treatment, alpha = 0.8), fatten = 0.5, alpha = 0.8, lwd = 0.6)+
  geom_jitter(position = position_jitter(width = .25), size = 0.6)+
  facet_wrap(facets = "age", nrow = 1)+
  theme_classic()+
  scale_fill_viridis_d()+
  theme(axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.title.x = element_blank(), axis.text = element_text(colour = "black", size = 20), strip.background = element_rect(size = 2), axis.title = element_text(size = 24), strip.text = element_text(size = 20))+
  guides(fill = "none", alpha = "none")
p.beta

p.eK= ggplot(com, aes(y = Ek, x = treatment)) +
  geom_boxplot(outlier.shape = NA, aes(fill= treatment, alpha = 0.8), fatten = 0.5, alpha = 0.8, lwd = 0.6)+
  geom_jitter(position = position_jitter(width = .25), size = 0.6)+
  facet_wrap(facets = "age", nrow = 1)+
  theme_classic()+
  scale_fill_viridis_d()+
  labs(y= ~Ek ~(μmol ~m^-2 ~s^-1))+
  scale_y_continuous(limits = c(0, 400, 100))+
  theme(axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.title.x = element_blank(), axis.text = element_text(colour = "black", size = 20), strip.background = element_rect(size = 2), axis.title = element_text(size = 24), strip.text = element_text(size = 20))+
  guides(fill = "none", alpha = "none")
p.eK

p.fvfm= ggplot(com, aes(y = fvfm, x = treatment)) +
  geom_boxplot(outlier.shape = NA, aes(fill= treatment, alpha = 0.8), fatten = 0.5, alpha = 0.8, lwd = 0.6)+
  geom_jitter(position = position_jitter(width = .25), size = 0.6)+
  facet_wrap(facets = "age", nrow = 1)+
  theme_classic()+
  scale_fill_viridis_d()+
  labs(y= ~F[V] ~F[M])+
  theme(axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.title.x = element_blank(), axis.text = element_text(colour = "black", size = 20), strip.background = element_rect(size = 2), axis.title = element_text(size = 24), strip.text = element_text(size = 20))+
  guides(fill = "none", alpha = "none")
p.fvfm


library(cowplot)

rlc.plots <- plot_grid(p.etr, p.alpha, p.eK, p.fvfm, labels = c('A', 'B', 'C', 'D'),label_x = -0.1,
                       label_y = 0.985, label_size = 40,ncol = 1, align = "v", byrow = F, hjust =2)

rlc.plots
ggsave("PAM_exsitu_presentations_StyloSpat.jpeg", plot = rlc.plots, width = 30, height = 40,dpi=600, 
       units = "cm")
ggsave("PAM_exsitu_presentations_StyloSpat.pdf", plot = rlc.plots, width = 30, height = 40,dpi=600, 
       units = "cm")
