# Settlement and Mortality InSitu 2021&2022 translocation expt.

library(dplyr)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(rstatix)
rm

setwd("/Volumes/CBP_Students/Jessica Bellworthy/Jessica Bellworthy/NSF-BSF Genomic/R/Github")
d = read.csv("set_mort_insitu_2021&2022.csv")
View(d)


d$age = as.factor(d$age)
d$treatment = as.factor(d$treatment)
d$treatment = factor(d$treatment, levels = c("SS","SD","DS","DD"))
d$chamber = as.factor(d$chamber)
d$swim = as.numeric(d$swim)
d$set = as.numeric(d$set)
d$dead = as.numeric(d$dead)


Sum_all <- d %>% 
  group_by(treatment, age, year) %>% 
  summarise_each(funs(median(., na.rm=TRUE), mean(., na.rm = TRUE), max(., na.rm = TRUE), sd(., na.rm = TRUE),
                      se = sd(., na.rm=TRUE)/sqrt(sum(!is.na(.)))), swim:pc_dead)

View(Sum_all)
write.csv(file="Set_Mort_Averages_InSitu_2021&2022.csv", Sum_all)


mytheme = theme_classic()+
  theme(axis.title.x = element_blank(), axis.text = element_text(colour = "black", size = 10), axis.title = element_text(size = 12))



# DATA EXPLORATION
# are there differences between the years or can data be combined?
# settlement by day and treatment and year
ggplot(d, aes(x= treatment, y = pc_set, fill = treatment))+
  geom_boxplot(aes(fill = treatment), outlier.shape = NA, fatten = 0.5, alpha = 0.8, lwd = 0.2)+
  geom_jitter(size = 0.2, width = 0.3)+
  facet_grid(age~year)+
  mytheme +
  guides(fill= "none")+
  labs(y = "Mean settlement %")+
  scale_y_continuous(limits = c(0, 100))+
  scale_fill_manual(values = c("#FDE725FF", "#238A8DFF", "#55C667FF",  "#404788FF"))

ggplot(d, aes(x= treatment, y = pc_dead, fill = treatment))+
  geom_boxplot(aes(fill = treatment), outlier.shape = NA, fatten = 0.5, alpha = 0.8, lwd = 0.2)+
  geom_jitter(size = 0.2, width = 0.3)+
  facet_grid(age~year)+
  mytheme +
  guides(fill= "none")+
  labs(y = "Mean mortality %")+
  scale_y_continuous(limits = c(0, 100))+
  scale_fill_manual(values = c("#FDE725FF", "#238A8DFF", "#55C667FF",  "#404788FF"))

# For both settlement and mortality data, 2021 data visually overlaps with that of 2022.

# Test this visual overlap statistically
# Assumption testing for in situ data by year
# Subset data by group
SD_8_21 = d %>% subset(age == '8' & treatment == "SD" & year == "2021")
SS_8_21 = d%>% subset(age == '8' & treatment == "SS" & year == "2021")
SD_60_21 = d%>% subset(age == '60' & treatment == "SD" & year == "2021")
SS_60_21 = d %>% subset(age == '60' & treatment == "SS" & year == "2021")

SD_8_22 = d %>% subset(age == '8' & treatment == "SD" & year == "2022")
SS_8_22 = d%>% subset(age == '8' & treatment == "SS" & year == "2022")
SD_60_22 = d%>% subset(age == '60' & treatment == "SD" & year == "2022")
SS_60_22 = d %>% subset(age == '60' & treatment == "SS" & year == "2022")


# 8 days - assumption OK, unless written otherwise
shapiro.test(SD_8_21$pc_set) 
shapiro.test(SS_8_21$pc_set) 
shapiro.test(SD_8_22$pc_set) 
shapiro.test(SS_8_22$pc_set) 

shapiro.test(SD_8_21$pc_dead) 
shapiro.test(SS_8_21$pc_dead) 
shapiro.test(SD_8_22$pc_dead) 
shapiro.test(SS_8_22$pc_dead) 


# 60 days - assumption OK, unless written otherwise
shapiro.test(SD_60_21$pc_set) 
shapiro.test(SS_60_21$pc_set) 
shapiro.test(SD_60_22$pc_set) #p = 0.044
shapiro.test(SS_60_22$pc_set) 

shapiro.test(SD_60_21$pc_dead) 
shapiro.test(SS_60_21$pc_dead) 
shapiro.test(SD_60_22$pc_dead) #p = 0.044
shapiro.test(SS_60_22$pc_dead) 


## Assumption testing for in situ data both years by day
# Subset data by group
SD_8 = d %>% subset(age == '8' & treatment == "SD")
SS_8 = d%>% subset(age == '8' & treatment == "SS")
SD_60 = d%>% subset(age == '60' & treatment == "SD")
SS_60 = d %>% subset(age == '60' & treatment == "SS")


# 8 days - assumption OK, unless written otherwise
shapiro.test(SD_8$pc_set) 
shapiro.test(SS_8$pc_set) 
levene_test(pc_set~treatment,d=(subset(d, age == '8')))

shapiro.test(SD_8$pc_dead) 
shapiro.test(SS_8$pc_dead) 
levene_test(pc_dead~treatment,d=(subset(d, age == '8')))


# 60 days - assumption OK, unless written otherwise
shapiro.test(SD_60$pc_set) # p = 0.016
shapiro.test(SS_60$pc_set) 
levene_test(pc_set~treatment,d=(subset(d, age == '60')))

shapiro.test(SD_60$pc_dead) # p = 0.016
shapiro.test(SS_60$pc_dead) 
levene_test(pc_dead~treatment,d=(subset(d, age == '60')))

# STATS for differences between years
# For settlement data there are no significant differences between 2021 & 2022
SD_8 %>%
  group_by(age) %>%
  wilcox_test(pc_set ~ year) %>%
  add_significance() %>%
  p_round(digits = 3) %>%
  add_xy_position(x = "treatment")

SS_8 %>%
  group_by(age) %>%
  wilcox_test(pc_set ~ year) %>%
  add_significance() %>%
  p_round(digits = 3) %>%
  add_xy_position(x = "treatment")

SD_60 %>%
  group_by(age) %>%
  wilcox_test(pc_set ~ year) %>%
  add_significance() %>%
  p_round(digits = 3) %>%
  add_xy_position(x = "treatment")

SS_60 %>%
  group_by(age) %>%
  wilcox_test(pc_set ~ year) %>%
  add_significance() %>%
  p_round(digits = 3) %>%
  add_xy_position(x = "treatment")


# STATS for differences between years
# Mortality all OK - no differences between years
SD_8 %>%
  group_by(age) %>%
  wilcox_test(pc_dead ~ year) %>%
  add_significance() %>%
  p_round(digits = 3) %>%
  add_xy_position(x = "treatment")

SS_8 %>%
  group_by(age) %>%
  wilcox_test(pc_dead ~ year) %>%
  add_significance() %>%
  p_round(digits = 3) %>%
  add_xy_position(x = "treatment")

SD_60 %>%
  group_by(age) %>%
  wilcox_test(pc_dead ~ year) %>%
  add_significance() %>%
  p_round(digits = 3) %>%
  add_xy_position(x = "treatment")

SS_60 %>%
  group_by(age) %>%
  wilcox_test(pc_dead ~ year) %>%
  add_significance() %>%
  p_round(digits = 3) %>%
  add_xy_position(x = "treatment")

###### So years can be grouped together   ######



#Stats when grouping data from both years together
stat_dead = d %>%
  group_by(age) %>%
  wilcox_test(pc_dead ~ treatment) %>%
  add_significance() %>%
  p_round(digits = 3) %>%
  add_xy_position(x = "treatment")
stat_dead 
# use p adjusted value for multiple comparisons - not significant in any case

stat_set = d %>%
  group_by(age) %>%
  wilcox_test(pc_set ~ treatment) %>%
  add_significance() %>%
  p_round(digits = 3) %>%
  add_xy_position(x = "treatment")
stat_set
# use p adjusted value for multiple comparisons - not significant in any case



# Mean settlement and mortality by day and treatment and year
p.set = ggplot(d, aes(x= treatment, y = pc_set))+
  geom_boxplot(aes(fill = treatment), outlier.shape = NA, fatten = 0.5, alpha = 0.8, lwd = 0.2)+
  geom_jitter(size = 0.2, width = 0.3)+
  facet_wrap(facets = "age")+
  mytheme +
  guides(fill= "none")+
  stat_pvalue_manual(stat_set, label = "{p}", tip.length = 0.005, hide.ns = TRUE, size = 2.5)+
  labs(y = "Mean settlement %")+
  scale_y_continuous(limits = c(-2, 101))+
  scale_fill_manual(values = c("#FDE725FF", "#238A8DFF", "#55C667FF",  "#404788FF"))
p.set

p.dead = ggplot(d, aes(x= treatment, y = pc_dead))+
  geom_boxplot(aes(fill = treatment), outlier.shape = NA, fatten = 0.5, alpha = 0.8, lwd = 0.2)+
  geom_jitter(size = 0.2, width = 0.3)+
  facet_wrap(facets = "age")+
  mytheme +
  guides(fill= "none")+
  stat_pvalue_manual(stat_dead, label = "{p}", tip.length = 0.005, hide.ns = TRUE, size = 2.5)+
  labs(y = "Mean mortality %")+
  scale_y_continuous(limits = c(-2, 101))+
  scale_fill_manual(values = c("#FDE725FF", "#238A8DFF", "#55C667FF",  "#404788FF"))
p.dead



#####################
# Kaplan Meier Curves

library(purrr)
library(tidyr)
library(tidyverse)
library(survminer)
library(survival)
library(ggsci)
library(splitstackshape)
library(MuMIn)


################### Non filtered data, all data ############### 
library(splitstackshape)
library(survival)
library(survminer)


setwd("/Volumes/CBP_Students/Jessica Bellworthy/Jessica Bellworthy/NSF-BSF Genomic/R")
k = read.csv("KMinsituAll.csv")
View(k)


Sum_all <- k %>% 
  group_by(treatment, age_days, dead_status) %>% 
  summarise_each(funs(sum), count)
View(Sum_all)


K = expandRows(k, "count")  #Expands row by count and Removes rows where count = 0 (n = 8 rows removed)
View(K)

K$treatment = factor(K$treatment, levels = c("SS","SD","DS","DD"))

### Fit KM curves ###

# Settlement
fitK = survfit(Surv(age_days, set_status)~ treatment, data = K)
summary(fitK)

names(fitK$strata) <- gsub("treatment=", "", names(fitK$strata))

setplot <- ggsurvplot(fitK, size = 2,
                      pval = TRUE, pval.method = TRUE,
                      pval.coord = c(0.4, 0.4),
                      pval.method.coord = c(0.4, 0.48),
                      fun = "event",
                      ggtheme = theme_classic()+
                        theme(axis.ticks = element_line(size = 2), axis.text = element_text(colour = "black", size = 16),
                              axis.line = element_line(size = 1), axis.title = element_text(size = 22), legend.background = element_blank(), legend.text = element_text(size = 16)),
                      palette = c("#FDE725FF", "#238A8DFF", "#55C667FF",  "#404788FF"),
                      legend = c(0.15, 0.9), legend.title = "", xlab="Days", ylab = "Settlement probability", legend.labs = c("SS","SD","DS","DD"))

setplot


# Settlement - Pairwise comparison between treatments, fixed effects only, p value adj. method BH
res <- pairwise_survdiff(Surv(age_days, set_status) ~ treatment,
                         data= K)
res   
# Two shallow treatments are significantly more likely to settle than the deep treatments p <0.0001 (Log-rank test)



# Mortality
fitK = survfit(Surv(age_days, dead_status)~ treatment, data = K)
summary(fitK)

names(fitK$strata) <- gsub("treatment=", "", names(fitK$strata))

deadplot <- ggsurvplot(fitK, size = 2,
                       pval = TRUE, pval.method = TRUE,
                       pval.coord = c(0.15, 0.20),
                       pval.method.coord = c(0.25, 0.30),
                      palette = c("#FDE725FF", "#238A8DFF", "#55C667FF",  "#404788FF"),
                      ggtheme = theme_classic()+
                        theme(axis.ticks = element_line(size = 2), axis.text = element_text(colour = "black", size = 16),
                              axis.line = element_line(size = 1), axis.title = element_text(size = 22), legend.text = element_text(size = 16), legend.background = element_blank()),
                      legend = c(0.15, 0.58), legend.title = "", xlab="Days", ylab = "Survival probability", legend.labs = c("SS","SD","DS","DD"))

deadplot


#Mortality Pairwise comparison between treatments, fixed effects only, p value adj. method BH
res <- pairwise_survdiff(Surv(age_days, dead_status) ~ treatment,
                         data= K)
res  
# All treatments significantly different from each other, except SS and SD



# Print/ export KM plots from plot viewer window

# Combine plots (will not print together with the KM graphs. These need saving separately)
p <- plot_grid(deadplot, setplot, p.dead, p.set, labels = c('A', 'B', 'C', 'D'),
                         label_y = 0.985, label_size = 11,ncol = 2, align = "h", byrow = F)

ggsave("Set.Dead_InSitu_2021&2022_final.jpeg", plot = p, width = 14, height = 10,dpi=300, 
       units = "cm")
ggsave("Set.Dead_InSitu_2021&2022_final.pdf", plot = p, width = 14, height = 10,dpi=300, 
       units = "cm")


# For combining only the box plot graphs
p <- plot_grid(p.dead, p.set, labels = c('A', 'B'),
               label_y = 0.985, label_size = 11, ncol = 1, align = "h", byrow = F)
ggsave("Set.Dead_InSitu_2021&2022_jitter.jpeg", plot = p, width = 10, height = 14,dpi=600, 
       units = "cm")
