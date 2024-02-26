library(lmerTest)
library(lme4)
library(ggplot2)
library(dplyr)
library("FSA")
library(cowplot)
library(rstatix)
library(MuMIn)
library(ggpubr)
library(magick)
library(nlme)
library(MASS)
library(car)
library(predictmeans)
library(ggpubr)

setwd("/Volumes/CBP_Students/Jessica Bellworthy/Jessica Bellworthy/NSF-BSF Genomic/R/Github")
d = read.csv("physioSpats_insitu.2022.csv")


d$age = as.factor(d$age)
d$treatment = as.factor(d$treatment)
d$treatment = factor(d$treatment, levels = c("SS","SD","DS","DD"))
d$replicate = as.factor(d$replicate)
str(d)

Sum_all <- d %>% 
  group_by(treatment, age) %>% 
  summarise_each(funs(median(., na.rm=TRUE), mean(., na.rm = TRUE), sd(., na.rm = TRUE),
                      se = sd(., na.rm=TRUE)/sqrt(sum(!is.na(.)))), ng.chl.spat:ng.chl.ug.protein)

View(Sum_all)
write.csv(file="SpatPhysio_Averages_InSitu_2022.csv", Sum_all)


########                                       ############
## Assumption testing for in situ data by year
#Subset data by group
SD_8 = d %>% subset(age == '8' & treatment == "SD")
SS_8 = d%>% subset(age == '8' & treatment == "SS")
SD_60 = d%>% subset(age == '60' & treatment == "SD")
SS_60 = d %>% subset(age == '60' & treatment == "SS")


# 8 days - assumption OK, unless written otherwise
shapiro.test(SD_8$ng.chl.spat) 
shapiro.test(SS_8$ng.chl.spat) 
shapiro.test(SD_8$ng.chl.spat) 
shapiro.test(SS_8$ng.chl.spat) 
leveneTest(ng.chl.spat~treatment,d=(subset(d, age == '8')))

shapiro.test(SD_8$ug.protein.spat) 
shapiro.test(SS_8$ug.protein.spat) 
shapiro.test(SD_8$ug.protein.spat) 
shapiro.test(SS_8$ug.protein.spat) 
leveneTest(ug.protein.spat~treatment,d=(subset(d, age == '8')))

shapiro.test(SD_8$cells.spat) 
shapiro.test(SS_8$cells.spat) 
shapiro.test(SD_8$cells.spat) 
shapiro.test(SS_8$cells.spat) 
leveneTest(cells.spat~treatment,d=(subset(d, age == '8')))

# 60 days - assumption OK, unless written otherwise
shapiro.test(SD_60$ng.chl.spat) 
shapiro.test(SS_60$ng.chl.spat) 
shapiro.test(SD_60$ng.chl.spat) 
shapiro.test(SS_60$ng.chl.spat) 
leveneTest(ng.chl.spat~treatment,d=(subset(d, age == '60')))

shapiro.test(SD_60$ug.protein.spat) 
shapiro.test(SS_60$ug.protein.spat) 
shapiro.test(SD_60$ug.protein.spat) 
shapiro.test(SS_60$ug.protein.spat) 
leveneTest(ug.protein.spat~treatment,d=(subset(d, age == '60')))

shapiro.test(SD_60$cells.spat) 
shapiro.test(SS_60$cells.spat) 
shapiro.test(SD_60$cells.spat) 
shapiro.test(SS_60$cells.spat) 
leveneTest(cells.spat~treatment,d=(subset(d, age == '60')))


# STATS for graphics
stat.chl <- d %>%
  group_by(age) %>%
  t_test(ng.chl.spat ~ treatment) %>%
  add_significance() %>%
  p_round(digits = 3) %>%
  add_xy_position(x = "treatment")
stat.chl

stat.prot <- d %>%
  group_by(age) %>%
  t_test(ug.protein.spat ~ treatment) %>%
  add_significance() %>%
  p_round(digits = 3) %>%
  add_xy_position(x = "treatment")
stat.prot

stat.cells <- d %>%
  group_by(age) %>%
  t_test(cells.spat ~ treatment) %>%
  add_significance() %>%
  p_round(digits = 3) %>%
  add_xy_position(x = "treatment")
stat.cells



##### GRAPHICS ####
mytheme = theme_classic()+
  theme(axis.title.x = element_blank(), axis.text = element_text(colour = "black", size = 8), axis.title = element_text(size = 10))

chl= ggplot(d, aes(y = ng.chl.spat, x = treatment))+
geom_boxplot(outlier.shape = NA, aes(fill= treatment, alpha = 0.8), fatten = 0.5, alpha = 0.7, lwd = 0.2)+
 geom_jitter(position = position_jitter(width = .25), size = 0.7, fill = "black")+
facet_wrap(facets = "age", nrow = 1)+
 mytheme+
  labs(y= ~ng ~chlorophyll ~ per ~spat)+
  scale_fill_manual(values = c("#f0f921", "#FCA510"))+
scale_y_continuous(limits = c(0, 700, 200))+
  theme(axis.title.x = element_blank(), axis.text = element_text(colour = "black"), strip.background = element_rect(size = 0.5))+
  stat_pvalue_manual(stat.chl, label = "{p.signif} {p}", tip.length = 0.005, hide.ns = FALSE, size = 2.5)+
  guides(fill = "none", alpha = "none")
chl


prot= ggplot(d, aes(y = ug.protein.spat, x = treatment))+
  geom_boxplot(outlier.shape = NA, aes(fill= treatment, alpha = 0.8), fatten = 0.5, alpha = 0.7, lwd = 0.2)+
  geom_jitter(position = position_jitter(width = .25), size = 0.7, fill = "black")+
  facet_wrap(facets = "age", nrow = 1)+
  mytheme+
  labs(y= ~μg ~protein ~ per ~spat)+
  scale_fill_manual(values = c("#f0f921", "#FCA510"))+
  scale_y_continuous(limits = c(0, 270, 50))+
  theme(axis.title.x = element_blank(), axis.text = element_text(colour = "black"), strip.background = element_rect(size = 0.5))+
  stat_pvalue_manual(stat.prot, label = "{p.signif} {p}", tip.length = 0.005, hide.ns = FALSE, size = 2.5)+
  guides(fill = "none", alpha = "none")
prot



cell= ggplot(d, aes(y = cells.spat, x = treatment))+
  geom_boxplot(outlier.shape = NA, aes(fill= treatment, alpha = 0.8), fatten = 0.5, alpha = 0.7, lwd = 0.2)+
  geom_jitter(position = position_jitter(width = .25), size = 0.7, fill = "black")+
  facet_wrap(facets = "age", nrow = 1)+
  mytheme+
  labs(y= ~symbiont ~cells ~per ~spat)+
  scale_fill_manual(values = c("#f0f921", "#FCA510"))+
  scale_y_continuous(limits = c(5000, 23000, 5000))+
  theme(axis.title.x = element_blank(), axis.text = element_text(colour = "black"), strip.background = element_rect(size = 0.5))+
  stat_pvalue_manual(stat.cells, label = "{p.signif} {p}", tip.length = 0.005, hide.ns = FALSE, size = 2.5)+
  guides(fill = "none", alpha = "none")
cell


# combine plots to grid
physio.plots <- plot_grid(chl, cell,prot,labels = c('A', 'B', 'C'),label_x = 0.028,
                      label_y = 0.985, label_size = 14,ncol = 1, align = "v", byrow = F, hjust =2)

physio.plots
ggsave("physio_insituStyloSpat_2022_reclour.jpeg", plot = physio.plots, width = 12, height = 16,dpi=300, 
       units = "cm")

ggsave("physio_insituStyloSpat_2022_reclour.pdf", plot = physio.plots, width = 12, height = 16,dpi=300, 
       units = "cm")

