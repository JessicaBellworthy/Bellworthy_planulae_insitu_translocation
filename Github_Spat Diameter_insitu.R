# Depth Gradient Translocation Expt., Spat Size Analysis, In situ
library(dplyr)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(car)
library(rstatix)



setwd("/Volumes/CBP_Students/Jessica Bellworthy/Jessica Bellworthy/NSF-BSF Genomic/R/Github")
data = read.csv("SpatSize_insitu_2021&2022.csv")
View(data)

data$age_days = as.factor(data$age)
data$treatment = as.factor(data$treatment)
data$chamber = as.factor(data$chamber)
data$year = as.factor(data$year)
data$treatment = factor(data$treatment, levels=c("SS", "SD"))
str(data)

Sum_all <- data %>% 
  group_by(treatment, age_days) %>% 
  summarise_each(funs(mean(., na.rm=TRUE), sd = sd, n = sum(!is.na(.)),
                      se = sd(., na.rm=TRUE)/sqrt(sum(!is.na(.)))), diameter)

View(Sum_all)
write.csv(file="SpatSizeAverages.csv", Sum_all)


mytheme = theme_classic()+
  theme(axis.title.x = element_blank(), axis.text = element_text(colour = "black", size = 6), axis.title = element_text(size = 9))


### Assumption testing for in situ all data by day
#Subset data by group
head(data)
SD_8 = data %>% subset(age_days == 8 & treatment == "SD")
SS_8 = data  %>% subset(age_days == 8 & treatment == "SS")
SD_60 = data  %>% subset(age_days == 60 & treatment == "SD")
SS_60 = data  %>% subset(age_days == 60 & treatment == "SS")

# 8 days - assumption OK, unless written otherwise
shapiro.test(SD_8$diameter) 
shapiro.test(SS_8$diameter) 
leveneTest(diameter~treatment,d=(subset(data , age_days == 8))) # not homogeneous

# 60 days - assumption OK, unless written otherwise
shapiro.test(SD_60$diameter) 
shapiro.test(SS_60$diameter) 
leveneTest(diameter~treatment,d=(subset(data , age_days == 60))) # not homogeneous


# Is there a difference between 2021 & 2022 data?
data%>%
  group_by(treatment, age_days) %>%
  wilcox_test(diameter ~ year) %>%
  add_significance() %>%
  p_round(digits = 3) %>%
  add_xy_position(x = "treatment")

# Years are not significantly different from each other.
# There is however a difference at day 8 SD between years (stat = 120, p = 0.01)
# Combine data but mark years with different symbols



# Is there an increase in size between 8 and 60 days in situ?
data%>%
  group_by(treatment) %>%
  wilcox_test(diameter ~ age_days) %>%
  add_significance() %>%
  p_round(digits = 3) %>%
  add_xy_position(x = "treatment")

# There is a significant increase in size between day 8 and day 60 in the SD treatment only (stat = 678, p = 0.000526)


# Is there a difference in spat size between treatments at 8 and 60 days in situ?
# STATS for graphics - not homogenous variances so non-parametric Wilcoxon test used
stat.size <- data%>%
  group_by(age_days) %>%
  wilcox_test(diameter ~ treatment) %>%
  add_significance() %>%
  p_round(digits = 3) %>%
  add_xy_position(x = "treatment")
stat.size

p5 = ggplot(data, aes(y = (diameter), x = treatment)) +
  geom_boxplot(aes(fill = treatment), outlier.shape = "", fatten = 0.5, alpha = 0.7, lwd = 0.2)+
  geom_jitter(aes(shape= year), position = position_jitter(width = .25), fill = "black", size = 0.4)+
  facet_wrap(facets = "age_days", nrow = 1)+
  scale_fill_manual(values = c('#FDE725FF','#31688EFF'))+
  scale_shape_manual(values = c(19, 17))+
  labs(y= ~Spat ~Diameter ~(mm))+
  mytheme+
  theme(axis.title.x = element_blank(), axis.text = element_text(colour = "black"))+
  stat_pvalue_manual(stat.size, label = "{p.signif} {p}", tip.length = 0.005, hide.ns = FALSE, size = 2.5)+
  guides(fill = "none", alpha = "none", shape = "none")+
  scale_y_continuous(limits = c(0, 3.5), breaks = seq(0,3.5,0.5))
p5


ggsave("SpatDiameter_insituStyloSpat2021&2022.jpeg", plot = p5, width = 8, height = 6, dpi=300, 
       units = "cm")
ggsave("SpatDiameter_insituStyloSpat2021&2022.pdf", plot = p5, width = 8, height = 6,dpi=300, 
       units = "cm")


