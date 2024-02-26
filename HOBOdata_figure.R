library(ggplot2)
library(dplyr)
library(plotrix)
library(cowplot)
library(scales) # to access breaks/formatting functions



# 2021 in situ April - June
setwd("/Volumes/CBP_Students/Jessica Bellworthy/Jessica Bellworthy/NSF-BSF Genomic/R/HOBO data")

d<-read.csv("HOBO data NSF BSF Genomic 2021 R.csv", strip.white=T)
d$Depth <- factor(d$Depth, levels=c("Shallow","Deep"))
d$Date<-as.Date(d$Date, "%m/%d/%Y")
d$Date_Time<- as.POSIXct(paste(d$Date, d$Time), "%Y-%m-%d %H:%M", tz = "GMT", usetz = FALSE)

View(d)
str(d)

sum(is.na(d))

temp<-d%>%
  group_by(Depth,Date)%>%
  summarize_at(vars(Temp), list(mean = mean, sd = sd,  min = min, max = max), na.rm=T)
View(temp)
write.csv(temp, file = "temperature_dailymean_2021.csv")

count = d%>%
  group_by(Depth,Date)%>% summarise(n = n())
View(count)

count = count %>% filter(n != 288)

# Change Date to date format
temp$Date<- as.POSIXct(temp$Date)


temp.plot.21 = 
  ggplot(d, aes(y=Temp, x=Date_Time, color=Depth))+ 
  geom_point(size=0.05, alpha = 0.1)+
  geom_line(data=temp, aes(y=mean, x=Date),size=1.25,alpha=0.8)+
  scale_fill_manual("Depth", values=c("Shallow"="red", "Deep"= "dark blue"))+
  scale_color_manual("Depth", values=c("Shallow"= "red","Deep"="dark blue"))+
  scale_y_continuous(limits = c(20, 30), breaks = seq(20, 30, 2)) +
  labs(y = "Temperature (°C)")+
  theme_classic()+
  theme(axis.text.x=element_text(vjust=0.5,size=12),#angling the labels on the x-axis
        panel.background= element_rect(fill=NA, color='black'),#this is making the black box around the graph
        strip.background = element_blank(), 
        strip.text = element_blank(),
        legend.title = element_blank(),
        legend.text =element_text(size = 12),
        legend.background = element_blank(),
        legend.position=c(0.85,0.2),
        axis.text.y = element_text(vjust=0.5, size=12), #making the axis text larger 
        axis.title.x = element_blank(),#making the axis title larger 
        axis.title.y = element_text(size=12))#making the axis title larger 

temp.plot.21

light<-d%>%
  group_by(Depth,Date)%>%
  summarize_at(vars(Light), list(mean = mean, sd = sd, min = min, max = max), na.rm=T)
View(light)
write.csv(light, file = "light_dailymean_2021.csv")

light$Date<- as.POSIXct(light$Date)

light.plot.21 = ggplot(d, aes(y=Light/10000, x=Date_Time, color=Depth))+ 
  #geom_point(size=0.05, alpha = 0.1)+
  geom_line(data=light, aes(y=mean/10000, x=Date),size=1.25,alpha=0.8)+
  scale_fill_manual("Depth", values=c("Shallow"="red", "Deep"= "dark blue"))+
  scale_color_manual("Depth", values=c("Shallow"= "red","Deep"="dark blue"))+
  scale_y_continuous(limits = c(0, 1.5), breaks = seq(0, 1.5, 0.5)) +
  labs(y = "Light Intensity (lux x10000)")+
  theme_classic()+
  theme(axis.text.x=element_text(vjust=0.5,size=12),#angling the labels on the x-axis
        panel.background= element_rect(fill=NA, color='black'),#this is making the black box around the graph
        strip.background = element_blank(), 
        strip.text = element_blank(),
        legend.title = element_blank(),
        legend.text =element_text(size = 12),
        legend.background = element_blank(),
        legend.position=c(0.85,0.9),
        axis.text.y = element_text(vjust=0.5, size=12), #making the axis text larger 
        axis.title.x = element_blank(),#making the axis title larger 
        axis.title.y = element_text(size=12))#making the axis title larger 

light.plot.21

HOBO.plots <- plot_grid(temp.plot.21, light.plot.21, labels = c('A', 'B'),
                        label_x = 0.1, label_y = 0.985, label_size = 14,ncol = 1, align = "h", byrow = F, hjust =2)
HOBO.plots
ggsave("HOBOplots_insitu_2021.jpeg", plot = HOBO.plots, width = 10, height = 18,dpi=600, 
       units = "cm")



# 2022 in situ March - July
setwd("/Volumes/CBP_Students/Jessica Bellworthy/Jessica Bellworthy/NSF-BSF Genomic/R/HOBO data")

d<-read.csv("HOBO data NSF BSF Genomic 2022 R.csv", strip.white=T)
d$Depth <- factor(d$Depth, levels=c("Shallow","Deep"))
d$Date<-as.Date(d$Date, "%d/%m/%Y")
d$Date_Time<- as.POSIXct(paste(d$Date, d$Time), "%Y-%m-%d %H:%M", tz = "GMT", usetz = FALSE)

View(d)
str(d)

temp<-d%>%
  group_by(Depth,Date)%>%
  summarize_at(vars(Temp), list(mean = mean, sd = sd,  min = min, max = max), na.rm=T)
View(temp)
write.csv(temp, file = "temperature_dailymean_2022.csv")

count = d%>%
  group_by(Depth,Date)%>% summarise(n = n())
count = count %>% filter(n != 288)

temp$Date<- as.POSIXct(temp$Date)

temp.plot.22 = 
  ggplot(d, aes(y=Temp, x=Date_Time, color=Depth))+ 
  geom_point(size=0.05, alpha = 0.1)+
  geom_line(data=temp, aes(y=mean, x=Date),size=1.25,alpha=0.8)+
  scale_fill_manual("Depth", values=c("Shallow"="red", "Deep"= "dark blue"))+
  scale_color_manual("Depth", values=c("Shallow"= "red","Deep"="dark blue"))+
  scale_y_continuous(limits = c(20, 30), breaks = seq(20, 30, 2)) +
  labs(y = "Temperature (°C)")+
  theme_classic()+
  theme(axis.text.x=element_text(vjust=0.5,size=12),#angling the labels on the x-axis
        panel.background= element_rect(fill=NA, color='black'),#this is making the black box around the graph
        strip.background = element_blank(), 
        strip.text = element_blank(),
        legend.title = element_blank(),
        legend.text =element_text(size = 12),
        legend.background = element_blank(),
        legend.position=c(0.85,0.2),
        axis.text.y = element_text(vjust=0.5, size=12), #making the axis text larger 
        axis.title.x = element_blank(),#making the axis title larger 
        axis.title.y = element_text(size=12))#making the axis title larger 

temp.plot.22

light<-d%>%
  group_by(Depth,Date)%>%
  summarize_at(vars(Light), list(mean = mean, sd = sd,  min = min, max = max), na.rm=T)
View(light)
write.csv(light, file = "light_dailymean_2022.csv")

light$Date<- as.POSIXct(light$Date)

light.plot.22 = ggplot(d, aes(y=Light/10000, x=Date_Time, color=Depth))+ 
  # geom_point(size=0.05, alpha = 0.1)+
  geom_line(data=light, aes(y=mean/10000, x=Date),size=1.25,alpha=0.8)+
  scale_fill_manual("Depth", values=c("Shallow"="red", "Deep"= "dark blue"))+
  scale_color_manual("Depth", values=c("Shallow"= "red","Deep"="dark blue"))+
  scale_y_continuous(limits = c(0, 1.5), breaks = seq(0, 1.5, 0.5)) +
  labs(y = "Light Intensity (lux x10000)")+
  theme_classic()+
  theme(axis.text.x=element_text(vjust=0.5,size=12),#angling the labels on the x-axis
        panel.background= element_rect(fill=NA, color='black'),#this is making the black box around the graph
        strip.background = element_blank(), 
        strip.text = element_blank(),
        legend.title = element_blank(),
        legend.text =element_text(size = 12),
        legend.background = element_blank(),
        legend.position=c(0.85,0.9),
        axis.text.y = element_text(vjust=0.5, size=12), #making the axis text larger 
        axis.title.x = element_blank(),#making the axis title larger 
        axis.title.y = element_text(size=12))#making the axis title larger 

light.plot.22


HOBO.plots <- plot_grid(temp.plot.22, light.plot.22, labels = c('A', 'B'),
                       label_x = 0.1, label_y = 0.985, label_size = 14, ncol = 1, align = "h", byrow = F, hjust =2)
HOBO.plots
ggsave("HOBOplots_insitu_2022.jpeg", plot = HOBO.plots, width = 10, height = 18,dpi=600, 
       units = "cm")



# Combine all plots to one figure

HOBO.plots <- plot_grid(temp.plot.21, light.plot.21, temp.plot.22, light.plot.22, labels = c('A', 'C', 'B', 'D'),
                        label_x = 0.08, label_y = 0.985, label_size = 14,ncol = 2, align = "h", byrow = F, hjust =2)
HOBO.plots
ggsave("HOBOplots_insitu_2021&2022.jpeg", plot = HOBO.plots, width = 20, height = 18,dpi=600, 
       units = "cm")

