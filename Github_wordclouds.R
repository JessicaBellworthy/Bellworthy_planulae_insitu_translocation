#Making work clouds for FE by life stage

library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(ggwordcloud)


# Making word clouds from my data
# First step: Quality control and merge redundant terms in the csv


#Try with ggplot to make the word clouds - more customization
#can also combine all data rows together and facet_wrap by age/ stage
setwd("/Volumes/CBP_Students/Jessica Bellworthy/HIVE - zoom and scripts/output/working data/SemanticSimilarity")

#### ADULTS
adults = read.csv("adultTerms.csv")
adultUp = filter(adults, direction == 'upregulated')
View(adultUp)
adultDown = filter(adults, direction == 'downregulated')
View(adultDown)
# removes regulated terms

tb = table(adultUp$term_name)
tb1 = data.frame(word = names(tb), freq =as.numeric(tb))
tb1 = as.data.frame(tb1)
tb1$direction <- "green"
View(tb1)

tb2 = table(adultDown$term_name)
tb2 = data.frame(word = names(tb2), freq =as.numeric(tb2))
tb2 = as.data.frame(tb2)
tb2$direction <- "red"
View(tb2)

tb3 = rbind(tb1, tb2) # bind together up and down regulated
View(tb3)
str(tb3)

set.seed(1)# for reproducibility 
data = tb3


data <- data %>%
  filter(freq > 1) %>% # keep terms only represented more than once
  mutate(angle = 90 * sample(c(0, 1), n(), replace = TRUE, prob = c(50, 50))) # give an angle to some terms which will be rotated in the wordcloud
View(data)

#Plot
ggplot(data, aes(label = word, size = freq, angle = angle, color = direction)) +
  geom_text_wordcloud(area_corr = TRUE) +
  scale_color_discrete(direction = -1)+ 
  scale_size_area(max_size = 10) +
  theme_minimal()

###########################################################


#### PLANULAE
planulae = read.csv("planulaeTerms.csv")
pUp = filter(planulae, direction == 'upregulated')
View(pUp)
pDown = filter(planulae, direction == 'downregulated')
View(pDown)

tb = table(pUp$term_name)
tb1 = data.frame(word = names(tb), freq =as.numeric(tb))
tb1 = as.data.frame(tb1)
tb1$direction <- "green"
View(tb1)

tb2 = table(pDown$term_name)
tb2 = data.frame(word = names(tb2), freq =as.numeric(tb2))
tb2 = as.data.frame(tb2)
tb2$direction <- "red"
View(tb2)

tb3 = rbind(tb1, tb2) # bind together up and down regulated
View(tb3)
str(tb3)

set.seed(1)# for reproducibility 
data = tb3


data <- data %>%
  filter(freq > 4) %>% # keep terms only represented more than once
  mutate(angle = 90 * sample(c(0, 1), n(), replace = TRUE, prob = c(20, 80))) # give an angle to some terms which will be rotated in the wordcloud
View(data)

#Plot
ggplot(data, aes(label = word, size = freq, angle = angle, color = direction)) +
  geom_text_wordcloud(area_corr = TRUE) +
  scale_color_discrete(direction = -1)+ 
  scale_size_area(max_size = 25) +
  theme_minimal()




##############################################################
#### 8 days
d8 = read.csv("8dayTerms.csv") # no upregulated terms - no need to subset

tb2 = table(d8$term_name)
tb2 = data.frame(word = names(tb2), freq =as.numeric(tb2))
tb2 = as.data.frame(tb2)
tb2$direction <- "red"
View(tb2)

set.seed(1)# for reproducibility 
data = tb2

data <- data %>%
 # filter(freq > 4) %>% # keep terms only represented more than once
  mutate(angle = 90 * sample(c(0, 1), n(), replace = TRUE, prob = c(50, 50))) # give an angle to some terms which will be rotated in the wordcloud
View(data)

nrow(data)

#Plot
ggplot(data, aes(label = word, size = freq, 
                 #angle = angle, # dont rotate with only five words - looks odd
                 color = direction)) +
  geom_text_wordcloud(area_corr = TRUE) +
  scale_color_discrete(direction = -1)+ 
  scale_size_area(max_size = 10) +
  theme_minimal()


#### 60 days
d60 = read.csv("60dayTerms.csv")
d60Up = filter(d60, direction == 'upregulated')
View(d60Up)
d60Down = filter(d60, direction == 'downregulated')
View(adultDown)

tb = table(d60Up$term_name)
tb1 = data.frame(word = names(tb), freq =as.numeric(tb))
tb1 = as.data.frame(tb1)
tb1$direction <- "green"
View(tb1)

tb2 = table(d60Down$term_name)
tb2 = data.frame(word = names(tb2), freq =as.numeric(tb2))
tb2 = as.data.frame(tb2)
tb2$direction <- "red"
View(tb2)

tb3 = rbind(tb1, tb2) # bind together up and down regulated
View(tb3)
str(tb3)

set.seed(1)# for reproducibility 
data = tb3


data <- data %>%
  filter(freq > 4) %>% # keep terms only represented more than once
  mutate(angle = 90 * sample(c(0, 1), n(), replace = TRUE, prob = c(80, 20))) # give an angle to some terms which will be rotated in the wordcloud
View(data)

nrow(data)

#Plot
ggplot(data, aes(label = word, size = freq, angle = angle, color = direction)) +
  geom_text_wordcloud(area_corr = TRUE) +
  scale_color_discrete(direction = -1)+ 
  scale_size_area(max_size = 10) +
  theme_minimal()


# Save/ Print ggplot clouds from the Plot viewer