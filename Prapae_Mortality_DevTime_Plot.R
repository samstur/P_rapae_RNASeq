

setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/INTBio grant work/Writing/Input Files/Input files for Dev Mortality Script")
data <- read.csv(file="cat_viral_records.csv")
head(data)
data$rearing_temp <- substr(data$Treatment, 1, 2) # this takes the first three letters of each sampleName to give the treatment name without the rep number
data$SIB <- as.character(data$SIB)

## subset only the caterpillars who were tested and had a clear result. Takes sample size from 120 to 108
data_tested <- subset(data,data$inf_status == "positive" | data$inf_status == "negative")


library("lme4")
library("car")
library(predictmeans)
library(gridExtra)
library(ggplot2)

View(data_tested)
## statistical analysis of dev time in caterpillar which were tested for the virus
data_tested$days_to_fourth.log <- log(data_tested$days_to_fourth)
data_tested$days_to_fourth.sqrt <- sqrt(data_tested$days_to_fourth)
model1 <- lmer(days_to_fourth.sqrt ~ rearing_temp+inf_status + (1|SIB),
               data = data_tested)
summary(model1)
Anova(model1)
df.residual(model1) #107
residplot(model1) 

library(ggbeeswarm)
lb1 <- c("DevTemp***")
lb2 <- c("Infection^ns")
dev_plot<-ggplot(data,aes(x=rearing_temp,y=days_to_fourth))+
  geom_quasirandom(size=0.5,alpha=0.3)+
  stat_summary(fun.data=mean_se,size=0.5)+
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               fun.args = list(mult = 2), width = 0.15,size=0.6)+
  theme_bw()+
  scale_x_discrete(limits=c("16","30"),
                   labels=c("16±5","30±5"))+
  theme(text = element_text(face="plain"))+
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 11),
        plot.tag = element_text(size= 14,face="bold"))+
  xlab("Developmental Temperature (°C)")+
  annotate("text", x=1.9, y=18, label=lb1, parse=F, size = 3,hjust = 0)+
  annotate("text", x=1.9, y=17.25, label=lb2, parse=T, size = 3,hjust = 0)+
  labs(tag="B")+
  ylab("Development Time (days)")
dev_plot


############# now for death data ######################

setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/INTBio grant work/Writing/Input Files/Input files for Dev Mortality Script")
datadead <- read.csv(file="fullcatdata.csv")
datadead$SIB <- as.character(datadead$SIB)
head(datadead)
## get rid of blank rows
datadead <- datadead[c(1:187),]

## retain only individuals that reached 4th instar or died 
datadead <- subset(datadead,datadead$fate == "4th" | datadead$fate == "pmd" )

# this takes the first three letters of each sampleName 
## to give the treatment name without the rep number
datadead$rearing_temp <- substr(datadead$Treatment, 1, 2) 

## add binomial death variable
datadead$death <- with(datadead, ifelse(fate == "4th", 0, 1))

library(ggplot2)
lb1_dead <- c("DevTemp^ns")
death_plot<-ggplot(datadead,aes(x=rearing_temp,y=death))+
  geom_quasirandom(size=0.5,alpha=0.3)+
  stat_summary(fun.data=mean_se,size=0.6)+
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               fun.args = list(mult = 2), width = 0.15,size=0.5)+
  theme_bw()+
  scale_x_discrete(limits=c("16","30"),
                   labels=c("16±5","30±5"))+
  theme(text = element_text(face="plain"))+
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 11),
        plot.tag = element_text(size=14,face="bold"))+
  xlab("Developmental Temperature (°C)")+
  annotate("text", x=1.9, y=0.95, label=lb1_dead, parse=T, size = 3,hjust = 0)+
  ylab("Incidence of Larval Mortality")+
  labs(tag="A")+
  coord_cartesian(ylim=c(0:1))
death_plot

setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/INTBio grant work/Writing/Final Scripts/Final Figures")
ggsave("devdeath_plot_testedcats_smallpoints.png",plot=grid.arrange(death_plot,dev_plot, 
                                                   nrow=1,ncol=2,widths=c(0.5,0.5),
                                                   heights=c(1)),
       dpi=600,units='in',width=6,height=4.5)

View(datadead)
modeldead <- glmer(death ~ rearing_temp + (1|SIB),
               data = datadead,family = binomial)
summary(modeldead)
Anova(modeldead)
df.residual(modeldead) #167
residplot(modeldead) 

######################### get means and SEs for each group
library(dplyr)
install.packages("plotrix")
library("plotrix")

df_dev <- data %>% 
  group_by(rearing_temp) %>% 
  summarise(mean = mean(days_to_fourth),
            std = 2*(std.error(days_to_fourth)))
df_dev

table(data$rearing_temp)

df_dead <- datadead %>% 
  group_by(rearing_temp) %>% 
  summarise(mean = mean(death),
            std = 2*(std.error(death)))
df_dead

## get totals of infected caterpillars for each treatment and test for a difference in inf across rearing temperatures
data_tested$inf_status_yn <- with(data_tested, ifelse(inf_status == "negative", 0, 1))
data_tested %>% count(rearing_temp, inf_status, sort = TRUE)
modelinf <- glmer(inf_status_yn ~ rearing_temp + (1|SIB),
                   data = data_tested,family = binomial)
summary(modelinf)
Anova(modelinf)
df.residual(modelinf) #109
residplot(modelinf) 


table(datadead$rearing_temp)
########### combine death data and viral infection data for a final data file to upload with paper
head(data_tested)
data_tested$unique_ID <- paste(data_tested$SIB,"-",data_tested$ID)
head(datadead)
datadead$unique_ID <- paste(datadead$SIB,"-",datadead$ID)
View(data_tested)
data_tested_subset <- data_tested[,c(14,19)]

## need to add testing data to all the datadead rows (or unknown if not tested)

data_final <- merge(datadead,data_tested_subset,by="unique_ID",all=T)
data_final$Treatment_new <- ifelse(data_final$Treatment == "30C", '30-30',
                                   ifelse(data_final$Treatment == "30F", '30-16',
                                          ifelse(data_final$Treatment == "16C", '16-16',
                                                 ifelse(data_final$Treatment == "16F", '16-30',NA))))
data_final <- data_final[,c(1,2,18,3,4,5,16,10,14,17)]
View(data_final)
data_final$unique_ID <- as.character(data_final$unique_ID)

setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/INTBio grant work/Writing/Final Scripts/Script_outputfiles/Prapae_Mortality_DevTime_Plot_outputfiles")
write.csv(data_final,"lifehistory_data_final.csv",row.names = F)

data <- read.csv("lifehistory_data_final.csv")
View(data)

library(dplyr)
table <- data %>% count(SIB, Treatment_new)
View(table)
setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/INTBio grant work/Writing/Final Scripts/Script_outputfiles/Prapae_Mortality_DevTime_Plot_outputfiles")
write.csv(table,"SIB_numbers_table.csv",row.names = F)



