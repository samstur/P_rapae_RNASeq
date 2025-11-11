

setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/INTBio grant work/Writing/Final Scripts/Script_outputfiles/Prapae_Mortality_DevTime_Plot_outputfiles")
data <- read.csv(data_final,"lifehistory_data_final.csv",row.names = F)

library("lme4")
library("car")
library(predictmeans)
library(gridExtra)
library(ggplot2)

data$rearing_temp <- substr(data$Treatment, 1, 2) 
data$testing_temp <- substr(data$Treatment, 4, 5)
data$SIB <- as.factor(data$SIB)

## statistical analysis of dev time in caterpillar which were tested for the virus
data$days_to_fourth.log <- log(data$days_to_fourth)
data$days_to_fourth.sqrt <- sqrt(data$days_to_fourth)
model1 <- lmer(days_to_fourth.sqrt ~ rearing_temp+inf_status + (1|SIB),
               data = data)
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

library(ggplot2)
lb1_dead <- c("DevTemp^ns")
death_plot<-ggplot(data,aes(x=rearing_temp,y=death))+
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
               data = data,family = binomial)
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
data$inf_status_yn <- with(data, ifelse(inf_status == "negative", 0, 1))
data %>% count(rearing_temp, inf_status, sort = TRUE)
modelinf <- glmer(inf_status_yn ~ rearing_temp + (1|SIB),
                   data = data,family = binomial)
summary(modelinf)
Anova(modelinf)
df.residual(modelinf) #109
residplot(modelinf) 

## get counts of sib families across treatments
library(dplyr)
table <- data %>% count(SIB, Treatment)
View(table)
setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/INTBio grant work/Writing/Final Scripts/Script_outputfiles/Prapae_Mortality_DevTime_Plot_outputfiles")
write.csv(table,"SIB_numbers_table.csv",row.names = F)



