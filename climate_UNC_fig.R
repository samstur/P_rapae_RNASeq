

## want to look at the climate experienced near Chapel Hill where the population was collected 
# read in list of site positions in lat and long


install.packages("daymetr")
library("daymetr")


## pull down climate data from Perry-Winkle Farm in Chapel Hill where the caterpillars were collected
data1<-download_daymet(lat = 35.866688,
                       lon = -79.203601,
                       start = 1985,
                       end = 2024,
                       internal = TRUE)
temp_data<-data.frame(data1$data)
head(temp_data)
temp_data<- temp_data[c(1,2,7,8)]
colnames(temp_data)[3] ="tmax"
colnames(temp_data)[4] ="tmin"
View(temp_data)

days<-unique(temp_data$yday)
days
df<-data.frame(yday="fake",tmin_mean=NA,tmin_min=NA,tmax_mean=NA,tmax_max=NA)
for (day in days){
  x <- day
  tmin_min <- min(subset(temp_data,temp_data$yday==x)$tmin, na.rm = TRUE)
  tmin_mean <- mean(subset(temp_data,temp_data$yday==x)$tmin, na.rm = TRUE)
  tmax_max <- max(subset(temp_data,temp_data$yday==x)$tmax, na.rm = TRUE)
  tmax_mean <- mean(subset(temp_data,temp_data$yday==x)$tmax, na.rm = TRUE)
  newdata<-data.frame(x,tmin_mean,tmin_min,tmax_mean,tmax_max)
  names(newdata)<- c("yday","tmin_mean","tmin_min","tmax_mean","tmax_max")
  df <- rbind(df, newdata)
  }
View(df)
# remove first dummy row
df<-df[-1,]

breakfunc <- function(x) {
  origin <- as.Date("2021-01-01")
  days <- origin + x
  origin <- as.POSIXlt(origin)
  dayseq <- as.POSIXlt(seq(days[1], days[2], by = "day"))
  with(dayseq, yday[mday == 1] + 365*(year[mday == 1] - origin$year ))
}

labelfunc <- function(x) {
  origin <- as.Date("2021-01-01")
  format(origin + x, format = "%b")
}

## plot
library(ggplot2)
plot<-ggplot(data=df,aes(x=as.numeric(df$yday),y=df$tmin_min))+
  geom_line(linetype = "dashed",size=0.5, colour = "darkblue")+
  #scale_color_manual(limits = c("South","North"),
  #                   labels = c("Southern Temperate", "Northern Temperate"),
  #                   values = c("dark red", "dark blue"))+
  geom_line(data=df,aes(x=as.numeric(df$yday),y=df$tmin_mean),
            color="darkblue",
            size = 1,
            linetype="solid")+
  geom_line(data=df,aes(x=as.numeric(df$yday),y=df$tmax_mean),
            color="darkred",
            size = 1,
            linetype="solid")+
  geom_line(data=df,aes(x=as.numeric(df$yday),y=df$tmax_max),
            color="darkred",
            size = 0.5,
            linetype="dashed")+
  theme_bw()+
  xlab("Date")+
  ylab("Temperature (°C)")+
  theme(axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 15),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size=15))+
  scale_y_continuous(breaks=c(-20,-10,0,10,20,30,40))+
  scale_x_continuous(breaks=breakfunc, labels=labelfunc)
plot

setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/INTBio grant work/RNA-Seq workflow")
ggsave("climate_plot_UNC.png",plot=plot,
       dpi=600,units='in',width=8,height=5)

## convert yday to month and day columns so that I can get monthly average tmin and tmax

# change doy to dates using mutate()
df <- df %>% mutate(date_= as.Date(as.numeric(yday)-1, origin="25-01-01"), 
             month= strftime(date_, "%m"), 
             day=strftime(date_,"%d")) 
monthlyaverages_min <- df %>%
  group_by(month) %>%
  dplyr::summarize(tmin_average = mean(tmin_mean, na.rm=TRUE))

monthlyaverages_max <- df %>%
  group_by(month) %>%
  dplyr::summarize(tmax_average = mean(tmax_mean, na.rm=TRUE))

monthlyaverages_min_record <- df %>%
  group_by(month) %>%
  dplyr::summarize(tmin_record = mean(tmin_min, na.rm=TRUE))

monthlyaverages_max_record <- df %>%
  group_by(month) %>%
  dplyr::summarize(tmax_record = mean(tmax_max, na.rm=TRUE))

monthlyaverages <- merge(monthlyaverages_min,
                         monthlyaverages_max,by="month")
monthlyaverages <- merge(monthlyaverages, 
                         monthlyaverages_min_record,by="month")
monthlyaverages <- merge(monthlyaverages, 
                         monthlyaverages_max_record,by="month")               

## save this climate data averaged over months  
setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/INTBio grant work/Writing/Final Scripts/Script_outputfiles/")
write.csv(monthlyaverages,"UNC_climate_data.csv",row.names = F)

