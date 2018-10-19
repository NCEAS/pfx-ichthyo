library(tidyverse)


#read in data from seldovia- raw data, not de-trended
seld.mml=read.csv('data/CO-OPS__9455500__ml.csv')

#read in data from sandpoint
sand.mml=read.csv('data/CO-OPS__9459450__ml.csv')

seld.filter<- seld.mml %>%
    select(Year, Month, MSL) %>%
    filter(Year<2014 & Year>1980 & Month<7) %>%
    group_by(Year)%>%
    mutate(avg=mean(MSL)) %>%
    ungroup() %>%
    spread(Month, MSL) 

matplot(seld.filter[,-1])
names(seld.filter)=c('year','avg', 'jan.MSL', 'feb.MSL', 'mar.MSL','apr.MSL', 'may.MSL','jun.MSL')

##pull the sandpoint station tide gauge data
sand.filter<- sand.mml %>%
  select(Year, Month, MSL) %>%
  filter(Year<2014 & Year>1980 & Month<7) %>%
  group_by(Year)%>%
  mutate(avg=mean(MSL)) %>%
  ungroup() %>%
  spread(Month, MSL) 

names(sand.filter)=c('year','avg', 'jan.MSL', 'feb.MSL', 'mar.MSL','apr.MSL', 'may.MSL','jun.MSL')

#read in detrended data with seasonal cycle removed
seld.anom<-read.csv('data/9455500_intannvar.csv', head=T)

seld.anom.filter<- seld.anom %>%
  select(Year, Month, Interannual_Variation) %>%
  filter(Year<2014 & Year>1980 & Month<7) %>%
  group_by(Year)%>%
  mutate(avg.anom=mean(Interannual_Variation)) %>%
  ungroup() %>%
  spread(Month, Interannual_Variation) 
names(seld.anom.filter)=c('Year','avg.anom', 'jan.anom', 'feb.anom', 'mar.anom','apr.anom', 'may.anom','jun.anom')

#pull in sandpoint anomaly data
sand.anom<-read.csv('data/9459450_intannvar.csv', head=T)
sand.anom.filter<- sand.anom %>%
  select(Year, Month, Interannual_Variation) %>%
  filter(Year<2014 & Year>1980 & Month<7) %>%
  group_by(Year)%>%
  mutate(avg.anom=mean(Interannual_Variation)) %>%
  ungroup() %>%
  spread(Month, Interannual_Variation) 
names(sand.anom.filter)=c('Year','avg.anom', 'jan.anom', 'feb.anom', 'mar.anom','apr.anom', 'may.anom','jun.anom')






#read in responses
d <- readRDS("correlation-matrix.rds")
#d=cbind(d,seld.filter[,2:8])
#d=cbind(d,sand.anom.filter[,2:8])
d=cbind(d,sand.filter[,2:8])

#predictors <- select(d, c(ENSO:sst, avg, V2:V5))
responses <- select(d, trend1:IchSW)

predictors<-select(d, avg:jun.MSL)

#names(predictors)=c('avg', 'jan.MSL', 'feb.MSL', 'mar.MSL','apr.MSL', 'may.MSL','jun.MSL')


names(predictors)=c('ENSO','NPGO','NPI','PDO','Upwell','SST','MSL','POLL','PCOD','POP','ARROW')

names(responses)=c('Trend 1','Trend 2', 'Sp Rich','Shannon')

cor(cbind(d,seld.filter), use="complete.obs")
#covariates
reshape()
for (i in 1:6) {
  seld.month= seld.filter %>%
    filter(Month==i)
  cor()
}