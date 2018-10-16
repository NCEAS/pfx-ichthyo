library(dplyr)
library(tidyr)

#read in data from seldovia- raw data, not de-trended
seld.mml=read.csv('data/CO-OPS__9455500__ml.csv')

seld.filter<- seld.mml %>%
    select(Year, Month, MSL) %>%
    filter(Year<2014 & Year>1980 & Month<7) %>%
    group_by(Year)%>%
    mutate(avg=mean(MSL)) %>%
    ungroup() %>%
    spread(Month, MSL) 

matplot(seld.filter[,-1])

#read in responses
d <- readRDS("correlation-matrix.rds")


cor(cbind(d,seld.filter), use="complete.obs")
#covariates
reshape()
for (i in 1:6) {
  seld.month= seld.filter %>%
    filter(Month==i)
  cor()
}