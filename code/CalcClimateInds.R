##This file calculates the 6 month climate indices

library(httr)
library(plyr)
library(dplyr)
library(XML)
library(curl)
library(rvest)
library(tidyr)
library(stringr)

## Steps for data cleaning: 
## 1) read in data
## 2) format to annual estimates (2 column dataframe with cols=Year,spEstimate)


#############    
###  Multivariate ENSO Index (MEI):
URL_enso <- "http://www.esrl.noaa.gov/psd/enso/mei/table.html"
enso_pre <- xpathSApply(xmlParse(content(GET(URL_enso))),"/html/body/pre", xmlValue)
enso_cols <- scan(textConnection(enso_pre), skip=10, nlines=1, what=character()) # get header row
enso <- read.csv(file=textConnection(enso_pre), skip=11, stringsAsFactors=F, sep="\t",
                 header=FALSE, col.names=enso_cols)
enso_df <- enso[1:66,]  # removes the text at bottom of file
#
ENSO_spring <- enso_df %>%
  rename(Year=YEAR) %>% # rename data columns
  filter(Year %in% c(1981:2013)) %>% # selects years 1975 - 2015
  select(Year:JUNJUL)%>%
  gather(Months, ENSO, -Year) %>% # reshapes data to be column-wise
  filter(!is.na(ENSO)) %>% # remove NA values
  group_by(Year) %>%
  summarise(ENSO_anul_mn=mean(ENSO)) %>% # get annual means
  ungroup()  #

###NPGO
# North Pacific Gyre Oscillation Index (NPGO):
URL_npgo <- "http://www.o3d.org/npgo/npgo.php"
npgo_pre <- xpathSApply(xmlParse(content(GET(URL_npgo))),"/html/body/pre", xmlValue)
npgo_cols <- scan(textConnection(npgo_pre), skip=25, nlines=1, what=character())# Get header row
npgo_cols <- npgo_cols[2:4] # select column names
npgo_df <- read.csv(file=textConnection(npgo_pre), skip=26, stringsAsFactors=F, sep="",
                    header=FALSE, col.names=npgo_cols, strip.white=TRUE)

npgo_spring <- npgo_df %>%
  rename(Year=YEAR) %>% # rename data columns
  filter(Year %in% c(1981:2013) & MONTH %in% c(1:6)) %>% # selects years 1975 - 2015
  group_by(Year) %>%
  summarise(NPGO_anul_mn=mean(NPGO)) %>% # get annual means
  ungroup()  

###### North Pacific Index (NPI) for sea level pressure:
URL_npi <- "https://climatedataguide.ucar.edu/sites/default/files/climate_index_files/npindex_monthly.ascii"
npiGet <- GET(URL_npi)
npi1 <- content(npiGet, as='text')
npi <- read.table(file=textConnection(npi1),stringsAsFactors=F, sep=" ", header=TRUE, fill=TRUE)
npi[1:50,]
#
NPI <- npi %>%
  rename(YearMon=X, SeaLevelPressure_hPa=and) %>% # rename columns with data
  select(YearMon, SeaLevelPressure_hPa) %>% # remove columns without data
  mutate(Year=substring(YearMon,1,4),   # creates Year column
         Month=substring(YearMon,5,6)) %>%  # creates Month column
  filter(Year %in% c(1981:2013)) %>% # selects years 1975 - 2015
  filter(!is.na(SeaLevelPressure_hPa),
         SeaLevelPressure_hPa != -999.00) %>% # remove NA values, and -999 values which are NAs
  filter(Month %in% c("01","02","03","04","05","06")) %>%
  group_by(Year) %>%
  summarise(SeaLevelPressure_Spring_hPa=mean(SeaLevelPressure_hPa)) %>%  # winter means
  #summarise(SeaLevelPressure_mean_hPa=mean(SeaLevelPressure_hPa)) %>% # get annual means
  ungroup()
#

#############
###  Pacific Decadal Oscillation Index (PDO):
URL_pdo <- "http://jisao.washington.edu/pdo/PDO.latest"
pdo_raw <- read_html(URL_pdo)
pdo_pre <- pdo_raw %>%
  html_node("p") %>%
  html_text()
pdo_cols <- scan(textConnection(pdo_pre), skip=31, nlines=1, what=character())# Get header row
pdo_df <- read.table(file=textConnection(pdo_pre), skip=32, nrows=117, stringsAsFactors=F, sep="",
                     header=FALSE, col.names=pdo_cols, strip.white=TRUE, fill=TRUE)
pdo_df$YEAR <- substr(pdo_df$YEAR, 1, 4)  # removes asterisks from years 2002-2015

pdo_spr <- pdo_df %>%
  rename(Year=YEAR) %>% # rename data columns
  filter(Year %in% c(1981:2013)) %>% # selects years 1974 - 2015
  select(Year, JAN:JUN) %>% # select only Dec, Jan, Feb data
  gather(Month, PDO, -Year) %>% # reshapes data to be column-wise
  mutate(Year = as.integer(Year)) %>%
  group_by(Year) %>%
  summarise(PDO_spring = mean(as.numeric(as.character(PDO)), na.rm = TRUE)) %>% # get winter means
  ungroup()


###Upwelling
up.anom=read.table('data/MonthlyUpwellingIndexAnom.txt', header=T) ##reads in monthly upwelling index, 
##data from http://www.pfeg.noaa.gov/products/PFELData/upwell/monthly/upanoms.mon

up.anom$spr.avg=rowMeans(up.anom[,4:9]) ## Jan to June values only

up.spr=up.anom$spr.avg[up.anom$YEAR>1980 & up.anom$YEAR<2014]
#plot(seq(1981,2013),up.anom$spr.avg[up.anom$YEAR>1980], type='o')


upwell=read.csv('data/ClimateCovars.csv', header=T)
#points(seq(1981,2013),upwell$Upwelling[upwell$YEAR>1974], col='red', type='o')

ccspr=data.frame(cbind(Year=seq(1981,2013), ENSO=ENSO_spring$ENSO_anul_mn, NPGO=npgo_spring$NPGO_anul_mn, 
                 NPI=NPI$SeaLevelPressure_Spring_hPa, PDO=pdo_spr$PDO_spring, Upwelling=up.spr))

write.csv(ccspr,'data/ClimateCovarsSpr.csv')
