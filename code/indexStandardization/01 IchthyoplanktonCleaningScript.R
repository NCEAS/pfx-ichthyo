#######################################################################
#####  EcoFOCI Shelikof Strait Ichthyoplankton Cleaning Script    #####
#####        Script by Colette Ward (ward at nceas.ucsb.edu)      #####
#####             and                                             #####
#######################################################################

## load packages
library(dplyr)
library(httr)

#Load data (use the following lines starting with # once the dataset is posted on AOOS)
#URL_Ich <- "url here"
#IchGet <- GET(URL_Ich)
#IchA <- content(IchGet, as='text')
#Ich <- read.csv(file=textConnection(IchA), stringsAsFactors=F, na.strings = c("NA", " ", ""))

# for now, file is in PFx shared folder -> Ichthyoplankton on Google Drive
setwd("/users/kristinmarshall/dropbox/Ichtyo_IndexStandardization/") # set to whatever working directory you're using

Ich <- read.csv('IchTime Series_for binning.csv', header=T, stringsAsFactors=F, na.strings = c("NA", " ", ""))
head(Ich)
str(Ich)

# basic cleaning & processing:
Ich1 = Ich %>%
  filter(GeographicArea == "GOA") %>%
  mutate(gmtDate=as.Date(GMTDate, "%d-%B-%y")) %>%   # split GMTDate into day, month, year
  mutate(Year=strsplit(as.character(gmtDate),split="-") %>% #
           sapply(function(x) x[1])) %>%
  mutate(year=as.numeric(Year)) %>%
  mutate(Month=strsplit(as.character(gmtDate),split="-") %>%
           sapply(function(x) x[2])) %>%
  mutate(month=as.numeric(Month)) %>%
  mutate(Day=strsplit(as.character(gmtDate),split="-") %>%
           sapply(function(x) x[3])) %>%
  mutate(day=as.numeric(Day)) %>%
  select(-DispVol, -CatchPer1000m3, -IchStageName, -IchCode, -IchStageCode, -NumberCaught, -NumberMeasuredOrStaged, 
         -GMTDate, -gmtDate, -Year, -Month, -Day, -LatitudeHemisphere, -LongitudeHemisphere)
head(Ich1)



# visualize sample timing
IchMonths=Ich1 %>%
  mutate(myp=paste(month,year,Purpose)) %>%
  filter(!duplicated(myp)) %>%
  select(month,year,Purpose)
plot(IchMonths$month ~ IchMonths$year, pch=16) # consider removing March & April samples in 1985


# Binning:
unique(sort(Ich1$IchName)) 
# Janet advises that some taxa should be binned prior to analyses.
# Her advice for binning is in "Final IchName list for binning_Sept 28 2015.xlsx" in the shared Ichthyoplankton folder on Google Drive

# I'm not sure how you want the spreadsheet to be organized for analyses so it's probably best you take over the scripting from here
# I'd been thinking about doing something like this:
# create new empty BinnedTaxa column in Ich1
# then create binning instructions; example for Anoplarchus spp.:
Ich1$IchName = as.character(Ich1$IchName)
Ich1$IchName[which(Ich1$IchName %in% c("Anoplarchus insignis", "Anoplarchus purpurescens", "Anoplarchus spp."))] = "Anoplarchus spp."
Ich1$IchName[which(Ich1$IchName %in% c("Artedius fenestralis"," Artedius harringtoni","Artedius spp."))] = "Artedius spp."
Ich1$IchName[which(Ich1$IchName %in% c("Atherestes spp.","Atherestes stomias"))] = "Atherestes spp."
Ich1$IchName[which(Ich1$IchName %in% c("Bathymaster spp.","Bathymasteridae","Ronquilus jordani"))] = "Bathymaster spp."
Ich1$IchName[which(Ich1$IchName %in% c("Gymnocanthus spp."))] = "Gymnocanthus spp."
Ich1$IchName[which(Ich1$IchName %in% c("Icelinus borealis","Icelinus spp.","Icelus spp.","Icelinus spp.-Icelus spp."))] = "Icelinus spp."
Ich1$IchName[which(Ich1$IchName %in% c("Leptagonus frenatus","Leptagonus spp."))] = "Leptagonus spp."
Ich1$IchName[which(Ich1$IchName %in% c("Myctophidae","Myctophidae B"))] = "Myctophidae spp."
Ich1$IchName[which(Ich1$IchName %in% c("Myoxocephalus B","Myoxocephalus G","Myoxocephalus polyacanthocephalus"))] = "Myoxocephalus spp."
Ich1$IchName[which(Ich1$IchName %in% c("Nautichthys oculofasciatus","Nautichthys pribilovius","Nautichthys robustus"))] = "Nautichthys spp."
Ich1$IchName[which(Ich1$IchName %in% c("Pholis laeta"))] = "Pholis spp."
Ich1$IchName[which(Ich1$IchName %in% c("Psychrolutes sigalutes"))] = "Psychrolutes spp."
Ich1$IchName[which(Ich1$IchName %in% c("Radulinus asprellus","Radulinus boleoides"))] = "Radulinus spp."
Ich1$IchName[which(Ich1$IchName %in% c("Triglops forficatus","Triglops macellus","Triglops pingeli","Triglops scepticus"))] = "Triglops spp."

# Output data frame to file
write.csv(Ich1, "IchTime Series_binned.csv", row.names=F)

