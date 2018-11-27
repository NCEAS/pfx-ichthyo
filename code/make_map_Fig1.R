
library(knitr)
library(PBSmapping)
library(INLA)
library(dplyr)
library(ggplot2)

setwd('/Users/kristin.marshall/downloads/')
#Load in the data -- this is the larvae only data that Janet's cleaned up. 
ichthyo = read.csv("mergedData_final_binned.csv")
ichthyo = ichthyo[-which(is.na(ichthyo$Latitude+ichthyo$Longitude)==T),]

# lump species
ichthyo$IchName = as.character(ichthyo$IchName)
ichthyo$IchName[ichthyo$IchName%in%c("Bathymasteridae")] = "Bathymaster spp."
ichthyo$IchName[ichthyo$IchName%in%c("Hexagrammos decagrammus","Hexagrammos lagocephalus")] = "Hexagrammidate spp."
ichthyo$IchName[ichthyo$IchName%in%c("Myctophidae")] = "Myctophidae spp."
ichthyo$IchName[ichthyo$IchName%in%c("Myoxocephalus scorpius")] = "Myoxocephalus spp."
ichthyo$IchName[ichthyo$IchName%in%c("Gadus chalcogrammus")] = "Theragra chalcogramma"

# For species that were binned, record the total density for that group (instead of multiple observations)
ich_new = group_by(ichthyo, IchName, Latitude, Longitude, Year, Haul) %>%
  summarize(CatchPer10m_2 = sum(CatchPer10m_2, na.rm=T)) %>% 
  as.data.frame()
ichthyo = ich_new

this.spp = "Theragra chalcogramma"
  
  subset = ichthyo[which(ichthyo$IchName==this.spp),]
  df = data.frame("Y"=subset$Latitude,"X"=subset$Longitude,"Z"=subset$CatchPer10m_2,"T"=subset$Year)
  # Extract year from the date variable

  # Convert date to integer, 1-31
  df$T = df$T - 1981 + 1
  
  #Convert to UTM
  attr(df, "projection") <-"LL"
  attr(df, "zone") <-7
  dfUTM = convUL(df)

  #Remove stations from N and E side of Kodiak -- tend to be sparsely sampled in time. Focus on core 90% -
  indx = which(dfUTM$X > quantile(dfUTM$X,0.05) & dfUTM$X < quantile(dfUTM$X,0.95) & dfUTM$Y > quantile(dfUTM$Y,0.05) & dfUTM$Y < quantile(dfUTM$Y,0.95))
  dfUTM$incl = "Include"
  dfUTM$incl[indx] = "Not Include"

  gridproj = expand.grid("Y" = seq(6100,6500,10), "X" = seq(10,570,10))

  data(nepacLLhigh)
  nepacLLhigh_utm = convUL(nepacLLhigh)

  indx = which(df$X > quantile(df$X,0.05) & df$X < quantile(df$X,0.95) & df$Y > quantile(df$Y,0.05) & df$Y < quantile(df$Y,0.95))
  df$incl = "Include"
  df$incl[indx] = "Not Include"

g = ggplot(data=dfUTM, aes(X, Y)) + coord_cartesian(ylim = c(6000,6900), xlim=c(min(dfUTM$X),max(dfUTM$X))) + geom_point(aes(colour = factor(incl)), size = 0.3) +  
 scale_color_manual(values=c("grey70","dodgerblue")) + 
 geom_polygon(data = nepacLLhigh_utm, aes(X, Y, group = PID), fill=grey(0.3)) + xlab("Longitude") + ylab("Latitude") + theme_minimal() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position="none") + xlab("Longitude") + ylab("Latitude")
print(g)


g = ggplot(data=df, aes(X, Y)) + coord_cartesian(ylim = c(52, 60), xlim=c(min(df$X),max(df$X))) + geom_point(aes(colour = factor(incl)), size = 0.2) +  
 scale_color_manual(values=c("grey80","#5e3c99")) + 
 geom_polygon(data = nepacLLhigh, aes(X, Y, group = PID), fill=grey(0.4)) + xlab("Longitude") + ylab("Latitude") + theme_minimal() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position="none") + xlab("Longitude") + ylab("Latitude")
print(g)

pdf("map.pdf", width=90/25.4, height=90/25.4)
print(g)
dev.off()

