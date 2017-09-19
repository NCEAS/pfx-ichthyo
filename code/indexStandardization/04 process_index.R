# Start by loading necessary packages  

# Try readxl package
if("readxl" %in% rownames(installed.packages()) == FALSE) {
  devtools::install_github("hadley/readxl")}
require(readxl)
library(MARSS)
library(knitr)
library(PBSmapping)
library(INLA)
library(dplyr)
library(vegan)

#Load in the data -- this is the larvae only data that Janet's cleaned up. 
ichthyo = read.csv("raw_data/mergedData_final_binned.csv")
ichthyo = ichthyo[-which(is.na(ichthyo$Latitude+ichthyo$Longitude)==T),]

# lump species
ichthyo$IchName = as.character(ichthyo$IchName)
ichthyo$IchName[ichthyo$IchName%in%c("Bathymasteridae")] = "Bathymaster spp."
ichthyo$IchName[ichthyo$IchName%in%c("Hexagrammos decagrammus","Hexagrammos lagocephalus")] = "Hexagrammos spp."
ichthyo$IchName[ichthyo$IchName%in%c("Myctophidae")] = "Myctophidae spp."
ichthyo$IchName[ichthyo$IchName%in%c("Myoxocephalus scorpius")] = "Myoxocephalus spp."
ichthyo$IchName[ichthyo$IchName%in%c("Gadus chalcogrammus")] = "Theragra chalcogramma"

#name updates march 20 2017
ichthyo$IchName[ichthyo$IchName%in%c("Theragra chalcogramma")]="Gadus chalcogrammus"
ichthyo$IchName[ichthyo$IchName%in%c("Ammodytes hexapterus")]="Ammodytes personatus"

# For species that were binned, record the total density for that group (instead of multiple observations)
ich_new = group_by(ichthyo, IchName, Latitude, Longitude, Year, Haul) %>%
  summarize(CatchPer10m_2 = sum(CatchPer10m_2, na.rm=T)) %>% 
  as.data.frame()
ichthyo = ich_new

#Summarize the most commonly occuring species
t = sort(table(ichthyo[which(ichthyo$CatchPer10m_2 > 0),]$IchName))
t = rev(t) # reverse sort
t = t[1:40] # only use top 25 spp
t = t/sum(t) # normalize
print(cumsum(t))

#Start with analyzing the 25 most common species, representing 90% of the total occurrences...
occAll = matrix(NA, length(seq(1981,2013)), length(t)) # stores log density estimates for DFA
#names(dfAll) = names(cumsum(rev(t))[1:25])
colnames(occAll) = names(cumsum(rev(t))[1:length(t)])

pdf("Ichthyo_occurrence.pdf")
for(ii in seq(1,length(t))) {
  this.spp = names(cumsum(t))[ii]
  
  # Calculate raw occurrence in data
  subset = ichthyo[which(ichthyo$IchName==this.spp),]
  df = data.frame("Y"=subset$Latitude,"X"=subset$Longitude,"Z"=subset$CatchPer10m_2,"T"=subset$Year)

  # Convert date to integer, 1-31
  df$T = df$T - min(df$T) + 1
  
  #Convert to UTM
  attr(df, "projection") <-"LL"
  attr(df, "zone") <-5
  dfUTM = convUL(df)
  
  #Remove stations from N and E side of Kodiak -- tend to be sparsely sampled in time. Focus on core 90% -
  indx = which(dfUTM$X > quantile(dfUTM$X,0.05) & dfUTM$X < quantile(dfUTM$X,0.95) & dfUTM$Y > quantile(dfUTM$Y,0.05) & dfUTM$Y < quantile(dfUTM$Y,0.95))
  dfUTM = dfUTM[indx,]
  dfUTM$Zbin = dfUTM$Z
  dfUTM$Zbin[which(dfUTM$Zbin>0)]=1
  rawOcc = aggregate(dfUTM$Zbin,by=list(dfUTM$T),mean)
  
  load(paste("binomial/",this.spp,"_binomial.Rdata",sep=""))
  projectedLatentGrid = plogis(projectedLatentGrid)
  print(dim(projectedLatentGrid)[3])
  # Calculate average occurrence, by year. Use MCMC samples to represent uncertainty
  meanOccurrence = apply(projectedLatentGrid,c(2,3),mean)
  occurrenceStats = (apply(meanOccurrence,2,quantile,c(0.025,0.5,0.975)))
  occAll[ ,ii] = occurrenceStats[2,]
  plot(seq(1981,2013), occurrenceStats[2,],type="b",lwd=2,ylab="Median occurrence (black=model,red=raw)",ylim=range(occurrenceStats),xlab="",main=this.spp)
  lines(seq(1981,2013),occurrenceStats[1,],col="grey30",lty=3)
  lines(seq(1981,2013),occurrenceStats[3,],col="grey30",lty=3)
  points(rawOcc$Group.1+1980, rawOcc$x, col="red",lwd=2)
  
  os = t(occurrenceStats)
  os = cbind(os,apply(meanOccurrence,2,mean))
  colnames(os)=c("lower95","median","upper95","mean")
  write.table(os,file=paste("binomial/",this.spp,"_occurrenceStats.csv",sep=""),sep=",",row.names=F)
  
}
dev.off()


write.csv(occAll, "occurrenceStatsForDFA.csv", row.names = T)
occAll = read.csv("occurrenceStatsForDFA.csv")[,-1]
speciesRichness = apply(occAll,1,sum,na.rm=T)
SRout=data.frame(year=seq(1981,2013), SR=speciesRichness)
write.csv(SRout,"OccSR.csv")

# Derive some stats related to richness and diversity
pdf("Expected species richness (core area).pdf")
plot(1981:2013, speciesRichness, xlab="", ylab="Expected species richness", lwd=2, type="b")
dev.off()

#pdf("Shannon diversity (core area).pdf")
#plot(1981:2013, -apply(occAll * log(occAll), 1, sum,na.rm=T), xlab="", ylab="Expected species", lwd=2, type="b")
#dev.off()

#pdf("Simpson diversity (core area).pdf")
#plot(1981:2013, 1/apply(occAll^2, 1, sum, na.rm=T), xlab="", ylab="Expected species", lwd=2, type="b")
#dev.off()

#occAll=read.csv("occurrenceStatsForDFA.csv", header=T)
IchSR=apply(occAll,1,sum,na.rm=T)


# Calculate within year variation by species
pdf("Ichthyo_variationInSpatialOccurrence.pdf")
for(ii in seq(1,length(t))) {
  this.spp = names(cumsum(t))[ii]
  
  load(paste("binomial/",this.spp,"_binomial.Rdata",sep=""))
  projectedLatentGrid = plogis(projectedLatentGrid)
  
  cv = apply(projectedLatentGrid,c(2,3),sd)/apply(projectedLatentGrid,c(2,3),mean)
  plot(1981:2013,apply(cv,2,mean),type="b",lwd=2, xlab="",ylab="CV of within-year spatial variation",main=this.spp)
}
dev.off()

#Start with analyzing the 25 most common species, representing 90% of the total occurrences...
dfAll = matrix(NA, length(seq(1981,2013)), length(t)) # stores log density estimates for DFA
#names(dfAll) = names(cumsum(rev(t))[1:25])
colnames(dfAll) = names(cumsum(t)[1:length(t)])

logdaycoefs = 0

pdf("Ichthyo_total density.pdf")
for(ii in seq(1,length(t))) {
  
  this.spp = names(cumsum(t))[ii]
  
  # Calculate raw occurrence in data
  subset = ichthyo[which(ichthyo$IchName==this.spp),]
  df = data.frame("Y"=subset$Latitude,"X"=subset$Longitude,"Z"=subset$CatchPer10m_2,"T"=subset$Year)
  
  # Convert date to integer, 1-31
  df$T = df$T - 1981 + 1
  
  #Convert to UTM
  attr(df, "projection") <-"LL"
  attr(df, "zone") <-5
  dfUTM = convUL(df)
  
  #Remove stations from N and E side of Kodiak -- tend to be sparsely sampled in time. Focus on core 90% -
  indx = which(dfUTM$X > quantile(dfUTM$X,0.05) & dfUTM$X < quantile(dfUTM$X,0.95) & dfUTM$Y > quantile(dfUTM$Y,0.05) & dfUTM$Y < quantile(dfUTM$Y,0.95))
  dfUTM = dfUTM[indx,]
  rawDens = aggregate(dfUTM$Z,by=list(dfUTM$T),mean)
  
  load(paste("binomial/",this.spp,"_binomial.Rdata",sep=""))
  projectedLatentGrid.bin = plogis(projectedLatentGrid)
  
  load(paste("positive/",this.spp,"_positive2.Rdata",sep=""))
  #logdaycoefs[ii] = (summary(inlaModel)$fixed[1,1])
  projectedLatentGrid.pos = exp(projectedLatentGrid)
  
  projectedLatentGrid.tot = projectedLatentGrid.bin * projectedLatentGrid.pos
  
  # Calculate average occurrence, by year. Use MCMC samples to represent uncertainty
  meanDens = apply(projectedLatentGrid.tot,c(2,3),mean)
  densityStats = (apply(meanDens,2,quantile,c(0.025,0.5,0.975), na.rm=T))
  dfAll[,ii] = log(densityStats[2,])
  logmeanDensity = log(densityStats)
  
  # This is for average density
  plot(seq(1981,2013), logmeanDensity[2,],type="b",lwd=2,ylab="Mean density (black=model,red=raw)",ylim=range(logmeanDensity[,-c(4,6,32)], na.rm=T),xlab="",main=paste(this.spp,"(log)"))
  lines(seq(1981,2013),logmeanDensity[1,],col="grey30",lty=3)
  lines(seq(1981,2013),logmeanDensity[3,],col="grey30",lty=3)
  points(rawDens$Group.1+1980, log(rawDens$x), col="red",lwd=2)
  
  # This is for total density
  plot(seq(1981,2013), densityStats[2,],type="b",lwd=2,ylab="Mean density (black=model,red=raw)",ylim=c(0, quantile(densityStats[,-c(4,6,32)],0.9,na.rm=T)),xlab="",main=paste(this.spp))
  lines(seq(1981,2013),densityStats[1,],col="grey30",lty=3)
  lines(seq(1981,2013),densityStats[3,],col="grey30",lty=3)
  points(rawDens$Group.1+1980, rawDens$x, col="red",lwd=2)
  
  os = t(densityStats)
  os = cbind(os,apply(meanDens,2,mean))
  colnames(os)=c("lower95","median","upper95","mean")
  write.table(os,file=paste("positive/",this.spp,"_densityStats.csv",sep=""),sep=",")
  
}
dev.off()

write.csv(exp(dfAll), "DensityForDFA.csv",row.names=F)
write.csv(dfAll, "logDensityForDFA.csv",row.names=F)


dfAll=read.csv("DensityForDFA.csv",header=T)
IchSW=diversity(dfAll, index="shannon")

pdf("Standardized density estimates.pdf")
matplot(1981:2013,apply(exp(dfAll),2,scale),type="l",ylab="Standardized density")
dev.off()

pdf("Standardized log-density estimates.pdf")
matplot(1981:2013,apply(dfAll,2,scale),type="l",ylab="Standardized log(density)")
dev.off()

