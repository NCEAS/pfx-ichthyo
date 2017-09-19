# Start by loading necessary packages  
# Try readxl package
if("readxl" %in% rownames(installed.packages()) == FALSE) {
  devtools::install_github("hadley/readxl")}
require(readxl)
library(knitr)
library(PBSmapping)
library(INLA)
library(date)
library(dplyr)

#Load in the data -- this is the larvae only data that Janet's cleaned up. 
ichthyo = read.csv("raw_data/mergedData_final_binned.csv")
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

#Summarize the most commonly occuring species
t = sort(table(ichthyo[which(ichthyo$CatchPer10m_2 > 0),]$IchName))
t = rev(t) # reverse sort
t = t[1:40] # only use top 40 spp
t = t/sum(t) # normalize
print(cumsum(t))

plot(cumsum(t), xlab="Number of species", ylab="Percent of all records (positive)")

#For each of the top 40 species in the analysis, we restricted the survey spatially to the only including data points within 90% of the range, defined by latitude and longitude (converted to UTM - zone 5). This exclusion removed sampling stations in regions that were not well - sampled; examples include the northeast side of Kodiak, which was not sampled at all in some years. For the remaining species, we separated the data into occurrence (presence-absence) and the conditional positive model (distribution when observed). We used latitude 
#and longitude to define the mesh, and included year as a fixed effect. We used a gamma distribution to model positive density, and for both the presence-absence and positive models, was allowed the spatial random field to be an autoregressive model. Once each of the models was fit, we used MCMC to jointly sample from the distributions of estimated fixed and random effects. We then created a 10x10 km grid within the survey region (bounded 10-570, 6100-6500) and projected the random effects estimated from the knot locations to the grid. We created total density estimates for each grid cell by multiplying probability of occurrence with conditional positive density, and summed the grid-specific densities across grid cells to generate a total abundance by species. This process was replicated for each MCMC draw, allowing us to generate 95% credible intervals, and other measures of uncertainty on the estimates of abundance for each species in each year. 

#Combined with our spatial These parameter 

#Start with analyzing the 40 most common species, representing 90% of the total occurrences...
for(ii in seq(1, length(t))) {

  this.spp = names(cumsum(t))[ii]

  subset = ichthyo[which(ichthyo$IchName==this.spp),]
  df = data.frame("Y"=subset$Latitude,"X"=subset$Longitude,"Z"=subset$CatchPer10m_2,"T"=subset$Year) # "Day"=subset$jday)

  # Convert date to integer, 1-31
  df$T = df$T - 1981 + 1
  #df$T = df$T - min(df$T) + 1
  
  #Convert to UTM
  attr(df, "projection") <-"LL"
  attr(df, "zone") <-5
  dfUTM = convUL(df)
  
  #Plot the raw data
  #par(mfrow=c(6,6),mai=c(0.3,0.3,0.1,0.1))
  #for(i in 1:31) {
  #plot(dfUTM$X[which(times==i)],dfUTM$Y[which(times==i)],cex=0.4, xlim = range(dfUTM$X), #ylim=range(dfUTM$Y),col=as.numeric(subset$Purpose[which(times==i)]))
  #}
  #dev.off()
  
  #Remove stations from N and E side of Kodiak -- tend to be sparsely sampled in time. Focus on core 90% -
  indx = which(dfUTM$X > quantile(dfUTM$X,0.05) & dfUTM$X < quantile(dfUTM$X,0.95) & dfUTM$Y > quantile(dfUTM$Y,0.05) & dfUTM$Y < quantile(dfUTM$Y,0.95))
  dfUTM = dfUTM[indx,]
  df = df[indx,]
  subset = subset[indx,]
  
  #Create grid locations for projection. 
  gridproj = expand.grid("Y" = seq(6100,6500,10), "X" = seq(10,570,10))
  # create distance matrix between data and projections
  allpts = rbind(gridproj,dfUTM[,1:2])
  alldist = as.matrix(dist(allpts))
  inclproj = rep(0, dim(gridproj)[1]) # should grid cell be included?
  for(i in 1:dim(gridproj)[1]) {
    indx = which(alldist[i,-c(1:dim(gridproj)[1])]<10)
    if(length(indx)>0) {
      nyears= length(unique(dfUTM$T[which(alldist[i,-c(1:dim(gridproj)[1])]<10)]))
      if(nyears >= 15) {inclproj[i]=1}
    }
    #if(min(alldist[i,-c(1:dim(gridproj)[1])]) < 10) {inclproj[i]=1}
  }
  gridproj$Include = inclproj
  
  #Create grid locations for projection.
  coords = unique(paste(dfUTM[,1],dfUTM[,2]))
  coords = unlist(strsplit(coords," "))
  coordsLat = as.numeric(coords[seq(1,length(coords),2)])
  coordsLon = as.numeric(coords[seq(2,length(coords),2)])
  x = cbind(coordsLon,coordsLat)
  library(INLA)
  bnd = inla.nonconvex.hull(x, convex=60)
  #  bnd = inla.nonconvex.hull(subcoords, convex=150)
  # increase cutoff to ~ 150 to create much coarser mesh
  #mesh1 = inla.mesh.2d(boundary=bnd,max.edge=c(60,1500),cutoff=150,offset=c(120,180))
  mesh1 = inla.mesh.2d(boundary=bnd,max.edge=c(100,100),cutoff=60,offset=c(50,50))

  plot(mesh1)
  summary(mesh1)
  points(x,col="blue",cex=0.5)
  points(gridproj[which(gridproj$Include==1),c(2,1)],col="red",cex=0.1)
  
  #Create grid locations for projection.
  # Make SPDE based on mesh
  spde=inla.spde2.matern(mesh1, alpha=3/2)
  n=nrow(dfUTM)
  
  # Pos model, turn pos values -> NA
  z = as.numeric(dfUTM$Z)
  z[which(z == 0)]=NA
  
  k=33 #max(dfUTM$T) # max number of years
  
  dat = data.frame(y=z, time=dfUTM$T, xcoo=dfUTM$X,ycoo=dfUTM$Y) #,logday=log(dfUTM$Day))
  
  # tack on year effects
  YEARS.lab = paste("Y",seq(1,k),sep="")
  dat[YEARS.lab] = 0
  # Make a design matrix where the first year is the intercept
  dat[,YEARS.lab[1]] = 1 # all other year coefs are offsets, like factors in glm()
  for(j in 1:length(YEARS.lab)){
    dat[dat$time == j,YEARS.lab[j]]	<-	1 
  }
  
  # construct index for iid / ar1 model    
  iset = inla.spde.make.index("i2D", n.spde=mesh1$n, n.group = k) 
  
  # Make the covariates
  X.1 = dat[,-c(1:4)]
  Covar.names <- colnames(X.1)
  XX.list <- as.list(X.1)
  effect.list <- list()    				
  #   effect.list[[1]] <- c(iset, list(Intercept=1))
  effect.list[[1]] <- c(iset)
  for (Z in 1:ncol(X.1)) effect.list[[Z+1]] <- XX.list[[Z]]
  names(effect.list) <- c("1", Covar.names)
  
  ### Make projection points stack.
  A <- inla.spde.make.A(mesh=mesh1, loc=cbind(dat$xcoo, dat$ycoo), group = dat$time)
  A.list = list()
  A.list[[1]] = A
  for (Z in 1:ncol(X.1)) A.list[[Z+1]] <- 1
  
  ### Make projection points stack.
  sdat <- inla.stack(tag='stdata', data=list(y=dat$y), A=A.list, effects=effect.list)
  
  formula = as.formula(paste0("y ~ -1 +",  paste(Covar.names, collapse="+"), "+ f(i2D, model=spde, group = i2D.group, control.group = list(model='ar1'))"))		# field evolves with AR1 by year
  
  inlaModel <- inla(formula, family = "gamma", data=inla.stack.data(sdat), control.predictor=list(compute=TRUE, A=inla.stack.A(sdat)), 
                    verbose = TRUE, debug=TRUE, keep=FALSE, control.compute = list(dic=TRUE, cpo=TRUE, config=TRUE), 
                    control.fixed = list(correlation.matrix=TRUE),control.inla = list(lincomb.derived.correlation.matrix=TRUE))

  # Generate MCMC samples
  nMCMC = 1000 # takes about 6 minutes to do the sampling
  inla.mcmc = inla.posterior.sample(nMCMC, inlaModel)
  
  # get indices of parameters in MCMC output
  re.indx = grep("i2D",rownames(inla.mcmc[[1]]$latent))# Get indx of random effects
  fe.indx = grep("Y1.1",rownames(inla.mcmc[[1]]$latent))[1]# Fixed effects stacked after random
  logday.indx = grep("logday.1",rownames(inla.mcmc[[1]]$latent))[1]# Fixed effects stacked after random
  
  projectedLatentGrid = array(0, dim = c(length(which(gridproj$Include==1)), nMCMC, k))
  # Loop over all years, and do MCMC projections for that year
  for(yr in 1:1) {
    #Grab random effects from this year
    indx = which(iset$i2D.group==yr)
    
    projMatrix <- inla.spde.make.A(mesh1, loc=as.matrix(gridproj[which(gridproj$Include==1),c(2,1)]))
    
    # Multiply this projection matrix x 
    for(n in 1:nMCMC) {
      # combine random effects, fixed effects
      projectedLatentGrid[,n,yr] = as.numeric(projMatrix%*%inla.mcmc[[n]]$latent[re.indx][indx]) + 
        inla.mcmc[[n]]$latent[fe.indx - 1 + yr] #*log(148)
    } # end mcmc loop
  } # end year loop
  
  # Loop over all years, and do MCMC projections for that year
  for(yr in 2:k) {
    #Grab random effects from this year
    indx = which(iset$i2D.group==yr)
    
    projMatrix <- inla.spde.make.A(mesh1, loc=as.matrix(gridproj[which(gridproj$Include==1),c(2,1)]))
    
    # Multiply this projection matrix x 
    for(n in 1:nMCMC) {
      # combine random effects, fixed effects
      projectedLatentGrid[,n,yr] = as.numeric(projMatrix%*%inla.mcmc[[n]]$latent[re.indx][indx]) + 
        inla.mcmc[[n]]$latent[fe.indx - 1 + 1] +   inla.mcmc[[n]]$latent[fe.indx - 1 + yr]
      #*log(148)
    } # end mcmc loop
  } # end year loop
  
  save(projectedLatentGrid, file=paste("positive/", this.spp,"_positive2.Rdata",sep=""))
  save(gridproj, file=paste("positive/", this.spp,"_gridpositive2.Rdata",sep=""))
}