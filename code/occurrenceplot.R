
occdat=read.csv('data/occurrenceStatsForDFA.csv',header=T)

par(mfrow=c(4,5), mar=c(0,0,0,0), oma=c(4,4,1,1))
par(mgp = c(.6, .6, 0))
par(tcl=.2)
for(i in 1:length(spp)){
  plot(years.plot,dat.z[i,],xlab="",ylab="Index", bty="L", pch=16, col="blue", type="b", ylim=c(min(dat.z, na.rm=T),max(dat.z, na.rm=T)), axes=F)
  
  if(i %in% c(16:20,36:40)) axis(1,at=seq(1990,2010, by=10), cex.axis=0.8)
  if(i %in% c(1,6,11,16, 21, 26, 31,36)) axis(2, at=seq(-2,5, by=1), cex.axis=0.8)
  box()
  mtext(spp[i], side=3, line=-1, cex=0.5)
}


mtext("Year", side=1, outer=TRUE, line=2, cex=0.7)
mtext("Standardized abundance", side=2, outer=TRUE, line=2, cex=0.7)