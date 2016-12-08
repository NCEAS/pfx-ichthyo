

up.ind=read.table('data/MonthlyUpwellingIndex.txt', header=T) ##reads in monthly upwelling index, 
##data from http://www.pfeg.noaa.gov/products/PFELData/upwell/monthly/upindex.mon

up.ind$spr.avg=rowMeans(up.ind[,4:8])

up.ind$spr.avg[up.ind$YEAR>1980]


up.anom=read.table('data/MonthlyUpwellingIndexAnom.txt', header=T) ##reads in monthly upwelling index, 
##data from http://www.pfeg.noaa.gov/products/PFELData/upwell/monthly/upanoms.mon

up.anom$spr.avg=rowMeans(up.anom[,4:9]) ##

up.anom$spr.avg[up.anom$YEAR>1980]
plot(seq(1981,2015),up.anom$spr.avg[up.anom$YEAR>1980], type='o')


upwell=read.csv('data/ClimateCovars.csv', header=T)
points(seq(1981,2015),upwell$Upwelling[upwell$Year>1980], col='red', type='o')
