
#########        ForestPlot funtion prints out graphs from Lignum simulation of a forest stand
#########        Call: ForestPlot(infile, aplot, pick)
#########        infile: name of input HDF5 file (if not in the current directory include path). The file
#########                must be in the HDF5 File Format (see below)
#########        aplot   Area of the simulated plot in m2
#########        pick    Every pick tree will appear in the height and diameter growth graphs.
#########                (If there are e.g. 600 trees on plot graphs become a mess if all are plotted)


#The Lignum forest simulation results are stored in a file in the HDF5 File Format
#(https://en.wikipedia.org/wiki/Hierarchical_Data_Format).

#ForestPlot assumes the input file contains the following datasets: 
#              /  ForestTreeData H5I_DATASET  FLOAT 51 x 599 x 21
#              /       StandData H5I_DATASET  FLOAT       19 x 21

#In this example 599 is the numberof trees on plot and 21 is the number times (years) at which
#results are stored

#In order to run ForestPlot in R a suitable package is needed. ForestPlot has been tested using
#Biomanager package (https://www.rdocumentation.org/packages/BiocManager/versions/1.30.17):
#install.packages("BiocManager")
#BiocManager::install("rhdf5")

#----------------------------------------------------------------------------
gini <- function(v){
v <- na.omit(v)
n <- length(v)
sum(outer(v, v, FUN=function(x,y){abs(x-y)}),na.rm=TRUE) / (2 * n * sum(v,na.rm=TRUE))
}

ForestPlot <- function(infile,aplot,pick) {


d <- H5Fopen(infile)

pdf_file <- paste(infile,".pdf",sep="")


pdf(pdf_file)

y <- d$StandData[1,]

#Height
plot(y,d$StandData[11,], type="l", ylim=c(0,1.2*d$StandData[11,length(y)]), lwd=2, xlab="time (y)", ylab="Tree height, nin, mean, max (m)", main="Mean, min and max stand height") #mean
points(y,d$StandData[12,], type="l", lwd=2, lty=2)   #min
points(y,d$StandData[13,], type="l",lwd=2, lty=2)   #max

#Base diameter
plot(y,100*d$StandData[5,], type="l", lwd=2, ylim=c(0,1.2*100*d$StandData[7,length(y)]),xlab="time (y)", ylab="Base diam, nin, mean, max (cm)", main="Mean, min and max diameter at base in the stand") #mean
points(y,100*d$StandData[6,], type="l",lwd=2, lty=2)   #min
points(y,100*d$StandData[7,], type="l",lwd=2, lty=2)   #max

# longest and shortest trees
h <- d$ForestTreeData[7,,length(y)]

mh <- max(d$ForestTreeData[7,,length(y)],na.rm=TRUE)
largest <- which(h>0.999*mh)[1]

mh <- min(d$ForestTreeData[7,,length(y)],na.rm=TRUE)
smallest <- which(h<1.001*mh)[1]

mh <- median(d$ForestTreeData[7,,length(y)],na.rm=TRUE)
med <- which(h<1.002*mh&h>0.98*mh)[1]

#Height
plot(y,d$ForestTreeData[7,largest,], type="l", ylim=c(0,1.2*d$ForestTreeData[7,largest,length(y)]), lty=1,xlab="time (y)", ylab="Tree height (m)",lwd=2, main=paste("Height growth of (at age ", as.character(length(y)-1),") shortest (red), median (green),\n tallest (blue) tree",sep=""), col="blue") #largest
points(y,d$ForestTreeData[7,med,], type="l",lwd=2, col="darkgreen")    #median
points(y,d$ForestTreeData[7,smallest,], type="l",lwd=2,col="red")    #smallest

#Base diameter
plot(y,100*d$ForestTreeData[8,largest,], type="l", ylim=c(0,1.5*100*d$ForestTreeData[8,largest,length(y)]), lty=1,xlab="time (y)", ylab="Diameter at base (cm)",lwd=2, main=paste("Diameter growth of (at age ", as.character(length(y)-1),") shortest (red), median (green),\n tallest (blue) tree",sep=""),col="blue") #largest
points(y,100*d$ForestTreeData[8,med,], type="l",lwd=2, col="darkgreen")    #median
points(y,100*d$ForestTreeData[8,smallest,], type="l",lwd=2,col="red")    #median

#Density
aplot1 <- aplot/1e4          #area in ha
plot(y,d$StandData[3,]/aplot1, type="l", ylim=c(0,1.1*d$StandData[3,1]/aplot1),lty=1, lwd=2, xlab="time (y)", ylab="No. trees / ha", main="Stand density")

#Self thinning plot
aplot1 <- aplot/1e4          #area in ha
plot(log(d$StandData[5,]),log(d$StandData[3,]/aplot1), xlim=c(log(0.001),log(0.5)),ylim=c(log(100),log(20000)),type="l", lty=1, lwd=2, xlab="log(mean base diameter)", ylab="log(No. trees / ha)", main="Self-thinning curve")


#Basal area
plot(y,d$StandData[14,]*1e4, ylim=c(0,50),type="l", lwd=2,xlab="time (y)", ylab= "m2/ha", main="Basal area")

#Basal area at crown base
plot(y,d$StandData[15,]*1e4, ylim=c(0,30),type="l", lwd=2,xlab="time (y)", ylab= "m2/ha", main="Basal area ar crown base")

#Stem volume
plot(y,d$StandData[16,]*1e4, ylim=c(0,500),type="l", lwd=2,xlab="time (y)", ylab= "m3/ha", main="Stem volume")


#LAI and specific leaf area
par(mar = c(5, 4, 4, 4) + 0.3)              # Additional space for second y-axis
plot(y,d$StandData[17,], ylim=c(0,15),type="l", lty=1, lwd=2,xlab="time (y)", ylab="All-sided LAI (m2/m2)",main="LAI = continuous,  specific LA = dashed") #LAI
par(new = TRUE)
plot(y,d$StandData[17,]/d$StandData[18,],type="l", lty=2,lwd=2, axes=FALSE,xlab = "", ylab = "", ylim=c(0,32))
axis(side = 4, at = pretty(c(0,32)))
mtext("Specific leaf area (m2/ kg C)", side = 4, line = 3)             # Add second axis label


#P / Af
plot(y,apply(d$ForestTreeData[17,,]/d$ForestTreeData[15,,],2,mean,na.rm=TRUE), ylim=c(0,0.2),type="l", lwd=2,xlab="time (y)", ylab= "P/Af min, mean, max (kg C / MJ PAR)", main="Photosynthetic rate per unit needle area") #   P/Af
points(y,apply(d$ForestTreeData[17,,]/d$ForestTreeData[15,,],2,min,na.rm=TRUE), type="l", lty=2, lwd=2)   #min
points(y,apply(d$ForestTreeData[17,,]/d$ForestTreeData[15,,],2,max,na.rm=TRUE), type="l", lty=2, lwd=2)   #max

#Lambda
plot(y,apply(d$ForestTreeData[51,,],2,mean,na.rm=TRUE), type="l", lty=2, lwd=2,xlab="time (y)", ylab= expression(paste(lambda," min, mean, max")),ylim=c(0,7), main = expression(paste(lambda," in Eq: New growth(",lambda,") = P - M")))
points(y,apply(d$ForestTreeData[51,,],2,min,na.rm=TRUE), type="l", lty=1, lwd=2)   #min
points(y,apply(d$ForestTreeData[51,,],2,max,na.rm=TRUE), type="l", lty=1, lwd=2)   #max

#Lambda largest, median shortest
plot(y,d$ForestTreeData[51,largest,], type="l", ylim=c(0,2), lty=1,xlab="time (y)", ylab="lambda",lwd=2, main=paste("Progression of lambda in shortest, median, tallest (at age ", as.character(length(y)-1),") tree"),col="blue") #largest
points(y,d$ForestTreeData[51,med,], type="l",lwd=2,col="green")    #median
points(y,d$ForestTreeData[51,smallest,], type="l",lwd=2,col="red")    #median



Ntrees <- d$StandData[3,1]

#Tree heights
plot(y,d$ForestTreeData[7,1,], ylim=c(0,30), type="l", main=paste("Individual tree heights\nevery ",as.character(pick),"th tree",sep=""), ylab="Tree height (m)",xlab="time (y)")
points(y,d$ForestTreeData[11,1,],type="l",lty=2)
for(i in 2:min(Ntrees/pick)) {
	points(y,d$ForestTreeData[7,i*pick,], type="l")
	points(y,d$ForestTreeData[11,i*pick,],type="l",lty=2)
}


#Tree diameters at base
plot(y,100*d$ForestTreeData[8,1,], ylim=c(0,30), type="l",main=paste("Individual tree diameter at base\nevery ",as.character(pick),"th tree",sep=""),xlab="time (y)", ylab="Tree diameter (cm)")
for(i in 2:min(Ntrees/pick)) {
	points(y,100*d$ForestTreeData[8,pick*i,], type="l")
}

#Crown ratio
plot(y,1-d$ForestTreeData[11,1,]/d$ForestTreeData[7,1,], ylim=c(0,1), type="l", main=paste("Crown ratios\nevery ",as.character(pick),"th tree",sep=""), ylab="Crown ratio",xlab="time (y)")
for(i in 2:min(Ntrees/pick)) {
	points(y,1-d$ForestTreeData[11,1,]/d$ForestTreeData[7,i*pick,], type="l")
}


#Height and diameter distributions
h <- d$ForestTreeData[7,,length(y)]
hist(h[h>0.9*d$StandData[12,length(y)]], main=paste("Height distribution at age ", as.character(max(y)),sep=""), xlab="Tree height (m)")
#
db <- d$ForestTreeData[13,,length(y)]
hist(db[h>0.9*d$StandData[12,length(y)]], main=paste("Distribution of diameter at base at age ", as.character(max(y)),sep=""), xlab="Diameter (cm)")
##

plot(y,apply(d$ForestTreeData[7,,],2,gini),ylim=c(0,1),type="l", lwd=2,xlab="time (y)", ylab= "Gini coefficient", main="Gini coeff., solid = height, dashed = base diameter") 
points(y,apply(d$ForestTreeData[8,,],2,gini),type="l", lwd=2, lty=2)

lc <- Lc(na.omit(d$ForestTreeData[7,,length(y)]))
plot(lc$p,lc$L, type="l", lwd=2, main="Lorenz curve of tree heights(solid) and\ndiameters (dashed), red = all equal", xlab="Cumulative % of trees from smallest to largest", ylab="Cumulative % total height or diameter")
lc <- Lc(na.omit(d$ForestTreeData[8,,length(y)]))
points(lc$p,lc$L, type="l",lwd=2,lty=2)
abline(0,1,col="red")



dev.off()

}


