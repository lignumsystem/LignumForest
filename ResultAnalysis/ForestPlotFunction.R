library(rhdf5)
library(ineq)
#########        ForestPlot function prints out graphs from Lignum simulation of a forest stand
#########        Call: ForestPlot(infile, aplot, pick, GYdata = "", center = "a")
#########        infile: name of input HDF5 file (if not in the current directory include path). The file
#########                must be in the HDF5 File Format (see below)
#########                Output is writen to a pdf file, name is infile.pdf
#########        aplot   Area of the simulated plot in m2
#########        pick    Every pick tree will appear in the height and diameter growth graphs.
#########                (If there are e.g. 600 trees on plot graphs become a mess if all are plotted)
#########        GYdata = ""       Path to the Growh & Yield data used for comparison
#########        center = "a"      If all trees (= default) or only in center of plot are used for evaluation
#########                          If center = "c", center of plot, any other character all trees, the extent
#########                          of the excluded area is read from the VoxelSpace of the hdf5 file.
#########                          In the case of center of plot the name of output file is infile_c.pdf
#########         Wfall = "n"      If plots Wf vs tree diameter and height are shown (Wfall = "Y")
#########         Gini = "n"       If plots of Gini coefficient and Lorennz curve are shown (Gini = "Y")



#The Lignum forest simulation results are stored in a file in the HDF5 File Format
#(https://en.wikipedia.org/wiki/Hierarchical_Data_Format).

#ForestPlot assumes the input file contains the following datasets: 
#              /  ForestTreeData H5I_DATASET  FLOAT 51 x 599 x 21
#              /       StandData H5I_DATASET  FLOAT       19 x 21
#              / CenterStandData H5I_DATASET  FLOAT       19 x 81

#In this example 599 is the numberof trees on plot and 21 is the number times (years) at which
#results are stored

#In order to run ForestPlot in R a suitable package is needed. ForestPlot has been tested using
#Biomanager package (https://www.rdocumentation.org/packages/BiocManager/versions/1.30.17):
#install.packages("BiocManager")
#BiocManager::install("rhdf5")

#For Lorenz curve (inequality) install the package "ineq"
#install.packages("ineq")
#
#To run ForestPlot load libraries first
#library("rhdf5")
#library("ineq")
#Then source "ForestPlotFunction.R"
#source("ForestPlotFunction.R")
#Finally run the ForestPlot function, e.g.:
#ForestPlot("../HDF5ForestData770_ajo157_3100.h5",400,2,"ResultAnalysis/","c")
#The figures will appear in "../HDF5ForestData770_ajo157_3100.h5.pdf"

#NOTE! Three data files are are used: Va27.txt, VVV-40.txt, ksto-mt.dat and Ilvessalo.txt
#(in LignumForest/Resultanalysis).
#Va27.txt = Varmola M (1987) Männyn viljelytaimikoiden kasvumalli. Lic. For & Agric.
#Thesis, Department of Forest Mensuration, University of Helsinki, 89 p.
#VVV-40.txt = Vuokila Y, Väliaho H (1980) Viljeltyjen havumetsiköiden kasvatusmallit.
#Communicationes Instituti Forestalis Fenniae 99, 271.
#ksto-mt.dat.  Finnish yiels tables compiled by Koivisto
#Ilvessalo.txt Ilvessalo G&Y tables
#----------------------------------------------------------------------------

gini <- function(v){
v <- na.omit(v)
n <- length(v)
sum(outer(v, v, FUN=function(x,y){abs(x-y)}),na.rm=TRUE) / (2 * n * sum(v,na.rm=TRUE))
}

in_center <- function(x,y, pX, pY, dx, dy) {
	x >= dx & x <= pX-dx & y >= dy & y <= pY-dy	
}

ForestPlot <- function(infile, pick, GYdata = "", center = "a", WFall = "n", Gini = "n") {


d <- H5Fopen(infile)

va27 <- read.table(paste(GYdata,"Va27.txt",sep=""),header=FALSE)
colnames(va27) <- c("a", "Hd",   "HgM",    "DgM",   "V",  "Hc", "G")

vv <- read.table(paste(GYdata,"VVV-40.txt",sep=""),header=FALSE)
colnames(vv) <- c("age",    "DBH",		"H",		"Hcb",		"Wf", "V")

ksto <- read.table(paste(GYdata,"ksto-mt.dat",sep=""), header=TRUE)

ftdata<-read.table(paste(GYdata,"Ilvessalo.txt",sep=""),header=TRUE,sep='')





#The size of the plot and is read from VoxelSpace.txt & information about
#plot boudary
	
# The x and y sizes of the stand (plot) are given on the first line of VoxelSpace.txt
# as two first items, they are read in to plot_x, plot_y
# The distance of center area from stand edge (x and y) is given on the third line
# of VoxelSpace.txt as: dist_x dist_y

#values are separated by spaces and lines by '\n' in VoxelSpace.txt

pos <- unlist(gregexpr(" ",d$VoxelSpace))   #First line, X Y Z values
plot_x <- as.numeric(substr(d$VoxelSpace,1,pos[1]-1))
plot_y <- as.numeric(substr(d$VoxelSpace,pos[1]+1,pos[2]-1))

pos <- unlist(gregexpr('\n',d$VoxelSpace))
line3 <- substr(d$VoxelSpace,pos[2]+1,pos[3]-1)  #line 3 between 2nd and 3rd '\n'
pos <- unlist(gregexpr(" ",line3))
dist_x <- as.numeric(substr(line3,1,pos[1]-1))   #line3 contains two numbers separated by a space
dist_y <- as.numeric(substr(line3,pos[1]+1,nchar(line3)))
	


#All trees on the plot or only in the center center of the plot
trees <- d$ForestTreeData

if(center != "c") {
	stand <- d$StandData
	pdf_file <- paste(infile,".pdf",sep="")
	print("All trees included")
	mtxt = "All trees included"
	mukana <- 1:length(trees[1,,1])
	aplot = plot_x * plot_y			#Plot area in m2
} else {
	stand <- d$CenterStandData
	pdf_file <- paste(infile,"_c.pdf",sep="")
	print("Center stand")
	mtxt = "Center stand"
	# If center stand, individual trees must be checked whether they belong to the center
	mukana <- which(in_center(trees[2,,1],trees[3,,1],plot_x,plot_y,dist_x,dist_y))
	aplot = (plot_x - 2*dist_x) * (plot_y - 2*dist_y)			#Plot area in m2
}

pdf(pdf_file)


y <- stand[1,]
ymax = max(y, na.rm=TRUE)


###Height
plot(y,stand[11,], type="l", ylim=c(0,40), lwd=2, xlab="time (y)", ylab="Tree height, nin, mean, max (m)", main=paste("Mean, min and max stand height\nIn the plots: ",mtxt,sep=""),sub=paste("File: ",infile,sep=""))##mean
legend('topleft',inset=0.05,c("Lignum mean","Lignum (min, max)"," Vuokila/Väliaho 1980"," Varmola M 1987","Koivisto: kasvu- ja tuotostaulukot"),col=c('black','black','blue','red','darkgreen'),
       lty=c(1,2,1,1,1),lwd=2)

points(y,stand[12,], type="l", lwd=2, lty=2)   #min
points(y,stand[13,], type="l",lwd=2, lty=2)   #max
points(va27$a,va27$HgM,type="l",lwd=3,col="red")
points(vv$age,vv$H,type="l",lwd=3,col="blue")
points(ksto$year,ksto$Hav,type="l",lwd=3,col="darkgreen")




# longest and shortest trees
h <- trees[7,mukana,ymax]

mh <- max(trees[7,mukana,ymax],na.rm=TRUE)
largest <- which(h>0.98*mh)[1]

mh <- min(trees[7,mukana,ymax],na.rm=TRUE)
smallest <- which(h<1.02*mh)[1]

mh <- median(trees[7,mukana,ymax],na.rm=TRUE)
med <- which(h<1.02*mh&h>0.98*mh)[1]


#DBH
plot(y,100*stand[8,], type="l", lwd=2, ylim=c(0,40),xlab="time (y)", ylab="Diameter, min, mean, max (cm)", main="Mean, min and max BH diameter in the stand")##mean
legend('topleft',inset=0.05,c("Lignum (mean)","Lignum (min,max)"," Vuokila/Väliaho 1980"," Varmola M 1987","Koivisto: kasvu- ja tuotostaulukot"),col=c('black','black','blue','red','darkgreen'),
       lty=c(1,2,1,1,1),lwd=2)

points(y,100*stand[9,], type="l",lwd=2, lty=2)   #min
points(y,100*stand[10,], type="l",lwd=2, lty=2)   #max
points(va27$a,va27$DgM,type="l",lwd=3,col="red")
points(vv$age,vv$DBH,type="l",lwd=3,col="blue")
points(ksto$year,ksto$Dbhav,type="l",lwd=3,col="darkgreen") 


#Height
plot(y,trees[7,mukana[largest],], type="l", ylim=c(0,1.2*trees[7,mukana[largest],ymax]), lty=1,xlab="time (y)", ylab="Tree height (m)",lwd=2, main=paste("Height growth of (at age ", as.character(ymax-1),") shortest (red), median (green),\n tallest (blue) tree",sep=""), col="blue") #largest
points(y,trees[7,mukana[med],], type="l",lwd=2, col="darkgreen")  #median
points(y,trees[7,mukana[smallest],], type="l",lwd=2,col="red")    #smallest


#Height vs breast height diameter
plot(100*stand[8,],stand[11,], type="l", xlim=c(0,30), ylim=c(0,30), lwd=2, xlab="BH diameter (cm)", ylab="Mean tree height (m)", main="Height vs breast height diameter\nGreen = Koivisto, Varmola, Vuok&V:o, dH=f(CR)*dD (blue)")
points(ksto$Dbhav,ksto$Hav,type="l",lwd=2,col="darkgreen")
points(va27$DgM,va27$HgM, type="l",lwd=2,col="darkgreen")
points(vv$DBH,vv$H, type="l",lwd=2,col="darkgreen")

#Height increment from diam. increment using a model component of Sievänen, R. 1993.
#A process-based model for dimensional growth of even-aged stands. Scand. J. For. Res. 8:28-48 that
#has been adjusted to Finnhs pine forests
Hinc <-(50+130*(1-apply(1-trees[11,mukana,]/trees[7,mukana,],2,mean,na.rm=TRUE)[1:ymax+1]))*diff(100*stand[8,])

points(100*stand[8,1:ymax],1+cumsum(Hinc)/100,type="l",lwd=2,col="blue")



###Density
aplot1 <- aplot/1e4          #area in ha

t1 <- stand[1,]
sdensity <- stand[3,]/aplot1
plot(t1,sdensity,type="l",lty=1, lwd=2, xlab="Time (y)", ylab="No. trees / ha",ylim=c(0,17000),
     main="Stand density",sub="Data: Ilvessalo")
legend('bottomleft',inset=0.05,c("Lignum","MT","VT","CT","Koivisto MT"),col=c('black','darkgreen','red','brown', 'blue', title="Forest type"),
       lty=1,lwd=2)
t2<-ftdata$Year
mt<-ftdata$MT
vt<-ftdata$VT
ct<-ftdata$CT
lines(t2,mt,col='darkgreen',lty=1,lwd=2)
lines(t2,vt,col='red',lty=1,lwd=2)
lines(t2,ct,col='brown',lty=1,lwd=2)
points(ksto$year,ksto$N,type="l",lwd=3,col="blue")


#self thinning plot
aplot1 <- aplot/1e4          	#area in ha
gt0 <- which(stand[8,] > 0)		#at the beginning Dbh = 0
plot(stand[8,gt0],stand[3,gt0]/aplot1,log="xy",ylim=c(500,20000),type="l", lty=1, lwd=2,
xlab="log(mean Dbh)", ylab="log(No. trees / ha)", main="Self-thinning curve")
legend('bottomleft',inset=0.05,cex=0.8,c("Lignum",
expression(paste(N == alpha*bar(D)^{-3/2}," (", N," = Density",", ",bar(D)," = RMS stand diameter",","," Reineke ",1933^(1),")"),"Koivisto: kasvu- ja tuotostaulukot")),col=c('black','red','darkgreen'),lty=1,lwd=2)
text(0.14,460,cex=0.6,"(1) Here in the context of Koivisto")
points(ksto$Dbhav/100,91000.0*(ksto$Dbhav)**(-3/2),type="l",lwd=2,col="red")
points(ksto$Dbhav/100,ksto$N,type="l",lwd=3,col="darkgreen")

#Basal area
plot(y,stand[14,]*1e4, ylim=c(0,80),type="l", lwd=2,xlab="time (y)", ylab= "m2/ha", main="Basal area")
legend('topleft',inset=0.05,c("Lignum","Varmola M 1987","Koivisto: kasvu- ja tuotostaulukot"),col=c('black','red','darkgreen'),
       lty=1,lwd=2)
points(va27$a,va27$G,type="l",lwd=3,col="red")
points(ksto$year,ksto$G,type="l",lwd=3,col="darkgreen")


#Basal area at crown base
plot(y,stand[15,]*1e4, ylim=c(0,60),type="l", lwd=2,xlab="time (y)", ylab= "m2/ha", main="Basal area at crown base")

#Stem volume
plot(y,stand[16,]*1e4, ylim=c(0,1000),type="l", lwd=2,xlab="time (y)", ylab= "m3/ha", main="Stem volume")
legend('topleft',inset=0.05,c("Lignum","Vuokila/Väliaho 1980","Varmola M 1987","Koivisto: kasvu- ja tuotostaulukot"),col=c('black',"blue",'red','darkgreen'),
       lty=1,lwd=2)
points(va27$a,va27$V,type="l",lwd=3,col="red")
points(vv$age,vv$V,type="l",lwd=3,col="blue")
points(ksto$year,ksto$V,type="l",lwd=3,col="darkgreen")


#LAI and specific leaf area
par(mar = c(5, 4, 4, 4) + 0.3)              # Additional space for second y-axis
plot(y,stand[17,], ylim=c(0,30),type="l", lty=1, lwd=2,xlab="time (y)", ylab="All-sided LAI (m2/m2)",main="Leaf Area Index and Specific Leaf Area") ###LAI
legend('bottomright',inset=0.05,c("LAI","SLA"),col='black',lty=c(1,2),lwd=2)
par(new = TRUE)
plot(y,stand[17,]/stand[18,],type="l", lty=2,lwd=2, axes=FALSE,xlab = "", ylab = "", ylim=c(0,32))
axis(side = 4, at = pretty(c(0,32)))
mtext("SLA (m2/ kg C)", side = 4, line = 3)             # Add second axis label


#Photosynthetic production and respiration
plot(y,1e-3*apply(trees[17,mukana,],2,sum,na.rm=TRUE)/aplot1,type="l", lwd=2,xlab="time (y)", ylab= "t C / ha", main="GPP and Respiration") ###   P/Af
legend('bottomright',inset=0.05,c("GPP","Respiration"),col='black',lty=c(1,2),lwd=2)    
points(y,1e-3*apply(trees[18,mukana,],2,sum,na.rm=TRUE)/aplot1, type="l", lty=2, lwd=2)   #min



######       ALL TREES

#Ntrees <- stand[3,1]
Ntrees <- length(mukana)


#Tree heights
plot(y,trees[7,mukana[1],], ylim=c(0,40), type="l", main=paste("Individual tree heights\nevery ",as.character(pick),"th tree",sep=""), ylab="Tree height (m)",xlab="time (y)")
points(y,trees[11,mukana[1],],type="l",lty=2)
for(i in 2:min(Ntrees/pick)) {
	points(y,trees[7,mukana[i*pick],], type="l")
	points(y,trees[11,mukana[i*pick],],type="l",lty=2)
}


#Tree diameters at BH
plot(y,100*trees[9,mukana[1],], ylim=c(0,30), type="l",main=paste("Individual tree diameter at BH\nevery ",as.character(pick),"th tree",sep=""),xlab="time (y)", ylab="Tree diameter (cm)")
for(i in 2:min(Ntrees/pick)) {
	points(y,100*trees[9,mukana[pick*i],], type="l")
}


#Tree diameters at Crown Base
plot(y,100*trees[10,mukana[1],], ylim=c(0,15), type="l",main=paste("Individual tree diameter at Crown Base\nevery ",as.character(pick),"th tree",sep=""),xlab="time (y)", ylab="Tree diameter (cm)")
for(i in 2:min(Ntrees/pick)) {
	points(y,100*trees[10,mukana[pick*i],], type="l")
}


#Crown ratio
plot(y,1-trees[11,mukana[1],]/trees[7,mukana[1],], ylim=c(0,1), type="l", main=paste("Crown ratios and mean\nevery ",as.character(pick),"th tree",sep=""), ylab="Crown ratio",xlab="time (y)",xlim=c(0,ymax))
legend('bottomleft',inset=0.05,c("Lignum trees","Lignum trees mean"),col=c('black','red'),lty=1,lwd=2)
for(i in 1:min(Ntrees/pick)) {
	points(y,1-trees[11,mukana[i*pick],]/trees[7,mukana[i*pick],], type="l")
}
points(y,apply(1-trees[11,mukana,]/trees[7,mukana,],2,mean,na.rm=TRUE),type="l",lwd=2,col="red")



###Height vs diameter
plot(100*trees[9,mukana[1],],trees[7,mukana[1],], ylim=c(0,30), xlim=c(0,30), type="l",main=paste("Height vs diameter BH\nevery ",as.character(pick),"th tree, dH=f(CR)*dD (blue)",sep=""),
     xlab="Tree diameter (cm)", ylab="Tree height (m)")
legend('bottomright',inset=0.05,c("Lignum trees","y=x"),col=c('black','red'),lty=1,lwd=2)
for(i in 2:min(Ntrees/pick)) {
	points(100*trees[9,mukana[i*pick],],trees[7,mukana[i*pick],], type="l")
}

#This is the same as in #Height vs breast height diameter at stand level
points(100*stand[8,1:ymax],1+cumsum(Hinc)/100,type="l",lwd=2,col="blue")


#Mean Branch legth
#bmax <- max(c(trees[48,mukana[largest],ymax],trees[48,mukana[med],ymax],trees[48,mukana[smallest],ymax]))
plot(y,trees[48,mukana[1],], type="l", ylim=c(0,4), lty=1,xlab="time (y)", ylab="D2 weighted mean length (m)",lwd=2, main=paste("Mean branch length every ",as.character(pick),"th tree",sep=""))
for(i in 2:min(Ntrees/pick)) {
	points(y,trees[48,mukana[i*pick],], type="l")
}


###Foliage mass vs cross-sectional area at crown base
plot((pi/4)*100^2*trees[10,mukana[1],]^2,2*trees[23,mukana[1],], ylim=c(0,15), xlim=c(0,200), type="l",
     main=paste("Foliage mass vs stem cross sectional area at crown base\nevery ",as.character(pick),"th tree",sep=""),
     xlab="Stem cross section area at crown base  (cm2)", ylab="Foliage mass (kg DM)")
legend('topleft',inset=0.05,c("Lignum trees","Lignum trees mean","y=0.055x"),col=c('black','red','blue'),lty=1,lwd=2)
for(i in 2:min(Ntrees/pick)) {
	points(100^2*trees[10,mukana[i*pick],]^2,2*trees[23,mukana[i*pick],], type="l")
}
points(apply((pi/4)*100^2*trees[10,,]^2,2,mean,na.rm=TRUE),apply(2*trees[23,,],2,mean,na.rm=TRUE),type="l",lwd=2,col="red")
abline(0,0.055,col="blue",lwd=2)

if(WFall == "Y") {
	###Foliage mass vs diameter at breast height
	plot(100*trees[9,mukana[1],],2*trees[23,mukana[1],], ylim=c(0,15), xlim=c(0,30), type="l",
     main=paste("Foliage mass vs DBH\nevery ",as.character(pick),"th tree",sep=""),
     xlab="Diameter at breast height (cm)", ylab="Foliage mass (kg DM)")

	for(i in 2:min(Ntrees/pick)) {
		points(100*trees[9,mukana[i*pick],],2*trees[23,mukana[i*pick],], type="l")
	}

	###Foliage mass vs tree height
	plot(trees[7,mukana[1],],2*trees[23,mukana[1],], ylim=c(0,15), xlim=c(0,30), type="l",
     main=paste("Foliage mass vs tree height\nevery ",as.character(pick),"th tree",sep=""),
     xlab="Tree height (m)", ylab="Foliage mass (kg DM)")
	for(i in 2:min(Ntrees/pick)) {
		points(trees[7,mukana[i*pick],],2*trees[23,mukana[i*pick],], type="l")
	}

}

#Height and diameter distributions
h <- trees[7,mukana,ymax]
hist(h[h>0.9*stand[12,ymax]], main=paste("Height distribution at age ", as.character(ymax),sep=""), xlab="Tree height (m)")
#
db <- trees[9,mukana,ymax]
hist(db[h>0.9*stand[12,ymax]], main=paste("Distribution of BH diameter at age ", as.character(max(y)),sep=""), xlab="Diameter (cm)")
##

if(Gini == "Y") {
	plot(y,apply(trees[7,mukana,],2,gini),ylim=c(0,1),type="l", lwd=2,xlab="time (y)", ylab= "Gini coefficient", 	main="Gini coeff., solid = height, dashed = base diameter") 
	points(y,apply(trees[8,mukana,],2,gini),type="l", lwd=2, lty=2)

	lc <- Lc(na.omit(trees[7,mukana,ymax]))
	plot(lc$p,lc$p-lc$L, type="l", lwd=2, main="Difference of Lorenz curve of tree heights (solid) and\ndiameters 	(dashed) to all equal", xlab="Cumulative % of trees from smallest to largest", ylab="Cumulative % total height 	or diameter",ylim=c(0,0.5))
	lc <- Lc(na.omit(trees[8,mukana,ymax]))
points(lc$p,lc$p-lc$L, type="l",lwd=2,lty=2)
}


#P / Af
plot(y,apply(trees[17,mukana,]/trees[15,mukana,],2,mean,na.rm=TRUE), ylim=c(0,0.2),type="l", lwd=2,xlab="time (y)", ylab= "P/Af min, mean, max (kg C / MJ PAR)", main="Photosynthetic rate per unit needle area\n(min, mean, max)") #   P/Af
points(y,apply(trees[17,mukana,]/trees[15,mukana,],2,min,na.rm=TRUE), type="l", lty=2, lwd=2)   #min
points(y,apply(trees[17,mukana,]/trees[15,mukana,],2,max,na.rm=TRUE), type="l", lty=2, lwd=2)   #max


plot(y,trees[35,mukana[1],]/(trees[32,mukana[1],]*trees[15,mukana[1],]), type="l",main=paste("Radiation capture efficiency, related to STAR of tree\nevery ",as.character(pick),"th tree",sep=""), xlab="time (y)", ylab="Qabs/(Af*QinTop)",ylim=c(0,0.2))
for(i in 2:min(Ntrees/pick)) {points(y,trees[35,mukana[pick*i],]/(trees[32,mukana[pick*i],]*trees[15,mukana[pick*i],]), type="l")}

#Lambda
#at time = 0, lambda does not have a value
plot(y[2:ymax],apply(trees[51,mukana,2:ymax],2,mean,na.rm=TRUE), type="l", lty=1, lwd=2,xlab="time (y)", ylab= expression(paste(lambda," min, mean, max")),ylim=c(0,7), main = expression(paste(lambda," in Eq: New growth(",lambda,") = P - M (min, mean, max)")))
points(y[2:ymax],apply(trees[51,mukana,2:ymax],2,min,na.rm=TRUE), type="l", lty=2, lwd=2)   #min
points(y[2:ymax],apply(trees[51,mukana,2:ymax],2,max,na.rm=TRUE), type="l", lty=3, lwd=2)   #max

#Lambda largest, median shortest
plot(y[2:ymax],trees[51,mukana[largest],2:ymax], type="l", ylim=c(0,2), lty=1,xlab="time (y)", ylab="lambda",lwd=2, main=paste("Progression of lambda in shortest, median, tallest (at age ", as.character(ymax-1),"\nblue = largest, green = meadian, red = smallest", sep=""),col="blue") #largest
points(y[2:ymax],trees[51,mukana[med],2:ymax], type="l",lwd=2,col="green")    #median
points(y[2:ymax],trees[51,mukana[smallest],2:ymax], type="l",lwd=2,col="red")    #median

dev.off()

}
