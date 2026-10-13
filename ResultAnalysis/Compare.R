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

in_center <- function(x,y, pX, pY, dx, dy) {
        x >= dx & x <= pX-dx & y >= dy & y <= pY-dy     
}

#files <- c('A-koe23.h5','N10000.h5','N7000.h5','N5000.h5','N3300.h5','N2500.h5')
#legends <- c('N15000','N10000','N7000','N5000','N3300','N2500')
#center = 'c'
#traits <- c("height","diameter","heightvsdiameter","selfthinning","basalarea","basalareaatcrownbase",
# "stemvolume","leafareaindex", "crownratio", "dominantheight","basalareavsdominantheight")

Compare <- function(traits, files, legends, center = "c", pdfd = "") {



d <- H5Fopen(files[1])
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
        print("All trees included")
        if(pdfd == "") {
		pdf_file = "compare_all.pdf"
		} else {
		pdf_file = paste(pdfd,"/","compare_all.pdf",sep="")
		}
	mtxt = "All trees included"
	aplot = plot_x * plot_y			#Plot area in m2
} else {
	print("Center stand")
	        if(pdfd == "") {
		pdf_file = "compare_c.pdf"
		} else {
		pdf_file = paste(pdfd,"/","compare_c.pdf",sep="")
		}
	mtxt = "Center stand"
	aplot = (plot_x - 2*dist_x) * (plot_y - 2*dist_y)	#Plot area in m2
}

H5Fclose(d)


pdf(pdf_file)

#
#First explanation of the graphs
#
colors <- c('black','red','darkgreen','blue', 'brown', 'cyan','orangered3','yellow3')

plot(0,0,axes=FALSE, ylab="",xlab="", type="l")

legend('top',"Explanation",box.lty=0)

legend('top', inset=0.1, legends, col= colors[1:length(legends)],lty=1,lwd=2,box.lty=0)

legend('center',"FILES:", box.lty=0)

legend('bottom', files, box.lty=0)


#
# 					Tree height
#

if(length(which(traits=="height") > 0)) {
	d <- H5Fopen(files[1])
	if(center != "c") {
        stand <- d$StandData
	} else {
        stand <- d$CenterStandData
	}
	y <- stand[1,]
	plot(y,stand[11,], type="l", ylim=c(0,30), lwd=2, xlab="time (y)", ylab="Tree height (m)", main= 	paste("Mean height ",mtxt,sep=""), col=colors[1])
	H5Fclose(d)

for(i in 2:length(files)) {
d <- H5Fopen(files[i])
if(center != "c") {
        stand <- d$StandData
} else {
        stand <- d$CenterStandData
}
y <- stand[1,]
points(y,stand[11,], type="l", lwd=2, col=colors[i])
H5Fclose(d)
}
}

#
# 					Breast height diameter
#
if(length(which(traits=="diameter") > 0)) {
d <- H5Fopen(files[1])
if(center != "c") {
        stand <- d$StandData
} else {
        stand <- d$CenterStandData
}
y <- stand[1,]
plot(y,100*stand[8,], type="l", ylim=c(0,40), lwd=2, xlab="time (y)", ylab="Diameter (m)", main= "Mean Breast height diameter", col=colors[1])
H5Fclose(d)


for(i in 2:length(files)) {
d <- H5Fopen(files[i])
if(center != "c") {
        stand <- d$StandData
} else {
        stand <- d$CenterStandData
}
y <- stand[1,]
points(y,100*stand[8,], type="l", lwd=2, col=colors[i])
H5Fclose(d)
}
}

#
# 					Height vs breast height diameter
#

if(length(which(traits=="heightvsdiameter") > 0)) {
d <- H5Fopen(files[1])
if(center != "c") {
        stand <- d$StandData
} else {
        stand <- d$CenterStandData
}
plot(100*stand[8,],stand[11,], type="l", ylim=c(0,30), xlim=c(0,30),lwd=2, xlab="BH diameter (cm)", ylab="Mean tree height (m)", main="Height vs breast height diameter", col=colors[1])
H5Fclose(d)

for(i in 2:length(files)) {
d <- H5Fopen(files[i])
if(center != "c") {
        stand <- d$StandData
} else {
        stand <- d$CenterStandData
}
y <- stand[1,]
points(100*stand[8,], stand[11,], type="l", lwd=2, col=colors[i])
H5Fclose(d)
}
}
#
# 					Self-thinning 
#

if(length(which(traits=="selfthinning") > 0)) {
aplot1 <- aplot /1e4     #plot area in ha
d <- H5Fopen(files[1])
if(center != "c") {
        stand <- d$StandData
} else {
        stand <- d$CenterStandData
}
plot(stand[16,]*1e4/(stand[3,]/aplot1),stand[3,]/aplot1, log="xy",ylim=c(500,20000),type="l", xlim=c(0.001,0.5),lwd=2, xlab="Mean tree stem volume (m3)", ylab="No. trees / ha", main="Self-thinning curve",col=colors[1])
points(seq(0.02,0.5,0.1),1100*(seq(0.02,0.5,0.1))**(-2/3),type="l",lwd=2,col="magenta3")
legend('topright',inset=0.01,cex=0.8,expression(paste(N == alpha*v^{-2/3}," (", N," = Density",", ",v," = mean V",")")),lty=1,lwd=2,col="magenta3",box.lty = 0)
H5Fclose(d)

for(i in 2:length(files)) {
d <- H5Fopen(files[i])
if(center != "c") {
        stand <- d$StandData
} else {
        stand <- d$CenterStandData
}
y <- stand[1,]
points(stand[16,]*1e4/(stand[3,]/aplot1),stand[3,]/aplot1, type="l", lwd=2, col=colors[i])
H5Fclose(d)
}
}

#
# 					Basal area
#

if(length(which(traits=="basalarea") > 0)) {
d <- H5Fopen(files[1])
if(center != "c") {
        stand <- d$StandData
} else {
        stand <- d$CenterStandData
}
y <- stand[1,]
plot(y,stand[14,]*1e4, type="l", ylim=c(0,80),lwd=2, xlab="time (y)", ylab="m2/ha", main="Basal area", col=colors[1])
H5Fclose(d)

for(i in 2:length(files)) {
d <- H5Fopen(files[i])
if(center != "c") {
        stand <- d$StandData
} else {
        stand <- d$CenterStandData
}
y <- stand[1,]
points(y,stand[14,]*1e4, type="l", lwd=2, col=colors[i])
H5Fclose(d)
}
}
#
# 					Basal area at crown base
#

if(length(which(traits=="basalareaatcrownbase") > 0)) {
d <- H5Fopen(files[1])
if(center != "c") {
        stand <- d$StandData
} else {
        stand <- d$CenterStandData
}
y <- stand[1,]
plot(y,stand[15,]*1e4, type="l", ylim=c(0,60),lwd=2, xlab="time (y)", ylab="m2/ha", main="Basal area at crown base", col=colors[1])
H5Fclose(d)

for(i in 2:length(files)) {
d <- H5Fopen(files[i])
if(center != "c") {
        stand <- d$StandData
} else {
        stand <- d$CenterStandData
}
y <- stand[1,]
points(y,stand[15,]*1e4, type="l", lwd=2, col=colors[i])
H5Fclose(d)
}
}

#
# 					Stem volume
#

if(length(which(traits=="stemvolume") > 0)) {
d <- H5Fopen(files[1])
if(center != "c") {
        stand <- d$StandData
} else {
        stand <- d$CenterStandData
}
y <- stand[1,]
plot(y,stand[16,]*1e4, type="l", ylim=c(0,1000),lwd=2, xlab="time (y)", ylab="m3/ha", main="Stem volume", col=colors[1])
H5Fclose(d)

for(i in 2:length(files)) {
d <- H5Fopen(files[i])
if(center != "c") {
        stand <- d$StandData
} else {
        stand <- d$CenterStandData
}
y <- stand[1,]
points(y,stand[16,]*1e4, type="l", lwd=2, col=colors[i])
H5Fclose(d)
}
}

#
# 					Stem volume
#

if(length(which(traits=="stemvolumegrowth") > 0)) {
d <- H5Fopen(files[1])
if(center != "c") {
        stand <- d$StandData
} else {
        stand <- d$CenterStandData
}
y <- stand[1,]
plot(y[1:(length(y)-1)],smooth(diff(stand[16,]))*1e4, type="l", ylim=c(-5,15),lwd=2, xlab="time (y)", ylab="m3/ha per year", main="Stem volume growth", col=colors[1])
H5Fclose(d)

for(i in 2:length(files)) {
d <- H5Fopen(files[i])
if(center != "c") {
        stand <- d$StandData
} else {
        stand <- d$CenterStandData
}
y <- stand[1,]
points(y[1:(length(y)-1)],smooth(diff(stand[16,]))*1e4, type="l", lwd=2, col=colors[i])
H5Fclose(d)
}
}


#
# 					Leaf Area Index LAI
#

if(length(which(traits=="leafareaindex") > 0)) {
d <- H5Fopen(files[1])
if(center != "c") {
        stand <- d$StandData
} else {
        stand <- d$CenterStandData
}
y <- stand[1,]
plot(y,stand[17,], type="l", ylim=c(0,30),lwd=2, xlab="time (y)", ylab="All-sided LAI (m2/m2)", main="Leaf Area Index", col=colors[1])
H5Fclose(d)

for(i in 2:length(files)) {
d <- H5Fopen(files[i])
if(center != "c") {
        stand <- d$StandData
} else {
        stand <- d$CenterStandData
}
y <- stand[1,]
points(y,stand[17,], type="l", lwd=2, col=colors[i])
H5Fclose(d)
}
}

#
# 					Crown ratio
#

if(length(which(traits=="crownratio") > 0)) {
d <- H5Fopen(files[1])
trees <- d$ForestTreeData
if(center != "c") {
        stand <- d$StandData
        mukana <- 1:length(trees[1,,1])
} else {
        stand <- d$CenterStandData
        mukana <- which(in_center(trees[2,,1],trees[3,,1],plot_x,plot_y,dist_x,dist_y))
}
y <- stand[1,]
plot(y,apply(1-trees[11,mukana,]/trees[7,mukana,],2,mean,na.rm=TRUE), type="l", ylim=c(0,1),lwd=2, xlab="time (y)", ylab="Crown ratio", main="Crown ratio", col=colors[1])
H5Fclose(d)

for(i in 2:length(files)) {
d <- H5Fopen(files[i])
trees <- d$ForestTreeData
if(center != "c") {
        stand <- d$StandData
        mukana <- 1:length(trees[1,,1])
} else {
        stand <- d$CenterStandData
        mukana <- which(in_center(trees[2,,1],trees[3,,1],plot_x,plot_y,dist_x,dist_y))
}
y <- stand[1,]
points(y,apply(1-trees[11,mukana,]/trees[7,mukana,],2,mean,na.rm=TRUE), type="l", lwd=2, col=colors[i])
H5Fclose(d)
}
}

#
# 					Dominant height
#

if(length(which(traits=="dominantheight") > 0)) {
d <- H5Fopen(files[1])
trees <- d$ForestTreeData
if(center != "c") {
        stand <- d$StandData
        mukana <- 1:length(trees[1,,1])
} else {
        stand <- d$CenterStandData
        mukana <- which(in_center(trees[2,,1],trees[3,,1],plot_x,plot_y,dist_x,dist_y))
}
#    How many trees correspond to 100 trees/ha
nhd <- as.integer(100*aplot/1e4 + 0.5)


y <- stand[1,]
hd <- 1:length(y)
for(i in 11:length(y)) {
	nt <- length(trees[9,mukana,i])
	o <- order(trees[9,mukana,i], na.last=FALSE)     #order by dbh in year i
	hd[i] <- mean(trees[7,mukana[o[(nt-nhd+1):nt]],i], na.rm=TRUE)   #mean height of nhd thickest trees
}

plot(y[11:length(y)],hd[11:length(y)], type="l", ylim=c(0,30),lwd=2, xlab="time (y)", ylab="Height (m)", main="Dominant height", col=colors[1],xlim = c(0,max(y)))
H5Fclose(d)

for(i in 2:length(files)) {
d <- H5Fopen(files[i])
trees <- d$ForestTreeData
if(center != "c") {
        stand <- d$StandData
        mukana <- 1:length(trees[1,,1])
} else {
        stand <- d$CenterStandData
        mukana <- which(in_center(trees[2,,1],trees[3,,1],plot_x,plot_y,dist_x,dist_y))
}

y <- stand[1,]
hd <- 1:length(y)
for(k in 11:length(y)) {
	nt <- length(trees[9,mukana,k])
	o <- order(trees[9,mukana,k], na.last = FALSE)     #order by dbh in year i
	hd[k] <- mean(trees[7,mukana[o[(nt-nhd+1):nt]],k], na.rm = TRUE)   #mean height of nhd thickest trees
}

points(y[11:length(y)],hd[11:length(y)], type="l", lwd=2, col=colors[i])
H5Fclose(d)
}
}


#
# 					Basal area vs dominant height
#

if(length(which(traits=="basalareavsdominantheight") > 0)) {
d <- H5Fopen(files[1])
trees <- d$ForestTreeData
if(center != "c") {
        stand <- d$StandData
        mukana <- 1:length(trees[1,,1])
} else {
        stand <- d$CenterStandData
        mukana <- which(in_center(trees[2,,1],trees[3,,1],plot_x,plot_y,dist_x,dist_y))
}
#    How many trees correspond to 100 trees/ha
nhd <- as.integer(100*aplot/1e4 + 0.5)


ba <- 1e4*stand[14,]        #Stand basal area
hd <- 1:length(ba)
for(i in 11:length(y)) {
	nt <- length(trees[9,mukana,i])
	o <- order(trees[9,mukana,i], na.last=FALSE)     #order by dbh in year i
	hd[i] <- mean(trees[7,mukana[o[(nt-nhd+1):nt]],i], na.rm=TRUE)   #mean height of nhd thickest trees
}

plot(hd[11:length(y)],ba[11:length(y)], type="l", xlim=c(0,30), ylim=c(0,40),lwd=2, xlab="Dominant height (m)", ylab="Basal area (m2/ha)", main="Basal area vs dominant height", col=colors[1])
H5Fclose(d)

if(length(files) > 1) {
for(i in 2:length(files)) {
d <- H5Fopen(files[i])
trees <- d$ForestTreeData
if(center != "c") {
        stand <- d$StandData
        mukana <- 1:length(trees[1,,1])
} else {
        stand <- d$CenterStandData
        mukana <- which(in_center(trees[2,,1],trees[3,,1],plot_x,plot_y,dist_x,dist_y))
}

ba <- 1e4*stand[14,]           #Stand basal area
hd <- 1:length(ba)
for(k in 11:length(y)) {
	nt <- length(trees[9,mukana,k])
	o <- order(trees[9,mukana,k], na.last = FALSE)     #order by dbh in year i
	hd[k] <- mean(trees[7,mukana[o[(nt-nhd+1):nt]],k], na.rm = TRUE)   #mean height of nhd thickest trees
}

points(hd[11:length(y)], ba[11:length(y)], type="l", lwd=2, col=colors[i])
H5Fclose(d)
}
}
}

dev.off()
}





