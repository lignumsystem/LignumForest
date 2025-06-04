library("rhdf5")

#########        PlotCheck function prints out graphs of tree height or diameter versus
#########        distance from center of the plot in the Quarz window
#########        Call: PlotCheck(infile, age, D)
#########        infile: name of input HDF5 file 
#########        age:    age of trees of the plot
#########        D:      If "D" plot Dbh, otherwise plot height, default is not "D"

PlotCheck <- function(infile, age, D ="H") {

d <- H5Fopen(infile)

stand <- d$StandData
Ntrees <- stand[3,1]

y <- stand[1,]
ymax = max(y, na.rm=TRUE)

pos <- unlist(gregexpr(" ",d$VoxelSpace))   #First line, X Y Z values
plot_x <- as.numeric(substr(d$VoxelSpace,1,pos[1]-1))
plot_y <- as.numeric(substr(d$VoxelSpace,pos[1]+1,pos[2]-1))

pos <- unlist(gregexpr('\n',d$VoxelSpace))
line3 <- substr(d$VoxelSpace,pos[2]+1,pos[3]-1)  #line 3 between 2nd and 3rd '\n'
pos <- unlist(gregexpr(" ",line3))
dist_x <- as.numeric(substr(line3,1,pos[1]-1))   #line3 contains two numbers separated by a space
dist_y <- as.numeric(substr(line3,pos[1]+1,nchar(line3)))

trees <- d$ForestTreeData

main = paste("Age ",as.character(age),"  /   line = least squares fit", sep="")

if(D != "D") {
plot(sqrt((trees[2,,age+1]-plot_x/2)^2+(trees[3,,age+1]-plot_y/2)^2),trees[7,,age+1], xlab="Distance from center (m)", ylab="Tree height (m)", main=main)

abline(lsfit(sqrt((trees[2,,age+1]-plot_x/2)^2+(trees[3,,age+1]-plot_y/2)^2),trees[7,,age+1]))
} else {
plot(sqrt((trees[2,,age+1]-plot_x/2)^2+(trees[3,,age+1]-plot_y/2)^2),100*trees[9,,age+1], xlab="Distance from center (m)", ylab="Dbh (cm)", main=main)

abline(lsfit(sqrt((trees[2,,age+1]-plot_x/2)^2+(trees[3,,age+1]-plot_y/2)^2),100*trees[9,,age+1]))
}

}