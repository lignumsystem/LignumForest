#### Usage:
#### source("/path/to/LignumForest/ResultAnalysis/ForestPlotAndConfig.R",chdir=TRUE)
#### ForestPlotAndConfig(infile,pick,GYdata="/path/to/LignumForest/ResultAnalysis/",outdir="/path/to/outdir")

library(rhdf5)
source("ForestPlotFunction.R")
source("ExtractSimulationConfig.R")

###Create pdf files with ForestPlot for both the forest plot and the center area
###Extract configuration files to Config subdirectory
ForestPlotAndConfig<-function(infile,pick, GYdata = "",outdir="."){
    h5f<-H5Fopen(infile)
    ### Retrieve border forest width
    s<-unlist(h5f$VoxelSpace)
    s<-strsplit(s,split="\\s+")
    s<-sapply(s,as.double)
    b1<-s[length(s)-1]
    b2<-s[length(s)]
    ### Retrieve the actual final forest plot area
    ls<-h5f$VoxelSpaceSizes
    m<-t(matrix(unlist(ls),nrow=12))
    dims<-dim(m)
    rows<-dims[1]
    cols<-dims[2]
    area<-m[rows,cols]
    width<-m[rows,cols-3]
    length<-m[rows,cols-2]
    width_c<-width-2.0*b1
    length_c<-length-2.0*b2
    area_c<-width_c*length_c
    cat("PlotWidth(X)",width,"PlotLength(Y)",length,"PlotArea",area,"\n")
    cat("BorderWidth(X)",b1,"BorderLength(Y)",b2,"\n")
    cat("CentreWidth(X)",width_c,"CentreLength(Y)",length_c,"CenterArea",area_c,"\n")
    h5closeAll()
    ForestPlot(basename(infile),area,pick,GYdata,center="a",outdir)
    ForestPlot(basename(infile),area_c,pick,GYdata,center="c",outdir)
    ExtractAllFiles(infile,paste(outdir,"/","Config",sep=""))
}
