#### Create pdf results from the HDF5 file and move the result files to outdir.
#### Note that the required packages (filesstrings,rhdf5 and ineq)  will be installed and the installation may
#### ask for the R package repository location (61 is Sweden).
#### Usage:
#### source("/path/to/LignumForest/ResultAnalysis/ForestPlotAndConfig.R",chdir=TRUE)
#### ForestPlotAndConfig(/path/to/infile,pick,GYdata="/path/to/LignumForest/ResultAnalysis/",outdir="/path/to/outdir")
if (!require("filesstrings",quietly=TRUE)){
    print("Installing package filesstrings")
    install.packages("filesstrings")
}
if (!require("rhdf5",quietly=TRUE)){
    print("Installing package rhdf5")
    install.packages("rhdf5")
}
if (!require("ineq",quietly=TRUE)){
    print("Installing package ineq")
    install.packages("ineq")
}
library(filesstrings)
library(rhdf5)
source("ForestPlotFunction.R")
source("ExtractSimulationConfig.R")

###ForestPlot results for both the forest plot ("a") and its center area ("c").
###Optionally:
###1) Move the two ForestPlot result files to the outdir. Assume file suffixes used in ForestPlot.
###2) Extract simulation configuration files to Config subdirectory in the outdir.
###3) Move simulation HDF5 files to outdir
###Parameters:
###infile: HDF5 file with simulation results and configuration
###pick: select every n:th tree for the result plots
###GYdata: "ResultAnalysis/" (default), Growth and Yield data files directory
###outdir: FALSE (default) or directory name. If directory name string move the pdf result files,
###        extract and move simulation configuration files and move HDF5 files to  outdir directory
ForestPlotAndConfig<-function(infile,pick=1, GYdata = "ResultAnalysis/",outdir=FALSE){
    ForestPlot(infile,pick,GYdata,center="a")
    ForestPlot(infile,pick,GYdata,center="c")
    if (is.character(outdir)){
        print("Move PDF files")
        ###ForestPlot file suffixes
        cat(paste(infile,'.pdf',sep=''),'->',paste(outdir,'/',infile,'.pdf'),"\n")
        cat(paste(infile,'_c.pdf',sep=''),'->',paste(outdir,'/',infile,'_c.pdf'),"\n")
        file.move(paste(infile,'.pdf',sep=''),outdir,overwrite=TRUE)
        file.move(paste(infile,'_c.pdf',sep=''),outdir,overwrite=TRUE)
        cat("Extract configuration",infile,'->',paste(outdir,'/','Config',sep=''),"\n")
        ExtractAllFiles(infile,paste(outdir,"/","Config",sep=''))
        print("Move HDF5 files")
        cat(infile,'->',outdir,"\n")
        ###File prefix for trees in XML format  
        cat(paste('TreesXML_',infile,sep=''),'->',outdir,"\n")
        file.move(infile,outdir,overwrite=TRUE)
        file.move(paste('TreesXML_',infile,sep=''),outdir,overwrite=TRUE)
    }
}
