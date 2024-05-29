## h5file: Lignum HDF5 data file after simulation
## densityfile: Ilvessalo.txt
## area: VoxelSpace area (see -voxelspace from command line) 
library(rhdf5)
PlotDensity<-function(h5file,densityfile,area){
    h5data<-H5Fopen(h5file)
    ftdata<-read.table(densityfile,header=TRUE,sep='')
     ## m^2->Hectars
    aplot<-area/1e4
    t1 <- h5data$StandData[1,]
    sdensity <- h5data$StandData[3,]/aplot
    plot(t1,sdensity,type="l",lty=1, lwd=2, xlab="Time (y)", ylab="No. trees / ha",ylim=c(0,17000),
         main="Stand density\n(Vuokila/Ilvessalo)")
    t2<-ftdata$Year
    mt<-ftdata$MT
    vt<-ftdata$VT
    ct<-ftdata$CT
    lines(t2,mt,col='green',lty=1,lwd=2)
    lines(t2,vt,col='red',lty=1,lwd=2)
    lines(t2,ct,col='brown',lty=1,lwd=2)
    legend('bottomleft',inset=0.05,c("Lignum","MT(ksto-mt.dat)","VT","CT"),col=c('black','green','red','brown'),lty=1,
           title="Forest type")
}
    
