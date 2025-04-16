library(rhdf5)


ExtractFunctions<-function(h5file,outdir){
    if (!dir.exists(outdir)){
        dir.create(outdir,recursive=TRUE)
    }
    data <- H5Fopen(h5file)

    namels<-names(data$AllFunctionFiles)
    for (name in namels){
        fndata<-data$AllFunctions[[name]]
        write.table(t(fndata),paste(outdir,'/',name,sep=""))
    }
}

ExtractFunctionFiles<-function(h5file,outdir){
    if (!dir.exists(outdir)){
        dir.create(outdir,recursive=TRUE)
    }
    data <- H5Fopen(h5file)

    namels<-names(data$AllFunctionFiles)
    for (name in namels){
        filedata<-data$AllFunctionFiles[[name]]
        writeLines(filedata,paste(outdir,'/',name,sep=""))
    }
}

ExtractParameterFiles<-function(h5file,outdir){
    if (!dir.exists(outdir)){
        dir.create(outdir,recursive=TRUE)
    }
    data <- H5Fopen(h5file)

    namels<-names(data$AllParameterFiles)
    for (name in namels){
        filedata<-data$AllParameterFiles[[name]]
        writeLines(filedata,paste(outdir,'/',name,sep=""))
    }
}

ExtractMetaFiles<-function(h5file,outdir){
    if (!dir.exists(outdir)){
        dir.create(outdir,recursive=TRUE)
    }
    data <- H5Fopen(h5file)

    namels<-names(data$AllMetaFiles)
    for (name in namels){
        filedata<-data$AllMetaFiles[[name]]
        writeLines(filedata,paste(outdir,'/',name,sep=""))
    }
}

ExtractCommandLine<-function(h5file,outdir){
    if (!dir.exists(outdir)){
        dir.create(outdir,recursive=TRUE)
    }
    data <- H5Fopen(h5file)

    commandline<-data$CommandLine
    writeLines(commandline,paste(outdir,'/',"commandline.sh",sep=""))
}

## Extract Meta files, parameter and  function files and command line
## from HDF5 result file h5file to outdir directory
ExtractAllFiles<-function(h5file,outdir){
    if (!dir.exists(outdir)){
        dir.create(outdir,recursive=TRUE)
    }
    ExtractCommandLine(h5file,outdir)
    ExtractMetaFiles(h5file,outdir)
    ExtractParameterFiles(h5file,outdir)
    ExtractFunctionFiles(h5file,outdir)
}
