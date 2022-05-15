#########        ExtractXML function reads Lignum trees stored in a HDF5 file and stores as XML files
#########        The function stores the XML files of shortest, median and longest (at final
#########        time) trees at given times.
#########        The names of the XML files are xmfile_TreeNN_year.xml
#########        Call: ForestPlot(datafile, xmlfile, years)
#########        datafile: name of input HDF5 file Containing tree information
#########        xmlfile:  HDF5 file containing the XML data
#########        years     vector of times for storing the XML files.


#The datafile must be stored iin the HDF5 File Format
#(https://en.wikipedia.org/wiki/Hierarchical_Data_Format).

#The datafile must contain the following dataset: 
#              /  ForestTreeData H5I_DATASET  FLOAT 51 x 599 x 21
#In this example 599 is the numberof trees on plot and 21 is the number times (years) at which
#results are stored

#The xmlfile must contain the XMLs of trees as follows:
#           group     name       otype dclass dim
#           /TreeXML/5  Tree_85 H5I_DATASET STRING   1
#This example shows the xml data of tree #85 at time (year) 5

#In order to run ForestPlot in R a suitable package is needed. ForestPlot has been tested using
#Biomanager package (https://www.rdocumentation.org/packages/BiocManager/versions/1.30.17):
#install.packages("BiocManager")
#BiocManager::install("rhdf5")



#----------------------------------------------------------------------------

ExtractXML <- function(datafile, xmlfile, years) {

# longest and shortest and median trees

d <- H5Fopen(datafile)

datayrs <- d$StandData[1,]         #Tree information is storad at datayrs years
maxdyrs <- datayrs[length(datayrs)]

ym <- max(years)

if(ym > maxdyrs) { print("years contains a too large year"); return(-1)}

k <- which(datayrs == ym)
if(length(k) == 0) {print("Invalid year in years"); return(-1)}

h <- d$ForestTreeData[7,,k]
mh <- max(d$ForestTreeData[7,,k],na.rm=TRUE)
largest <- which(h>0.999*mh)[1]              #Takes first index that fulfills the condition

mh <- min(d$ForestTreeData[7,,k],na.rm=TRUE)
smallest <- which(h<1.001*mh)[1]

mh <- median(d$ForestTreeData[7,,k],na.rm=TRUE)
med <- which(h<1.001*mh&h>0.999*mh)[1]

lid <- d$ForestTreeData[1,largest,1]        #Tree id
sid <- d$ForestTreeData[1,smallest,1]
mid <- d$ForestTreeData[1,med,1]


#Write the XML strings to files
x <- H5Fopen(xmlfile)

for(i in 1:length(years)) {
	t <- paste("TreeXML/",as.character(years[i]),"/Tree_",as.character(sid),sep="")
	if(H5Lexists(x,t)) {
		tree <- h5read(x,t)
		oname <- paste(xmlfile,"_Tree",as.character(sid),"_",as.character(years[i]),".xml",sep="")
		cat(tree,file=oname)
		}
	t <- paste("TreeXML/",as.character(years[i]),"/Tree_",as.character(mid),sep="")
	if(H5Lexists(x,t)) {
		tree <- h5read(x,t)
		oname <- paste(xmlfile,"_Tree",as.character(mid),"_",as.character(years[i]),".xml",sep="")
		cat(tree,file=oname)
		}
	t <- paste("TreeXML/",as.character(years[i]),"/Tree_",as.character(lid),sep="")
	if(H5Lexists(x,t)) {
		tree <- h5read(x,t)
		oname <- paste(xmlfile,"_Tree",as.character(lid),"_",as.character(years[i]),".xml",sep="")
		cat(tree,file=oname)
		}			
	}
}

