#########        ExtractXMLYear function reads Lignum trees stored in a HDF5 file and stores as XML files
#########        The function stores the XML files of shortest, median and longest trees at a given
#########        one year. NOTE that the  difference to function ExtractXML is that this does only
#########        for one year and evaluates shortest, median and longest trees this year.
#########        The names of the XML files are xmfile_TreeNN_Yyear.xml
#########        NOTE: year in the file name is written as Yyear to distinquish from output of ExtractXML
#########        Call: ForestPlot(datafile, xmlfile, year)
#########        datafile: name of input HDF5 file Containing tree information
#########        xmlfile:  HDF5 file containing the XML data
#########        year      the year of evaluation.


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
#library("rhdf5")



#----------------------------------------------------------------------------

ExtractXMLYear <- function(datafile, xmlfile, year) {

# longest and shortest and median trees

d <- H5Fopen(datafile)

datayrs <- d$StandData[1,]         
maxdyrs <- max(datayrs, na.rm=TRUE)

if(year > maxdyrs) { print("year is too large"); return(-1)}

#h <- d$ForestTreeData[7,,year+1]
mh <- max(d$ForestTreeData[7,,year+1],na.rm=TRUE)
largest <- which(h>0.999*mh)[1]              #Takes first index that fulfills the condition

mh <- min(d$ForestTreeData[7,,year+1],na.rm=TRUE)
smallest <- which(h<1.001*mh)[1]

mh <- median(d$ForestTreeData[7,,year+1],na.rm=TRUE)
med <- which(h<1.001*mh&h>0.99*mh)[1]

lid <- d$ForestTreeData[1,largest,1]        #Tree id
sid <- d$ForestTreeData[1,smallest,1]
mid <- d$ForestTreeData[1,med,1]


#Write the XML strings to files
x <- H5Fopen(xmlfile)

	t <- paste("TreeXML/",as.character(year),"/Tree_",as.character(sid),sep="")
	if(H5Lexists(x,t)) {
		tree <- h5read(x,t)
		oname <- paste(xmlfile,"_Tree",as.character(sid),"_Y",as.character(year),".xml",sep="")
		cat(tree,file=oname)
		}
	t <- paste("TreeXML/",as.character(year),"/Tree_",as.character(mid),sep="")
	if(H5Lexists(x,t)) {
		tree <- h5read(x,t)
		oname <- paste(xmlfile,"_Tree",as.character(mid),"_Y",as.character(year),".xml",sep="")
		cat(tree,file=oname)
		}
	t <- paste("TreeXML/",as.character(year),"/Tree_",as.character(lid),sep="")
	if(H5Lexists(x,t)) {
		tree <- h5read(x,t)
		oname <- paste(xmlfile,"_Tree",as.character(lid),"_Y",as.character(year),".xml",sep="")
		cat(tree,file=oname)
		}			
	}

