## @package generatelocations
## @brief Generate tree positions
##
## Generate tree positions and forest plot file that can be used in LignumForest.
## Given the forest plot size, border forest width and the forest density
## generate tree positions in a 2D regular rectangular lattice.
## Other geometries can be implemented as needed.
##
## Required software can be installed with MacPorts and `pip` in Python3
##
##       sudo port install python313
##
## Create Python virtual environment and activate it
##
##       /opt/local/bin/python3 -m venv scipyenv
##       scipyenv/bin/activate
##
## Install Python packages required
##
##       pip install --upgrade pip
##       pip install pandas
##       pip install matplotlib
##       pip install numpy
##       pip install scipy
## 
## For the usage type
## 
##       python3 generatelocations.py -h
## 
import argparse
import math
import csv
import pandas as pd
import matplotlib.pyplot as plt

if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Generate tree locations as a regular 2D square lattice")

    parser.add_argument('-s','--side',type=int,required=True,dest="s",help="Length (m) of one side of a square forest plot")
    parser.add_argument('-b','--border',type=int,required=True,dest="b",help="Border forest width (m) in the square forest plot")
    parser.add_argument('-d','--density',type=int,required=True,dest="d",help="Stand density trees/ha")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-p','--plot',action='store_true',help="Plot locations with matplotlib")
    group.add_argument('-o','--output', type=str,dest="o",help="Output file name")
    args = parser.parse_args()
        
    side = args.s
    border = args.b 
    density = args.d

    #Tree distance in a regular grid when density is given for 1 ha stand
    #One side of the 1ha square forest plot is 100m and the tree distance along that side
    tree_distance = 100.0/math.sqrt(density)

    #Create the list of locations in a regular lattice pattern
    locations_ls =[]
    x=0
    y=0
    while (x+tree_distance) < side:
        x = x+tree_distance
        y = 0
        while (y+tree_distance) < side:
            y = y+tree_distance
            locations_ls.append((x,y))
    df_locations = pd.DataFrame(locations_ls)
    #One can adjust marker size with the option 's'.
    #For example: plt.scatter(df[0],df[1],s=0.1)
    if args.plot:
        plt.scatter(df_locations[0],df_locations[1])
        plt.show()
        quit()
    
    if args.o:
        print("Writing tree location file:", args.o)
        file_name = args.o
        ls = []
        stand_comment_ls=["The upper right corner of stand (left corner =(0.0)), border_x, border_y"]
        stand_ls =[side,side,border,border]
        tree_locations_comment_ls =["Locations of trees, each (x,y) coordinate on its own line"]
        ls.append(stand_comment_ls)
        ls.append(stand_ls)
        ls.append(tree_locations_comment_ls)
        df_stand = pd.DataFrame(ls)
        df_stand=pd.concat([df_stand,df_locations])
        df_stand.to_csv(file_name,sep=' ',header=False,index=False)
        print("Done")
