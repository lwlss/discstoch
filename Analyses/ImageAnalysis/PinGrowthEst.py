############ GETTING THE TOTAL AREA OF EACH PIN OVER TIME #################

import muqfatc.imageanalysis as ia
import muqfatc.pingrowth as pinia
# The myqfatc package requires cv2 (version 3.0.0) to be intalled 
import sys, os
import numpy

# Find all available time courses in the current working directory
syspath = os.path.dirname(sys.argv[0]) 
fullpath = os.path.abspath(syspath) 
allobj=os.listdir(fullpath)
# Getting all the folder names
myfolders=[]
for f in allobj:
    if len(f)==6 and "." not in f and f[0]=="R" and f[3]=="C" and f!='R03C03':
        myfolders.append(f)

# Getting time and area for pin growth estimates for 90 observations
area,time=pinia.pintimecourse(myfolders,90,fullpath)

# Write Output to File
numpy.savetxt("PopulationArea.txt",area,delimiter="\t")
numpy.savetxt("PopulationTime.txt",time,delimiter="\t")
text_file = open("PopulationFolders.txt", "w")
for item in myfolders:
    text_file.write("%s\n" % item)
text_file.close()
