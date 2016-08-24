############ GETTING THE TOTAL AREA OF EACH COLONY OVER TIME #################

import muqfatc.imageanalysis as ia
import muqfatc.lineagegrowth as linia
# The myqfatc package requires cv2 (version 3.0.0) to be intalled 
import sys, os
import numpy

##import cv2
##from PIL import Image
##import matplotlib.pyplot as plt

# Final Photo
FINALPHOT=20 # must be a number divisible by 5 when save_pics = True
# Border size around the image
DX, DY = 25,25

# Filter for only meaningful clonal colony growth curves?
apply_filt=True
# Save blob images?
save_pics=True
# Show dilation erosion images?
show_im=False
# If save_pics==True, should growth curves be plotted on the log scale?
log=False

# Find all available time courses
syspath = os.path.dirname(sys.argv[0]) #current working directory
fullpath = os.path.abspath(syspath) 
allobj=os.listdir(fullpath) #list of all objects in the directory

# Getting all folder names

myfolders=[]
for f in allobj:
    if len(f)==6 and "." not in f and f[0]=="R" and f[3]=="C":
        myfolders.append(f)

# Generating a grey background of median pixels of all first images
bk=ia.makeBackground(myfolders)   

# Generate time course estimates (as output matrices saved in text files)
# and time course images (saved in respective output folders). 
linia.lineagetimecourse(myfolders,FINALPHOT,bk,DX,DY,fullpath,apply_filt=True,
                       save_pics=True,log=False,show_im=False)


