import os
import numpy
import sys
from PIL import Image
import cv2.cv as cv

import re
# Solution of http://stackoverflow.com/questions/12093940/reading-files-in-a-particular-order-in-python
numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts

syspath = os.path.dirname(sys.argv[0])
fullpath = os.path.abspath(syspath)

# Splitting out individual timecourse images for each blob
folder="R08C05"
#folder="R03C03"
FINALPHOT=35
outputdir=os.path.join(fullpath,"Blobs_"+folder)
if not os.path.exists(outputdir):
    os.makedirs(outputdir)

locs=numpy.loadtxt(folder+"_LOC.txt", dtype=numpy.int, delimiter="\t")

tifs=[os.path.join(fullpath,folder,f) for f in os.listdir(os.path.join(fullpath,folder)) if ".tif" in f]
tifs=sorted(tifs,key=numericalSort) #makes sure files are sorted numerically
for j in range(min(len(tifs),FINALPHOT)):
    print("Getting blobs for "+tifs[j])
    im=Image.open(tifs[j])
    imw,imh=im.size
    for i in range(len(locs)): 
        x,y,w,h=locs[i]
        blob=im.crop((x,y,x+w,y+h))
        blob.save(os.path.join(outputdir,'Blob{:04d}_Image{:04d}.jpg'.format(i,j)))
