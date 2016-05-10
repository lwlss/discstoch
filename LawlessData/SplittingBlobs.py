import os
import numpy
import sys
from PIL import Image

syspath = os.path.dirname(sys.argv[0])
fullpath = os.path.abspath(syspath)

# Splitting out individual timecourse images for each blob
folder="R08C05"
FINALPHOT=35
outputdir=os.path.join(fullpath,"Blobs_"+folder)
if not os.path.exists(outputdir):
    os.makedirs(outputdir)

locs=numpy.loadtxt(folder+"_LOC.txt",dtype=numpy.int,delimiter="\t")
tifs=[os.path.join(fullpath,folder,f) for f in os.listdir(os.path.join(fullpath,folder)) if ".tif" in f]
for j in range(min(len(tifs),FINALPHOT)):
    print("Getting blobs for "+tifs[j])
    im=Image.open(tifs[j])
    for i in range(len(locs)):
        x,y,w,h=locs[i]
        blob=im.crop((x-w/2,y-h/2,x+w/2,y+h/2))
        blob.save(os.path.join(outputdir,'Blob{:04d}_Image{:04d}.jpg'.format(i,j)))
