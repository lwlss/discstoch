#### Functions required for estimating population growth in each pin of muQFA
#### Requires OpenCV version 3.0.0 to be installed

import cv2
import numpy
import os
import muqfatc.imageanalysis as ia

def totalArea(image,col):
    '''Calculates total area of colonies darker than the input col
    in that image'''
    im=cv2.imread(image,3)
    lower=numpy.array((0, 0, 0),dtype="uint8")
    upper=numpy.array((col, col, col),dtype="uint8")
    mask=cv2.inRange(im,lower,upper)
    celldens=cv2.countNonZero(mask)
    return(celldens)

def pintimecourse(folders,final_im,path):
    """Iterates through tiff images up to the final image specified for each of
    the folders and returns area and and time matrices."""
    all_area=numpy.zeros((len(folders), final_im))
    all_time=numpy.zeros((len(folders), final_im))
    foldercount=0
    for f in folders:
        print(f)
        # Find all photos in that folder
        files=os.listdir(os.path.join(path,f))
        # Sort files numerically 
        files=sorted(files,key=ia.numericalSort)
        # Only look at tiff files in that folder 
        tiffs=[]
        tiffcount=0
        for filename in files:
            if ".tif" in filename and tiffcount<final_im:
                tiffs.append(filename)
                tiffcount+=1
        # Iterate through tiff files
        area=[]
        time=[]
        for filename in tiffs:
            # Get total area of yeast in each image
            impath=os.path.join(f,filename)
            # This colour or darker sets yeast colonies appart from background
            col=115
            tiffarea=totalArea(impath,col)
            # Get time at which image was generated 
            imtim=os.path.getmtime(impath)
            area.append(tiffarea)
            time.append(imtim)

        # Store data for all folders in area and time matrix
        time=[(i-time[0])/3600 for i in time]
        all_time[foldercount,:]=time
        all_area[foldercount,:]=area
        foldercount+=1
    return(all_area,all_time)
