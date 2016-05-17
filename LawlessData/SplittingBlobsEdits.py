import os
import numpy
import sys
from PIL import Image
import cv2.cv as cv
import matplotlib.pyplot as plt

import re

numbers = re.compile(r'(\d+)')

def numericalSort(value):
    '''To sort file names numerically'''
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts

def revert(x):
    '''For reverting the location names'''
    NewCol=int(re.findall('\d+',x)[1])
    Col=str(1-NewCol+24)
    if Col<10:
       Col='0'.join(Col)
    Newx=[x[:-2],Col]
    Newx=''.join(Newx)
    return(Newx)

syspath = os.path.dirname(sys.argv[0])
fullpath = os.path.abspath(syspath)

# Data 
folders=("R03C04","R03C05","R03C06","R03C07","R03C08","R03C09","R03C10",
         "R04C03","R04C04","R04C05","R04C06","R04C07","R04C08","R04C09","R04C10",
         "R05C03","R05C04","R05C05","R05C06","R05C07","R05C08","R05C09","R05C10"
         "R06C03","R06C04","R06C05","R06C06","R06C07","R06C08","R06C09","R06C10",
         "R07C03","R07C04","R07C05","R07C06","R07C07","R07C08","R07C09","R07C10",
         "R08C03","R08C04","R08C05","R08C06","R08C07","R08C08","R08C09","R08C10",
         "R09C03","R09C04","R09C05","R09C06","R09C07","R09C08","R09C09","R09C10",
         "R10C03","R10C04","R10C05","R10C06","R10C07","R10C08","R10C09","R10C10")
weird_folders=("R05C10", "R10C05")
#### look into these folders "R05C10", 'R10C05'
#### something doesn't quite match up; need to go back to the original script that generates the data
area=numpy.loadtxt("Lawless_area.txt", dtype=numpy.int)
data=numpy.loadtxt("Lawless_data.txt", dtype=numpy.str)
    
# Splitting out individual timecourse images for each blob
for i in range(len(folders)):
    folder = weird_folders[i]
    print(folder)
    revert_folder=revert(folder)
    print(revert_folder)
    folder_indices=numpy.where(data[:,1]=='"{}"'.format(revert_folder))
    FINALPHOT=35
    outputdir=os.path.join(fullpath,"Blobs_"+folder)
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
    locs=numpy.loadtxt(folder+"_LOC.txt", dtype=numpy.int, delimiter="\t")
    tifs=[os.path.join(fullpath,folder,f) for f in os.listdir(os.path.join(fullpath,folder)) if ".tif" in f]
    tifs=sorted(tifs,key=numericalSort)
    images=numpy.empty([FINALPHOT,len(locs)],dtype=object) # to store all images of blobs in a matrix
    for j in range(min(len(tifs),FINALPHOT)):
        print("Getting blobs for "+tifs[j])
        im=Image.open(tifs[j]) # plate image for each time point; there is 35 in total (entire well) 
        imw,imh=im.size
        for l in range(len(locs)): # iterating through all identified blobs
            # Blob Image 
            x,y,w,h=locs[l]
            blob=im.crop((x,y,x+w,y+h))
            if blob.size[0]<50:
                factor=50/blob.size[0]
                blob=blob.resize((factor*blob.size[0],factor*blob.size[1]))                
            images[j][l]=blob
            #blob.save(os.path.join(outputdir,'Blob{:04d}_Image{:04d}.jpg'.format(i,j)))
    for k in range(len(locs)):
        if folder in weird_folders:
            if k==0:
                print("Missing Growth Curve")
                blobs = images[:,k]
                width, height = blobs[-1].size
                width+=1
                height+=1
                tcimage=Image.new('RGB',(7*width,5*height))
                x=0
                for h in xrange(0, 5*height,height):
                    for w in xrange(0, 7*width,width):
                        tcimage.paste(blobs[x],(w,h))
                        x=x+1
                tcimage.save(os.path.join(outputdir,'Folder{}_Blob{:04d}.jpg'.format(folder,k)))
            else: 
                index=folder_indices[0][k-1]
                print(index)
                plt.figure(k-1)
                plt.plot(area[index,:],marker='o',ls='--')
                plt.ylabel('Area')
                plt.xlabel('Time')
                plt.title('Growth Curve {}'.format(index))
                plt.savefig('Folder{}_Blob{:04d}.jpg'.format(folder.replace('"', ""),k-1))
                plt.close()
                growth_curve=Image.open('Folder{}_Blob{:04d}.jpg'.format(folder.replace('"', ""),k-1))
                blobs = images[:,k]
                width, height = blobs[-1].size
                width+=1
                height+=1
                tcimage=Image.new('RGB',(7*width,5*height))
                x=0
                for h in xrange(0, 5*height,height):
                    for w in xrange(0, 7*width,width):
                        tcimage.paste(blobs[x],(w,h))
                        x=x+1
                tcimage.save(os.path.join(outputdir,'Folder{}_Blob{:04d}.jpg'.format(folder,k)))
                # Plotting the growth curve together with the time course
                gcw, gch = growth_curve.size
                final_image=Image.new('RGB', (gcw,gch+(5*height)), "white")
                final_image.paste(tcimage,((gcw-(7*width))/2,0))
                final_image.paste(growth_curve,(0,5*height))
                final_image.save(os.path.join(outputdir,'Folder{}_Blob{:04d}_TC.jpg'.format(folder,k)))
        else:
            # Growth Curve Image
            index=folder_indices[0][k]
            print(index)
            plt.figure(k)
            plt.plot(area[index,:],marker='o',ls='--')
            plt.ylabel('Area')
            plt.xlabel('Time')
            plt.title('Growth Curve {}'.format(index))
            plt.savefig('Folder{}_Blob{:04d}.jpg'.format(folder.replace('"', ""),k))
            plt.close()
            growth_curve=Image.open('Folder{}_Blob{:04d}.jpg'.format(folder.replace('"', ""),k))
            blobs = images[:,k]
            width, height = blobs[-1].size
            width+=1
            height+=1
            tcimage=Image.new('RGB',(7*width,5*height))
            x=0
            for h in xrange(0, 5*height,height):
                for w in xrange(0, 7*width,width):
                    tcimage.paste(blobs[x],(w,h))
                    x=x+1
            tcimage.save(os.path.join(outputdir,'Folder{}_Blob{:04d}.jpg'.format(folder,k)))
            # Plotting the growth curve together with the time course
            gcw, gch = growth_curve.size
            final_image=Image.new('RGB', (gcw,gch+(5*height)), "white")
            final_image.paste(tcimage,((gcw-(7*width))/2,0))
            final_image.paste(growth_curve,(0,5*height))
            final_image.save(os.path.join(outputdir,'Folder{}_Blob{:04d}_TC.jpg'.format(folder,k)))

               
##area=numpy.loadtxt("Lawless_area.txt", dtype=numpy.int)
##data=numpy.loadtxt("Lawless_data.txt", dtype=numpy.str)
###for i in range(numpy.shape(area)[0]):
##i=1
##strain=data[i,0]
##folder=data[i,1]
##folder_indices=numpy.where(data[:,1]==folder)
##blob_number=numpy.where(i == folder_indices[0])[0][0]
##plt.plot(area[i,:],marker='o',ls='--')
##plt.ylabel('Area')
##plt.xlabel('Time')
##plt.title('Folder {} Blob {:04d}'.format(folder,blob_number))
##plt.savefig('Folder{}_Blob{:04d}.jpg'.format(folder.replace('"', ""),blob_number))
