#### Functions to track individual clonal colony area estimates over time
#### in muQFA
#### Requires OpenCV version 3.0.0 and PIL to be installed 

import cv2
from PIL import Image
import math
import numpy
import os
import matplotlib.pyplot as plt
import muqfatc.imageanalysis as ia

def getBlobs(image,bk,showIms=False,DX=25,DY=25,np=False,apply_filt=True):
    '''Get masks representing microcolony sizes and positions'''
    # Open image as colour
    if np is False:
        im=cv2.imread(image,3)
        siz=im.shape 
    else:
        im=image
        siz=im.shape
    
    if showIms:
        img=Image.fromarray(im,'RGB')
        img.save("1ShowImages.png")

    # Erosion & Dilation Steps 
    if np is False:
        border_im=ia.makeBorder(im,bk,DX=DX,DY=DY)
        if showIms:
            img=Image.fromarray(border_im,'RGB')
            img.save("2ShowImages.png")
        # Convert image to gray scale
        im_gray=cv2.cvtColor(border_im,cv2.COLOR_BGR2GRAY) #im before
        # Generate Canny Edge Map
        canny_im=cv2.Canny(im_gray,10,50,3)
        if showIms:
            img=Image.fromarray(canny_im)
            img.save("3ShowImages.png")
        # Find contours and fill in the area
        kernel = numpy.ones((5,5),numpy.uint8)
        dilate_im=cv2.dilate(canny_im,kernel,iterations=5)
        if showIms:
            img=Image.fromarray(dilate_im)
            img.save("4ShowImages.png")
        erode_im=cv2.erode(dilate_im,kernel,iterations=5)
        if showIms:
            img=Image.fromarray(erode_im)
            img.save("5ShowImages.png")
    else:
        # Convert image to gray scale; use original image for canny edge map! 
        im_gray=cv2.cvtColor(im,cv2.COLOR_BGR2GRAY)
        # Generate Canny Edge Map
        canny_im=cv2.Canny(im_gray,10,50,3)
        if showIms:
            img=Image.fromarray(canny_im)
            img.save("3ShowImages.png")
        # Find contours and fill in the area
        kernel=numpy.matrix([[0,1,1,1,0],[1,1,1,1,1],[1,1,1,1,1],[1,1,1,1,1],[0,1,1,1,0]],numpy.uint8)
        dilate_im=cv2.dilate(canny_im,kernel,iterations=2)
        if showIms:
            img=Image.fromarray(dilate_im)
            img.save("4ShowImages.png")
        erode_im=cv2.erode(dilate_im,kernel,iterations=2)
        if showIms:
            img=Image.fromarray(erode_im)
            img.save("5ShowImages.png")

    # Get colony contours from image 
    final_im=erode_im.copy()
    (contour,_)=cv2.findContours(final_im,cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)
    if showIms:
        contour_im=numpy.zeros(siz,numpy.uint8)
        cv2.drawContours(contour_im,contour,-1,(200,0,0),1)
        img=Image.fromarray(contour_im)
        img.save("6ShowImages.png")
    FCA=[]
    FAA=[]
    # Get the area for for each contour and check for circularity
    for c in contour:
        area=cv2.contourArea(c)
        (x,y,),r=cv2.minEnclosingCircle(c)
        area2=math.pi*(r**2)
        if np is False and apply_filt is True:
            if area/area2 > (0.4+(area*0.0003)): #area/area2 > 0.7:
                FCA.append(c)
                FAA.append(area)
        else:
            if area/area2 > (0.32+(area*0.0003)) and apply_filt is True:
                #this is the most stringent cut-off I would set to get a
                #"very high-quality data set" as this eliminates all previously
                #identified issues; could change the 0.32 to 0.3 or eliminate the
                #area cut-off at this point entirely to get a greater number of
                #growth rates which do contain more experimental errors though! 
                FCA.append(c)
                FAA.append(area)
        if apply_filt is False:
                FCA.append(c)
                FAA.append(area)
    if showIms:
        final_contour_im=numpy.zeros(siz,numpy.uint8)
        cv2.drawContours(final_contour_im,FCA,-1,(0,200,0),1)
        img=Image.fromarray(final_contour_im)
        img.save("7ShowImages.png")
    return(FCA,FAA)

def lineagetimecourse (folders,finalim,bk,DX,DY,fullpath,apply_filt=True,
                       save_pics=True,log=False,show_im=False):
    '''This function generates blob and area time courses for all tiff images in
    the provided folder. No return value.'''
    # Iterating through all folders
    for f in folders:
        print(f)
        finalphoto="img_%09d__000.tif"%finalim
        impath=os.path.join(f,finalphoto)
        blbs,blbs_area=getBlobs(impath,bk,showIms=show_im,DX=DX,DY=DY,
                                np=False,apply_filt=apply_filt)
        # Paint all final blobs to an empty image
        if save_pics:
            black=numpy.zeros((bk.shape[0],bk.shape[1]),numpy.uint8)
            cv2.drawContours(black,blbs,-1,(255,255,255),-1)
            img=Image.fromarray(black)
            img.save("Selected_Colonies_{}.png".format(f))
        # Find all photos available (before "final" photo)
        files=os.listdir(os.path.join(fullpath,f))
        tiffs=[]
        tiffcount=0
        # Iterating through all the files in the folder
        files=sorted(files,key=ia.numericalSort) #to order files numerically 
        for filename in files:
            if ".tif" in filename and tiffcount<(finalim+1):
                tiffs.append(filename)
                tiffcount+=1
        # Folder to save the output
        outputdir=os.path.join(fullpath,"Blobs_{}_".format(finalim)+f)
        if not os.path.exists(outputdir):
            os.makedirs(outputdir)
        # Getting the ROI for all blobs in the cut-off image
        blbno=0
        counter=0
        for blb in blbs:
            x,y,w,h=cv2.boundingRect(blb)
            if blbs_area[blbno]<((bk.shape[0]*bk.shape[1])*0.000215): #300px
                x=int(round(x-((bk.shape[0]*bk.shape[1])*7.2e-06))) #10px
                y=int(round(y-((bk.shape[0]*bk.shape[1])*7.2e-06)))
                w=int(round(w+(2*((bk.shape[0]*bk.shape[1])*7.2e-06)))) #20px
                h=int(round(h+(2*((bk.shape[0]*bk.shape[1])*7.2e-06))))
            # Crop blob image and calculate the contours
            timecourse_area=[]
            timecourse_time=[]
            blob_images=numpy.empty([1,finalim],dtype=object)
            blob_bw_images=numpy.empty([1,finalim],dtype=object)
            for imno in range(0,finalim):
                imname=tiffs[imno]
                impath=os.path.join(fullpath,f,imname)
                imtim=os.path.getmtime(impath)
                currim=cv2.imread(impath,3)
                currim=ia.makeBorder(currim,bk,DX=DX,DY=DY)
                ROI=currim[y:y+h,x:x+w]
                ROI_img=Image.fromarray(ROI)
                colony,colony_area=getBlobs(ROI,bk,showIms=show_im,DX=DX,DY=DY,
                                            np=True,apply_filt=apply_filt)
                # Only save data for single clonal colonies 
                if len(colony) is 1 and apply_filt is True:
                    if imno is 0 and colony_area[0]>((bk.shape[0]*bk.shape[1])*0.000215):
                        #300px
                        print("Too big for a single colony")
                        pass
                    else: 
                        timecourse_area.append(colony_area[0])
                        timecourse_time.append(imtim)
                        blob_images[0,imno]=ROI_img
                if apply_filt is False:
                        if len(colony_area) is 0:
                            timecourse_area.append(0)
                        elif len(colony_area) is 1:
                            timecourse_area.append(colony_area[0])
                        else:
                            timecourse_area.append(sum(colony_area))
                        timecourse_time.append(imtim)
                        blob_images[0,imno]=ROI_img
            # Only save data for full time courses
            if len(timecourse_area) is finalim and apply_filt is True:
                if counter is 0:
                    folder_area=timecourse_area
                    folder_time=timecourse_time
                    folder_images=blob_images[0,:]
                    counter+=1
                else:
                    folder_area=numpy.vstack((folder_area,timecourse_area))
                    folder_time=numpy.vstack((folder_time,timecourse_time))
                    folder_images=numpy.vstack((folder_images,blob_images))
                if save_pics:
                # Time course image for each blob
                    col=finalim/5
                    timecourse_image=Image.new('RGB',(col*w,5*h))
                    blackwhite_image=Image.new('RGB',(col*w,5*h))
                    x=0
                    for i in range(0,5*h, h):
                        for j in range(0,col*w,w):
                            timecourse_image.paste(blob_images[0,x],(j,i))
                    timecourse_image.save(os.path.join(
                        outputdir,'Folder{}_Blob{:04d}_TimeCourse.jpg'.format(f,blbno)))
                    # Growth Curve for each blob
                    time=[(i-timecourse_time[0])/3600 for i in timecourse_time]
                    plt.figure(figsize=(6,4))
                    if log:
                        timecourse_area=[math.log(i) for i in timecourse_area]
                        plt.ylabel('log(Area) (log(px))')
                        output_name='Folder{}_Blob{:04d}_Log_'.format(f,blbno)
                    else:
                        plt.ylabel('Area (px)')
                        output_name='Folder{}_Blob{:04d}_'.format(f,blbno)
                    plt.plot(time,timecourse_area,marker='o',ls='--')
                    plt.xlabel('Time (h)')
                    plt.title('Growth Curve for {} Blob {}'.format(f,blbno))
                    plt.savefig(os.path.join(
                        outputdir,output_name+'GrowthCurve.jpg'))
                    plt.close()
                    growth_curve=Image.open(os.path.join(
                        outputdir,output_name+'GrowthCurve.jpg'.format(f,blbno)))
                    # Growth Curves + Time course Images for each blob
                    gcw,gch=growth_curve.size
                    width,height=timecourse_image.size
                    width+=1
                    height+=1
                    final_image=Image.new('RGB', (gcw,gch+(height)), "white")
                    final_image.paste(timecourse_image,((gcw-(width))/2,0))
                    final_image.paste(growth_curve,(0,height))
                    final_image.save(os.path.join(
                        outputdir,output_name+'TCGC.jpg'))
            # Time courses for when there is no filter
            if apply_filt is False: #and sum(timecourse_area)>0
                if counter is 0:
                    folder_area=timecourse_area
                    folder_time=timecourse_time
                    folder_images=blob_images[0,:]
                    counter+=1
                else:
                    folder_area=numpy.vstack((folder_area,timecourse_area))
                    folder_time=numpy.vstack((folder_time,timecourse_time))
                    folder_images=numpy.vstack((folder_images,blob_images))
                # Time course image for each blob
                if sav_pics: 
                    col=finalim/5
                    timecourse_image=Image.new('RGB',(col*w,5*h))
                    blackwhite_image=Image.new('RGB',(col*w,5*h))
                    x=0
                    for i in range(0,5*h, h):
                        for j in range(0,col*w,w):
                            timecourse_image.paste(blob_images[0,x],(j,i))
                            if display_image:
                                blackwhite_image.paste(blob_bw_images[0,x],(j,i))
                            x+=1
                    timecourse_image.save(os.path.join(
                        outputdir,'Folder{}_Blob{:04d}_TimeCourse.jpg'.format(f,blbno)))
                    if display_image:
                        blackwhite_image.save(os.path.join(
                            outputdir,'Folder{}_Blob{:04d}_TimeCourse_BW.jpg'.format(f,blbno)))  
                    # Growth Curve for each blob
                    log=False
                    time=[(i-timecourse_time[0])/3600 for i in timecourse_time]
                    plt.figure(figsize=(6,4))
                    if log:
                        timecourse_area=[math.log(i) for i in timecourse_area]
                        plt.ylabel('log(Area) (log(px))')
                        output_name='Folder{}_Blob{:04d}_Log_'.format(f,blbno)
                    else:
                        plt.ylabel('Area (px)')
                        output_name='Folder{}_Blob{:04d}_'.format(f,blbno)
                    plt.plot(time,timecourse_area,marker='o',ls='--')
                    plt.xlabel('Time (h)')
                    plt.title('Growth Curve for {} Blob {}'.format(f,blbno))
                    plt.savefig(os.path.join(
                        outputdir,output_name+'GrowthCurve.jpg'))
                    plt.close()
                    growth_curve=Image.open(os.path.join(
                        outputdir,output_name+'GrowthCurve.jpg'.format(f,blbno)))
                    # Growth Curves + Time course Images for each blob
                    gcw,gch=growth_curve.size
                    width,height=timecourse_image.size
                    width+=1
                    height+=1
                    final_image=Image.new('RGB', (gcw,gch+(height)), "white")
                    final_image.paste(timecourse_image,((gcw-(width))/2,0))
                    final_image.paste(growth_curve,(0,height))
                    final_image.save(os.path.join(
                        outputdir,output_name+'TCGC.jpg'))
            blbno+=1
        numpy.savetxt(f+"_AREA_{}.txt".format(finalim),folder_area,delimiter="\t")
        numpy.savetxt(f+"_TIME_{}.txt".format(finalim),folder_time,delimiter="\t")
     
