############### IMAGE ANALYSIS IN CV2 - RESTRUCTURED ###########################

import cv2
import PIL, sys, os, numpy, random, math, re
from PIL import Image
import matplotlib.pyplot as plt

############### Functions ######################################################

def numericalSort(value):
    '''To sort file names numerically'''
    numbers = re.compile(r'(\d+)')
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts

def makeBackground(folders):
    ''' Open stack of first images from each spot and get median pixel
    intensity to get a good background lighting map '''
    example=cv2.imread(os.path.join(folders[0],"img_000000000__000.tif"),3)
    images=numpy.zeros(example.shape+(len(folders),),numpy.uint8)
    for f in xrange(len(folders)):
        currim=cv2.imread(os.path.join(folders[f],"img_000000000__000.tif"),3)
        images[:,:,:,f]=currim
    backg=numpy.median(images,axis=3)
    new_image=numpy.zeros(example.shape,numpy.uint8)
    new_image[:]=backg
    return(new_image)

def getBlobs(image,bk,showIms=False,DX=25,DY=25,np=False):
    '''Get masks representing microcolony sizes and positions'''
    # Open image as colour
    if np is False:
        im=cv2.imread(image,3)
        siz=im.shape #row, column = height, width
    else:
        im=image
        siz=im.shape #row, column = height, width
        
    if showIms:
        img=Image.fromarray(im,'RGB')
        img.save("1ShowImages.png")

    if np is False:
        # Add border around image (based on background)
        border_im=cv2.resize(bk,(2*DY+siz[1],2*DX+siz[0]),
                             interpolation = cv2.INTER_CUBIC)
        new_siz=border_im.shape
        border_im[DY:new_siz[0]-DY,DX:new_siz[1]-DX]=im
        if showIms:
            img=Image.fromarray(border_im,'RGB')
            img.save("2ShowImages.png")
        # Convert image to gray scale; use original image for canny edge map! 
        im_gray=cv2.cvtColor(im,cv2.COLOR_BGR2GRAY)
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
    final_im=erode_im.copy()
    (contour,_)=cv2.findContours(final_im,cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)
    if showIms:
        contour_im=numpy.zeros(siz,numpy.uint8)
        cv2.drawContours(contour_im,contour,-1,(200,0,0),1)
        img=Image.fromarray(contour_im)
        img.save("6ShowImages.png")
    FCA=[]
    FAA=[]
    for c in contour:
        area=cv2.contourArea(c)
        (x,y,),r=cv2.minEnclosingCircle(c)
        area2=math.pi*(r**2)
        if np is False:
            if area/area2 > 0.65:
                FCA.append(c)
                FAA.append(area)
        else:
            if area/area2 > (0.32+(area*0.0003)):
                #this is the most stringent cut-off I would set to get a
                #"very high-quality data set" as this eliminates all previously
                #identified issues; could change the 0.32 to 0.3 or eliminate the
                #area cut-off at this point entirely to get a greater number of
                #growth rates which do contain more experimental errors though! 
                FCA.append(c)
                FAA.append(area)
    if showIms:
        final_contour_im=numpy.zeros(siz,numpy.uint8)
        cv2.drawContours(final_contour_im,FCA,-1,(0,200,0),1)
        img=Image.fromarray(final_contour_im)
        img.save("7ShowImages.png")
    return(FCA,FAA)

############### Main Script ####################################################

# Final Photo
FINALPHOT=35
# Border size around the image
DX, DY = 25,25

# Find all available time courses
syspath = os.path.dirname(sys.argv[0]) #current working directory
fullpath = os.path.abspath(syspath) #same as above
allobj=os.listdir(fullpath) #list of all objects in the directory

# Getting all folder names
folders=[]
for f in allobj:
    if len(f)==6 and "." not in f and f[0]=="R" and f[3]=="C":
        folders.append(f)

# Generating a grey background of median pixels of all first images
bk=makeBackground(folders) #random.sample(folders,30)

# Make new folder for writing output images
if not os.path.exists("OutputImages"):
    os.makedirs("OutputImages")

imdict={}
# Iterating through all folders
for f in folders:
    print(f)
    finalphoto="img_%09d__000.tif"%FINALPHOT
    impath=os.path.join(f,finalphoto)
    blbs,blbs_area=getBlobs(impath,bk,showIms=False,DX=DX,DY=DY,np=False)

    # Paint all final blobs to an empty image
    black=numpy.zeros((bk.shape[0],bk.shape[1]),numpy.uint8)
    cv2.drawContours(black,blbs,-1,(255,255,255),-1)
    img=Image.fromarray(black)
    img.save("Selected_Colonies_{}.png".format(f))

    # Make directory for output images
    if not os.path.exists(os.path.join(fullpath,"OutputImages",f)):
        os.makedirs(os.path.join(fullpath,"OutputImages",f))

    # Find all photos available (before "final" photo)
    imdict[f]=[]
    files=os.listdir(os.path.join(fullpath,f))
    #gets files with folder name in them 
    tiffs=[]
    tiffcount=0

    # Iterating through all the files in the folder
    files=sorted(files,key=numericalSort) #to order files numerically 
    for filename in files:
        if ".tif" in filename and tiffcount<36: #only for files up to the final photo
            tiffs.append(filename)
            tiffcount+=1

    # Folder to save the output
    outputdir=os.path.join(fullpath,"New_Blobs_"+f)
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
    
    # Getting the ROI for all blobs in the cut-off image
    blbno=0
    counter=0
    for blb in blbs:
        x,y,w,h=cv2.boundingRect(blb)
        if blbs_area[blbno]<300:
            x=x-10
            y=y-10
            w=w+10
            h=h+10
        
        # Crop blob image and calculate the contours (should be only one!)
        timecourse_area=[]
        timecourse_time=[]
        blob_images=numpy.empty([1,FINALPHOT],dtype=object)
        for imno in range(0,FINALPHOT):
            imname=tiffs[imno]
            impath=os.path.join(fullpath,f,imname)
            imtim=os.path.getmtime(impath)
            currim=cv2.imread(impath,3)
            ROI=currim[y:y+h,x:x+w]
            ROI_img=Image.fromarray(ROI)
            colony,colony_area=getBlobs(ROI,bk,showIms=False,DX=DX,DY=DY,np=True)
            
            # Only save data for single clonal colonies 
            if len(colony) is 1:
                if imno is 0 and colony_area[0]>300.0:
                    print("Too big for a single colony")
                    pass
                else: 
                    timecourse_area.append(colony_area[0])
                    timecourse_time.append(imtim)
                    blob_images[0,imno]=ROI_img
                
        # Only save data for full time courses
        if len(timecourse_area) is 35:
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
            timecourse_image=Image.new('RGB',(7*w,5*h))
            x=0
            for i in range(0,5*h, h):
                for j in range(0,7*w,w):
                    timecourse_image.paste(blob_images[0,x],(j,i))
                    x+=1
            timecourse_image.save(os.path.join(
                outputdir,'Folder{}_Blob{:04d}_TimeCourse.jpg'.format(f,blbno)))
            # Growth Curve for each blob
            plt.figure(figsize=(6,4))
            plt.plot(timecourse_area,marker='o',ls='--')
            plt.ylabel('Area')
            plt.xlabel('Time')
            plt.title('Growth Curve for {} Blob {}'.format(f,blbno))
            plt.savefig(os.path.join(
                outputdir,'Folder{}_Blob{:04d}_GrowthCurve.jpg'.format(f,blbno)))
            plt.close()
            growth_curve=Image.open(os.path.join(
                outputdir,'Folder{}_Blob{:04d}_GrowthCurve.jpg'.format(f,blbno)))
            # Growth Curves + Time course Images for each blob
            gcw,gch=growth_curve.size
            width,height=timecourse_image.size
            width+=1
            height+=1
            final_image=Image.new('RGB', (gcw,gch+(height)), "white")
            final_image.paste(timecourse_image,((gcw-(width))/2,0))
            final_image.paste(growth_curve,(0,height))
            final_image.save(os.path.join(
                outputdir,'Folder{}_Blob{:04d}_TCGC.jpg'.format(f,blbno)))
        blbno+=1


        
    

