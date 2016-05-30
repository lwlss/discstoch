############# IMAGE ANALYSIS IN CV2 ##############

import cv2
import cv2.cv as cv 
import PIL, sys, os, numpy, random, math, re
from PIL import Image

################# Functions ######################

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

def getBlobs(image,bk,showIms=False,DX=25,DY=25):
    '''Get masks representing microcolony sizes and positions'''
    # Open image as colour
    im=cv2.imread(image,3)
    siz=im.shape #row, column = height, width
    if showIms:
        img=Image.fromarray(im,'RGB')
        img.save("1ShowImages.png")
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
    dilate_im=cv2.dilate(canny_im,kernel,iterations=4)
    if showIms:
        img=Image.fromarray(dilate_im)
        img.save("4ShowImages.png")
    erode_im=cv2.erode(dilate_im,kernel,iterations=4)
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
    for c in contour:
        area=cv2.contourArea(c)
        (x,y,),r=cv2.minEnclosingCircle(c)
        area2=math.pi*(r**2)
        if area/area2 > 0.7:
            FCA.append(c)
    if showIms:
        final_contour_im=numpy.zeros(siz,numpy.uint8)
        cv2.drawContours(final_contour_im,FCA,-1,(0,200,0),1)
        img=Image.fromarray(final_contour_im)
        img.save("7ShowImages.png")
    return(FCA)

################# Main Script ####################

# Final Photo
FINALPHOT=35
# Border size around the image
DX, DY = 25,25

# Colour Specifications (#Not sure if these are useful, come back or delete)
colourbig=numpy.uint8([[[150,150,150]]])
#deprectated: fnt=cv.InitFont(cv.CV_FONT_HERSHEY_PLAIN,1.5,1.5)
colourex=numpy.uint8([[[255,255,255]]]) #white

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
for f in folders[0:3]: #change this back to folders at the end 
    print(f)
    finalphoto="img_%09d__000.tif"%FINALPHOT
    impath=os.path.join(f,finalphoto)
    blbs=getBlobs(impath,bk,showIms=True,DX=DX,DY=DY)
    Nblbs=len(blbs) #number of blobs in the final photo
    res=numpy.zeros((FINALPHOT,Nblbs+1),numpy.float)
    locs=numpy.zeros((Nblbs,4),numpy.int)

    # Make directory for output images
    if not os.path.exists(os.path.join(fullpath,"OutputImages",f)):
        os.makedirs(os.path.join(fullpath,"OutputImages",f))

    # Find all photos available (before "final" photo)
    imdict[f]=[]
    files=os.listdir(os.path.join(fullpath,f)) #gets files with folder name in them 
    tiffs=[]
    tiffcount=0

    # Iterating through all the files in the folder
    files=sorted(files,key=numericalSort) #to order files numerically 
    for filename in files:
        if ".tif" in filename and tiffcount<36: #only for files up to the final photo
            tiffs.append(filename)
            tiffcount+=1   
    for imno in xrange(0,FINALPHOT): #iterating through all image numbers less than 36
        imname=tiffs[imno]
        imtim=os.path.getmtime(os.path.join(fullpath,f,imname)) #gets time of when image was saved? 
        res[imno,0]=imtim #stores the time 
        # Get current blobs
        impath=os.path.join(fullpath,f,imname)
        cblbs=getBlobs(impath,bk,showIms=True,DX=DX,DY=DY)
        currim=cv2.imread(impath,3)
        #cv2.imshow("output",currim)
        #cv2.waitKey(0)
        #gets ALL current blobs for that image number in that folder

        # Paint all current blobs to an empty image
        black=numpy.zeros((bk.shape[0],bk.shape[1]),numpy.uint8)
        cv2.drawContours(black,cblbs,-1,(255,255,255),-1)

        # Cropping the final photos
        blbno=1
        for blb in blbs: #would have to iterate through cblbs here to crop all images
            x,y,w,h=cv2.boundingRect(blb)
            ROI=black[y:y+h,x:x+w]
            maskedarea=cv2.countNonZero(ROI)
            res[imno,blbno]=maskedarea
            blob=currim[y:y+h,x:x+w]
            #########CONTINUE HERE, NEED TO SAVE ALL BLOBS
            #########WOULD LIKE TO GET TIME COURSE IMAGES AND AREA GROWTH CURVE IN ONE GO
            locs[blbno-1]=[x,y,w,h]       
            blbno+=1
            
            

### Write results for folder to file
##numpy.savetxt(f+"_OUT.txt",res,delimiter="\t")
##    #numpy.savetxt(f+"_LOC.txt",locs,delimiter="\t")
##numpy.savetxt(f+"_LOC.txt",locs,delimiter="\t",fmt="%s")
##cv.SaveImage(f+"LOC.png",black)
                        
    
    #cv2.imshow("output",new_im)
    #cv2.waitKey(0)
    

cv2.destroyAllWindows()



########## Old Stuff (To Be Deleted at the End) #########

##def makeBackground(folders):
##    NoSpots=len(folders)
##    example=cv2.imread(os.path.join(folders[0],"img_000000000__000.tif"),3)
##    siz=example.shape
##    megarr=numpy.zeros((siz[0],siz[1],3,NoSpots),numpy.uint16)
##    spot=0
##    for f in folders:
##        currim=cv2.imread(os.path.join(f,"img_000000000__000.tif"),3)
##        megarr[:,:,:,spot]=numpy.asarray(currim)
##        spot+=1
##    bkg=numpy.array(numpy.round(numpy.median(megarr,axis=3)),dtype=numpy.uint16)
##    cv2.copyMakeBorder(bkg,25,25,25,25,cv2.BORDER_REPLICATE)
####    cv2.imshow("output",bkg)
####    cv2.waitKey(0)
##    return(bkg)
##
