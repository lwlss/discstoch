import cv2.cv as cv
import cv2
import PIL, sys, os, numpy, random, math
from PIL import Image

def showImage(im):
    '''Convert cv image to PIL image and display'''
    if cv.GetImage(im).depth!=8:
        tmp = cv.CreateImage(cv.GetSize(im), cv.IPL_DEPTH_8U, im.channels)
        cv.Convert(im,tmp)
    if im.channels==1:
        pim = Image.fromstring("L", cv.GetSize(im), im.tostring())
    elif im.channels==3:
        pim = Image.fromstring("RGB", cv.GetSize(im), im.tostring())
    else:
        print("Strange number of channels there!")
    pim.show()

def dist(a,b):
    '''Distance between two points'''
    return((a[0]-b[0])**2+(a[1]-b[1])**2)

def findNearest(point,plist):
    '''Find element of plist nearest to point'''
    dmin=999999999999999
    nearest=-1
    for p in plist:
        d=dist(p,point)
        if d<dmin:
            dmin=d
            nearest=p
    return(nearest)

def makeBackground(folders):
    ''' Open stack of first images from each spot and get median pixel intensity to get a good background lighting map '''
    NoSpots=len(folders)
    example=cv.LoadImage(os.path.join(folders[0],"img_000000000__000.tif"),3)
    siz=cv.GetSize(example)
    megarr=numpy.zeros((siz[1],siz[0],3,NoSpots),numpy.uint16)

    spot=0
    for f in folders:
        currim=cv.LoadImageM(os.path.join(f,"img_000000000__000.tif"),3)
        megarr[:,:,:,spot]=numpy.asarray(currim)
        spot+=1
    bkg=numpy.array(numpy.round(numpy.median(megarr,axis=3)),dtype=numpy.uint16)
    backg=cv.fromarray(bkg)

    blnk=cv.CreateMat(siz[1]+51,siz[0]+51,cv.CV_16UC3)
    bkgIPL=cv.GetImage(backg)
    cv.CopyMakeBorder(backg,blnk,(25,25),1)
    return(blnk)

def getBlobs(image,bk,showIms=False,DX=25,DY=25):
    '''Get masks representing microcolony sizes and positions'''
    imnum=0
    # Open image as colour
    im=cv.LoadImage(image,3)
    siz=cv.GetSize(im)
    tot=float(siz[0]*siz[1])
    if showIms:
        cv.SaveImage("1ShowImages%05d.png"%imnum,im) #input image 
        imnum+=1
    # Add border around image (based on background)
    im2=cv.GetImage(bk)
    if showIms:
        cv.SaveImage("2ShowImages%05d.png"%imnum,im2) #final image 
        imnum+=1
    cv.SetImageROI(im2,(DX,DY,siz[0],siz[1]))
    #temporarily crop the image according to background
    cv.Copy(im,im2) #makes a copy and renames it 
    cv.ResetImageROI(im2)
    if showIms:
        cv.SaveImage("3ShowImages%05d.png"%imnum,im2) #input image 
        imnum+=1
    gray=cv.CreateImage(cv.GetSize(im2),8,1)
    cv.CvtColor(im2,gray,cv.CV_RGB2GRAY)

    # Generate Canny edge map
    dst = cv.CreateImage(cv.GetSize(im2),8,1)
    cv.Canny(gray,dst,10,50,3)
    if showIms:
        cv.SaveImage("4ShowImages%05d.png"%imnum,dst) #detects contours in image 
        imnum+=1

    # Find objects and fill holes with Contour
    storage=cv.CreateMemStorage()    
    element_shape = cv.CV_SHAPE_RECT 
    #Changed the dilation so that now it actually fills in the circles fully 
    cv.Dilate(dst,dst,iterations=5)
    if showIms:
        cv.SaveImage("5ShowImages%05d.png"%imnum,dst)
        imnum+=1
    cv.Erode(dst,dst,iterations=5)
    if showIms:
       cv.SaveImage("6ShowImages%05d.png"%imnum,dst)
       imnum+=1
    clone=cv.CloneImage(dst)
    imcol=cv.CloneImage(im2)
    contour=cv.FindContours(clone, storage, cv.CV_RETR_TREE, cv.CV_CHAIN_APPROX_SIMPLE, (0, 0))
    FCA=[]
    while contour:
        area=cv.ContourArea(contour)
        peri=cv.ArcLength(contour)
        i,(x,y),radius=cv.MinEnclosingCircle(contour)
        circumference=2*math.pi*radius
        area2=math.pi*(radius**2)
        if (0<area/area2>0.80):
            FCA.append(contour)
        contour=contour.h_next()
    return(FCA)


def imageFromBlob(blobs,bk):
    tot=1392640
    clist=([255,255,0],[255,0,0],[0,255,0],[0,0,255],[255,0,255],[0,255,255])
    black=cv.CreateImage(cv.GetSize(bk),8,3)
    #colourex=cv.CV_RGB(255,255,255)
    for blob in blobs:
        area=cv.ContourArea(blob)
        if area>15.0:
            mom=cv.Moments(blob)
            cx=mom.m01/mom.m00
            cy=mom.m10/mom.m00
            col=random.sample(clist,1)[0]
            colourex=cv.CV_RGB(col[0],col[1],col[2])
            #colourex=cv.CV_RGB(random.randint(0,255),random.randint(0,255),random.randint(0,255))
            #colourex=cv.CV_RGB(255,255,255)
            colourfnt=cv.CV_RGB(255,255,255)
            colourbig=cv.CV_RGB(255,0,255)
            cv.DrawContours(black,blob,colourex,colourex,0,1,8)
            cv.DrawContours(black,blob,colourex,colourex,0,cv.CV_FILLED,8)
            #cv.Circle(black,(int(round(cy)),int(round(cx))),2,cv.CV_RGB(0,0,255),-1)
            cv.PutText(black,str(round(100.0*area/tot,3)),(int(round(cy+20)),int(round(cx+20))),fnt,colourfnt)
            #cv.DrawContours(black,blob,colourex,colourex,0,cv.CV_FILLED,8)
    return(black)

FINALPHOT=35
# Size of border (px) to put around image
DX,DY=25,25

bigfnt=cv.InitFont(cv.CV_FONT_HERSHEY_SIMPLEX,2.0,2.0) #not sure why; not used 
colourbig=cv.CV_RGB(150,150,150) #specifies the colour grey
# Find all available image timecourses
syspath = os.path.dirname(sys.argv[0]) #current working directory
fullpath = os.path.abspath(syspath) #same as above
allobj=os.listdir(fullpath) #list of all objects in the directory

#Getting all of the clonal colony folder names (well or pinid or whatever it's called) 
folders=[]
for f in allobj:
    if len(f)==6 and "." not in f and f[0]=="R" and f[3]=="C":
        folders.append(f)

# Generate a background lighting map from the median of a stack of the first image from each spot
backg=makeBackground(random.sample(folders,30)) #function call;
#use first image with no cells of 30 randomly sampled plates and use average as background colour 
bk = cv.CreateImage(cv.GetSize(backg), cv.IPL_DEPTH_8U, 3)
cv.Convert(backg,bk)
#showImage(bk)

# Make new folder for writing output images
if not os.path.exists("OutputImages"):
    os.makedirs("OutputImages")

imdict={} #dictionary to store images in...
fnt=cv.InitFont(cv.CV_FONT_HERSHEY_PLAIN,1.5,1.5) #inititalizes font for writing text 
bigfnt=cv.InitFont(cv.CV_FONT_HERSHEY_SIMPLEX,4.0,4.0) #again not sure why; not used 

#iterating through all the folders 
for f in folders:
    print(f)
    # Get blobs from final photo
    finalphoto="img_%09d__000.tif"%FINALPHOT
    blbs=getBlobs(os.path.join(fullpath,f,finalphoto),bk,showIms=False,DX=DX,DY=DY) #showIms=True,
    ##THIS IS THE FUNCTION TO EDIT 

    # Make numpy array for storing results
    NBlbs=len(blbs) #how many blobs in the final photo
    print(NBlbs)
    res=numpy.zeros((FINALPHOT,NBlbs+1),numpy.float)
    locs=numpy.zeros((NBlbs,4),numpy.int)

    # Make directory for output images
    if not os.path.exists(os.path.join(fullpath,"OutputImages",f)):
        os.makedirs(os.path.join(fullpath,"OutputImages",f))

    # Find all photos available (before "final" photo)
    imdict[f]=[] #f is the folder 
    files=os.listdir(os.path.join(fullpath,f)) #gets files with folder name in them 
    tiffs=[]
    tiffcount=0
    #iterating through all of the files of each folder 
    for filename in files:
        if ".tif" in filename and tiffcount<36: #only for files up to the final photo 
            tiffs.append(filename)
            tiffcount+=1
    for imno in xrange(0,FINALPHOT): #iterating through all image numbers less than 36
        imname=tiffs[imno]
        #print(imno)
        imtim=os.path.getmtime(os.path.join(fullpath,f,imname)) #gets time of when image was saved? 
        res[imno,0]=imtim #stores the time 
        # Get current blobs
        cblbs=getBlobs(os.path.join(fullpath,f,imname),bk,showIms=False)
        #gets ALL current blobs for that image number in that folder
        
        # Paint all current blobs to an empty image
        black=cv.CreateImage(cv.GetSize(bk),8,1)
        colourex=cv.CV_RGB(255,255,255) #white
        for cblb in cblbs:
            #draws a white contour around the blob
            cv.DrawContours(black,cblb,colourex,colourex,0,cv.CV_FILLED,8)
        blbno=1
        for blb in blbs: #iterating through final photos            
            ROI=cv.BoundingRect(blb) #cuts out the rectangle of the last blob
            cv.SetImageROI(black,ROI) #for background specifies where rectangle goes
            maskedarea=cv.CountNonZero(black) #Area
            #print maskedarea, blbno, f
            #if maskedarea is not 0:
            cv.ResetImageROI(black)
            res[imno,blbno]=maskedarea
            locs[blbno-1]=[ROI[0]-DX,ROI[1]-DY,ROI[2],ROI[3]]
            blbno+=1
        blbno=0
        for blb in blbs:
            ROI=cv.BoundingRect(blb)
            cv.PutText(black,str(blbno),(ROI[0],ROI[1]),fnt,colourbig)
            blbno+=1
            
            
    # Write results for folder to file
    numpy.savetxt(f+"_OUT.txt",res,delimiter="\t")
    #numpy.savetxt(f+"_LOC.txt",locs,delimiter="\t")
    numpy.savetxt(f+"_LOC.txt",locs,delimiter="\t",fmt="%s")
    cv.SaveImage(f+"LOC.png",black)
