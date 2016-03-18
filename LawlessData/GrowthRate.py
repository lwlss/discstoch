import cv, PIL, sys, os, numpy, random
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

# Find all available image timecourses
syspath = os.path.dirname(sys.argv[0])
fullpath = os.path.abspath(syspath)
allobj=os.listdir(fullpath)
folders=[]
for f in allobj:
    if len(f)==6 and "." not in f and f[0]=="R" and f[3]=="C":
        folders.append(f)

# Generate a background lighting map from the median of a stack of the first image from each spot
backg=makeBackground(random.sample(folders,30))
bk = cv.CreateImage(cv.GetSize(backg), cv.IPL_DEPTH_8U, 3)
cv.Convert(backg,bk)
#showImage(bk)

# Make new folder for writing output images
if not os.path.exists("OutputImages"):
    os.makedirs("OutputImages")

imdict={}
fnt=cv.InitFont(cv.CV_FONT_HERSHEY_PLAIN,1.0,1.0)
bigfnt=cv.InitFont(cv.CV_FONT_HERSHEY_SIMPLEX,4.0,4.0)

for f in folders[0:3]:
    print(f)
    if not os.path.exists(os.path.join(fullpath,"OutputImages",f)):
        os.makedirs(os.path.join(fullpath,"OutputImages",f))
    imdict[f]=[]
    files=os.listdir(os.path.join(fullpath,f))
    tiffs=[]
    tiffcount=0
    for filename in files:
        if ".tif" in filename:
            tiffs.append(filename)
    for imname in tiffs:
        imtim=os.path.getmtime(os.path.join(fullpath,f,imname))

        # Open image as colour
        im=cv.LoadImage(os.path.join(fullpath,f,imname),3)
        siz=cv.GetSize(im)
        tot=float(siz[0]*siz[1])

        # Paste original image into padded image
##        backgIPL=cv.GetImage(bk)
##        cv.SetImageROI(backgIPL,(25,25,siz[0],siz[1]))
##        cv.Copy(im,backgIPL)
##        cv.ResetImageROI(backgIPL)
        im2=cv.GetImage(bk)
        cv.SetImageROI(im2,(25,25,siz[0],siz[1]))
        cv.Copy(im,im2)
        cv.ResetImageROI(im2)
        #showImage(backgIPL)

        # Edge detection
        # Convert image to greyscale
        gray=cv.CreateImage(cv.GetSize(im2),8,1)
        cv.CvtColor(im2,gray,cv.CV_RGB2GRAY)

        # Generate Canny edge map
        dst = cv.CreateImage(cv.GetSize(im2),8,1)
        #cv.Canny(gray,dst,50,100)
        cv.Canny(gray,dst,10,50,3)
        #showImage(dst)

        # Adaptive thresholding of image
        # Convert image to greyscale
        #gray=cv.CreateImage(cv.GetSize(im),8,1)
        #cv.CvtColor(im,gray,cv.CV_RGB2GRAY)
        #dst = cv.CreateImage(cv.GetSize(gray), cv.IPL_DEPTH_8U, 1)
        #cv.AdaptiveThreshold(gray,dst,255,cv.CV_ADAPTIVE_THRESH_MEAN_C,cv.CV_THRESH_BINARY_INV,101,5)
        #showImage(dst)

        # Gaussian blur
        ##cv.Smooth(dst,dst,cv.CV_GAUSSIAN,55,0,0,0)
        ##showImage(dst)
        ##cv.Threshold(dst,dst,100,255,cv.CV_THRESH_BINARY)
        ##showImage(dst)

        # Find objects and fill holes with Contour
        storage=cv.CreateMemStorage()
        element_shape = cv.CV_SHAPE_RECT

        for x in xrange(1,10):
            cv.Dilate(dst,dst)
            cv.Erode(dst,dst)

        #showImage(dst)

        clone=cv.CloneImage(dst)
        if imname==tiffs[0]:
            firstim=cv.CloneImage(clone)
        imcol=cv.CloneImage(im2)
        contour=cv.FindContours(clone, storage, cv.CV_RETR_TREE, cv.CV_CHAIN_APPROX_SIMPLE, (0, 0))
        totalarea=0.0
        cellcount=0
        while contour:
            area=cv.ContourArea(contour)
            if area>5.0:
                mom=cv.Moments(contour)
                cx=mom.m01/mom.m00
                cy=mom.m10/mom.m00
                totalarea+=area
                cellcount+=1
                #colourex=cv.CV_RGB(random.randint(0,255),random.randint(0,255),random.randint(0,255))
                colourex=cv.CV_RGB(255,255,255)
                colourfnt=cv.CV_RGB(255,255,0)
                colourbig=cv.CV_RGB(255,0,255)
                cv.DrawContours(imcol,contour,colourex,colourex,0,1,8)
                cv.DrawContours(imcol,contour,colourex,colourex,0,cv.CV_FILLED,8)
                #cv.Circle(imcol,(int(round(cy)),int(round(cx))),2,cv.CV_RGB(0,0,255),-1)
                cv.PutText(imcol,str(round(100.0*area/tot,3)),(int(round(cy)),int(round(cx))),fnt,colourfnt)
            contour=contour.h_next()
        # Overlay original cell sizes in red
        contour=cv.FindContours(firstim, storage, cv.CV_RETR_TREE, cv.CV_CHAIN_APPROX_SIMPLE, (0, 0))
        while contour:
            area=cv.ContourArea(contour)
            if area>5.0:
                colourex=cv.CV_RGB(255,0,0)
                cv.DrawContours(imcol,contour,colourex,colourex,0,1,8)
            contour=contour.h_next()
        cv.PutText(imcol,"% Area: "+str(round(100.0*totalarea/tot,2)),(int(round(0.1*siz[0])),int(round(0.8*siz[1]))),bigfnt,colourbig)
        cv.PutText(imcol,"Count: "+str(cellcount),(int(round(0.1*siz[0])),int(round(0.9*siz[1]))),bigfnt,colourbig)
        #showImage(imcol)
        cv.SaveImage(os.path.join(fullpath,"OutputImages",f,"Mod%05d.png"%tiffcount),imcol)
        cv.SaveImage(os.path.join(fullpath,"OutputImages",f,"Orig%05d.png"%tiffcount),im)
        tiffcount+=1
        imdict[f].append([imtim,totalarea,cellcount])

# Write results to file
out=open("ReportAA.txt","w")
out.write("%s\t%s\t%s\t%s\n"%("Culture","Time","TotalArea","Number"))
for f in imdict.keys():
    res=imdict[f]
    for row in res:
        out.write("%s\t%f\t%f\t%d\n"%(f,row[0],row[1],row[2]))
out.close()




