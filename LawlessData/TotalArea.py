############ GETTING THE TOTAL AREA OF EACH PIN OVER TIME #################

import cv2
import PIL, sys, os, numpy, random, math, re
from PIL import Image
import matplotlib.pyplot as plt

############ Functions ###################

def numericalSort(value):
    '''To sort file names numerically'''
    #Function was taken from here:
    #http://stackoverflow.com/questions/4623446/how-do-you-sort-files-numerically 
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

def makeBorder(image,bk,DX=25,DY=25):
    '''Adds a border based on background around the image.'''
    siz=image.shape #row, column = height, width
    border_im=cv2.resize(bk,(2*DY+siz[1],2*DX+siz[0]),
                         interpolation = cv2.INTER_CUBIC)
    new_siz=border_im.shape
    border_im[DY:new_siz[0]-DY,DX:new_siz[1]-DX]=image
    return(border_im)

def totalArea(image,col):
    '''Calculates total area of colonies in that image'''
    im=cv2.imread(image,3)

    lower=numpy.array((0, 0, 0),dtype="uint8")
    upper=numpy.array((col, col, col),dtype="uint8")
    mask=cv2.inRange(im,lower,upper)
##    cv2.imshow("Range", mask)
##    cv2.waitKey(0)
##    cv2.imshow("Range", im)
##    cv2.waitKey(0)
##    cv2.destroyAllWindows()
    celldens=cv2.countNonZero(mask)
    return(celldens)

############ Main ########################

# Find all available time courses
syspath = os.path.dirname(sys.argv[0]) #current working directory
fullpath = os.path.abspath(syspath) #same as above
allobj=os.listdir(fullpath) #list of all objects in the directory

# Getting all folder names
folders=[]
for f in allobj:
    if len(f)==6 and "." not in f and f[0]=="R" and f[3]=="C" and f!='R03C03':
        folders.append(f)

# Make background
bk=makeBackground(folders)
DX=25
DY=25

##plt.figure()

all_area=numpy.zeros((len(folders), 90))
all_time=numpy.zeros((len(folders), 90))

foldercount=0
for f in folders:
    print(f)
    # Find all photos in that folder
    files=os.listdir(os.path.join(fullpath,f))

    # Sort files numerically 
    files=sorted(files,key=numericalSort)

    # Only look at tiff files in that folder 
    tiffs=[]
    tiffcount=0
    for filename in files:
        if ".tif" in filename and tiffcount<90:
            tiffs.append(filename)
            tiffcount+=1

    # Iterate through tiff files
    area=[]
    time=[]
    for filename in tiffs:
        # Get total area of cells in each image
        impath=os.path.join(f,filename)
        col=110
        imtim=os.path.getmtime(impath)
        imkernel=1
        tiffarea=totalArea(impath,col)
        area.append(tiffarea)
        time.append(imtim)

    time=[(i-time[0])/3600 for i in time]
    all_time[foldercount,:]=time
    all_area[foldercount,:]=area
    foldercount+=1

numpy.savetxt("PopulationArea.txt",all_area,delimiter="\t")
numpy.savetxt("PopulationTime.txt",all_time,delimiter="\t")

text_file = open("PopulationFolders.txt", "w")
for item in folders:
  text_file.write("%s\n" % item)
text_file.close()
