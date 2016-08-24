#### Generic functions required for Image Analysis
#### These functions are used both for estimating pin and colony growth
#### Requires OpenCV version 3.0.0 to be installed

import cv2
import numpy
import re, os

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

def red_yellow_colourmap(original_col,shade=10):
    '''Red colour specification is input and returns
    a red to yellow colour map'''
    col_list=[]
    col_list.append(original_col)
    curr_col=original_col
    z=max(original_col)/float(shade)
    for i in range(0,shade-1):
        new_col=(curr_col[0],curr_col[1]+z,curr_col[2])
        col_list.append(new_col)
        curr_col=new_col
    return(col_list)

def white_black_colourmap(original_col,shade=10):
    '''White colour specification is input and returns
    a white to black colour map'''
    col_list=[]
    col_list.append(original_col)
    curr_col=original_col
    z=max(original_col)/float(shade)
    for i in range(0,shade-1):
        new_col=(curr_col[0]-z,curr_col[1]-z,curr_col[2]-z)
        col_list.append(new_col)
        curr_col=new_col
    return(col_list)

def numericalSort(value):
    '''To sort file names numerically'''
    #Function was taken from here:
    #http://stackoverflow.com/questions/4623446/how-do-you-sort-files-numerically 
    numbers = re.compile(r'(\d+)')
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts
