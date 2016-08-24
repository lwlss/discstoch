# Colour-coding the image plate accoding to residual range
# and estimated growth rate

import muqfatc.resratemap as rrmap
import muqfatc.imageanalysis as ia
import os, sys
import numpy
from PIL import Image, ImageFont, ImageDraw

# Find all available time courses
syspath = os.path.dirname(sys.argv[0]) #current working directory
fullpath = os.path.abspath(syspath) 
allobj=os.listdir(fullpath) #list of all objects in the directory

FINALPHOT=20
myfolders=[]
for f in allobj:
    if len(f)==6 and "." not in f and f[0]=="R" and f[3]=="C":
        myfolders.append(f)

# Generate backgroun image
DX=25
DY=25
bk=ia.makeBackground(myfolders)

# Colour-coding pins accoding to residual range and estimated growth rate
filename="Lawless_ResRateRange_Filtered_Norm_{}.txt".format(FINALPHOT)
rrmap.makerateresim(myfolders,FINALPHOT,filename,bk,DX=25,DY=25)

#########################################################

# Generate colour coded residual and rate image for the WHOLE PLATE
# This requires folders for the whole plate with the top-left pin

# Convert microQFA pinid to QFA pinid (plate is upside down)
def pinid(pid):
    Col=int(pid[4:])
    NewCol=24-(Col-1)
    NewFolder=pid[:4]+str(NewCol)
    return(NewFolder)

data=numpy.loadtxt("Parser_Lawless_Namespace.txt",dtype=numpy.str)
folders=data[:,0]

example=Image.open("Residuals_Filtered_Norm_{}{}.png".format(FINALPHOT,eval(folders[0])))
w,h=example.size
col_chart1=Image.open("Colour_Chart_Res_Filtered_Norm_{}.png".format(FINALPHOT))
col_chart2=Image.open("Colour_Chart_Rate_Filtered_Norm_{}.png"/format(FINALPHOT))
col_w,col_h=col_chart1.size
plate=Image.new('RGB',(w*8+col_w+(9*10),h*8+(8*10)),"grey")
plate.paste(col_chart2,(0,0))
plate.paste(col_chart1,(0,col_h))

counter=0
for f in folders:
    f=eval(f)
    print(f)
    curr_im=Image.open("Residuals_Filtered_Norm_{}{}.png".format(FINALPHOT,f))
    true_f=pinid(f)
    col=int(true_f[4:])
    ro=int(true_f[1:3])
    plate.paste(curr_im,(w*(col-15)+(col-14)*10+col_w,h*(ro-3)+(ro-2)*10))
    draw=ImageDraw.Draw(plate)
    draw.text((w*(col-15)+(col-14)*10+col_w,h*(ro-3)+(ro-2)*10),f,(255,255,255))

plate.save("ResRateColMap_Filtered_Norm_{}.png".format(FINALPHOT))
plate.show()
