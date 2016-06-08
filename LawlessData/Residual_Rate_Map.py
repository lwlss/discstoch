# Generate colour coded residual and rate image
import numpy
from PIL import Image, ImageFont, ImageDraw
import sys, os

# Convert microQFA pinid to QFA pinid (plate is upside down)
def pinid(pid):
    Col=int(pid[4:])
    NewCol=24-(Col-1)
    NewFolder=pid[:4]+str(NewCol)
    return(NewFolder)

data=numpy.loadtxt("Parser_Lawless_Namespace.txt",dtype=numpy.str)
folders=data[:,0]

example=Image.open("Residuals_{}.png".format(eval(folders[0])))
w,h=example.size
col_chart1=Image.open("Colour_Chart_Res.png")
col_chart2=Image.open("Colour_Chart_Rate.png")
col_w,col_h=col_chart1.size
plate=Image.new('RGB',(w*8+col_w,h*8),"grey")
plate.paste(col_chart2,(0,0))
plate.paste(col_chart1,(0,col_h))

counter=0
for f in folders:
    f=eval(f)
    print(f)
    curr_im=Image.open("Residuals_{}.png".format(f))
    true_f=pinid(f)
    col=int(true_f[4:])
    ro=int(true_f[1:3])
    plate.paste(curr_im,(w*(col-15)+col_w,h*(ro-3)))
    draw=ImageDraw.Draw(plate)
    draw.text((w*(col-15)+col_w,h*(ro-3)),f,(255,255,255))

plate.show()
    
