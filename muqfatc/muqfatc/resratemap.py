import muqfatc.imageanalysis as ia
import muqfatc.lineagegrowth as linia

import cv2
import numpy
from PIL import Image
import os, sys


def residualrange (inputfile):
    '''Defines a residual range and colour map for the time course
    residuals specified by the input file'''
    res_rate_data=numpy.loadtxt(inputfile,dtype=numpy.str)
    dist=[float(eval(i)) for i in res_rate_data[:,0]]
    d_id=[eval(a)+eval(b) for a,b in zip(res_rate_data[:,2],
                                            res_rate_data[:,3])]
    d=[round(q,1) for q in dist]
    colours=ia.red_yellow_colourmap((250,0,0))
    c=colours[::-1] #so that high residuals are red
    r_r=[p/100.0 for p in range(int(round(min(d)*100)),
                                      int(round(max(d)*100)),
                                      (int(round(max(d)*100))-
                                       int(round(min(d)*100)))/10)]
    return(d,d_id,c,r_r)

def raterange (inputfile):
    '''Defines a rate range and colour map for the time course rates
    specified by the input file'''
    res_rate_data=numpy.loadtxt(inputfile,dtype=numpy.str)
    r=[float(eval(s)) for s in res_rate_data[:,1]]
    colours2=ia.white_black_colourmap((250,250,250))
    c2=colours2[::-1] #so that fast rates are white
    r=[round(m,1) for m in r]
    r_r=[l/100.0 for l in range(0,int(round(max(r)*100)),
                                        int(round(max(r)*10)))]
    return(r,c2,r_r)

def makecolchart (col,title,siz,r_r):
    '''Plots a legend for rate and residual colour ranges'''
    colour_chart=numpy.zeros((1000,110,3),numpy.uint8)
    centers=range(50,1000,100)
    for x in range(0,len(col)):
        cv2.circle(colour_chart,(25,centers[x]+20),20,col[x],-1)
        cv2.putText(colour_chart,str(r_r[x]), (55,centers[x]+10),
                    cv2.FONT_HERSHEY_SIMPLEX, 0.8, (255,255,255))
    cv2.putText(colour_chart,title,(siz,30),
                    cv2.FONT_HERSHEY_SIMPLEX, 0.5, (255,255,255))
    img=Image.fromarray(colour_chart)
    return(img)

def makerateresim (folders,finalim,filename,bk,DX,DY):
    dist,dist_id,colours,res_range=residualrange(filename)
    rates,colours2,rates_range=raterange(filename)
    img_res=makecolchart(colours,"Res. Dist.",10,res_range)
    img_res.save("Colour_Chart_Res_Filtered_Norm_{}.png".format(finalim))
    img_rate=makecolchart(colours2,"Growth Rates",0,rates_range)
    img_rate.save("Colour_Chart_Rate_Filtered_Norm_{}.png".format(finalim))
    # Iterating through all folders
    for f in folders:
        print(f)
        finalphoto="img_%09d__000.tif"%finalim
        impath=os.path.join(f,finalphoto)
        blbs,blbs_area=linia.getBlobs(impath,bk,showIms=False,DX=DX,DY=DY,
                                      np=False,apply_filt=False)
        res_black=numpy.zeros(bk.shape,numpy.uint8)
        # Getting the ROI for all blobs in the cut-off image
        namespace=[]
        blbno=0
        counter=0
        for blb in blbs:
            x,y,w,h=cv2.boundingRect(blb)
            if blbs_area[blbno]<((bk.shape[0]*bk.shape[1])*0.000215): #300px
                x=int(round(x-((bk.shape[0]*bk.shape[1])*7.2e-06))) #10px
                y=int(round(y-((bk.shape[0]*bk.shape[1])*7.2e-06)))
                w=int(round(w+(2*((bk.shape[0]*bk.shape[1])*7.2e-06)))) #20px
                h=int(round(h+(2*((bk.shape[0]*bk.shape[1])*7.2e-06))))

            # Colouring the blobs in the final image according to residuals
            curr_id=f+str(blbno)
            index=[index for index,item in enumerate(dist_id) if item==curr_id]
            if len(index)>0:
                curr_dist=dist[index[0]]
                # Colour range which the current residual distance closest to
                col_val=min(res_range, key=lambda x:abs(x-round(curr_dist,1)))
                col_index=[ind for ind,it in enumerate(res_range) if it==col_val]
                cent=(int(x+(w/2.0)),int(y+h/2.0))
                cv2.circle(res_black,cent,20,colours[col_index[0]],-1)
                curr_rate=rates[index[0]]
                col_val2=min(rates_range, key=lambda n:abs(n-round(curr_rate,1)))
                col_index2=[ind for ind,it in enumerate(rates_range) if it==col_val2]
                cent2=(int(x+(w/2.0)),int(y+h/2.0))
                cv2.circle(res_black,cent2,10,colours2[col_index2[0]],-1)
            blbno+=1

        # Save residual range coloured image
        img=Image.fromarray(res_black)
        img.save("Residuals_Filtered_Norm_{}{}.png".format(finalim,f))
        # Save the blobs according to folder name
        numpy.savetxt(f+"_BLBNO_{}.txt".format(finalim),namespace,delimiter="\t")
