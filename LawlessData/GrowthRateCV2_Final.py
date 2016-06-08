############### IMAGE ANALYSIS IN CV2 - RESTRUCTURED ###########################

import cv2
import PIL, sys, os, numpy, random, math, re
from PIL import Image
import matplotlib.pyplot as plt

############### Functions ######################################################

def numericalSort(value):
    '''To sort file names numerically'''
    #Function was taken from here:
    #http://stackoverflow.com/questions/4623446/how-do-you-sort-files-numerically
    #Needs proper referencing or rewrite! 
    numbers = re.compile(r'(\d+)')
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts

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

    if np is False:
        border_im=makeBorder(im,bk,DX=DX,DY=DY)
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

############### Main Script ####################################################

# Final Photo
FINALPHOT=35 #35; must be a number divisible by 5
# Border size around the image
DX, DY = 25,25

#Filter for only meaningful clonal colony growth curves?
apply_filt=False
#Colour-code final residuals according to residual range? 
res=True
rates=True
if res:
    res_rate_data=numpy.loadtxt("Lawless_ResidualRange.txt",dtype=numpy.str)
    dist=[float(eval(i)) for i in res_rate_data[:,0]]
    dist_id=[eval(a)+eval(b) for a,b in zip(res_rate_data[:,2],res_rate_data[:,3])]
    dist=[round(q,1) for q in dist]
    colours=red_yellow_colourmap((250,0,0))
    colours=colours[::-1] #so that high residuals are red
    res_range=[p/10.0 for p in range(0,int(round(max(dist)))*10,int(round(max(dist))))]
    if rates:
        rates=[float(eval(s)) for s in res_rate_data[:,1]]
        colours2=white_black_colourmap((250,250,250))
        colours2=colours2[::-1] #so that fast rates are white
        rate=[round(m,1) for m in dist]
        rates_range=[l/10.0 for l in range(0,int(round(max(rates)))*10,int(round(max(rates))))]

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

# Making a colour chart for Residuals and Rates
if res:
    colour_chart_res=numpy.zeros((1000,110,3),numpy.uint8)
    centers=range(50,1000,100)
    for x in range(0,len(colours)):
        cv2.circle(colour_chart_res,(25,centers[x]+20),20,colours[x],-1)
        cv2.putText(colour_chart_res,str(res_range[x]), (55,centers[x]+10),
                    cv2.FONT_HERSHEY_SIMPLEX, 0.8, (255,255,255))
    cv2.putText(colour_chart_res,"Res. Dist.",(10,30),
                    cv2.FONT_HERSHEY_SIMPLEX, 0.5, (255,255,255))
    img_res=Image.fromarray(colour_chart_res)
    img_res.save("Colour_Chart_Res.png")
    img_res.show()
if rates:
    colour_chart_rate=numpy.zeros((1000,110,3),numpy.uint8)
    centers=range(50,1000,100)
    for z in range(0,len(colours2)):
        cv2.circle(colour_chart_rate,(25,centers[z]+20),10,colours2[z],-1)
        cv2.putText(colour_chart_rate,str(rates_range[z]), (55,centers[z]+10),
                    cv2.FONT_HERSHEY_SIMPLEX, 0.8, (255,255,255))
    cv2.putText(colour_chart_rate,"Growth Rates",(0,30),
                    cv2.FONT_HERSHEY_SIMPLEX, 0.5, (255,255,255))
    img_rate=Image.fromarray(colour_chart_rate)
    img_rate.save("Colour_Chart_Rate.png")
    img_rate.show()

# Iterating through all folders
for f in folders:
    print(f)
    finalphoto="img_%09d__000.tif"%FINALPHOT
    impath=os.path.join(f,finalphoto)
    blbs,blbs_area=getBlobs(impath,bk,showIms=False,DX=DX,DY=DY,
                            np=False,apply_filt=apply_filt)

    # Paint all final blobs to an empty image
    black=numpy.zeros((bk.shape[0],bk.shape[1]),numpy.uint8)
    cv2.drawContours(black,blbs,-1,(255,255,255),-1)
    img=Image.fromarray(black)
    img.save("Selected_Colonies_{}.png".format(f))
    if res:
        res_black=numpy.zeros(bk.shape,numpy.uint8)

    # Find all photos available (before "final" photo)
    files=os.listdir(os.path.join(fullpath,f))
    #gets files with folder name in them 
    tiffs=[]
    tiffcount=0

    # Iterating through all the files in the folder
    files=sorted(files,key=numericalSort) #to order files numerically 
    for filename in files:
        if ".tif" in filename and tiffcount<(FINALPHOT+1):
            tiffs.append(filename)
            tiffcount+=1

    # Folder to save the output
    outputdir=os.path.join(fullpath,"Unfiltered_Blobs_{}_".format(FINALPHOT)+f)
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
    
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
        if res:
            curr_id=f+str(blbno)
            index=[index for index,item in enumerate(dist_id) if item==curr_id]
            if len(index)>0:
                curr_dist=dist[index[0]]
                # Colour range which the current residual distance closest to
                col_val=min(res_range, key=lambda x:abs(x-round(curr_dist,1)))
                col_index=[ind for ind,it in enumerate(res_range) if it==col_val]
                cent=(int(x+(w/2.0)),int(y+h/2.0))
                cv2.circle(res_black,cent,20,colours[col_index[0]],-1)
                if rates:
                    curr_rate=dist[index[0]]
                    col_val2=min(rates_range, key=lambda x:abs(x-round(curr_rate,1)))
                    col_index2=[ind for ind,it in enumerate(rates_range) if it==col_val2]
                    cent2=(int(x+(w/2.0)),int(y+h/2.0))
                    cv2.circle(res_black,cent2,10,colours2[col_index2[0]],-1)

                # Mark growth rate as text 
                #cv2.drawContours(res_black,blb,-1,colours[col_index[0]],thickness=-1)
                ##the above code doesn't save them filled in, not sure why???
                
        
        # Crop blob image and calculate the contours
        timecourse_area=[]
        timecourse_time=[]
        blob_images=numpy.empty([1,FINALPHOT],dtype=object)
        blob_bw_images=numpy.empty([1,FINALPHOT],dtype=object)
        for imno in range(0,FINALPHOT):
            imname=tiffs[imno]
            impath=os.path.join(fullpath,f,imname)
            imtim=os.path.getmtime(impath)
            currim=cv2.imread(impath,3)
            currim=makeBorder(currim,bk,DX=DX,DY=DY)
            ROI=currim[y:y+h,x:x+w]
            ROI_img=Image.fromarray(ROI)
            display_image=False
            colony,colony_area=getBlobs(ROI,bk,showIms=display_image,DX=DX,DY=DY,
                                        np=True,apply_filt=apply_filt)
            
            # Only save data for single clonal colonies 
            if len(colony) is 1 and apply_filt is True:
                if imno is 0 and colony_area[0]>((bk.shape[0]*bk.shape[1])*0.000215): #300px
                    print("Too big for a single colony")
                    pass
                else: 
                    timecourse_area.append(colony_area[0])
                    timecourse_time.append(imtim)
                    blob_images[0,imno]=ROI_img
                    if display_image:
                        blob_bw_images[0,imno]=Image.open("5ShowImages.png")

            if apply_filt is False:
                    if len(colony_area) is 0:
                        timecourse_area.append(0)
                    elif len(colony_area) is 1:
                        timecourse_area.append(colony_area[0])
                    else:
                        timecourse_area.append(sum(colony_area))
                    timecourse_time.append(imtim)
                    blob_images[0,imno]=ROI_img
                    if display_image:
                        blob_bw_images[0,imno]=Image.open("5ShowImages.png")
                
        # Only save data for full time courses
        if len(timecourse_area) is FINALPHOT and apply_filt is True:
            if counter is 0:
                folder_area=timecourse_area
                folder_time=timecourse_time
                folder_images=blob_images[0,:]
                counter+=1
            else:
                folder_area=numpy.vstack((folder_area,timecourse_area))
                folder_time=numpy.vstack((folder_time,timecourse_time))
                folder_images=numpy.vstack((folder_images,blob_images))

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
        col=FINALPHOT/5
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
        namespace.append(blbno)
            
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

    # Save residual range coloured image
    if res:
        img=Image.fromarray(res_black)
        img.save("Residuals_{}.png".format(f))
    
    # Save the area and time according to folder name
    numpy.savetxt(f+"_UNFILTERED_AREA.txt",folder_area,delimiter="\t")
    numpy.savetxt(f+"_UNFILTERED_TIME.txt",folder_time,delimiter="\t")
    numpy.savetxt(f+"_UNFILTERED_BLBNO.txt",namespace,delimiter="\t")
