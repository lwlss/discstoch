#PARSER FOR LAWLESS DATA 
#To bring all three data sets into the same format 
#so that a single script can be used to extract growth curves 

# Convert microQFA pinid to QFA pinid (plate is upside down)
pinid=function(pid){
  Col=as.numeric(substr(pid,5,6))
  NewCol=24-(Col-1)
  return(paste(substr(pid,1,4),sprintf("%02d",NewCol),sep=""))
}

# Read in library description file describing which strain at each spot
makeORFGetter=function(libraries,lib="LiqCurve",plate=1){
  # Open the library descriptions
  libs=read.delim(libraries,sep="\t",header=TRUE,stringsAsFactors=FALSE)
  libs$ORF=toupper(libs$ORF)
  # What library are we using? e.g. library="AndrewsOE_v2" library="SDL_v2"
  liblist=unique(libs$Library)
  NLIB=length(liblist)
  libdict=1:NLIB
  names(libdict)=liblist
  #libs=libs[libs$Library==library,]
  # Make an array for storing info about spots
  NROW=max(libs$Row)
  NCOL=max(libs$Column)
  NPLATE=max(libs$Plate)
  SPOTARRAY=array("missing",dim=c(NLIB,NPLATE,NROW,NCOL))
  # Fill the spot array object
  for (x in 1:length(libs$Plate)) SPOTARRAY[libdict[[libs[x,"Library"]]],libs[x,"Plate"],libs[x,"Row"],libs[x,"Column"]]=libs[x,"ORF"]
  getORF<-function(row,col){
    SPOTARRAY[libdict[[lib]],plate,row,col]}
}

extractStuff=function(fname){
  strname=pinid(substr(fname,1,6)) #file name 
  row=as.numeric(substr(strname,2,3)) 
  col=as.numeric(substr(strname,5,6)) 
  strgene=getORF(row,col) #uses extracted numbers to get the strain_replicate name 
  print(strgene)
  microdat=read.delim(fname,header=FALSE,sep="\t") #includes area and time (time is the first column)
  numcolonies=dim(microdat)[2]
  
  #Column1 is the time
  tim=(microdat[,1]-microdat[1,1])/(3600) #start time at zero and convert from seconds to hours 

  #Calculating the starting area and the growth rate using NLS
  AList=c()
  rList=c()
  for(x in 2:numcolonies){
    dat=microdat[,x] #transposes data; single lineage data is now horizontal 
    dat=sapply(dat,max,1) #turns all zeros in the list into ones... why??? probably because cells can arise from nothing; need to have at least one 
    # if(length(unique(dat))>1){
    #   inocguess=max(min(dat),1) #minimum data point; but not smaller than 1 
    #   # Fit exponential model
    #   result=nls(y~A+r*x,data=data.frame(x=tim,y=log(dat)),start=list(A=log(inocguess),r=0.1)) #not accounting for outliers 
    #   res=result$m$getPars()
    #   A=exp(as.numeric(res[1])) #converting back from logs 
    #   r=max(0,as.numeric(res[2])) #excludes negative growth rates and turns them to zero; negative growth rate has no biological meaning 
    # } else{A=1;r=0}
    # AList=c(AList,A)
    # rList=c(rList,r)
    gene=substr(strgene,1,4)
    output1=data.frame(gene, strname, strgene) #Genotype, ClonalColony, Identifier 
    output2=data.frame(t(dat))
    output3=data.frame(t(tim))
    write.table(output1,file="Lawless_data.txt",append=TRUE,col.names=FALSE,row.names=FALSE)
    write.table(output2,file="Lawless_area.txt",append=TRUE,col.names=FALSE,row.names=FALSE)
    write.table(output3,file="Lawless_time.txt",append=TRUE,col.names=FALSE,row.names=FALSE)
  }
  
  # sorter=order(rList,decreasing=TRUE) #decreasing growth rate indices in list 
  # rsort=rList[sorter][1:(length(sorter)*0.10)] #extracting the 10% highest growth rates 
  # QFAsim=mean(rsort) #Not sure what this is telling me and what it is useful for...????? 
  # 
  # return(data.frame(Pin=strname,A=AList,r=rList,QFA=rep(QFAsim,length(AList)),Strain=rep(strgene,length(AList))))
  #This stuff is not actually need at this point in time, but probably useful for later
  #probably best to recalculate growth rates for all three data sets later so that the same calculations have been used 
}

getORF=makeORFGetter("LibraryDescription.txt") #Not sure I quite understand this yet! 

# Fit exponential model to all microcolonies
flist=list.files(pattern="*_OUT.txt")
#dat=extractStuff(flist[1])
for (f in flist[2:length(flist)]){
  params=extractStuff(f)
  #dat=rbind(dat,params)
}
