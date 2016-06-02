#Parser for Lawless data set upon improved image analysis

names=read.table("Parser_Lawless_Namespace.txt",header=FALSE)
names=names[-1,] #no data for R03C03; plates all blurry
names=as.matrix(names)

folders=names[,1]

for (i in 1:length(folders)){
  area=read.table(paste(folders[i],"_AREA.txt",sep=""),header=FALSE)
  time=read.table(paste(folders[i],"_TIME.txt",sep=""),header=FALSE)
  data=matrix(0,nrow=dim(area)[1],ncol=3)
  #Start time at zero and convert from seconds to hours
  for (j in 1:dim(time)[1]){
    time[j,]=(time[j,]-time[j,1])/3600
    #Save folder data 
    data[j,]=names[which(names[,1]==folders[i]),]
  }
  #Store Data 
  data[,c(1,2)]=data[,c(2,1)]
  #names(data)=c("genotype","clonalcolony","identifier")
  write.table(area,"Lawless_area_revised.txt",col.names=FALSE,row.names=FALSE,append=TRUE)
  write.table(time,"Lawless_time_revised.txt",col.names=FALSE,row.names=FALSE,append=TRUE)
  write.table(data,"Lawless_data_revised.txt",col.names=FALSE,row.names=FALSE,append=TRUE)
}
