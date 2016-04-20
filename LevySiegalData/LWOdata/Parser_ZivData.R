#PARSER FOR ZIV DATA 
#To bring all three data sets into the same format 
#so that a single script can be used to extract growth curves 

data=read.table("LWOdata.txt", header=TRUE)
area=read.table("LWOarea.txt", header=TRUE)

#Images were taken every hour for 24h 
tim=matrix(rep(seq(1,23),dim(area)[1]),ncol=23,nrow=dim(area)[1],byrow=TRUE)
tim=as.data.frame(tim)

SubData=data.frame(data$Strain, data$well, data$colony) #Genotype, ClonalColony, Identifier 

write.table(area,"Ziv_area.txt",col.names=FALSE,row.names=FALSE)
write.table(tim,"Ziv_time.txt",col.names=FALSE,row.names=FALSE)
write.table(SubData,"Ziv_data.txt",col.names=FALSE,row.names=FALSE)
