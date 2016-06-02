#Backwards Parser to easily access data associated with folder names

data=read.table("Lawless_data.txt")

# Revert pinid naming convention 
pinid=function(pid){
  Col=as.numeric(substr(pid,5,6))
  OldCol=1-Col+24
  return(paste(substr(pid,1,4),sprintf("%02d",OldCol),sep=""))
}

folders=unique(data$V2)
data_namespace=matrix(0,nrow=length(folders),ncol=3)
for (i in 1:length(folders)){
  indices=which(data$V2==folders[i])
  strain=unique(data[indices,]$V1)
  identifier=unique(data[indices,]$V3)
  data_namespace[i,]=c(pinid(folders[i]),toString(strain),toString(identifier))
}

write.table(data_namespace,"Parser_Lawless_Namespace.txt",row.names=FALSE,col.names=FALSE)