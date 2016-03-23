#Exploring the data by Ziv et al. 2013
#setwd("~/GitHub/discstoch/LevySiegalData/LWOdata")

data=read.table("LWOdata.txt", header=TRUE)
area=read.table("LWOarea.txt", header=TRUE)

dim(area)
dim(data)

head(area)
head(data)

names(area) #time points; see paper for what the intervals are... 
names(data) #is GR the growth rate? How has this been calculated? see paper... 

plates=unique(data$plate)
strains=unique(data$Strain)
wells=unique(data$well)

paste("Number of plates is", length(plates))
paste("Number of Strains is",length(strains))
paste("Number of wells is", length(wells))

for (i in 1:length(strains))
{
  plate_data<-subset(data,plate==plates[i],select=c(Strain,well))
  print(paste("Number of strains in plate", i, "is", length(unique(plate_data$Strain))))
  print(paste("Number of wells in plate", i, "is", length(unique(plate_data$well))))
}

#Overlayed growth curves

pdf(height = 16, width = 16, file = "ClonalGrowthCurves.pdf")
par(mfrow=c(3,3))

for(j in 1:length(wells))
{
  indices=which(data$well == wells[j])
  l_area=as.matrix(log(area[indices,]))
  l_area[which(l_area<0)]<-0 
  
  #finding the associated strain and plate
  w_data=subset(data,well==wells[j], select=c(Strain,plate))
  paint=c("black","darkblue","darkred","darkgreen")
  
  #plotting clonal growth curves so that line colour indicates strain
  #black: FY4-5; blue: OakBC248; red: WineBC241; green: HybBC252
  plot(j,type='n',xlim=c(0,dim(l_area)[2]), ylim=c(0,max(l_area)), xlab="Time (h)", ylab="Area (px)",
       main=paste("Log-linear Growth Curves"), cex.main=2, cex.lab=1.5, cex.axis=1.5)
  
  for (i in 1:dim(l_area)[1])
  {
    lines(l_area[i,], col=paint[which(strains==unique(w_data$Strain))])
  }
  #putting associated well and plate number in the legend
  legend("topleft",legend=c(wells[j],unique(w_data$plate)),cex=1.5)
  
}

dev.off()


