#Exploring the Data by Levy et al. 2012 
#setwd("~/GitHub/discstoch/LevySiegalData/Lawless")

library(grDevices)

############Exploring###############

load("d.110601.Rfile") #rep1
rep1=d
load("d.110612.Rfile") #rep3
rep3=d
load("d.110607.Rfile") #rep2
rep2=d
load("d.110625.Rfile")
rep5=d
load("d.110706.Rfile") #rep4
rep4=d

replicates=list(rep1,rep2,rep3,rep4,rep5)

###############Plotting###############

#Isogenic growth curves within the same replicates
rep=rep4

#conditions in the replicate 
conditions=as.list(rep$condition.names)
condition_names=unique(conditions)

#wells for condition
strain_wells=names(conditions[which(conditions==condition_names[[10]])])


extract_well <- function(Nwell,rep,sw){
  w=rep$well.list[[sw[Nwell]]]
  area=rep$areas[w,]
  time=as.numeric(rep$times[w[Nwell],]) # Assume that all timepoints for one well are the same
  area[area==0]=NA
  return(list("area"=area,"time"=time,"w"=w))
}

#focusing on one well for now
data=extract_well(1,rep4,strain_wells)
area=data$area
time=data$time
w=data$w


# Too many growth curves to see clearly
Nsample=100

output=par(mfrow=c(1,2))

#Plotting growth curves for one particular clonal population
plot(1,type='n', xlim=c(0,max(time)), ylim=range(area,na.rm=TRUE), xlab="Time (h)", ylab="Microcolony Area (px)",
     main="CTF4_4_1")

k=c()

for (i in sample(1:dim(area)[1],Nsample,replace=FALSE)) #sampling 20 growth curves of 1536 possible ones in that well 
{
  # Try transparency to see pattern in hundreds of growth curves
  lines(time,area[i,],col=rgb(0,0,0,0.3),lwd=2) #black transparent 
  
  #Linear Regression to find Growth rate 
  fit<-lm(log(area[i,]) ~ time)
  rate=fit$coefficient[[2]]
  k=c(k,rate)
}

#Growth rate Distribution
lo=trunc(min(k)*10)/10-0.1 #rounding down
hi=trunc(max(k)*10)/10+0.1 #rounding up 
cells=seq(lo,hi,0.01)
hist(k,breaks=cells,xlab="r (1/h)", main=paste(Nsample," microcolonies"))

par(output)




#Exploring the data by Ziv et al. 2013
#setwd("~/GitHub/discstoch/LevySiegalData/LWOdata")

library(grDevices)

data=read.table("LWOdata.txt", header=TRUE)
area=read.table("LWOarea.txt", header=TRUE)

plates=unique(data$plate)
strains=unique(data$Strain)
wells=unique(data$well)


#Overlayed growth curves

#Convert zero observations to NaN (Not a Number)
area[area==0]=NA

# Too many growth curves to see clearly
Nsample=80

indices=which(data$well == wells[196])
j_area=as.matrix(area[indices,])
  
#cannot take sample larger than the population when replace=FALSE
if (dim(j_area)[196] < Nsample) {
  Nsample=dim(j_area)[1]
}
  
#finding the associated strain and plate
w_data=subset(data,well==wells[196], select=c(Strain,plate))

output=par(mfrow=c(1,2))
  
#plotting clonal growth curves so that line colour indicates strain
plot(1,type='n',xlim=c(0,dim(j_area)[2]), ylim=range(j_area,na.rm=TRUE), xlab="Time (h)", ylab="Microcolony Area (px)",
     main=paste("Y4-5_196_110316"))

k<-c()
time=seq(1,23,1)

for (i in sample(1:dim(j_area)[1],Nsample,replace=FALSE))
{
  lines(j_area[i,],col=rgb(0,0,0,0.3),lwd=2)
  
  #Linear Regression to find Growth rate 
  fit<-lm(log(j_area[i,]) ~ time)
  rate=fit$coefficient[[2]]
  k=c(k,rate)
}

#Growth rate Distribution
lo=trunc(min(k)*10)/10-0.1 #rounding down
hi=trunc(max(k)*10)/10+0.1 #rounding up 
cells=seq(lo,hi,0.01)
hist(k,breaks=cells,xlab="r (1/h)", main=paste(Nsample," microcolonies"))

par(output)




