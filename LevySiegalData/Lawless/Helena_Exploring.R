#Exploring the Data by Levy et al. 2012 
#setwd("~/GitHub/discstoch/LevySiegalData/Lawless")

library(grDevices)

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

#Checking the sort of information that each replicate contains
#Do all replicates contain the same sort of information? 

for (i in 1:length(replicates))
{
  n=length(replicates[[i]])
  print(n)
  print(names(replicates[[i]]))
  
  for (k in 1:n){
    print(names(replicates[[i]][k]))
    print(names(replicates[[i]][[k]]))
    print(" ")
  }
  
}

#Isogenic growth curves within the same replicates
rep=rep4

#conditions in the replicate 
conditions=as.list(rep$condition.names)
condition_names=unique(conditions)

#wells for condition
strain_wells=names(conditions[which(conditions==condition_names[[10]])])

#focusing on one well for now
w=rep$well.list[[strain_wells[1]]]
area=rep$areas[w,]
time=as.numeric(rep$times[w[1],]) # Assume that all timepoints for one well are the same

# Convert zero observations to NaN (Not a Number)
area[area==0]=NA

#converting area to logs 
l_area=log(area)

# Too many growth curves to see clearly
Nsample=20

output=par(mfrow=c(2,1))
#Plotting growth curves for clonal population
plot(NULL, xlim=c(0,max(time)), ylim=range(area,na.rm=TRUE), xlab="Time (h)", ylab="Area (px)", log="y",
     main=paste("Log-Linear Growth Curves for Replicate 4 , Strain", condition_names[[10]], ", Well", strain_wells[1]))
for (i in sample(1:dim(area)[1],Nsample,replace=FALSE))
{
  # Try transparency to see pattern in hundreds of growth curves
  lines(time,area[i,],col=rgb(0,0,0,0.3),lwd=2)
}

#looking at the entire strain within that replicate
plot(NULL, xlim=c(0,max(time)), ylim=range(area,na.rm=TRUE), xlab="Time (h)", ylab="Area (px)", log="y",
     main=paste("Log-Linear Growth Curves for Replicate 4 , Strain", condition_names[[10]]))
rect(par("usr")[1],par("usr")[3],par("usr")[2],2*max(area,na.rm=TRUE),col = "lightgrey")
Nstrain_wells=length(strain_wells)
cols=adjustcolor(rainbow(Nstrain_wells),0.1)
for (k in 1:Nstrain_wells)
{
  w=rep$well.list[[strain_wells[k]]]
  area=rep$areas[w,]
  time=as.numeric(rep$times[w[1],])
  area[area==0]=NA
  l_area=log(area)
  for (i in sample(1:dim(area)[1],Nsample,replace=FALSE))
  {
    lines(time,area[i,],lwd=2,col=cols[k])
  }
}
par(output)
