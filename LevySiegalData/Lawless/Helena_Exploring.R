#Exploring the Data by Levy et al. 2012 
#setwd("~/GitHub/discstoch/LevySiegalData/Lawless")

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
time=rep$times[w,]

#converting area to logs and getting rid of -Inf values 
l_area=log(area)
l_area[which(l_area<0)]<-0

output=par(mfrow=c(2,1))
#Plotting growth curves for clonal population
plot(1, type='n', xlim=c(0,max(time)), ylim=c(0,max(l_area)), xlab="Time (h)", ylab="Area (px)", 
     main=paste("Log-Linear Growth Curves for Replicate 4 , Strain", condition_names[[10]], ", Well", strain_wells[1]))
for (i in 1:dim(time)[1])
{
  lines(time[i,],l_area[i,])
}

#looking at the entire strain within that replicate
plot(2, type='n', xlim=c(0,max(time)), ylim=c(0,10), xlab="Time (h)", ylab="Area (px)", 
     main=paste("Log-Linear Growth Curves for Replicate 4 , Strain", condition_names[[10]]))
for (k in 1:length(strain_wells))
{
  w=rep$well.list[[strain_wells[k]]]
  area=rep$areas[w,]
  time=rep$times[w,]
  l_area=log(area)
  l_area[which(l_area<0)]<-0
  for (i in 1:dim(time)[1])
  {
    lines(time[i,],l_area[i,])
  }
}
par(output)
