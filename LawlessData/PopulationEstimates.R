# Population Area Estimates

library(segmented)
library(detstocgrowth)

area=read.table("PopulationArea.txt",header=FALSE)
times=read.table("PopulationTime.txt",header=FALSE)
folders=read.table("PopulationFolders.txt",header=FALSE)
folders=as.vector(t(folders))
data=read.table("~/BayesianInference/Lawless_data_shortTC.txt",header=FALSE)
names(data)=c("genotype","clonalcolony","identifier","blobnumber")
strains=c()
for (i in folders){
  strains=c(strains,as.vector(unique(data[which(data$clonalcolony==i),]$genotype)))
}
strain="HTZ1"
sf=folders[which(strains==strain)]
sa=area[which(strains==strain),]
st=times[which(strains==strain),]

plot(NULL,xlim=c(0,28),ylim=c(9.5,14),xlab="Time (h)", ylab="log(Area)",
     main=paste("Population Estimates for",strain),cex.lab=1.4)

bp1=c()
bp2=c()
s1=c()
s2=c()
s3=c()

for (i in 1:dim(sa)[1]){
  a=log(as.numeric(sa[i,]))
  t=as.numeric(st[i,])
  linmod=lm(a~t)
  segmod=segmented(linmod,seg.Z =~t,psi=c(15,25))
  points(t,a,pch=16,col=adjustcolor("red",0.2))
  plot(segmod,add=T)
  bp1=c(bp1,segmod$psi[1,2])
  bp2=c(bp2,segmod$psi[2,2])
  s1=c(s1,slope(segmod)$t[1,1])
  s2=c(s2,slope(segmod)$t[2,1])
  s3=c(s3,slope(segmod)$t[3,1])
}

print(strain)
print(mean(bp1))
print(mean(bp2))
print(mean(s1)); print(mean(s2)); print(mean(s3))
print(mean(as.numeric(sa[,1]))/92.5)

params=c("1:100","101:200","201:300","301:400","401:500","501:600","601:700","701:800","801:900",
         "901:1000","1001:1100","1101:1200","1201:1300","1301:1400","1401:1500","1501:1600",
         "1601:1700","1701:1800","1801:1900")
total_params=matrix(0,ncol=2,nrow=length(params)*100+46)
for (i in 1:length(params)){
  #print((1:100)+(100*(i-1)))
  print(params[i])
  total_params[(1:100)+(100*(i-1)),]=as.matrix(read.table(paste("~/BayesianInference/Lawless_Bayes_parameters_BayesDetExp_",params[i],".txt",sep=""),header=TRUE))
}
total_params[1901:1946,]=as.matrix(read.table("~/BayesianInference/Lawless_Bayes_parameters_BayesDetExp_1901:1946.txt",header=TRUE))
total_rates_bayes=total_params[,2]
rates=total_rates_bayes

lineagearea=read.table("Lawless_area_shortTC.txt",header=FALSE)
lineagetime=read.table("Lawless_time_shortTC.txt",header=FALSE)

sortedfolders=as.vector(unique(data$clonalcolony))

xax=13
# pdf(file="PinGrowthRates.pdf",title="Observed and Simulated Population Growth",width=20,height=28)
# par(mfrow=c(8,8))
# plot(NULL)
for (f in sortedfolders[18]){
  print(f)
  # True Population Observation
  plot(NULL,xlim=c(0,19),ylim=c(9.5,xax),xlab="Time (h)", ylab="No. of Cells",axes=F,
       main=paste("Observed and Simulated Population \nGrowth in Pin",f),cex.lab=1.4)
  sa=area[which(folders==f),1:60]
  st=times[which(folders==f),1:60]
  a=log(as.numeric(sa))
  t=as.numeric(st)
  linmod=lm(a~t)
  segmod=segmented(linmod,seg.Z =~t,psi=c(15))
  points(t,a,pch=16,col=adjustcolor("red",0.7))
  plot(segmod,add=T)
  # Population Observation from Single Lineages 
  N=as.numeric(sa[1]/92.5)
  #N=exp(a[1])
  time=seq(0,19,1)
  iterations=1000
  folder=list()
  folder$area=lineagearea[which(data$clonalcolony==f),]
  folder$times=lineagetime[which(data$clonalcolony==f),]
  folder_rates=rates[which(data$clonalcolony==f)]
  # plot(1, type='n',xlim=c(0,20), ylim=c(N,10^5),xlab="Time (h)",
  #      ylab="No. of Cells",log='y',cex.lab=1.4, main=paste("Population Simulations Using \nSingle Linege Growth Rates \nObtained from",f))
  for (i in 1:iterations){
    popsim=pop_sim_dat(folder,folder_rates,N,time,1)
    popdata=popsim$PopData
    poprate=popsim$SimRates
    lines(time,log(popdata*92.5),col=adjustcolor("darkgoldenrod1",0.02))
  }
  box()
  axis(side=2,at=seq(10,xax),labels=c(round(exp(seq(10,xax))/92.5)))
  axis(side=1)
  legend("topleft",legend=c("Observed","Fit","Simulated"),pch=c(16,NA,NA),lty=c(NA,1,1),col=c("red","black","darkgoldenrod1"),lwd=c(NA,2,3))
}
# dev.off()
