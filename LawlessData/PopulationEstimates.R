# Population Area Estimates

library(segmented)

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
strain="HIS3"
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
