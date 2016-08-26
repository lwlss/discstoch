# Getting the precision vaue 

area=read.table("Lawless_area_shortTC.txt",header=FALSE)
times=read.table("Lawless_time_shortTC.txt",header=FALSE)

minsd1=100
minsd2=100
rss1=c()
se1=c()
rss2=c()
se2=c()
orires=matrix(0,ncol=length(times[1,]),nrow=dim(area)[1])
logres=matrix(0,ncol=length(times[1,]),nrow=dim(area)[1])
for (i in 1:dim(area)[1]){
  mod1=nls(y~A*exp(r*x),data=data.frame(x=as.numeric(times[i,]),y=as.numeric(area[i,])),start=list(A=area[i,1],r=0.2))
  mod2=lm(y~x,data=data.frame(x=as.numeric(times[i,]),y=as.numeric(log(area[i,]))))
  se1=c(se1,summary(mod1)$parameters[2,2])
  se2=c(se2,summary(mod2)$coefficients[2,2])
  rss1=c(rss1,sum(residuals(mod1)^2))
  rss2=c(rss2,sum(exp(residuals(mod2))^2))
  if(minsd1>sd(residuals(mod1))){
    minsd1=sd(residuals(mod1))
  }
  if(minsd2>sd(exp(residuals(mod2)))){
    minsd2=sd(exp(residuals(mod2)))
  }
  orires[i,]=residuals(mod1)
  logres[i,]=residuals(mod2)
}
print(minsd1)
print(minsd2)
print(range(se1))
print(range(se2))
print(range(rss1))
print(range(rss2))


op=par(mfrow=c(2,2))
plot(NULL,ylim=range(abs(orires)),xlim=c(range(area)),
     xlab="Colony Size (Area)",ylab="Absolute Residual Size",main="HIS3; Model Fit On The Original Scale",cex.lab=1.4)
for (i in which(data$genotype=="HIS3")){
  points(area[i,],abs(orires[i,]),col=adjustcolor("red",0.3))
}
plot(NULL,ylim=range(abs(logres)),xlim=c(range(area)),
     xlab="Colony Size (Area)",ylab="Absolute Residual Size",main="HIS3; Model Fit On The Log Scale",cex.lab=1.4)
for (i in which(data$genotype=="HIS3")){
  points(area[i,],abs(logres[i,]),col=adjustcolor("red",0.3))
}
plot(NULL,ylim=range(abs(orires)),xlim=c(range(area)),
     xlab="Colony Size (Area)",ylab="Absolute Residual Size",main="HTZ1; Model Fit On The Original Scale",cex.lab=1.4)
for (i in which(data$genotype=="HTZ1")){
  points(area[i,],abs(orires[i,]),col=adjustcolor("red",0.3))
}
plot(NULL,ylim=range(abs(logres)),xlim=c(range(area)),
     xlab="Colony Size (Area)",ylab="Absolute Residual Size",main="HTZ1; Model Fit On The Log Scale",cex.lab=1.4)
for (i in which(data$genotype=="HTZ1")){
  points(area[i,],abs(logres[i,]),col=adjustcolor("red",0.3))
}
par(op)