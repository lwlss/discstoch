# Growth Rate Distributions
library(data.table)
library(detstocgrowth)

dataset<-function(x){
  if (x == "Lawless"){
    # DataSet1: Lawless
    area=fread("Lawless_area_shortTC.txt",header=FALSE)
    times=fread("Lawless_time_shortTC.txt",header=FALSE)
    data=fread("Lawless_data_shortTC.txt",header=FALSE) #3rd column (Identifier) => strain_parentcolony
    names(data)=c("genotype","clonalcolony","identifier","blobnumber")
    return(list("area"=area,"data"=data,"times"=times,"residuals"=residuals))
  }
  else if (x == "Levy"){
    # DataSet2: Levy
    area=fread("Levy_area_filtered.txt",header=FALSE)
    times=fread("Levy_times_filtered.txt",header=FALSE)
    data=fread("Levy_data_filtered.txt",header=TRUE) #3rd column (Identifier) => replicate
    return(list("area"=area,"data"=data,"times"=times,"residuals"=residuals))
  }
  else if (x == "Ziv"){
    # DataSet3: Ziv
    area=fread("Ziv_area_filtered1.txt",header=FALSE)
    times=fread("Ziv_times_filtered1.txt",header=FALSE)
    data=fread("Ziv_data_filtered1.txt",header=TRUE) #3rd column (Identifier) => colony
    return(list("area"=area,"data"=data,"times"=times,"residuals"=residuals))
  }
  else {print("Not a valid dataset")}
}

# Choosing a data set
datsetname="Levy"
x=dataset(datsetname)
area=x$area
times=x$times
data=x$data

# Data
strain_names=unique(data$genotype)
print(strain_names)
pickstrain=strain_names[14] #choose strain here!
strain=detstocgrowth::subset_strain(data,area,times,pickstrain)

#pdf(height=5, width=10, file=paste("GrowthRateDist_",strain$name,".pdf",sep=""))

resid1=matrix(0,nrow=dim(area)[1],ncol=dim(area)[2])
resid2=matrix(0,nrow=dim(area)[1],ncol=dim(area)[2])
lm_r=c()
lm_i=c()
e_nlm_r=c()
e_nlm_A=c()
for (i in 1:dim(strain$area)[1]){
  op=par(mfrow=c(1,2))
  # Log-linear model using lm
  fit=detstocgrowth::LM_growthrate(strain$area[i,],strain$times[i,])$fit
  k=detstocgrowth::LM_growthrate(strain$area[i,],strain$times[i,])$rate
  int=detstocgrowth::LM_growthrate(strain$area[i,],strain$times[i,])$int
  resid1[i,]=abs(residuals(fit))
  lm_r=c(lm_r,k)
  lm_i=c(lm_i,int)
  # Exponential model using nlm
  result=nls(y~A*exp(r*x),data=data.frame(x=as.numeric(strain$times[i,]),
                                          y=as.numeric(strain$area[i,])),
             start=list(A=strain$area[i,]$V1,r=0.1))
  res=result$m$getPars()
  resid2[i,]=abs(residuals(result))
  A=as.numeric(res[1])
  r=max(0,as.numeric(res[2]))
  e_nlm_r=c(e_nlm_r,r)
  e_nlm_A=c(e_nlm_A,A)
  # Plotting on the log-scale
  # plot(as.numeric(strain$times[i,]),log(as.numeric(strain$area[i,])),
  #      ylab="log(Area)",xlab="Time (h)",main=paste("Growth Curve",i))
  # lines(as.numeric(strain$times[i,]),predict(fit),col="blue",lwd=3)
  # lines(as.numeric(strain$times[i,]),log(predict(result)),col="red",lwd=3)
  # legend("topleft",legend=c(paste("log-linear, r=",signif(k,2),sep=""),
  #                           paste("exponential, r=",signif(r,2),sep="")),
  #        lty=c(1,1), lwd=c(3,3),col=c("blue","red"))
  # # Plotting on the original scale
  # plot(as.numeric(strain$times[i,]),as.numeric(strain$area[i,]),
  #      ylab="Area",xlab="Time (h)",main=paste("Growth Curve",i))
  # lines(as.numeric(strain$times[i,]),exp(predict(fit)),col="blue",lwd=3)
  # lines(as.numeric(strain$times[i,]),predict(result),col="red",lwd=3)
  # legend("topleft",legend=c(paste("log-linear, r=",signif(k,2),sep=""),
  #                           paste("exponential, r=",signif(r,2),sep="")),
  #        lty=c(1,1), lwd=c(3,3),col=c("blue","red"))
}
par(op)

lm_r[which(lm_r<0)]=0

op=par(mfrow=c(2,1))
plot(NULL,xlim=c(0,max(area)),ylim=c(0,max(resid1)),main=paste("Residual Magnitude of Log. Lin. Model vs Colony Size for",strain$name),
     xlab="Colony Size (Area)", ylab="Absolute Residual Size",cex.lab=1.4,cex.main=1.2)
for (i in 1:dim(resid1)[1]){
  points(as.numeric(strain$area[i,]),as.numeric(resid1[i,]),col=adjustcolor("red",0.3),pch=19)
}
plot(NULL,xlim=c(0,max(area)),ylim=c(0,max(resid2)),main=paste("Residual Magnitude of Exp. Model vs Colony Size for",strain$name),
     xlab="Colony Size (Area)", ylab="Absolute Residual Size",cex.lab=1.4,cex.main=1.2)
for (i in 1:dim(resid2)[1]){
  points(as.numeric(strain$area[i,]),as.numeric(resid2[i,]),col=adjustcolor("red",0.3),pch=19)
}
par(op)

# Plotting the distributions
op=par(mfrow=c(1,2))
detstocgrowth::histo(lm_r,paste("Log.Lin.Mod.:",strain$name),c=seq(0,0.6,0.01),SE=FALSE)
box()
detstocgrowth::histo(e_nlm_r,paste("Exp.Mod.:",strain$name),c=seq(0,0.6,0.01),SE=FALSE)
box()
par(op)
op=par(mfrow=c(1,2))
detstocgrowth::histo(lm_i,paste("Log.Lin.Mod.:",strain$name),c=seq(0,6,0.25),SE=FALSE,int=TRUE)
box()
detstocgrowth::histo(log(e_nlm_A),paste("Exp.Mod.:",strain$name),c=seq(0,6,0.25),SE=FALSE,int=TRUE)
box()
par(op)

#dev.off()

# Comparing the fastest strains
top=10
lmf=unique(sort(lm_r,decreasing=TRUE))[1:top]
nlmf=unique(sort(e_nlm_r,decreasing=TRUE))[1:top]
both=intersect(as.numeric(lapply(lmf,function(x) which(x==lm_r))),as.numeric(lapply(nlmf,function(x) which(x==e_nlm_r))))
print(length(both))
print(both)

strain$indices[both]
strain$indices[as.numeric(lapply(lmf,function(x) which(x==lm_r)))]

