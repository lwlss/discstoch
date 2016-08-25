# Comparing Bayesian Stochastic, Bayesian Deterministic

library(detstocgrowth)
library(data.table)

dataset<-function(x){
  if (x == "Lawless"){
    # DataSet1: Lawless
    area=fread("Lawless_area_shortTC.txt",header=FALSE)
    times=fread("Lawless_time_shortTC.txt",header=FALSE)
    data=fread("Lawless_data_shortTC.txt",header=FALSE) #3rd column (Identifier) => strain_parentcolony
    return(list("area"=area,"data"=data,"times"=times))
  }
  else if (x == "Levy"){
    # DataSet2: Levy
    area=fread("Levy_area_filtered1.txt",header=FALSE)
    times=fread("Levy_times_filtered1.txt",header=FALSE)
    data=fread("Levy_data_filtered1.txt",header=FALSE) #3rd column (Identifier) => replicate
    return(list("area"=area,"data"=data,"times"=times))
  }
  else if (x == "Ziv"){
    # DataSet3: Ziv
    area=fread("Ziv_area_filtered1.txt",header=FALSE)
    times=fread("Ziv_times_filtered1.txt",header=FALSE)
    data=fread("Ziv_data_filtered1.txt",header=FALSE) #3rd column (Identifier) => colony
    return(list("area"=area,"data"=data,"times"=times))
  }
  else {print("Not a valid dataset")}
}

# Choosing a data set
datsetname="Lawless"

# Extracting data for the growth curve(s) on which to do inference
gc=252
x=dataset(datsetname)
area=x$area[gc,]
times=x$times[gc,]
data=x$data[gc,]

# Comparing the Deterministic and the Stochastic Inference
Det=read.table("GC252DetBayesExp.txt",header=TRUE)
Stoch=read.table("GC252HybBayesLog.txt",header=TRUE)
Stochest=mean(Stoch$r)
Detest=mean(Det$r)
Detlogest=LM_growthrate(area,times)$fit

#png(filename="CompStochDetN10000")
par(mfrow=c(1,1))
plot(as.numeric(times),as.numeric(area/min(area)),xlab="Time (h)",ylab="Cell Count",
     main="Comparing Stochastic and Deterministic Inference",cex.lab=1.4,pch=4)
for (i in 1:10000){
  s=simCellsStoch(50,sample(Stoch$r,1),1)
  lines(s$t,s$c,col=adjustcolor("darkgreen",0.01))
  d=detLog(1000000,sample(Det$r,1),mean(Det$x0)/min(area),c(as.numeric(times),7))
  lines(c(as.numeric(times),7),d,col=adjustcolor("red",0.01))
}
dlog=detLog(1000000,Detlogest$coefficients[2],Detlogest$coefficients[1],as.numeric(times))
lines(as.numeric(times),dlog,col="blue",lwd=2)
lines(as.numeric(times),as.numeric(area/min(area)),cex.lab=1.4,pch=4,type='p')
lines(as.numeric(times),as.numeric(area/min(area)),cex.lab=1.4,type='l')
legend("topleft",c("Bayes. Stoch.","Bayes. Det. (exp.)","Freq. Det. (log-lin.)"),lty=c(1,1,1),col=c("darkgreen","red","blue"),lwd=c(2,2,2),bty='n')
#dev.off()
