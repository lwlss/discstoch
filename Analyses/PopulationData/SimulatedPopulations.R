# This script is to analyze Population Level Simulations
# and produces all of the main Figures submitted in the Dissertation

library(data.table)

############################################## Data #################################################

dataset<-function(x){
  if (x == "Lawless"){
    # DataSet1: Lawless
    area=fread("~/discstoch/Analyses/LineageData/Lawless_area_shortTC.txt",header=FALSE)
    times=fread("~/discstoch/Analyses/LineageData/Lawless_time_shortTC.txt",header=FALSE)
    data=fread("~/discstoch/Analyses/LineageData/Lawless_data_shortTC.txt",header=FALSE) #3rd column (Identifier) => strain_parentcolony
    names(data)=c("genotype","clonalcolony","identifier","blobnumber")
    return(list("area"=area,"data"=data,"times"=times))
  }
  else if (x == "Levy"){
    # DataSet2: Levy
    area=fread("Levy_area_NF.txt",header=TRUE)
    times=fread("Levy_time_NF.txt",header=TRUE)
    data=fread("Levy_data_NF.txt",header=TRUE) #3rd column (Identifier) => replicate
    #info=read.table("Levy_GrowthRateInfo.txt",header=TRUE,row.names=1)
    return(list("area"=area,"data"=data,"times"=times))
  }
  else if (x == "Ziv"){
    # DataSet3: Ziv
    area=fread("Ziv_area_filtered1.txt",header=FALSE)
    times=fread("Ziv_times_filtered1.txt",header=FALSE)
    data=fread("Ziv_data_filtered1.txt",header=TRUE) #3rd column (Identifier) => colony
    info=read.table("Ziv_GrowthRateInfo.txt",header=TRUE,row.names=1)
    return(list("area"=area,"data"=data,"times"=times,"info"=info))
  }
  else {print("Not a valid dataset")}
}

# Choosing a data set
datsetname="Lawless"
x=dataset(datsetname)
area=x$area
times=x$times
data=x$data

# Deterministic Bayesian Parameters (currently logistic -> change to exponential ones)
params=c("1:100","101:200","201:300","301:400","401:500","501:600","601:700","701:800","801:900",
         "901:1000","1001:1100","1101:1200","1201:1300","1301:1400","1401:1500","1501:1600",
         "1601:1700","1701:1800","1801:1900")

total_params=matrix(0,ncol=2,nrow=length(params)*100+46)
for (i in 1:length(params)){
  total_params[(1:100)+(100*(i-1)),]=as.matrix(read.table(paste("~/discstoch/Analyses/LineageData/Lawless_Bayes_parameters_BayesDetExp_",params[i],".txt",sep=""),header=TRUE))
}
total_params[1901:1946,]=as.matrix(read.table("~/discstoch/Analyses/LineageData/Lawless_Bayes_parameters_BayesDetExp_1901:1946.txt",header=TRUE))
total_rates_bayes=total_params[,2]
rates=total_rates_bayes

strain_names=unique(data$genotype)
pickstrain="HIS3"
strain=detstocgrowth::subset_strain(data,area,times,pickstrain)
strain_rates=rates[which(data$genotype==pickstrain)]


################################################ Colour Map ########################################################

# Producing a colour map according to which the rest of the simulation are going to be colour-coded
par(mfrow=c(1,1))
growthraterange=seq(0,0.5,0.001)
yellowredmap(growthraterange)


######################################################### FishPlots ################################################

# Sampling strain rates and plotting the population compositions using the fishplot package
set.seed(1)
N=10000
sample_rates=sample(strain_rates,N,replace=TRUE)
time=seq(1,48)
pop_comoposition(sample_rates,time,growthraterange)


################################################## Percentile Lineage Contributions ################################

t=seq(1,48,1)
iter=1000
StartPops=c(50,100,500,1000,5000,10000)
#StartPops=c(50,100,500)
for (N in StartPops){
  newinfo=pop_composition_dat(t,iter,N,strain_rates)
  if (N==StartPops[1]){
    info=newinfo
  } else{
    info=data.frame(info,newinfo)
  }
}

# Plotting the Number of Lineages which make up more than 5% over time
colours=rep(c(rainbow(6)))
par(mfrow=c(1,2))
info=as.matrix(info)
plot(NULL,xlim=range(t),ylim=c(0,10),main=paste("Lineages which make up more \nthan 5% of the", strain$name, "Population"),
     xlab="Time (h)", ylab="Number of Lineages",cex.lab=1.4)
for(i in c(1,3,5,7,9)){
  lines(t,info[,i],lty=1,col=colours[i],lwd=3)
}
plot(NULL,xlim=range(t),ylim=c(0,100),main=paste("Dominance ratios for the fastest lineage \n in the",strain$name,"Population"),
     xlab="Time (h)", ylab="Percent",cex.lab=1.4)
for(i in c(2,4,6,8,10)){
  lines(t,info[,i],lty=1,col=colours[i-1],lwd=3)
}
legend("topright",legend=c("N=50","N=100","N=500","N=1000","N=5000","N=10000"),lwd=rep(3,6),col=colours[1:6])
#write.table(info,file=paste("MorePlotsData_",strain$name,".txt",sep=""),col.names=TRUE,row.names=FALSE)


################################################## Population Simulations ##########################################


iterations=100
yl=10^22
gr=seq(0,0.5,0.005)
N=1000
t=seq(1,48)
simpopdat=pop_sim_plot(strain,strain_rates,gr,yl,iterations,t)

print(pickstrain)
print(mean(simpopdat$bp_pr))
print(min(simpopdat$bp_pr))
print(time[mean(simpopdat$bp_loc)])
print(mean(simpopdat$pop_rates))
print(mean(simpopdat$pop_rates[-which(is.na(simpopdat$pop_rates_lag))]))
print(mean(simpopdat$pop_rates_lag,na.rm=TRUE))
tt=t.test(simpopdat$pop_rates,simpopdat$pop_rates_lag,conf.level=0.99,paired=TRUE)
print(tt$p.value)

h=hist(strain_rates,breaks=gr,plot=FALSE)
hcol=colours[sapply(h$mids, function(x) which(abs(gr_range-x)==min(abs(gr_range-x)))[1])]
plot(h,col=hcol,cex.lab=1.4)
abline(v=mean(strain_rates),col="red",lwd=3)
abline(v=mean(simpopdat$pop_rates),col="black",lwd=3)


##################################################### Lag Phase Estimates ##########################################

N=seq(100,10000,100)
time=seq(0,48,1)
iterations=100
strainlag=lagduration(strain=0,strain_rates,t,N,iterations)

# Plotting lag duration for the six different distributions
par(mfrow=c(1,1))
plot(c(N),c(strainlag$m),ylim=c(15,50),xlab="Population Start Size (No. of Cells)",
     ylab="Apparent Lag Phase",type='l',cex.lab=1.4,col="cyan",lwd=2,
     main="Break Point Estimates over Time")
lines(c(N),strainlag$uc,lty=2,col=adjustcolor("red",0.6))
lines(c(N),strainlag$lc,lty=2,col=adjustcolor("red",0.6))


################################################## Mean and Variance Plots #########################################

N=10000
iterations=100
t=seq(0,48,1)
init_pop=seq(1,N,1)

sizerates=meanvarsim(strain,strain_rates,N,iterations,t)

variance=c()
mean=c()
meanlag=c()
variancelag=c()
for (i in 1:length(seq(1,N,1))){
  variance=c(variance,var(sizerates$all_simpoprate[[i]],na.rm=TRUE))
  mean=c(mean,mean(sizerates$all_simpoprate[[i]],na.rm=TRUE))
  variancelag=c(variancelag,var(sizerates$all_simpopratelag[[i]],na.rm=TRUE))
  meanlag=c(meanlag,mean(sizerates$all_simpopratelag[[i]],na.rm=TRUE))
}
true_mean=mean(strain_rates)
true_variance=var(strain_rates)

#To get rid of -INF values in variance when taking the log
variance[which(variance==0)]=NA

#http://www.r-bloggers.com/r-single-plot-with-two-different-y-axes/
#png(paste(paste("Final_PopSimDat",datsetname,strain$name,signif(true_mean,3),signif(true_variance,3),sep="_"),".png",sep=""))
par(mar=c(5,5,2,5))
plot(NULL,xlim=range(init_pop),ylim=c(0,0.65),ylab="Growth Rate (1/h)",xlab="Initial Population Size (No. of cells)",cex.lab=1.4,cex.main=1.2)
for (i in 1:length(all_simpoprate)){
  points(rep(init_pop[i],iterations),all_simpoprate[[i]],col=adjustcolor("grey",0.7))
}
lines(init_pop,mean,col=adjustcolor("darkgreen",0.7),lwd=2)
lines(init_pop,meanlag,col=adjustcolor("green",0.7),lwd=2)
abline(h=true_mean,col=adjustcolor("red",0.7),lwd=3)
par(new=T)
d=data.frame(x=init_pop,v=variance,m=mean)
with(d,plot(x,v,type='l',col=adjustcolor("purple",0.7),ylim=c(0,max(variance,na.rm=TRUE)),axes=F,xlab=NA,ylab=NA,cex=1.2),lwd=2)
axis(side=4)
mtext(side=4,line=3,"Variance",cex=1.4)
lines(init_pop,variancelag,col=adjustcolor("blue",0.7),lwd=2)
legend(200,0.0085,
       legend=c("Mean Lag Growth","Mean Exponential Growth","True Mean","Lag Variance","Exponential Variance","Pop'n Simulations"),
       lty=c(1,1,1,1,1,0),pch=c(NA,NA,NA,NA,NA,1),lwd=c(3,3,3,3,3,3),col=c("green","darkgreen","red","blue","purple","grey"),cex=0.9)
#dev.off()
