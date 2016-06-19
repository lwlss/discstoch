# Analysis Script Using the Package
# Frequentist, Deterministic, Exponential Modelling

library(data.table)
library(detstocgrowth)
library(fishplot)
library(MASS)

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
datsetname="Lawless"
x=dataset(datsetname)
area=x$area
times=x$times
data=x$data

# Choosing a strain and extracting the data for it
strain_names=unique(data$genotype)
print(strain_names)
pickstrain=strain_names[1] #choose strain here!
strain=detstocgrowth::subset_strain(data,area,times,pickstrain)
detstocgrowth::plot_growth(strain$area,strain$times,strain$name,Nsample=dim(strain$area)[1],title=TRUE,hist=TRUE)

#Calculating the estimated growth rates for all growth curves of the strain
strain_rates=c()
strain_int=c()
dist=c()
for (i in 1:dim(strain$area)[1]){
  k=detstocgrowth::LM_growthrate(strain$area[i,],strain$times[i,])$rate
  intercept=detstocgrowth::LM_growthrate(as.numeric(strain$area[i,]),as.numeric(strain$times[i,]))$int
  fit=detstocgrowth::LM_growthrate(as.numeric(strain$area[i,]),as.numeric(strain$times[i,]))$fit
  res=residuals(fit)
  dist=c(dist,range(res)[1]-range(res)[2])
  strain_rates=c(strain_rates,k)
  strain_int=c(strain_int,intercept)
}

#Setting growth rates <0 equal to 0
strain_rates[which(strain_rates<0)]=0

# FishPlot using data simulated from the growth rates over a longer time course and a
# starting population of a single cell
sample_rates=sample(strain_rates,8)
sample_rates=c(max(strain_rates),sample_rates,0)
# Simulating growth curves using the above rates
t=seq(1,48,1)
sim_area=matrix(0,ncol=length(t),nrow=length(sample_rates))
for (i in 1:length(sample_rates)){
  sim_area[i,]=round(detstocgrowth::simExponential(sample_rates[i],t))
}
pop_data=colSums(sim_area)
frac.table=matrix(0,nrow=dim(sim_area)[1],ncol=length(t))
for (i in 1:dim(sim_area)[1]){
  frac.table[i,]=as.numeric(sim_area[i,])/as.numeric(pop_data)
}
frac.table=frac.table*100
frac.table=trunc(frac.table) #need to round this down otherwise creatFishObjct() gives error
parents=rep(0,dim(sim_area)[1])
fish = fishplot::createFishObject(frac.table,parents,timepoints=t,col=rainbow(length(sample_rates)))
fish = fishplot::layoutClones(fish)
fishplot::fishPlot(fish,shape="spline",title.btm=paste(strain$name),
         cex.title=0.5, vlines=c(0,max(t)),
         vlab=c("1h","48h"))
title(paste("FishPlot for",datsetname,strain$name))

# Using 100 growth curves for extrapolation
total=dim(strain$area)[1]
N=100
time=seq(0,35,0.5)

# Refitting the exponential model to the newly simulated data to obtain
# population growth rate
iterations=20
popsim=pop_sim_dat(strain,strain_rates,N,time,iterations)
popdata=popsim$PopData

pop_rates=c()
for (i in 1:dim(popdata)[1]){
  k=detstocgrowth::LM_growthrate(as.numeric(popdata[i,]),time)$rate
  pop_rates=c(pop_rates,k)
}

op=par(mfrow=c(2,1))
detstocgrowth::pop_plot_growth(popdata,time,strain$name,title=TRUE,hist=FALSE)
detstocgrowth::histo(strain_rates,strain$name,SE=TRUE)
#abline(v=mean(strain_rates),col="red",lwd=3)
abline(v=mean(pop_rates),col="blue",lwd=3)
legend("topright",legend=c("True Mean","Pop'n Mean"),col=c("red","blue"),lty=c(1,1))
par(op)

#Note indices for all iterations!
popdata_indices=popsim$indices

# Fitting a piece-wise regression to the Population Simulations
piece_op=detstocgrowth::piecewise_pop_rate(popdata,time)

# Analysing how much each growth rate contributes to the pop'n growth during
# the exponential phase
# contribution=fast_rate_contribution(piece_op$k2,strain_rates[popdata_indices])
# filename=paste("Contribution_",datsetname,strain$name,"_0.1N_c",
#                signif(contribution, digits=3),".png",sep="")

# # png(filename)
# op=par(mfrow=c(3,1))
# detstocgrowth::plot_growth(strain$area[popdata_indices,],strain$times[popdata_indices,],
#             strain$name,Nsample=100)
# detstocgrowth::pop_plot_growth(popdata,time,strain$name,hist=FALSE)
# detstocgrowth::strain_pop_hist(strain$name,strain_rates[popdata_indices],piece_op$k1,piece_op$k2)
# par(op)
# # dev.off()

# Changes in variance & value of growth rate estimates with an inceasing
# size of the starting population
N=dim(strain$area)[1]
iterations=50
total_popsim=detstocgrowth::pop_sim_dat(strain,strain_rates,N,time,iterations)
simdata=total_popsim$SimData
poprates=total_popsim$SimRates

#Sum up one hundred more rows each time
#N=100
init_pop=seq(1,N,1)
all_simpoprate=list()
all_meansimrates=c()
for (i in 1:length(init_pop)){
  init_simmeanrates=c()
  init_simpoprate=c()
  for (j in 1:iterations){
    if(i==1){
      simitdata=simdata[[j]][1,]
    }else{simitdata=colSums(simdata[[j]][1:init_pop[i],])}
    simpoprate=detstocgrowth::LM_growthrate(simitdata,time)$rates #Pop'n Rates
    init_simpoprate=c(init_simpoprate,simpoprate)
    simrates=poprates[j:N*j] #Orignial Single Lineage Rates
    simmmeanrates=mean(simrates) #Mean Original Rates
    init_simmeanrates=c(init_simmeanrates,simmmeanrates)
  }
  all_simpoprate[[i]]=init_simpoprate
  all_meansimrates=c(all_meansimrates,init_simmeanrates)
}
#Use all_simpoprate to calculate the mean and variance for each element in the list
#See whether variance is decreasing with increasing pop'n size
#See whether different between pop'n mean rate and single lineage mean rate differs
variance=c()
mean=c()
for (i in 1:length(init_pop)){
  variance=c(variance,var(all_simpoprate[[i]]))
  mean=c(mean,mean(all_simpoprate[[i]]))
}
true_mean=mean(poprates)
true_variance=var(poprates)

#To get rid of -INF values in variance when taking the log
variance[which(variance==0)]=NA

#http://www.r-bloggers.com/r-single-plot-with-two-different-y-axes/
#png(paste(paste("PopSimDat",datsetname,strain$name,signif(true_mean,3),signif(true_variance,3),"IT1",sep="_"),".png",sep=""))
par(mar=c(5,5,2,5))
plot(NULL,xlim=range(init_pop),ylim=range(poprates),ylab="Growth Rate",xlab="Initial Population Size")
for (i in 1:length(all_simpoprate)){
  points(rep(init_pop[i],iterations),all_simpoprate[[i]],col="grey")
}
lines(init_pop,mean,col="darkgreen")
abline(h=true_mean,col="red")
par(new=T)
d=data.frame(x=init_pop,v=log(variance),m=mean)
with(d,plot(x,v,type='p',col="cornflowerblue",axes=F,xlab=NA,ylab=NA,cex=1.2))
loess_fit<-loess(v~x,data=d)
init_pop=init_pop[which(variance>0)]
lines(init_pop,predict(loess_fit),col="darkblue",lwd=1,lty=1)
axis(side=4)
mtext(side=4,line=3,"log(Variance)")
legend("bottom",legend=c("Mean","True Mean","Variance","Pop'n Simulations"),lty=c(1,1,0,0),pch=c(NA,NA,16,16),col=c("forestgreen","red","cornflowerblue","grey"))
title(main=paste("Simulated Population Growth Rate for", datsetname, strain$name))
#dev.off()

####### Aside (One-off Calculations)

# #Calculating the estimate growth rates for all growth curves
# all_rates=c()
# all_int=c()
# all_dist=c()
# for (i in 1:dim(area)[1]){
#   k=LM_growthrate(area[i,],times[i,])$rate
#   intercept=LM_growthrate(as.numeric(area[i,]),as.numeric(times[i,]))$int
#   fit=LM_growthrate(as.numeric(area[i,]),as.numeric(times[i,]))$fit
#   res=studres(fit)
#   all_dist=c(all_dist,range(res)[1]-range(res)[2])
#   all_rates=c(all_rates,k)
#   all_int=c(all_int,intercept)
# }

# #Saving the range of the residuals according to blob number
# res_data=cbind(dist,strain$data$clonalcolony,strain$data$blobnumber)
#
# #Get rid of NA distances; this means not enough time points
# all_rates=all_rates[-which(is.na(all_dist))]
# all_data=data[-which(is.na(all_dist)),]
# all_dist=all_dist[-which(is.na(all_dist))]
# all_rates[which(all_rates<0)]=0
# #Make distance positive
# all_dist=abs(all_dist)
# #Save useful info
# all_res_data=cbind(all_dist,all_rates,data$clonalcolony,data$blobnumber)
# write.table(all_res_data,"Lawless_ResRateRange_Filtered_Norm_20.txt",col.names=FALSE,row.names=FALSE)
#all_rates[which(all_rates<0)]=0

# # Calibration Curve
# cells=breaks=seq(0,350,10)
# col_range=seq(min(strain_rates),max(strain_rates),(max(strain_rates)-min(strain_rates))/5)
# colours=rainbow(5)
# strain$area=as.matrix(strain$area)
# hist(strain$area[,1],main=paste("Starting Cell Size for",strain$name),xlab="Area at t=0",breaks=cells)
# hist(strain$area[which(strain_rates==0),1],col="grey",add=T,breaks=cells)
# hist(strain$area[,1],main=paste("Starting Cell Size for",strain$name),xlab="Area at t=0",breaks=cells)
# hist(strain$area[which(strain_rates==0),1],col="grey",add=T,breaks=cells)
# for (i in 1:5){
#   hist(strain$area[which(strain_rates>col_range[i] & strain_rates<=col_range[i+1]),1],col=adjustcolor(colours[i],0.4),breaks=cells,add=T)
# }
# legend("topright",legend=round(col_range,2),pch=rep(0,length(col_range)),col=c("grey",colours))
#
# mean(strain$area[which(strain_rates<=col_range[1]),1])
# median(strain$area[which(strain_rates<=col_range[1]),1])
#
# area=as.matrix(area)
# mean(area[which(all_rates==0),1])
# median(area[which(all_rates==0),1]) #Say this is equal to one cell
#
# area_cell=median(area[which(all_rates==0),1])
# calibrated_area=t(apply(area,1, function(x) x/area_cell))
#
# #Getting the blob number of the three fastest rates
# strain_rates_sorted=sort(strain_rates,decreasing=TRUE)
# max_rates=strain_rates_sorted[1:3]
# max_indices=which(strain_rates==max_rates[1]| strain_rates==max_rates[2]|strain_rates==max_rates[3])
# print(strain$data[max_indices,])

# #Overlay the 2 histograms to show that non-dividing cell start size distribution matches that of the entire pop'n
# hist(exp(strain_int),main=paste("Simulated Starting Size for",strain$name),xlab="Area (px)",breaks=cells)
# hist(exp(strain_int[which(strain_rates==0)]),col="grey",add=T,breaks=cells)
# hist(exp(strain_int),main=paste("Simulated Starting Size for",strain$name),xlab="Area (px)",breaks=cells)
# hist(exp(strain_int[which(strain_rates==0)]),col="grey",add=T,breaks=cells)
# for (i in 1:5){
#   hist(exp(strain_int[which(strain_rates>col_range[i] & strain_rates<=col_range[i+1])]),col=adjustcolor(colours[i],0.4),breaks=cells,add=T)
# }
# legend("topright",legend=round(col_range,2),pch=rep(0,length(col_range)),col=c("grey",colours))

# #Fish Plots (Could write this as a function vis_het_dat())
# sample=8
# indices1=sample(dim(strain$area)[1],sample)
# indices2=sample(max_indices,1)
# indices3=sample(which(strain_rates==0),1)
# indices=c(indices1,indices2,indices3)
# time_indices=seq(1,dim(strain$area)[2],1)
# timepoints=times[1,time_indices]
# area_m=as.matrix(strain$area)
# area_test=area_m[indices,time_indices]
# pop_data=colSums(area_test)
# frac.table1=matrix(0,nrow=dim(area_test)[1],ncol=length(timepoints))
# for (i in 1:dim(area_test)[1]){
#   frac.table1[i,]=as.numeric(area_test[i,])/as.numeric(pop_data)
# }
# frac.table1=frac.table1*100
# parents=rep(0,dim(area_test)[1])
# fish = createFishObject(frac.table1,parents,timepoints=timepoints,col=rainbow(sample+2))
# fish = layoutClones(fish)
# fishPlot(fish,shape="spline",title.btm="Sample1",
#          cex.title=0.5, vlines=c(0,150),
#          vlab=c("day 0","day 150"))
# title(paste("FishPlot for",datsetname,strain$name))
