# Analysis Script Using the Package
# Frequentist, Deterministic, Exponential Modelling
setwd("~/Documents/MSc/discstoch-master/FormattedData/FinalizedScripts")

library(data.table)
library(detstocgrowth)
library(fishplot)

dataset<-function(x){
  if (x == "Lawless"){
    # DataSet1: Lawless
    area=fread("Lawless_area_unfiltered.txt",header=FALSE)
    times=fread("Lawless_time_unfiltered.txt",header=FALSE)
    data=fread("Lawless_data_unfiltered.txt",header=FALSE) #3rd column (Identifier) => strain_parentcolony
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
pickstrain=strain_names[2]#choose strain here!
strain=subset_strain(data,area,times,pickstrain)
#plot_growth(strain$area,strain$times,strain$name,Nsample=dim(strain$area)[1],hist=TRUE)

#Exclude time points of all zeros from the analyses
zero_indices=which(rowSums(strain$area)==0)
strain$area=strain$area[-zero_indices,]
strain$times=strain$times[-zero_indices,]
strain$data=strain$data[-zero_indices,]

#Calculating the estimated growth rates for all growth curves of the strain
strain_rates=c()
strain_int=c()
dist=c()
for (i in 1:dim(strain$area)[1]){
  k=LM_growthrate(strain$area[i,],strain$times[i,])$rate
  intercept=LM_growthrate(as.numeric(strain$area[i,]),as.numeric(strain$times[i,]))$int
  fit=LM_growthrate(as.numeric(strain$area[i,]),as.numeric(strain$times[i,]))$fit
  res=residuals(fit)
  dist=c(dist,range(res)[1]-range(res)[2])
  strain_rates=c(strain_rates,k)
  strain_int=c(strain_int,intercept)
}

#Setting growth rates <0 equal to 0
strain_rates[which(strain_rates<0)]=0

#Getting the blob number of the three fastest rates
strain_rates_sorted=sort(strain_rates,decreasing=TRUE)
max_rates=strain_rates_sorted[1:3]
max_indices=which(strain_rates==max_rates[1]| strain_rates==max_rates[2]|strain_rates==max_rates[3])
print(strain$data[max_indices,])

#Fish Plots (Could write this as a function vis_het_dat())
sample=8
indices1=sample(dim(strain$area)[1],sample)
indices2=sample(max_indices,1)
indices3=sample(which(strain_rates==0),1)
indices=c(indices1,indices2,indices3)
time_indices=seq(1,dim(strain$area)[2],1)
timepoints=times[1,time_indices]
area_m=as.matrix(strain$area)
area_test=area_m[indices,time_indices]
pop_data=colSums(area_test)
frac.table=matrix(0,nrow=dim(area_test)[1],ncol=length(timepoints))
for (i in 1:dim(area_test)[1]){
  frac.table[i,]=as.numeric(area_test[i,])/as.numeric(pop_data)
}
frac.table=frac.table*100
parents=rep(0,dim(area_test)[1])
fish = createFishObject(frac.table,parents,timepoints=timepoints,col=rainbow(sample+2))
fish = layoutClones(fish)
fishPlot(fish,shape="spline",title.btm="Sample1",
         cex.title=0.5, vlines=c(0,150),
         vlab=c("day 0","day 150"))
title(paste("FishPlot for",datsetname,strain$name))

# Using 10% of the total data for extrapolation
total=dim(strain$area)[1]
N=round(total*0.1)
time=seq(0,35,0.5)

# Refitting the exponential model to the newly simulated data to obtain
# population growth rate
iterations=20
popsim=pop_sim_dat(strain,strain_rates,N,time,iterations)
popdata=popsim$PopData

pop_rates=c()
for (i in 1:dim(popdata)[1]){
  k=LM_growthrate(as.numeric(popdata[i,]),time)$rate
  pop_rates=c(pop_rates,k)
}

op=par(mfrow=c(2,1))
pop_plot_growth(popdata,time,strain$name,hist=FALSE)
histo(strain_rates,strain$name,SE=FALSE)
abline(v=mean(strain_rates),col="red",lwd=3)
abline(v=mean(pop_rates),col="blue",lwd=3)
legend("topleft",legend=c("True Mean","Pop'n Mean"),col=c("red","blue"),lty=c(1,1))
par(op)

#Note indices for all iterations!
popdata_indices=popsim$indices

# Fitting a piece-wise regression to the Population Simulations
piece_op=piecewise_pop_rate(popdata,time)

# Analysing how much each growth rate contributes to the pop'n growth during
# the exponential phase
# contribution=fast_rate_contribution(piece_op$k2,strain_rates[popdata_indices])
# filename=paste("Contribution_",datsetname,strain$name,"_0.1N_c",
#                signif(contribution, digits=3),".png",sep="")
# png(filename)
op=par(mfrow=c(3,1))
plot_growth(strain$area[popdata_indices,],strain$times[popdata_indices,],
            strain$name,Nsample=100)
pop_plot_growth(popdata,time,strain$name,hist=FALSE)
strain_pop_hist(strain$name,strain_rates[popdata_indices],piece_op$k1,piece_op$k2)
par(op)
# dev.off()

# #Overlay the 2 histograms to show that non-dividing cell start size distribution matches that of the entire pop'n
# hist(exp(strain_int),main="Starting Size for HIS3",xlab="Area (px)",breaks=seq(0,300,25))
# hist(exp(strain_int[which(strain_rates>0 & strain_rates<0.05)]),breaks=seq(0,300,25),col="lightblue",add=T)
# #Not really a convincing argument at the moment....


# Changes in variance & value of growth rate estimates with an inceasing
# size of the starting population
N=dim(strain$area)[1]
iterations=50
total_popsim=pop_sim_dat(strain,strain_rates,N,time,iterations)
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
    simpoprate=LM_growthrate(simitdata,time)$rates #Pop'n Rates
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
png(paste(paste("PopSimDat",datsetname,strain$name,signif(true_mean,3),signif(true_variance,3),"IT1",sep="_"),".png",sep=""))
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
dev.off()
