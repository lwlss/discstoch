# Single lineages data for population simulations
# Using frequentist deterministic modelling
# In accordance with FreqDetMod.R

#====================== Parsing the Data ======================================

# Subsetting the data according to strain (genotype)
#'Subsetting the formatted data according to strain.
#'
#'@param d - Data
#'@param a - Area
#'@param t - Times
#'@param strain - Strain Name
#'@return A list of all strain information: its \code{area, times, data, name} and \code{indices} which can
#'accessed using the \code{$} operator.
#'@export
subset_strain<-function(d,a,t,strain){
  s_name=toString(strain)
  indices=which(d$genotype == strain)
  s_area=a[indices,]; s_times=t[indices,]; s_data=d[indices,]
  return(list("area"=s_area,"times"=s_times,"data"=s_data,
              "name"=s_name, "indices"=indices))
}

# Simulating multiple population growth curves
#'@export
pop_sim_dat<-function(strain,k,Nsample,t,it){
  #where s is a single strain name, N is the sample size,
  #and it is the number of iterations, t is the time
  PopData=matrix(0,nrow=it,ncol=length(t))
  total_indices=c()
  total_simdata=list()
  total_rates=c()
  for (i in 1:it){
    #cannot take sample larger than the population when replace=F
    if (dim(strain$area)[1] < Nsample) {Nsample=dim(strain$area)[1]}
    indices=sample(1:dim(strain$area)[1],Nsample,replace=FALSE)
    rates=k[indices]
    expdata=strain$area[indices,]
    simdata=matrix(0,nrow=length(rates),ncol=length(t))
    for (j in 1:length(rates)){
      simdata[j,]=simExponential(rates[j],t)
    }
    total_simdata[[i]]=simdata
    PopData[i,]=colSums(simdata)
    total_indices=c(total_indices,indices)
    total_rates=c(total_rates,rates)
  }
  return(list("PopData"=PopData,"SimData"=total_simdata,"indices"=total_indices,"SimRates"=total_rates))
}

# Simulating multiple population growth curves from the tail distribution
#'@export
pop_sim_dat_tail<-function(strain,k,Nsample,t,it){
  PopData_Tail=matrix(0,nrow=it,ncol=length(t))
  cutoff=mean(k)+(sd(k))
  indices=which(k>=cutoff)
  #cannot take sample larger than the population when replace=F
  if (length(indices) < Nsample) {Nsample=length(indices)}
  for (i in 1:it){
    pickindices=sample(indices,Nsample,replace=FALSE)
    rates=strain$data$rate[pickindices]
    simdata=matrix(0,nrow=length(rates),ncol=length(t))
    for (j in 1:length(rates)){
      simdata[j,]=simExponential(rates[j],t)
    }
    PopData_Tail[i,]=colSums(simdata)
  }
  return(list("PopData_Tail"=PopData_Tail,"indices_tail"=total_indices))
}

# Returns the rates of a piece-wise linear fit
#'@export
piecewise_pop_rate<-function(a,t){
  k1=c()
  k2=c()
  pb=c()
  intval=c()
  for (i in 1:dim(a)[1]){
    rate1=piecewise_LM(log(a[i,]),t)$slope1
    rate2=piecewise_LM(log(a[i,]),t)$slope2
    pbreak=piecewise_LM(log(a[i,]),t)$pb
    intercept=piecewise_LM(log(a[i,]),t)$intercept
    k1=c(k1,rate1)
    k2=c(k2,rate2)
    pb=c(pb,pbreak)
    intval=c(intval,intercept)
  }
  return(list("k1"=k1,"k2"=k2,"pb"=pb, "intval"=intval))
}

# Calculates the contribution of a fasted sampled strain to the pop'n growth
#'@export
fast_rate_contribution<-function(k_p,k_i){
  r_i_max=max(strain$data$rate[popdata_indices])
  contribution=k_p/r_i_max
  return(contribution)
}

#====================== Plotting the Data =====================================

#Getting breaks for a histogram
#'@export
hist_cells<-function(dat,int){
  lo=trunc(min(dat)*10)/10-0.1 #rounding down
  hi=trunc(max(dat)*10)/10+0.1 #rounding up
  cells=seq(lo,hi,int)
  return(cells)
}

# Histogram of Growth Rate
#'@export
histo<-function(k,s,SE=TRUE){
  #cells=hist_cells(k,0.01)
  cells=seq(0,0.6,0.01)
  hist(k,breaks=cells,xlab="r",
       main=paste(s,"; N=",length(k),sep=""),
       cex.lab=1.4, cex.main=1.2,col="lightblue")
  if (SE==TRUE){
    abline(v=mean(k), col="red")
    #Standard Error
    abline(v=(mean(k)+(sd(k)/sqrt(length(k)))),col="red",lty=2)
    abline(v=(mean(k)-(sd(k)/sqrt(length(k)))),col="red",lty=2)
    abline(v=range(k)[1],col="lightblue")
    abline(v=range(k)[2],col="lightblue")
  }
}

# Plotting growth curves for sampled strain on log scale
#'Plotting growth curves on the log scale.'
#'
#'@param a - Area
#'@param t - Time
#'@param s - Identifier
#'@param Nsample - No. of growth curves sampled
#'@param title - Adds \code{"Identifier; Nsample"} as title
#'@param hist - Adds a histogram of the estimated growth rates for all areas supplied using \code{LM_growthrate}
#'@return A time versus area plot of sampled growth curves and (if \code{hist=TRUE}) a histogram of all growth rates.
#'@export
plot_growth<-function(a,t,s,Nsample=100,title=TRUE,hist=FALSE){
  #where area (a), times (t) and name (s) are required as input variables
  #cannot take sample larger than the population when replace=F
  if (dim(a)[1] < Nsample) {Nsample=dim(a)[1]}
  indices=sample(1:dim(a)[1],Nsample,replace=FALSE)
  if (hist==FALSE){
    plot(1,type='n', xlim=range(t[indices,],na.rm=TRUE),
         ylim=range(a[indices,],na.rm=TRUE), log='y',ylab="",xlab="")
    if (title==TRUE){
      title(main=paste(s,"; N=",Nsample,sep=""),xlab="Time (h)",
            ylab="Microcolony Area (px)",cex.lab=1.2)
    }
    for (i in indices){
      lines(as.numeric(t[i,]),as.numeric(a[i,]),col=rgb(0.3,0.3,0.3,0.3),lwd=2)
    }
    return(indices)
  }else{
    op=par(mfrow=c(2,1))
    k=c()
    plot(1,type='n', xlim=range(t[indices,],na.rm=TRUE),
         ylim=range(a[indices,],na.rm=TRUE),xlab="",
         ylab="", log='y')
    if (title==TRUE){
      title(main=paste(s,"; N=",Nsample,sep=""),xlab="Time (h)",
            ylab="Microcolony Area (px)",cex.lab=1.2)
    }
    for (i in indices){
      lines(as.numeric(t[i,]),as.numeric(a[i,]),col=rgb(0.3,0.3,0.3,0.3),lwd=2)
      rate=LM_growthrate(as.numeric(a[i,]),as.numeric(t[i,]))$rates
      k=c(k,rate)
    }
    k[which(k<0)]=0
    histo(k,s,SE=TRUE)
    par(op)
    return(indices)
  }
}

# Plotting multiple population growth according to the exponential model
#'Plotting simulated population growth curves on the log scale.'
#'
#'@param a - Area
#'@param t - Time
#'@param s - Identifier
#'@param title - Adds \code{"Identifier; Nsample"} as title
#'@param hist - Adds a histogram of the estimated growth rates for all areas supplied using \code{LM_growthrate}
#'@return A time versus area plot of simulated growth curves and (if \code{hist=TRUE}) a histogram of all growth rates.
#'@export
pop_plot_growth<-function(a,t,s,title=TRUE,hist=FALSE){
  if (hist==FALSE){
    plot(1,type='n', xlim=range(t), ylim=range(a),xlab="",
         ylab="",log='y')
    if (title==TRUE){
      title(main=paste(s, "; N=",a[1,1],"; Iterations=", dim(a)[1], sep=""),xlab="Time (h)",
            ylab="log(Area) (px)",cex.lab=1.2)
    }
    for (i in 1:dim(a)[1]){
      lines(t,a[i,],col=rgb(0.3,0.3,0.3,0.3),lwd=2)
    }
  }else{
    op=par(mfrow=c(2,1))
    plot(1,type='n', xlim=range(t), ylim=range(a),xlab="Time",
         ylab="log(Area)",log='y')
    if (title==TRUE){
      title(main=paste(s, "; N=",a[1,1],"; Iterations=", dim(a)[1], sep=""),xlab="Time (h)",
            ylab="log(Area) (px)",cex.lab=1.2)
    }
    k=c()
    for (i in 1:dim(a)[1]){
      lines(t,a[i,],col=rgb(0.3,0.3,0.3,0.3),lwd=2)
      rate=LM_growthrate(a[i,],t)$rates
      k=c(k,rate)
    }
    histo(k,s,SE=TRUE)
    par(op)
    return(k)
  }
}

# Plotting the histogram of the strain growth rate wih the simulated population
# growth rate
#'@export
strain_pop_hist<-function(s,k_s,k_pop1,k_pop2=0){
  histo(k_s,s,SE=FALSE)
  abline(v=mean(k_pop1), col="red")
  abline(v=(mean(k_pop1)+(sd(k_pop1)/sqrt(length(k_pop1)))),col="red",lty=2)
  abline(v=(mean(k_pop1)-(sd(k_pop1)/sqrt(length(k_pop1)))),col="red",lty=2)
  if(k_pop2!=0){
    abline(v=mean(k_pop2), col="green")
    abline(v=(mean(k_pop2)+(sd(k_pop2)/sqrt(length(k_pop2)))),col="green",lty=2)
    abline(v=(mean(k_pop2)-(sd(k_pop2)/sqrt(length(k_pop2)))),col="green",lty=2)
  }
}
