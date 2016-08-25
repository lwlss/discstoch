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

# Subsetting the data according to the identifier
#'Subsetting the formatted data according to strain.
#'
#'@param d - Data
#'@param a - Area
#'@param t - Times
#'@param identifier - Identifier Name
#'@return A list of all identifier information: its \code{area, times, data, name} and \code{indices} which can
#'accessed using the \code{$} operator.
#'@export
subset_identifier<-function(d,a,t,identifier){
  i_name=toString(identifier)
  indices=which(d$identifier == identifier)
  i_area=a[indices,]; i_times=t[indices,]; i_data=d[indices,]
  return(list("area"=i_area,"times"=i_times,"data"=i_data,
              "name"=i_name, "indices"=indices))
}

#'Subsetting the formatted data according to pin.
#'
#'@param d - Data
#'@param a - Area
#'@param t - Times
#'@param pin - Pin ID
#'@return A list of all pin information: its \code{area, times, data, name} and \code{indices} which can
#'accessed using the \code{$} operator.
#'@export
subset_pin<-function(d,a,t,pin){
  c_name=toString(pin)
  indices=which(d$clonalcolony == pin)
  c_area=a[indices,]; c_times=t[indices,]; c_data=d[indices,]
  return(list("area"=c_area,"times"=c_times,"data"=c_data,
              "name"=c_name, "indices"=indices))
}

#'Simulating multiple population growth curves.
#'
#'@param strain - Strain data as returned by the \code{subset_strain} function or else \code{strain=0}
#'@param k - Growth Rates
#'@param Nsamples - Simulated population Start Size
#'@param t - Simulated Time Series
#'@param it - Number of Iterations
#'@return A list of the simulated population data \code{PopData}, all simulated strain data
#'\code{total_simdata}, all sampled lineage growth rate \code{SimRates}, and their corresponding
#'growth curve indices \code{total_indices}.
#'@export
pop_sim_dat<-function(strain=0,k,Nsample,t,it){
  #where s is a single strain name, N is the sample size,
  #and it is the number of iterations, t is the time
  if (strain==0){
    PopData=matrix(0,nrow=it,ncol=length(t))
    total_indices=c()
    total_simdata=list()
    total_rates=c()
    for (i in 1:it){
      indices=sample(1:length(k),N,replace=TRUE)
      rates=k[indices]
      simdata=matrix(0,nrow=length(rates),ncol=length(t))
      for (j in 1:length(rates)){
        simdata[j,]=simExponential(rates[j],t)
      }
      total_simdata[[i]]=simdata
      PopData[i,]=colSums(simdata)
      total_indices=c(total_indices,indices)
      total_rates=c(total_rates,rates)
    }
    return(list("PopData"=PopData,"SimData"=total_simdata,
                "indices"=total_indices,"SimRates"=total_rates))
  }
  else{
    PopData=matrix(0,nrow=it,ncol=length(t))
    total_indices=c()
    total_simdata=list()
    total_rates=c()
    for (i in 1:it){
      #cannot take sample larger than the population when replace=F
      #if (dim(strain$area)[1] < Nsample) {Nsample=dim(strain$area)[1]}
      indices=sample(1:dim(strain$area)[1],Nsample,replace=TRUE)
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
    return(list("PopData"=PopData,"SimData"=total_simdata,
                "indices"=total_indices,"SimRates"=total_rates))
  }
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

#' Yellow to Red Map for Growth Rate
#'
#'@param gr - Range of Growth Rates
#'@return A yellow to red colour map which provides a measure for growth rate
#'@export
yellowredmap<-function(gr){
  # adapted from http://stackoverflow.com/questions/9314658/colorbar-from-custom-colorramppalette
  y=gr
  cols<-colorRampPalette(c("lightyellow","yellow","orange","red"))(n=length(gr))
  z=matrix(1:length(cols),nrow=1)
  x=1
  image(x,y,z,col=cols,axes=FALSE,xlab="",ylab="Growth Rate (1/h)",cex.lab=1.3)
  axis(2)
  box()
}

#' FishPlots for Population Simulations
#'
#'@param r - Sampled Growth Rates
#'@param t - Time sequence over which to simulate
#'@param gr - Range of Growth Rates
#'@return A FishPlot for the fastest growing strains
#'@export
pop_comoposition<-function(r,t,gr){
  sim_area=matrix(0,ncol=length(t),nrow=length(r))
  for (i in 1:length(r)){
    sim_area[i,]=round(detstocgrowth::simExponential(r[i],t))
  }
  pop_data=colSums(sim_area)
  frac.table=matrix(0,nrow=dim(sim_area)[1],ncol=length(t))
  for (i in 1:dim(sim_area)[1]){
    frac.table[i,]=as.numeric(sim_area[i,])/as.numeric(pop_data)
  }
  frac.table=frac.table*100
  frac.table=floor(frac.table*100)/100
  for (i in 1:dim(frac.table)[1]){
    if (sum(which(frac.table[i,]==0))>0){
      frac.table[i,(which(frac.table[i,]==0)[1]):(dim(frac.table)[2])]=0
    }
  }
  if(length(r)==10000){
    nogrowth=apply(frac.table,1,function(x) (sum(x)==0.01))
    frac.table=frac.table[-(which(nogrowth==TRUE)),]
    r=r[-(which(nogrowth==TRUE))]
  }
  parents=rep(0,dim(frac.table)[1])
  gr_range=gr
  cols<-colorRampPalette(c("lightyellow","yellow","orange","red"))(n = length(gr))
  ids=sapply(r,function(x) which(abs(gr_range-x)==min(abs(gr_range-x))))
  col_list=cols[ids]
  fish = fishplot::createFishObject(frac.table,parents,timepoints=t,col=col_list)
  fish = fishplot::layoutClones(fish)
  fishplot::fishPlot(fish,shape="spline",title.btm=" ",vlines=c(5,max(t)),
                     vlab=c("1h",paste(max(t),"h",sep="")),border=0.01,cex.vlab=1.3)
}

#' Calculates Population Compositions as fed into the Fishplots using \code{pop_comoposition}
#'
#'@param t - Time sequence over which to simulate
#'@param iter - Number of Iterations
#'@param N0 - Population Start Size
#'@param strain_rates - List of Growth Rates
#'@return The number of lineages which make up more than five percent and the percentage with which
#'the fastest lineage dominates the population are returned in a data frame.
#'@export
pop_composition_dat<-function(t,iter,N0,strain_rates){
  ds=matrix(0,nrow=iter,ncol=length(t))
  fs=matrix(0,nrow=iter,ncol=length(t))
  for(j in 1:iter){
    sample_rates=sample(strain_rates,N0,replace=TRUE)
    fid=c(which(sample_rates==max(sample_rates)))
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
    ds[j,]=apply(frac.table,2,function(x) sum(x>5))
    if(length(fid)>1){
      fs[j,]=colSums(frac.table[fid,])
    }else{fs[j,]=frac.table[fid,]}
  }
  p5=apply(ds,2,mean)
  dr=apply(fs,2,mean)
  name1=paste("5P_",strain$name,"_",N0,sep="")
  name2=paste("DR_",strain$name,"_",N0,sep="")
  newinfo=data.frame(p5,dr)
  names(newinfo)=c(name1,name2)
  return(newinfo)
}

#' Simulated Population Growth Curves
#'
#'@param strain - Strain Data from which to simulate
#'@param strain_rates - Strain rates from which to simulate
#'@param gr - Growth Rate range applied to colour-scale
#'@param yl - Maxmimum limits on the y axis
#'@param it - Number of iterations
#'@param t - Time sequence over which to simulate
#'@return The number of lineages which make up more than five percent and the percentage with which
#'the fastest lineage dominates the population are returned in a data frame.
#'@export
pop_sim_plot<-function(strain,strain_rates,gr,yl,it,t){
  plot(1,type='n', xlim=range(t), ylim=c(100,yl),xlab="Time (h)",ylab="No. of Cells",log='y',cex.lab=1.4)
  colours<-colorRampPalette(c("lightyellow","yellow","orange","red"))(n = length(gr))
  gr_range=gr
  bp_pr=c()
  bp_loc=c()
  pop_rates=c()
  pop_rates_lag=c()
  for (i in 1:it){
    popsim=pop_sim_dat(strain=0,strain_rates,N,t,1)
    popdata=popsim$PopData
    poprate=popsim$SimRates
    id=which(abs(gr_range-max(poprate))==min(abs(gr_range-max(poprate))))
    lines(t,popdata,col=colours[id],lwd=2)
    op_bcp=bcp::bcp(as.numeric(log(popdata)),t)
    max_prob=max(op_bcp[8]$posterior.prob,na.rm=TRUE)
    breakpoint_location=which(op_bcp[8]$posterior.prob==max_prob)[1]
    bp_pr=c(bp_pr,max_prob)
    bp_loc=c(bp_loc,breakpoint_location)
    if (max_prob==1){
      area=as.numeric(log(popdata))
      t=t
      linmod=lm(area~t)
      segmod=try(segmented::segmented(linmod,seg.Z =~t,psi=breakpoint_location),silent=TRUE)
      if (typeof(segmod)[1]=="list"){
        k=segmented::slope(segmod)$t[,1][2]
        klag=segmented::slope(segmod)$t[,1][1]
        abline(v=segmod$psi[2],col=adjustcolor("black",0.2),lty=1)
      }
      else{
        k=detstocgrowth::LM_growthrate(as.numeric(popdata),t)$rate
        klag=NA
      }
    }else{
      k=detstocgrowth::LM_growthrate(as.numeric(popdata),t)$rate
      klag=NA
    }
    pop_rates=c(pop_rates,k)
    pop_rates_lag=c(pop_rates_lag,klag)
  }
  return(list("bp_pr"=bp_pr,"bp_loc"=bp_loc,"pop_rates"=pop_rates,"pop_rates_lag"=pop_rates_lag))
}

#' Population Growth Rates with Increasing Start Sizes
#'
#'@param strain - Strain Data from which to simulate
#'@param strain_rates - Strain rates from which to simulate
#'@param N - Largest population size
#'@param iterations - Number of iterations
#'@param t - Time sequence over which to simulate
#'@return Returns two lists with the lag and exponential growth rate .
#'@export
meanvarsim<-function(strain,strain_rates,N,iterations,t){
  init_pop=seq(1,N,1)
  all_simpoprate=list()
  all_simpopratelag=list()
  for (i in 1:length(init_pop)){
    init_simpoprate=c()
    init_simpopratelag=c()
    for (j in 1:iterations){
      total_popsim=detstocgrowth::pop_sim_dat(strain=0,strain_rates,i,t,1)
      simdata=total_popsim$PopData
      simitdata=simdata
      op_bcp=bcp::bcp(as.numeric(log(simitdata)),t)
      max_prob=max(op_bcp[8]$posterior.prob,na.rm=TRUE)
      breakpoint_location=which(op_bcp[8]$posterior.prob==max_prob)[1]
      if (max_prob>0.5){
        area=as.numeric(log(simitdata))
        t=t
        linmod=lm(area~t)
        segmod=try(segmented::segmented(linmod,seg.Z =~t,psi=breakpoint_location),silent=TRUE)
        if (typeof(segmod)[1]=="list"){
          k=segmented::slope(segmod)$t[,1][2]
          klag=segmented::slope(segmod)$t[,1][1]
        }
        else{
          k=detstocgrowth::LM_growthrate(as.numeric(simitdata),t)$rate
          klag=NA
        }
      }else{
        k=detstocgrowth::LM_growthrate(as.numeric(simitdata),t)$rate
        klag=NA
      }
      init_simpoprate=c(init_simpoprate,k)
      init_simpopratelag=c(init_simpopratelag,klag)
    }
    all_simpoprate[[i]]=as.numeric(init_simpoprate)
    all_simpopratelag[[i]]=as.numeric(init_simpopratelag)
  }
  return(list("all_simpoprate"=all_simpoprate,"all_simpopratelag"=all_simpopratelag))
}
