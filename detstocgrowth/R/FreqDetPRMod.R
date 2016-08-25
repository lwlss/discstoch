# Single lineages data for population simulations
# Using frequentist deterministic modelling
# In accordance with FreqDetMod.R

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

#'Calculating the duration of the lag phase for various population start sizes.
#'
#'@param strain - Strain data as returned by the \code{subset_strain} function or else \code{strain=0}
#'@param k - Growth Rates
#'@param t - Simulated Time Series
#'@param N - Vector of population start sizes for which to estimate lag duration
#'@param it - Number of Iterations
#'@return Returns a list with the mean lag lengths \code{m}, the upper confidence intervals \code{cu},
#'the lower confidence intervals \code{cl}.
#'@export
lagduration<-function(strain=0,k,t,N,it){
  time=t
  r=k
  iters=it
  meanlaglen=c()
  upperconf=c()
  lowerconf=c()
  varlaglen=c()
  for(j in N){
    laglen=c()
    print(j)
    for (i in 1:iterations){
      popsim=pop_sim_dat(strain,r,N,time,1)
      popdata=popsim$PopData
      op_bcp=bcp::bcp(as.numeric(log(popdata)),time)
      max_prob=max(op_bcp[8]$posterior.prob,na.rm=TRUE)
      breakpoint_location=which(op_bcp[8]$posterior.prob==max_prob)[1]
      if (max_prob>0.5){
        area=as.numeric(log(popdata))
        t=time
        linmod=lm(area~t)
        segmod=try(segmented::segmented(linmod,seg.Z =~t,psi=breakpoint_location),silent=TRUE)
        if (typeof(segmod)[1]=="list"){
          laglen=c(laglen,segmod$psi[2])
        }
      }
    }
    meanlaglen=c(meanlaglen,mean(laglen))
    test=t.test(laglen)
    lowerconf=c(lowerconf,test$conf.int[1])
    upperconf=c(upperconf,test$conf.int[2])
    varlaglen=c(varlaglen,var(laglen))
  }
  return(list("m"=meanlaglen,"lc"=lowerconf,"uc"=upperconf))
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

# Histogram of Growth Rate (/ Intercepts)
#'@export
histo<-function(k,s,c,SE=TRUE,int=FALSE){
  if (sum(c)!=0){
    cells=c #  cells=seq(0,0.6,0.01)
  }else{cells=hist_cells(k,0.01)}
  if(int==FALSE){
    hist(k,breaks=cells,xlab="r (1/h)",
         main=paste(s,"; N=",length(k),sep=""),
         cex.lab=1.4, cex.main=1.2,col="lightblue")
  }else{
    hist(k,breaks=cells,xlab="Intercept (log(Area))",
         main=paste(s,"; N=",length(k),sep=""),
         cex.lab=1.4, cex.main=1.2,col="lightblue")
  }
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
    histo(k,s,c=0,SE=TRUE)
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
    histo(k,s,c=0,SE=TRUE)
    par(op)
    return(k)
  }
}

# Plotting the histogram of the strain growth rate wih the simulated population
# growth rate
#'@export
strain_pop_hist<-function(s,k_s,k_pop1,k_pop2=0){
  histo(k_s,s,c=0,SE=FALSE)
  abline(v=mean(k_pop1), col="red")
  abline(v=(mean(k_pop1)+(sd(k_pop1)/sqrt(length(k_pop1)))),col="red",lty=2)
  abline(v=(mean(k_pop1)-(sd(k_pop1)/sqrt(length(k_pop1)))),col="red",lty=2)
  if(k_pop2!=0){
    abline(v=mean(k_pop2), col="green")
    abline(v=(mean(k_pop2)+(sd(k_pop2)/sqrt(length(k_pop2)))),col="green",lty=2)
    abline(v=(mean(k_pop2)-(sd(k_pop2)/sqrt(length(k_pop2)))),col="green",lty=2)
  }
}

#'Simulating multiple population growth curves.
#'
#'@param strain - Strain data as returned by the \code{subset_strain} function or else \code{strain=0}
#'@param k - Growth Rates
#'@param t - Simulated Time Series
#'@param it - Number of Iterations
#'@param yl - Limits in y=axis
#'@param col - Vector of colours
#'@param colrange - Vector of k values which are to match the colour of vectors provide
#'@return Plots the population data as simulated by \code{pop_sim_dat}.
#'@export
plot_pop_sim_dat=function(strain=0,k,t,it,yl,col,colrange){
  plot(1,type='n', xlim=range(time), ylim=yl,xlab="Time (h)",
       ylab="No. of Cells",log='y',cex.lab=1.5)
  colours=col
  gr_range=colrange
  r=k
  time=t
  iterations=it
  lagr=c()
  expr=c()
  for (i in 1:iterations){
    popdata=pop_sim_dat(strain,r,N,time,1)$PopData
    op_bcp=bcp::bcp(log(as.numeric(popdata)),time)
    max_prob=max(op_bcp$posterior.prob,na.rm=TRUE)
    breakpoint_location=which(op_bcp$posterior.prob==max_prob)[1]
    if(max_prob>0.5){
      a=as.numeric(log(popdata))
      t=time
      linmod=lm(a~t)
      segmod=try(segmented::segmented(linmod,seg.Z =~t,psi=breakpoint_location),silent=TRUE)
      if (typeof(segmod)[1]=="list"){
        abline(v=segmod$psi[2],col=adjustcolor("black",0.1),lty=1)
        lagr=c(lagr,slope(segmod)$t[1])
        expr=c(expr,slope(segmod)$t[2])
        k=slope(segmod)$t[2]
      }else{k=detstocgrowth::LM_growthrate(as.numeric(popdata),time)$rate}
    }else{k=detstocgrowth::LM_growthrate(as.numeric(popdata),time)$rate}
    id=which(abs(gr_range-k)==min(abs(gr_range-k)))
    lines(time,popdata,col=adjustcolor(colours[id],0.2),lwd=2)
  }
  tp=t.test(lagr,expr,conf.level=0.99,paired=TRUE)
  print(tp$p.value)
  legend("bottomright",legend=c(paste("Rate 1 =",signif(mean(lagr),2)),paste("Rate 2 =",signif(mean(expr),2))),bty='n')
  return(list("pairedt"=tp,"lagr"=lagr,"expr"=expr))
}

#'Pin Population Estimates
#'
#'@param sa - Pin Area Estimates
#'@param st - Pin Time Points
#'@param strainname - Plot title
#'@return Plots the pin population data with a line of best fit assuming a lag, esponential and stationary phase.
#'The two break points and the three consecuive slope estimates of the fit are returned.
#'@export
pin_pop_plot<-function(strainname,sa,st){
  plot(NULL,xlim=c(0,28),ylim=c(10,15),xlab="Time (h)", ylab="log(Area)",
       main=paste("Population Estimates for",strainname),cex.lab=1.4)
  bp1=c()
  bp2=c()
  s1=c()
  s2=c()
  s3=c()
  for (i in 1:dim(sa)[1]){
    a=log(as.numeric(sa[i,]))
    t=as.numeric(st[i,])
    linmod=lm(a~t)
    segmod=segmented::segmented(linmod,seg.Z =~t,psi=c(15,25))
    points(t,a,pch=16,col=adjustcolor("red",0.2))
    plot(segmod,add=T)
    bp1=c(bp1,segmod$psi[1,2])
    bp2=c(bp2,segmod$psi[2,2])
    s1=c(s1,segmented::slope(segmod)$t[1,1])
    s2=c(s2,segmented::slope(segmod)$t[2,1])
    s3=c(s3,segmented::slope(segmod)$t[3,1])
  }
  return(list("bp2"=bp1,"bp2"=bp2,"s1"=s1,"s2"=s2,"s3"=s3))
}

#'Comparin Pin Population Estimates with Pin Population Simulations.
#'
#'@param area - Strain Pin Area Estimates
#'@param times - Strain Pin Time Points
#'@param sortedfolders - Pin Folder Names
#'@param obs - Number of Observations to Plot and Simulate
#'@param conv - Cell to Area Conversion Factor
#'@param iterations - Number of Population Growth Curves to simulate
#'@param leg - if TRUE legend is added to the plot
#'@return Plots the pin population data with a line of best fit assuming a lag and the simulated population data overlayed.
#'@export
PlotPinObsSim<-function(area,times,sortedfolders,obs,conv,iterations,leg=FALSE){
  for (f in sortedfolders){
    print(f)
    plot(NULL,xlim=c(0,obs),ylim=c(10,20),xlab="Time (h)", ylab="No. of Cells",axes=F,
         main=paste("Observed and Simulated Population \nGrowth in Pin",f),cex.lab=1.4)
    sa=area[which(folders==f),]
    st=times[which(folders==f),]
    a=log(as.numeric(sa))
    t=as.numeric(st)
    linmod=lm(a~t)
    segmod=segmented::segmented(linmod,seg.Z =~t,psi=c(15,25))
    points(t,a,pch=16,col=adjustcolor("red",0.7),lwd=2)
    plot(segmod,add=T)
    N=as.numeric(sa[1]/conv)
    time=seq(0,obs,1)
    iterations=it
    folder=list()
    folder$area=lineagearea[which(data$clonalcolony==f),]
    folder$times=lineagetime[which(data$clonalcolony==f),]
    folder_rates=rates[which(data$clonalcolony==f)]
    lagr=c()
    expr=c()
    for (i in 1:iterations){
      popsim=pop_sim_dat(strain=0,folder_rates,N,time,1)
      popdata=popsim$PopData
      poprate=popsim$SimRates
      lines(time,log(popdata*conv),col=adjustcolor("cyan2",0.1))
      op_bcp=bcp::bcp(as.numeric(log(popdata)),time)
      max_prob=max(op_bcp[8]$posterior.prob,na.rm=TRUE)
      breakpoint_location=which(op_bcp[8]$posterior.prob==max_prob)[1]
      if (max_prob>0.5){
        m=as.numeric(log(popdata))
        z=time
        linmod_pop=lm(m~z)
        segmod_pop=try(segmented::segmented(linmod_pop,seg.Z =~z,psi=6),silent=FALSE)
        if (typeof(segmod_pop)[1]=="list"){
          abline(v=segmod_pop$psi[2],col=adjustcolor("blue",0.1))
          lagr=c(lagr,segmented::slope(segmod_pop)$z[1])
          expr=c(expr,segmented::slope(segmod_pop)$z[2])
        }
      }
    }
    abline(v=segmod$psi[1,2],col="red",lwd=2)
    box()
    axis(side=2,at=seq(10,20),labels=c(round(exp(seq(10,20))/conv)))
    axis(side=1)
    if(leg==TRUE){
      legend("topleft",legend=c("Observed","Fit","Simulated"),pch=c(16,NA,NA),lty=c(NA,1,1),col=c("red","black","cyan2"),lwd=c(NA,2,3),bty='n')
      legend("right",legend=c("Obs. Lag", "Sim. Lag"),col=c("red","blue"),lty=c(1,1),lwd=c(2,2),bty='n')
      legend("bottomright",legend=c(paste("Rate 1 =",signif(mean(lagr),2)),paste("Rate 2 =",signif(mean(expr),2))),bty='n')
    } else{}
    tp=t.test(lagr,expr,conf.level=0.99,paired=TRUE)
    return(list("tp"=tp,"lagr"=lagr,"expr"=expr))
  }
}
