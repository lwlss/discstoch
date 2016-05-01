#Deterministic Modelling 
#Fitting and Analysing the exponential model 
#Log-linear Regression Analysis 

library(data.table)
library(segmented)
library(formattable)

# Choosing a data set to work with
dataset<-function(x){
  if (x == "Lawless"){
    # DataSet1: Lawless
    area=fread("Lawless_area.txt",header=FALSE)
    times=fread("Lawless_time.txt",header=FALSE)
    data=fread("Lawless_data.txt",header=TRUE) #3rd column (Identifier) => strain_parentcolony 
    residuals=fread("Lawless_residuals.txt",header=FALSE)
    return(list("area"=area,"data"=data,"times"=times,"residuals"=residuals))
  }
  else if (x == "Levy"){
    # DataSet2: Levy
    area=fread("Levy_area.txt",header=FALSE)
    times=fread("Levy_time.txt",header=FALSE)
    data=fread("Levy_data.txt",header=TRUE) #3rd column (Identifier) => replicate
    residuals=fread("Levy_residuals.txt",header=FALSE)
    return(list("area"=area,"data"=data,"times"=times,"residuals"=residuals))
  }
  else if (x == "Ziv"){
    # DataSet3: Ziv
    area=fread("Ziv_area.txt",header=FALSE)
    times=fread("Ziv_time.txt",header=FALSE)
    data=fread("Ziv_data.txt",header=TRUE) #3rd column (Identifier) => colony
    residuals=fread("Ziv_residuals.txt",header=FALSE)
    return(list("area"=area,"data"=data,"times"=times,"residuals"=residuals))
  }
  else {print("Not a valid dataset")}
}

# Subsetting the data according to strain (genotype)
subset_strain<-function(d,a,t,strain){ 
  s_name=toString(strain)
  indices=which(d$genotype == strain)
  s_area=a[indices,]; s_times=t[indices,]; s_data=d[indices,]
  return(list("area"=s_area,"times"=s_times,"data"=s_data, "name"=s_name, "indices"=indices))
}
subset_strain_normal<-function(d,a,t,strain,datset){
  s_name=toString(strain)
  import=paste(datset,"_residuals_normtests.txt",sep="")
  residuals=fread(import,header=TRUE)
  indices=which(d$genotype==strain)
  indices=which(residuals$Lilliefors.H0[indices]==TRUE) #only pick data for normally distributed residuals 
  s_area=a[indices,]; s_times=t[indices,]; s_data=d[indices,]
  return(list("area"=s_area,"times"=s_times,"data"=s_data, "name"=s_name, "indices"=indices))
}

#plotting growth cruves for sampled strain growth curves on log scale
plot_growth<-function(a,t,s,Nsample=100){ 
  #where area (a), times (t) and name (s) are required as input variables
  if (dim(a)[1] < Nsample) {Nsample=dim(a)[1]} #cannot take sample larger than the population when replace=F
  indices=sample(1:dim(a)[1],Nsample,replace=FALSE)
  k=c()
  plot(1,type='n', xlim=range(t[indices,],na.rm=TRUE), ylim=range(a[indices,],na.rm=TRUE),xlab="Time (h)", 
       ylab="Microcolony Area (px)",main=paste(s,"; N=",Nsample,sep=""),cex.lab=1.2, log='y')
  for (i in indices){
    lines(as.numeric(t[i,]),as.numeric(a[i,]),col=rgb(0.3,0.3,0.3,0.3),lwd=2)
  }
  return(list("rates"=k,"ind"=indices))
}

# Simulating exponential growth  
simExponential<-function(r,t,N0=1){
  N=N0*exp(r*t)
  return(N)
}

# Simulating multiple population growth curves 
pop_sim_dat<-function(strain,Nsample,t,it){
  #where s is a single strain name, N is the sample size, 
  #and it is the number of iterations, t is the time 
  PopData=matrix(0,nrow=it,ncol=length(t))
  total_indices=c()
  for (i in 1:it){
    if (dim(strain$area)[1] < Nsample) {Nsample=dim(strain$area)[1]} #cannot take sample larger than the population when replace=F
    indices=sample(1:dim(strain$area)[1],Nsample,replace=FALSE)
    rates=strain$data$rate[indices]
    expdata=strain$area[indices,]
    simdata=matrix(0,nrow=length(rates),ncol=length(t))
    for (k in 1:length(rates)){
      simdata[k,]=simExponential(rates[k],t)
    }
    PopData[i,]=colSums(simdata)
    total_indices=c(total_indices,indices)
  }
  #plot_res_9(dim(strain$area)[2],strain$name,total_indices,residuals)
  return(PopData)
}

# Simulating multiple population growth cruves from the tail distribution 
pop_sim_dat_tail<-function(strain,Nsample,t,it){
  PopData_Tail=matrix(0,nrow=it,ncol=length(t))
  cutoff=mean(strain$data$rate)+(sd(strain$data$rate))
  indices=which(strain$data$rate>=cutoff)
  if (length(indices) < Nsample) {Nsample=length(indices)} #cannot take sample larger than the population when replace=F
  for (i in 1:it){
    pickindices=sample(indices,Nsample,replace=FALSE)
    rates=strain$data$rate[pickindices]
    simdata=matrix(0,nrow=length(rates),ncol=length(t))
    for (k in 1:length(rates)){
      simdata[k,]=simExponential(rates[k],t)
    }
    PopData_Tail[i,]=colSums(simdata)
  }
  #plot_res_9(dim(strain$area)[2],strain$name,total_indices,residuals)
  return(PopData_Tail)
}

# Scatterplot of the color-coded residuals of 80 samples divided over 16 plots 
plot_res_16<-function(d2,name,indices,residuals){
  sample_indices=sample(indices,80,replace=FALSE)
  op=par(mfrow=c(4,4))
  for (i in 1:16){
    plot(NULL,main=paste("Residuals for ", toString(name), " (", toString(i), ")", sep=""),ylim=c(-3,3),xlim=c(1,d2),xlab="Index", ylab="Residual")
    rect(par("usr")[1],par("usr")[3],par("usr")[2],2*max(area,na.rm=TRUE),col = "lightgrey")
    cols=adjustcolor(rainbow(5),0.9)
    indices=sample_indices[seq(1,5)+(5*(i-1))]
    res=residuals[indices,]
    for (i in 1:length(indices)){
      lines(seq(1,d2),res[i,],type="p",col=cols[i]) 
    }
    abline(0,0,lty=2)
  }
  par(op)
} 
plot_res_9<-function(d2,name,indices,residuals){
  sample_indices=sample(indices,45,replace=FALSE)
  op=par(mfrow=c(3,3))
  for (i in 1:9){
    plot(NULL,main=paste("Residuals for ", toString(name), " (", toString(i), ")", sep=""),ylim=c(-1,1),xlim=c(1,d2),xlab="Index", ylab="Residual")
    rect(par("usr")[1],par("usr")[3],par("usr")[2],2*max(area,na.rm=TRUE),col = "lightgrey")
    cols=adjustcolor(rainbow(5),0.9)
    indices=sample_indices[seq(1,5)+(5*(i-1))]
    res=residuals[indices,]
    for (i in 1:length(indices)){
      lines(seq(1,d2),res[i,],type="p",col=cols[i]) 
    }
    abline(0,0,lty=2)
  }
  par(op)
}
plot_res_4<-function(d2,name,indices,residuals){
  sample_indices=sample(indices,20,replace=FALSE)
  op=par(mfrow=c(2,2))
  for (i in 1:4){
    plot(NULL,main=paste("Residuals for ", toString(name), " (", toString(i), ")", sep=""),ylim=c(-3,3),xlim=c(1,d2),xlab="Index", ylab="Residual")
    rect(par("usr")[1],par("usr")[3],par("usr")[2],2*max(area,na.rm=TRUE),col = "lightgrey")
    cols=adjustcolor(rainbow(5),0.9)
    indices=sample_indices[seq(1,5)+(5*(i-1))]
    res=residuals[indices,]
    for (i in 1:length(indices)){
      lines(seq(1,d2),res[i,],type="p",col=cols[i]) 
    }
    abline(0,0,lty=2)
  }
  par(op)
}

#Calculating the growth rate (LS for exponential model)
LM_growthrate<-function(ai,ti){
  fit<-lm(log(ai)~ti) #this is applying the exponential model 
  rate=fit$coefficient[[2]]
  intercept=fit$coefficient[[1]]
  return(list("rates"=rate,"fit"=fit,"int"=intercept))
}

# Plotting multiple population growth and calculating rate constants according to exponential model
pop_plot<-function(a,t,N,res=TRUE){ 
  if(res==TRUE){
    op=par(mfrow=c(3,1))
    plot(1,type='n', xlim=range(t), ylim=range(a),xlab="Time", 
         ylab="log(Area)",main="Population Simulations",cex.lab=1.4,cex.main=1.2,log='y')
    k=c()
    residuals=matrix(NA,ncol=length(time),nrow=dim(a)[1])
    for (i in 1:dim(a)[1]){
      lines(t,a[i,],col=rgb(0.3,0.3,0.3,0.3),lwd=2)
      rate=LM_growthrate(a[i,],t)$rates
      residuals[,]=resid(LM_growthrate(a[i,],t)$fit)
      k=c(k,rate)
    }
    plot(NULL,main="Residuals for Population Simulations",ylim=c(-3,3),xlim=c(1,dim(a)[2]),xlab="Index", ylab="Residual",cex.lab=1.4,cex.main=1.2)
    rect(par("usr")[1],par("usr")[3],par("usr")[2],2*max(area,na.rm=TRUE),col = "lightgrey")
    cols=adjustcolor(rainbow(dim(a)[1]),0.9)
    for (i in 1:dim(a)[1]){
      lines(seq(1,dim(a)[2]),residuals[i,],type="p",col=cols[i]) 
    }
    abline(0,0,lty=2)
    box()
    qqnorm(residuals,col=cols,cex.lab=1.4,cex.main=1.2)
    abline(0,1)
    par(op)
  }
  else{
    plot(1,type='n', xlim=range(t), ylim=range(a),xlab="Time", 
         ylab="log(Area)",main=paste("Population Simulations; N=",a[1,1],"; Iterations=", dim(a)[1], sep=""),cex.lab=1.4,cex.main=1.2,log='y')
    k=c()
    for (i in 1:dim(a)[1]){
      lines(t,a[i,],col=rgb(0.3,0.3,0.3,0.3),lwd=2)
      rate=LM_growthrate(a[i,],t)$rates
      k=c(k,rate)
    }
  }
  
  return(k)
}

piecewise_LM<-function(ai,t){
  dat=data.frame(x=t,y=ai)
  lin_mod=lm(y~x,data=dat)
  #fitNA<-segmented(lin_mod,seg.Z=~x,psi=NA,data=dat) 
  fit<-segmented(lin_mod,seg.Z=~x,psi=10,data=dat,control=seg.control(display=FALSE),model=TRUE)
  residuals=fit$residuals
  pb=fit$psi[2]
  intercept=exp(fit$model$y[which(fit$model$x==round(pb))]) #### converting intercept values from log(area) to area 
  ###### CONTINUE HERE 
  slope1=slope(fit)$x[1]
  slope2=slope(fit)$x[2]
  return(list("slope1"=slope1,"slope2"=slope2,"res"=residuals,"pb"=pb,"fit"=fit, "intercept"=intercept))
}
piecewise_pop_plot<-function(a,t, plot=TRUE){
  if(plot==TRUE){
    op=layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
    plot(1,type='n', xlim=range(t), ylim=range(log(a)),xlab="Time", 
         ylab="log(Area)",main="Population Simulations",cex.lab=1.4,cex.main=1.4)
    k1=c()
    k2=c()
    pb=c()
    intval=c()
    residuals=matrix(NA,ncol=length(time),nrow=dim(a)[1])
    for (i in 1:dim(a)[1]){
      lines(t,log(a[i,]),col=rgb(0.3,0.3,0.3,0.3),lwd=2)
      rate1=piecewise_LM(log(a[i,]),t)$slope1
      rate2=piecewise_LM(log(a[i,]),t)$slope2
      residuals[,]=piecewise_LM(log(a[i,]),t)$res
      pbreak=piecewise_LM(log(a[i,]),t)$pb
      intercept=piecewise_LM(log(a[i,]),t)$intercept
      plot(piecewise_LM(log(a[i,]),t)$fit,add=T,col=adjustcolor("blue",0.3),lty=2) #maybe take this out.. not sure yet 
      k1=c(k1,rate1)
      k2=c(k2,rate2)
      pb=c(pb,pbreak)
      intval=c(intval,intercept)
    }
    plot(NULL,main="Residuals for Population Simulations",ylim=c(-0.5,0.5),xlim=c(1,dim(a)[2]),xlab="Index", ylab="Residual",cex.lab=1.4,cex.main=1.2)
    rect(par("usr")[1],par("usr")[3],par("usr")[2],2*max(area,na.rm=TRUE),col = "lightgrey")
    cols=adjustcolor(rainbow(dim(a)[1]),0.9)
    for (i in 1:dim(a)[1]){
      lines(seq(1,dim(a)[2]),residuals[i,],type="p",col=cols[i])
    }
    abline(0,0,lty=2)
    box()
    qqnorm(residuals,col=cols,cex.lab=1.4,cex.main=1.2)
    abline(0,1)
    par(op)
  }
    else{
      k1=c()
      k2=c()
      pb=c()
      intval=c()
      residuals=matrix(NA,ncol=length(time),nrow=dim(a)[1])
      for (i in 1:dim(a)[1]){
        rate1=piecewise_LM(log(a[i,]),t)$slope1
        rate2=piecewise_LM(log(a[i,]),t)$slope2
        residuals[,]=piecewise_LM(log(a[i,]),t)$res
        pbreak=piecewise_LM(log(a[i,]),t)$pb
        intercept=piecewise_LM(log(a[i,]),t)$intercept
        k1=c(k1,rate1)
        k2=c(k2,rate2)
        pb=c(pb,pbreak)
        intval=c(intval,intercept)
      }
    }
  return(list("k1"=k1,"k2"=k2,"pb"=pb, "intval"=intval))
}
piecewise_pop_hist<-function(a,t,op,s){
  # Plotting Population Simulations 
  outp=par(mfrow=c(2,1))
  plot(1,type='n', xlim=range(t), ylim=range(log(a)),xlab="Time",
       ylab="log(Area)",main="Population Simulations",cex.lab=1.4,cex.main=1.4)
  for (i in 1:dim(a)[1]){
    lines(t,log(a[i,]),col=rgb(0.3,0.3,0.3,0.3),lwd=2)
    abline(v=piecewise_LM(log(a[i,]),t)$pb,col="red",lty=2)
    #plot(piecewise_LM(log(a[i,]),t)$fit,add=T,col=adjustcolor("blue",0.3),lty=2) 
  }
  #Plotting a histogram for all single lineage data + population mean rate on a histogram
  cells=hist_cells(s$data$rate, 0.01)
  hist(s$data$rate, breaks=cells, xlab="r (1/h)", main=paste(toString(s$name),"; N=", toString(total), sep=""),cex.lab=1.4, cex.main=1.2,col="lightblue")
  # abline(v=mean(pop_rates),col="blue",lwd=2)
  # abline(v=(mean(pop_rates)-(sd(pop_rates)/sqrt(length(pop_rates)))),col="blue",lty=2)
  # abline(v=(mean(pop_rates)+(sd(pop_rates)/sqrt(length(pop_rates)))),col="blue",lty=2)
  abline(v=mean(op$k2),col=adjustcolor("blue",0.9),lwd=2)
  abline(v=(mean(op$k2)-(sd(op$k2)/sqrt(length(op$k2)))),col="blue",lty=2)
  abline(v=(mean(op$k2)+(sd(op$k2)/sqrt(length(op$k2)))),col="blue",lty=2)
  abline(v=mean(op$k1),col="blue",lwd=2)
  abline(v=(mean(op$k1)-(sd(op$k1)/sqrt(length(op$k1)))),col="blue",lty=2)
  abline(v=(mean(op$k1)+(sd(op$k1)/sqrt(length(op$k1)))),col="blue",lty=2)
  par(outp) 
} #could put this function into the previous one 

#Getting breaks for a histogram 
hist_cells<-function(dat,int){
  lo=trunc(min(dat)*10)/10-0.1 #rounding down
  hi=trunc(max(dat)*10)/10+0.1 #rounding up 
  cells=seq(lo,hi,int)
  return(cells)
}
# Histogram only
histo<-function(k,s){
  cells=hist_cells(k,0.01)
  hist(k,breaks=cells,xlab="r (1/h)", main=paste(s,"; N=",length(k),sep=""),cex.lab=1.4, cex.main=1.2,col="lightblue")
  abline(v=mean(strain$data$rate), col="red")
  abline(v=(mean(strain$data$rate)+(sd(strain$data$rate))),col="red",lty=2)
}

#####################Main############################
datsetname="Levy"

# Choosing a data set 
x=dataset(datsetname)
area=x$area
times=x$times
data=x$data
residuals=x$residuals
area[area==0]=NA

strain_names=unique(data$genotype)
pickstrain=strain_names[26] #choose strain here!
#strain=subset_strain(data,area,times,pickstrain)
strain=subset_strain_normal(data,area,times,pickstrain,datsetname)

total=dim(strain$area)[1]
N=round(total*0.2)#using 20% of the total data for extrapolation
time=seq(0,35,0.5)

# Refitting the exponential model to the newly simulated data to obtain r_p
iterations=20
popdata=pop_sim_dat(strain,N,time,iterations)
#popdata_tail=pop_sim_dat_tail(strain,N,time,iterations)
#op_rates=pop_plot(popdata,time)
#piece_op_tail=piecewise_pop_plot(popdata_tail,time,plot=FALSE)
#piece_op=piecewise_pop_plot(popdata,time)

filename=paste("NormRes_",datsetname,strain$name,"_PopSim.png")
png(filename)
par(mfrow=c(2,1))
plot_growth(strain$area,strain$times,strain$name)
pop_plot(popdata,time,res=FALSE)
par(op)
dev.off()

# cutoff=mean(strain$data$rate)+(sd(strain$data$rate))
# sm1=mean(piece_op$k1); sm2=mean(piece_op$k2); 
# fm1=mean(piece_op_tail$k1); fm2=mean(piece_op_tail$k2)
# NCut=length(which(strain$data$rate>=cutoff))
# 
# filename=paste("PiecewieseFastRatePopSim_",datsetname,strain$name,"_N",N,"_Ncut",NCut,"_s1m",sm1,"_s2m",sm2,"_f1m",fm1,"_f2m",fm2,".png",sep="")
# png(filename)
# op=par(mfrow=c(2,1))
# pop_rates_tail=pop_plot(popdata_tail,time,res=FALSE)
# x=histo(strain$data$rate,strain$name)
# abline(v=mean(piece_op$k1),col="green")
# abline(v=mean(piece_op$k2),col="green")
# abline(v=mean(piece_op_tail$k2),col="blue")
# abline(v=mean(piece_op_tail$k1),col="blue")
# par(op)
# dev.off()

#piece_op=piecewiese_pop_plot(popdata_tail,time)

#piecewise_pop_hist(popdata_tail,time,piece_op,strain)


# DF<-data.frame("Calculations"=c("BreakPoint.x", "Intercept.y", "Slope1", "Slope2"),
#                "Mean"=c(toString(signif(mean(piece_op$pb),digits=3)),toString(signif(mean(piece_op$intval),digits=3)),toString(signif(mean(piece_op$k1),digits=3)),toString(signif(mean(piece_op$k2),digits=3))),
#                "StandardDeviation" =c(toString(signif(sd(piece_op$pb),digits=3)),toString(signif(sd(piece_op$intval),digits=3)),toString(signif(sd(piece_op$k1),digits=3)),toString(signif(sd(piece_op$k2),digits=3))),
#                "LNsampled"=c(rep(N,4)),
#                "PNSimulated"=c(rep(iterations,4)),
#                "Strain"=c(rep(strain$name,4)),
#                "Dataset"=c(rep(datsetname,4)))
# 
#  DF<-data.frame(c("BreakPoint.x", "Intercept.y", "Slope1", "Slope2"),
#                 c(toString(signif(mean(piece_op$pb),digits=3)),toString(signif(mean(piece_op$intval),digits=3)),toString(signif(mean(piece_op$k1),digits=3)),toString(signif(mean(piece_op$k2),digits=3))),
#                 c(toString(signif(sd(piece_op$pb),digits=3)),toString(signif(sd(piece_op$intval),digits=3)),toString(signif(sd(piece_op$k1),digits=3)),toString(signif(sd(piece_op$k2),digits=3))),
#                 c(rep(N,4)),
#                 c(rep(iterations,4)),
#                 c(rep(strain$name,4)),
#                 c(rep(datsetname,4)))
# 
# #formattable(DF)
# write.table(DF,"PopSimStats.txt",col.names=FALSE,row.names=FALSE,append=TRUE)
