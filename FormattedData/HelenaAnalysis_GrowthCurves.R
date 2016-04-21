# Extracting Growth Cruves from the Formatted Data Sets: Lawless, Levy & Ziv 

#####################Functions############################

# Choosing a data set to work with
dataset<-function(x){
  if (x == "Lawless"){
    # DataSet1: Lawless
    area=read.table("Lawless_area.txt",header=FALSE)
    times=read.table("Lawless_time.txt",header=FALSE)
    data=read.table("Lawless_data.txt",header=FALSE) #3rd column (Identifier) => strain_parentcolony 
    return(list("area"=area,"data"=data,"times"=times))
  }
  else if (x == "Levy"){
    # DataSet2: Levy
    area=read.table("Levy_area.txt",header=FALSE)
    times=read.table("Levy_time.txt",header=FALSE)
    data=read.table("Levy_data.txt",header=FALSE) #3rd column (Identifier) => replicate
    return(list("area"=area,"data"=data,"times"=times))
  }
  else if (x == "Ziv"){
    # DataSet3: Ziv
    area=read.table("Ziv_area.txt",header=FALSE)
    times=read.table("Ziv_time.txt",header=FALSE)
    data=read.table("Ziv_data.txt",header=FALSE) #3rd column (Identifier) => colony
    return(list("area"=area,"data"=data,"times"=times))
  }
  else {print("Not a valid dataset")}
}

# Subsetting the data according to strain (genotype)
subset_strain<-function(d,a,t,strain){ 
  s_name=toString(strain)
  indices=which(d$genotype == strain)
  s_area=a[indices,]; s_times=t[indices,]; s_data=d[indices,]
  return(list("area"=s_area,"times"=s_times,"data"=s_data, "name"=s_name))
}

# Subsetting the data according to an identifier 
subset_identifier<-function(d,a,t,identifier){
  i_name=toString(identifier)
  indices=which(d$identifier == identifier)
  i_area=a[indices,]; i_times=t[indices,]; i_data=d[indices,]
  return(list("area"=i_area,"times"=i_times,"data"=i_data, "name"=i_name))
}

# Subsetting the data according to a clonal colony
subset_colony<-function(d,a,t,colony){
  c_name=toString(colony)
  indices=which(d$clonalcolony == colony)
  c_area=a[indices,]; c_times=t[indices,]; c_data=d[indices,]
  return(list("area"=c_area,"times"=c_times,"data"=c_data, "name"=c_name))
}

# Subsetting the data according to a specific clonal colony, identifier and/ or genotype 
subset_3vars<-function(d,a,t,gen=total,cc=total,id=total){
  if (gen == "total"){
    indices1=seq(1:dim(a)[1])
  }else{indices1=which(d$genotype == gen)}
  if (cc == "total"){
    indices2=seq(1:dim(a)[1])
  }else{indices2=which(d$clonalcolony == cc)}
  if (id == "total"){
    indices3=seq(1:dim(a)[1])
  }else{indices3=which(d$identifier == id)}
  indices=intersect(indices1,indices2)
  indices=intersect(indices,indices3)
  if (length(indices) == 0){
    print("Input Error: Invalid combination")
  }
  else{
    f_name=paste("Genotype: ", toString(gen), " Colony: ", toString(cc), " Identifier: ", toString(id))
    f_area=a[indices,]; f_times=t[indices,]; f_data=d[indices,]
    return(list("area"=f_area,"times"=f_times,"data"=f_data, "name"=f_name))
  }
}

#Calculating the growth rate (LS for exponential model)
LM_growthrate<-function(ai,ti){
  fit<-lm(log(ai)~ti)
  rate=fit$coefficient[[2]]
  return(rate)
}

#Getting breaks for a histogram 
hist_cells<-function(dat,int){
  lo=trunc(min(dat)*10)/10-0.1 #rounding down
  hi=trunc(max(dat)*10)/10+0.1 #rounding up 
  cells=seq(lo,hi,int)
  return(cells)
}

# Plotting growth curve
plot_growth<-function(a,t,s,Nsample=100,plot=TRUE){ 
  if (dim(a)[1] < Nsample) {Nsample=dim(a)[1]} #cannot take sample larger than the population when replace=F
  #where area (a), times (t) and name (s) are required as input variables
  indices=sample(1:dim(a)[1],Nsample,replace=FALSE)
  k=c()
  if(plot==TRUE){
    op=par(mfrow=c(2,1))
    plot(1,type='n', xlim=range(t[indices,],na.rm=TRUE), ylim=range(a[indices,],na.rm=TRUE),xlab="Time (h)", 
         ylab="Microcolony Area (px)",main=paste(s),cex.lab=1.2)
    for (i in indices){
      lines(as.numeric(t[i,]),as.numeric(a[i,]),col=rgb(0.3,0.3,0.3,0.3),lwd=2)
      rate=LM_growthrate(as.numeric(area[i,]),as.numeric(times[i,]))
      k=c(k,rate)
    }
    #Growth rate Distribution
    #k=k[k>=0 & k<=0.5]
    #lo=0
    #hi=0.5
    # cells=seq(lo,hi,0.01)
    cells=hist_cells(k,0.01)
    hist(k,breaks=cells,xlab="r (1/h)", main=paste(Nsample," microcolonies"),cex.lab=1.2)
    par(op)
  }
  else {
    for (i in indices){
      rate=LM_growthrate(as.numeric(area[i,]),as.numeric(times[i,]))
      k=c(k,rate)
    }
  }
  return(list("rates"=k,"ind"=indices))
}

# Histogram only
histo<-function(a,t,s){
  k=c()
  for (i in 1:dim(a)[1]){
    rate=LM_growthrate(as.numeric(area[i,]),as.numeric(times[i,]))
    k=c(k,rate)
    op=par(mfrow=c(1,1))
    cells=hist_cells(k,0.05)
    hist(k,breaks=cells,xlab="r (1/h)", main=paste(s),cex.lab=1.4, cex.main=1.2,col="lightblue")
    par(op)
  }
}

# Plotting growth curves (colour-coding the identifiers)
plot_growth_colour<-function(a,t,d,s,Nsample=100){
  if (dim(a)[1] < Nsample) {Nsample=dim(a)[1]} #cannot take sample larger than the population when replace=F
  #where area (a), times (t), data (d), and name (s) are required as input variables
  indices=sample(1:dim(a)[1],Nsample,replace=FALSE)
  ids=d$identifier[indices]
  id_names=unique(ids)
  cols=adjustcolor(rainbow(length(id_names)),0.3)
  op=par(mfrow=c(2,1))
  plot(1,type='n', xlim=range(t[indices,],na.rm=TRUE), ylim=range(a[indices,],na.rm=TRUE),xlab="Time (h)", 
       ylab="Microcolony Area (px)",main=paste(s),cex.lab=1.2)
  rect(par("usr")[1],par("usr")[3],par("usr")[2],2*max(area,na.rm=TRUE),col = "lightgrey")
  k=c()
  for (i in indices){
    lines(as.numeric(t[i,]),as.numeric(a[i,]),col=cols[which(id_names==d[i,]$identifier)], lwd=2)
    rate=LM_growthrate(as.numeric(area[i,]),as.numeric(times[i,]))
    k=c(k,rate)
  }
  legend("topleft",legend=id_names,pch=15,col=cols)
  #Growth rate Distribution
  cells=hist_cells(k,0.01)
  hist(k,breaks=cells,xlab="r (1/h)",main=paste(Nsample," microcolonies"),cex.lab=1.2)
  par(op)
}

# Finding the range of all elements in a list
list_range<-function(l){
  #where l is the list of which to find the range 
  ma=c()
  mi=c()
  for (i in 1:length(l)){
    ma=c(ma,max(l[[i]],na.rm=TRUE))
    mi=c(mi,min(l[[i]],na.rm=TRUE))
  }
  return(list("ma"=max(ma),"mi"=min(mi)))
}

# Overlaying multiple strains/identifiers/clonal colonies in a plot 
plot_growth_overlay<-function(al,tl,sl,Nsample=100){
  maxval=c()
  indices=list()
  for(l in 1:length(al)){
    if (dim(al[[l]])[1] < Nsample) {Nsample=dim(al[l])[1]} #cannot take sample larger than the population when replace=F
    indices[[l]]=sample(1:dim(al[[l]])[1],Nsample,replace=FALSE)
    maxval=c(maxval,max(al[[l]][indices[[l]],],na.rm=TRUE))
  }
  op=par(mfrow=c(2,1))
  # plot(1,type='n', xlim=c(list_range(tl)$mi,list_range(tl)$ma), ylim=c(list_range(al)$mi,list_range(al)$ma),xlab="Time (h)", 
  #      ylab="Microcolony Area (px)",main="Overlaid Growth Curves",cex.lab=1.2)
  plot(1,type='n', xlim=c(list_range(tl)$mi,list_range(tl)$ma), ylim=c(0,max(maxval)),xlab="Time (h)",
       ylab="Microcolony Area (px)",main="Overlaid Growth Curves",cex.lab=1.2)
  rect(par("usr")[1],par("usr")[3],par("usr")[2],2*max(area,na.rm=TRUE),col = "lightgrey")
  cols=adjustcolor(rainbow(length(al)),0.3)
  rates=list()
  names=c()
  histcounts=c()
  for (l in 1:length(al)){
    # if (dim(al[[l]])[1] < Nsample) {Nsample=dim(al[l])[1]} #cannot take sample larger than the population when replace=F
    # indices=sample(1:dim(al[[l]])[1],Nsample,replace=FALSE)
    ind=indices[[l]]
    names=c(names,sl[[l]])
    k=c()
    for (i in ind){
      lines(as.numeric(tl[[l]][i,]),as.numeric(al[[l]][i,]),col=cols[l],lwd=2)
      rate=LM_growthrate(as.numeric(al[[l]][i,]),as.numeric(tl[[l]][i,]))
      k=c(k,rate)
    }
    rates[[l]]=k
    cells=hist_cells(k,0.01)
    histcounts=c(histcounts,max(hist(k,breaks=cells,plot=FALSE)$counts))
  }
  legend("topleft",legend=names,pch=15,col=cols)
  #Growth rate Distribution
  cells=hist_cells(c(list_range(rates)$mi,list_range(rates)$ma),0.01)
  hist(0,breaks=cells,xlab="r (1/h)", main=paste(Nsample," microcolonies"),cex.lab=1.2, ylim=c(0,max(histcounts)))
  rect(par("usr")[1],par("usr")[3],par("usr")[2],2*max(area,na.rm=TRUE),col = "lightgrey")
  for (l in 1:length(al)){
    hist(rates[[l]],breaks=cells,col=cols[l],add=T)
  }
  box()
  legend("topleft",legend=names,pch=15,col=cols)
  par(op)
}

# Simulating exponential growth  
simExponential<-function(r,t,N0=1){
  N=N0*exp(r*t)
  return(N)
}

#Function to plot a single estimated population growth curve
single_pop_plot<-function(Narea,Ntime,title,log=TRUE){
  if (log == TRUE){
    plot(Ntime,Narea,main=title,log="y",xlab="Time",ylab="log(Area)",cex.lab=1.4,cex.main=1.2)
  }
  else{
    plot(Ntime,Narea,main=title,xlab="Time ",ylab="Area",cex.lab=1.4,cex.main=1.2)
  }
}

# Simulating multiple population growth curves 
pop_sim_dat<-function(strain,N,t,it){
  #where s is a single strain name, N is the sample size, 
  #and it is the number of iterations
  PopData=matrix(0,nrow=it,ncol=length(t))
  for (i in 1:it){
    truevals=plot_growth(strain$area,strain$times,strain$name,Nsample=N,plot=FALSE)
    rates=truevals$rates
    expdata=strain$area[truevals$ind,]
    exptime=strain$times[truevals$ind,]
    simdata=matrix(0,nrow=length(rates),ncol=length(t))
    for (k in 1:length(rates)){
      simdata[k,]=simExponential(rates[k],t)
    }
    Npop_sim=colSums(simdata)
    PopData[i,]=Npop_sim
  }
  return(PopData)
}

# Plotting multiple population growth curves
pop_plot<-function(a,t){ 
  op=par(mfrow=c(1,1))
  plot(1,type='n', xlim=range(t), ylim=range(a),xlab="Time", 
       ylab="log(Area)",main="Population Simulations",cex.lab=1.4,cex.main=1.2,log='y')
  k=c()
  for (i in 1:dim(a)[1]){
    lines(t,a[i,],col=rgb(0.3,0.3,0.3,0.3),lwd=2)
    rate=LM_growthrate(a[i,],t)
    k=c(k,rate)
  }
  par(op)
  return(k)
}

#####################Main############################

# Choosing a data set 
x=dataset("Lawless")
area=x$area
times=x$times
data=x$data
colnames(data)=c("genotype","clonalcolony","identifier")
area[area==0]=NA

strain_names=unique(data$genotype)
pickstrain=strain_names[2] #choose strain here!
strain=subset_strain(data,area,times,pickstrain) 
total=length(which(data$genotype==pickstrain))
N=round(total*0.2)#using 20% of the total data for extrapolation
time=seq(0,35,0.5)

# Refitting the exponential model to the newly simulated data to obtain r_p
data=pop_sim_dat(strain,N,time,20)
pop_rates=pop_plot(data,time)
mean(pop_rates); sd(pop_rates)

histo(strain$area,strain$times,pickstrain)
abline(v=mean(pop_rates),col="blue",lwd=2)
abline(v=(mean(pop_rates)-sd(pop_rates)),col="blue",lty=2)
abline(v=(mean(pop_rates)+sd(pop_rates)),col="blue",lty=2)

# Single population plots (from before)

truevals=plot_growth(strain$area,strain$times,strain$name,Nsample=N)
rates=truevals$rates
expdata=strain$area[truevals$ind,]
exptime=strain$times[truevals$ind,]

# Entirely new data
time=seq(0,35,0.5)
simdata=matrix(0,nrow=length(rates),ncol=length(time))
for (k in 1:length(rates)){
  simdata[k,]=simExponential(rates[k],time)
}

# Using exisiting growth data and simulating forwards (This would probably make much more sense to do stochastically)
lasttime=exptime[1,dim(exptime)[2]] #simulating forward from the last measured time interval (assuming they're all roughly the same)
forwardtime=seq(lasttime+0.3,3*lasttime,0.3) #simulating forwards twice the length of the original experiment
forwarddata=matrix(0,nrow=length(rates),ncol=length(forwardtime))
newN0=expdata[,dim(expdata)[2]] #getting the last measured area as a starting pop'n for forward simulation
for (k in 1:length(rates)){
  forwarddata[k,]=simExponential(rates[k],forwardtime,N0=newN0[k])
}
totaldata=cbind(expdata,forwarddata)
totaltime=cbind(t(as.numeric(exptime[1,])),t(forwardtime))
#this is assuming that all the time points are exactly the same (which isn't entirely true but close enough...)

# Population growth derived from microcolony growth rates
Npop_sim=colSums(simdata)
Npop_forwardsim=as.numeric(colSums(totaldata))
#10 fold difference in starting value extrapolates exponentially

# Plotting the population growth
op=par(mfrow=c(2,1))
single_pop_plot(Npop_sim,time,title="Simulated Population Growth",log=TRUE)
single_pop_plot(Npop_forwardsim,totaltime,title="Forward Simulated Population Growth",log=TRUE)
par(op)
#according to the exponential model this should be a straight line!
