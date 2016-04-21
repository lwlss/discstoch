# Extracting Growth Cruves from the Formatted Data Sets: Lawless, Levy & Ziv 

# setwd("~/discstoch-master/FormattedData")

library(data.table)

#####################Functions############################

# Choosing a data set to work with
dataset<-function(x){
  print(c("Lawless","Levy","Ziv"))
  x=readline(prompt="Which dataset? ")
  if (x == "Lawless"){
    # DataSet1: Lawless
    print("Loading dataset...")
    area=fread("Lawless_area.txt",header=FALSE)
    times=fread("Lawless_time.txt",header=FALSE)
    data=fread("Lawless_data.txt",header=FALSE) #3rd column (Identifier) => strain_parentcolony 
    return(list("area"=area,"data"=data,"times"=times))
  }
  else if (x == "Levy"){
    # DataSet2: Levy
    print("Loading dataset...")
    area=fread("Levy_area.txt",header=FALSE)
    times=fread("Levy_time.txt",header=FALSE)
    data=fread("Levy_data.txt",header=FALSE) #3rd column (Identifier) => replicate
    return(list("area"=area,"data"=data,"times"=times))
  }
  else if (x == "Ziv"){
    # DataSet3: Ziv
    print("Loading dataset...")
    area=fread("Ziv_area.txt",header=FALSE)
    times=fread("Ziv_time.txt",header=FALSE)
    data=fread("Ziv_data.txt",header=FALSE) #3rd column (Identifier) => colony
    return(list("area"=area,"data"=data,"times"=times))
  }
  else {print("Not a valid dataset. Retry..."); dataset()}
}

# Subsetting the data according to strain (genotype)
subset_strain<-function(d,a,t){ 
  #where d is the formatted data;
  #a and t are the respective area and time 
  strain_names=unique(d$genotype)
  print(strain_names)
  strain=readline(prompt="Extract data for strain: ")
  if(length(intersect(strain,strain_names))!=1){
    print("Invalids strain name. Retry...")
    subset_strain(d,a,t)
  }
  else{
    s_name=toString(strain)
    indices=which(d$genotype == strain)
    s_area=a[indices,]; s_times=t[indices,]; s_data=d[indices,]
    return(list("area"=s_area,"times"=s_times,"data"=s_data, "name"=s_name))
  }
}

# Subsetting the data according to an identifier 
subset_identifier<-function(d,a,t){
  #where d is the formatted data;
  #a and t are the respective area and time 
  identifier_names=unique(d$identifier)
  print(identifier_names)
  identifier=readline(prompt="Extract data for identifier: ")
  if(length(intersect(identifier,identifier_names))!=1){
    print("Invalid identifier. Retry...")
    subset_identifier(d,a,t)
  }
  else{
    i_name=toString(identifier)
    indices=which(d$identifier == identifier)
    i_area=a[indices,]; i_times=t[indices,]; i_data=d[indices,]
    return(list("area"=i_area,"times"=i_times,"data"=i_data, "name"=i_name))
  }
}

# Subsetting the data according to a clonal colony
subset_colony<-function(d,a,t){
  #where d is the formatted data;
  #a and t are the respective area and time 
  colony_names=unique(d$clonalcolony)
  print(colony_names)
  colony=readline(prompt="Extract data for clonal colony: ")
  if(length(intersect(colony,colony_names))!=1){
    print("Invalid clonal colony. Retry...")
    subset_colony(d,a,t)
  }
  else{
    c_name=toString(colony)
    indices=which(d$clonalcolony == colony)
    c_area=a[indices,]; c_times=t[indices,]; c_data=d[indices,]
    return(list("area"=c_area,"times"=c_times,"data"=c_data, "name"=c_name)) 
  }
}

# Subsetting the data according to a specific clonal colony, identifier and/ or genotype 
subset_3vars<-function(d,a,t){
  #where d is the formatted data;
  #a and t are the respective area and time;
  gnames=unique(d$genotype)
  print(gnames)
  gen=readline(prompt="Which genotype? (Or enter 'total'.) ")
  if (length(intersect(gnames,gen))!=1 & length(intersect(gen,"total"))!=1){
    print("Invalid genotype. Retry...")
    subset_3vars(d,a,t)
  }
  else{
    cnames=unique(d$clonalcolony)
    print(cnames)
    cc=readline(prompt="Which clonal colony? (Or enter 'total'.) ")
    if (length(intersect(cnames,cc))!=1 & length(intersect(cc,"total"))!=1){
      print("Invalid clonal colony Retry...")
      subset_3vars(d,a,t)
    }
    else{
      inames=unique(d$identifier)
      print(inames)
      id=readline(prompt="Which identifier? (Or enter 'total'.) ")
      if (length(intersect(inames,id))!=1 & length(intersect(id,"total"))!=1){
        print("Invalid clonal colony Retry...")
        subset_3vars(d,a,t)
      }
      else{
        #gen, cc, and id are the genotype, clonalcolony, and id respectively
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
    }
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
plot_growth<-function(a,t,s,Nsample=100){ 
  if (dim(a)[1] < Nsample) {Nsample=dim(a)[1]} #cannot take sample larger than the population when replace=F
  #where area (a), times (t) and name (s) are required as input variables
  indices=sample(1:dim(a)[1],Nsample,replace=FALSE)
  op=par(mfrow=c(2,1))
  plot(1,type='n', xlim=range(t[indices,],na.rm=TRUE), ylim=range(a[indices,],na.rm=TRUE),xlab="Time (h)", 
       ylab="Microcolony Area (px)",main=paste(s),cex.lab=1.2)
  k=c()
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

# Plotting growth curves (colour-coding the identifiers)
plot_growth_colour<-function(a,t,d,s,Nsample=100){
  if (dim(a)[1] < Nsample) {Nsample=dim(a)[1]} #cannot take sample larger than the population when replace=F
  #where area (a), times (t), data (d), and name (s) are required as input variables
  indices=sample(1:dim(a)[1],Nsample,replace=FALSE)
  ids=d$identifier[indices]
  id_names=unique(ids)
  cols=adjustcolor(rainbow(length(id_names)),0.5)
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
  legend("topleft",legend=id_names,pch=15,col=cols,horiz=TRUE)
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
  plot(1,type='n', xlim=c(list_range(tl)$mi,list_range(tl)$ma), ylim=c(0,max(maxval)),xlab="Time (h)",
       ylab="Microcolony Area (px)",main="Overlaid Growth Curves",cex.lab=1.2)
  rect(par("usr")[1],par("usr")[3],par("usr")[2],2*max(area,na.rm=TRUE),col = "lightgrey")
  cols=adjustcolor(rainbow(length(al)),0.3)
  rates=list()
  names=c()
  histcounts=c()
  for (l in 1:length(al)){
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



#####################Main############################

# Choosing a Data set
x=dataset()
area=x$area
times=x$times
data=x$data
# All formated *_data.txt files should be in the form: Genotype, Clonal Colony (well), Identifier 
# NB: No headers on any of the files.
colnames(data)=c("genotype","clonalcolony","identifier")
# Getting rid of zero areas
# NB in Lawless's data zero's were converted to ones... 
area[area==0]=NA

# Extracting data for one strain 
strain=subset_strain(data,area,times)
# Plotting growth curves for that strain
plot_growth(strain$area,strain$times,strain$name, Nsample=50)
plot_growth_colour(strain$area,strain$times,strain$data,strain$name,Nsample=5)

# Extracting data for one identifier 
identifier=subset_identifier(data,area,times)
# Plotting growth curves for that identifier
plot_growth(identifier$area,identifier$times,identifier$name, Nsample=50)

# Extracting data for one clonal colony
clonalcolony=subset_colony(data,area,times) #no Nsample here (not necessary)
# Plotting growth curves for that clonal colony 
plot_growth(clonalcolony$area,clonalcolony$times,clonalcolony$name, Nsample=50)

# Subsetting data according to a specific genotype, clonal colony and/or identifier 
subdata=subset_3vars(data,area,times)
# Plotting the corresponding growth curves
plot_growth(subdata$area,subdata$times,subdata$name)

# Overlaying growth curves from various strains/identifiers/clonalcolonies
# NB: code works for more than 2 strains/identifiers/clonalcolonies
strain1=subset_strain(data,area,times)
strain2=subset_strain(data,area,times)
strain3=subset_strain(data,area,times)
strain4=subset_strain(data,area,times)
a4=list(strain1$area,strain2$area,strain3$area,strain4$area)
t4=list(strain1$times,strain2$times,strain3$times,strain4$times)
s4=list(strain1$name,strain2$name,strain3$name,strain4$name)
plot_growth_overlay(a4,t4,s4)
# identifier1=subset_identifier(data,area,times)
# identifier2=subset_identifier(data,area,times)
# identifier3=subset_identifier(data,area,times)
# identifier4=subset_identifier(data,area,times)
# a4=list(identifier1$area,identifier2$area,identifier3$area,identifier4$area)
# t4=list(identifier1$times,identifier2$times,identifier3$times,identifier4$times)
# s4=list(identifier1$name,identifier2$name,identifier3$name,identifier4$name)
# plot_growth_overlay(a4,t4,s4)
