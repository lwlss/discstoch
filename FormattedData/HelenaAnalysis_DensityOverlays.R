library(data.table)
library(segmented)

dataset<-function(x){
  if (x == "Lawless"){
    # DataSet1: Lawless
    area=fread("Lawless_area_shortTC.txt",header=FALSE)
    times=fread("Lawless_time_shortTC.txt",header=FALSE)
    data=fread("Lawless_data_shortTC.txt",header=FALSE) #3rd column (Identifier) => strain_parentcolony
    info=read.table("Lawless_GrowthRateInfo.txt",header=TRUE,row.names=1)
    names(data)=c("genotype","clonalcolony","identifier","blobnumber")
    return(list("area"=area,"data"=data,"times"=times,"info"=info))
  }
  else if (x == "Levy"){
    # DataSet2: Levy
    area=fread("Levy_area_filtered.txt",header=FALSE)
    times=fread("Levy_times_filtered.txt",header=FALSE)
    data=fread("Levy_data_filtered.txt",header=TRUE) #3rd column (Identifier) => replicate
    info=read.table("Levy_GrowthRateInfo.txt",header=TRUE,row.names=1)
    return(list("area"=area,"data"=data,"times"=times,"info"=info))
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
datsetname="Ziv"
x=dataset(datsetname)
area=x$area
times=x$times
data=x$data
info=x$info

rates=c()
for (i in 1:dim(area)[1]){
  k=detstocgrowth::LM_growthrate(area[i,],times[i,])$rate
  rates=c(rates,k)
}

colours=rev(rainbow(length(unique(data$genotype))))

plot(NULL,xlim=c(0,max(density(rates)$x)),ylim=c(0,11),
     xlab="Growth Rate (1/h)", ylab="Density", cex.lab=1.4,
     main="Growth Rate Density for Each Genotype")
for (i in 1:length(unique(data$genotype))){
  strain_rates=rates[which(data$genotype==unique(data$genotype)[i])]
  lines(density(strain_rates),col=colours[i],lwd=4)
  polygon(density(strain_rates)$x,density(strain_rates)$y,col=adjustcolor(colours[i],0.1))
}
legend("topleft",legend=(unique(data$genotype)),col=colours,
       lty=rep(1,(length(unique(data$genotype)))),
       lwd=rep(2,(length(unique(data$genotype)))))
