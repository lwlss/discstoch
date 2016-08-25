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
    area=fread("Levy_area.txt",header=FALSE)
    times=fread("Levy_time.txt",header=FALSE)
    data=fread("Levy_data.txt",header=FALSE) #3rd column (Identifier) => replicate
    info=read.table("Levy_GrowthRateInfo.txt",header=TRUE,row.names=1)
    names(data)=c("genotype","clonalcolony","identifier")
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
datsetname="Lawless"
x=dataset(datsetname)
area=x$area
times=x$times
data=x$data
info=x$info

if (datsetname=="Lawless"){
  rates=c()
  for (i in 1:dim(area)[1]){
    k=EXP_growthrate(area[i,],times[i])$rate
    rates=c(rates,k)
  }
} else{
  LevyP=c("YME1","PET9","YFR054C","YHR095W","SNF6","RAD50","NOT5","HTZ1")
  if(datsetname=="Levy"){
    area=area[which(data$genotype %in% LevyP),]
    times=times[which(data$genotype %in% LevyP),]
    data=data[which(data$genotype %in% LevyP),]
    info=info[which(data$genotype %in% LevyP),]
  }else{}
  for(i in 1:dim(area)[1]){
    area[i,which(area[i,]==0)]=NA
  }
  rates=c()
  for (i in 1:dim(area)[1]){
    k=detstocgrowth::LM_growthrate(area[i,],times[i,])$rate
    rates=c(rates,k)
  }
}

# Setting negative growth rates equal to zero.
rates[which(rates<0)]=0

#write.table(data.frame(Rate=t(rates),Data=data),file="LevyPublishedRates.txt",row.names=FALSE,col.names=TRUE)
op=par(mfrow=c(1,1))
DensityOverlay(rates,data,0.5,datsetname)
par(op)




# # Deterministic Bayesian Parameters (currently logistic -> change to exponential ones)
# params=c("1:100","101:200","201:300","301:400","401:500","501:600","601:700","701:800","801:900",
#          "901:1000","1001:1100","1101:1200","1201:1300","1301:1400","1401:1500","1501:1600",
#          "1601:1700","1701:1800","1801:1900")
#
# total_params=matrix(0,ncol=2,nrow=length(params)*100+46)
# for (i in 1:length(params)){
#   #print((1:100)+(100*(i-1)))
#   print(params[i])
#   total_params[(1:100)+(100*(i-1)),]=as.matrix(read.table(paste("Lawless_Bayes_parameters_BayesDetExp_",params[i],".txt",sep=""),header=TRUE))
# }
# total_params[1901:1946,]=as.matrix(read.table("Lawless_Bayes_parameters_BayesDetExp_1901:1946.txt",header=TRUE))
#
# # Intercept, Rate, CarryingCapacity
# total_rates_bayes=total_params[,2]
# rates=total_rates_bayes
#
# colours=c("darkgreen","orange")
#
# for (i in 1:length(unique(data$genotype))){
#   strain_rates=rates[which(data$genotype==unique(data$genotype)[i])]
#   lines(density(strain_rates),col=colours[i],lwd=4)
#   polygon(density(strain_rates)$x,density(strain_rates)$y,col=adjustcolor(colours[i],0.1))
# }
# legend("topright",legend=c("HIS3 Freq.","HIS3 Bayes.","HTZ1 Freq.","HTZ1 Bayes."),col=c("blue","red","yellow","orange"),
#        lty=rep(1,2*(length(unique(data$genotype)))),
#        lwd=rep(3,2*(length(unique(data$genotype)))),cex=0.9)
