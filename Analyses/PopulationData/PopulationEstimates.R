# Pin Population Estimates over Time

library(detstocgrowth)

# Getting the population data as produced by the Image Analysis tool muqfatc
area=read.table("~/discstoch/Analyses/ImageAnalysis/PopulationArea.txt",header=FALSE)
times=read.table("~/discstoch/Analyses/ImageAnalysis/PopulationTime.txt",header=FALSE)
folders=read.table("~/discstoch/Analyses/ImageAnalysis/PopulationFolders.txt",header=FALSE)
folders=as.vector(t(folders))
data=read.table("~/discstoch/Analyses/LineageData/Lawless_data_shortTC.txt",header=FALSE)
names(data)=c("genotype","clonalcolony","identifier","blobnumber")

strains=c()
for (i in folders){
  strains=c(strains,as.vector(unique(data[which(data$clonalcolony==i),]$genotype)))
}

# Sub-dividng the data according to strain
strain="HTZ1"
sf=folders[which(strains==strain)]
sa=area[which(strains==strain),]
st=times[which(strains==strain),]

# Plotting the Pin Population Estimates
PinDat=pin_pop_plot(strain,sa,st)
print(mean(PinDat$bp1))
print(mean(PinDat$bp2))
print(mean(PinDat$s1)); print(mean(PinDat$s2)); print(mean(PinDat$s3))

# Comparin Pin Estimates to the Population Simulation from Single Lineage Data

# Getting the Bayesian Parameter Estimates
params=c("1:100","101:200","201:300","301:400","401:500","501:600","601:700","701:800","801:900",
         "901:1000","1001:1100","1101:1200","1201:1300","1301:1400","1401:1500","1501:1600",
         "1601:1700","1701:1800","1801:1900")
total_params=matrix(0,ncol=2,nrow=length(params)*100+46)
for (i in 1:length(params)){
  total_params[(1:100)+(100*(i-1)),]=as.matrix(read.table(paste("~/discstoch/Analyses/LineageData/Lawless_Bayes_parameters_BayesDetExp_",params[i],".txt",sep=""),header=TRUE))
}
total_params[1901:1946,]=as.matrix(read.table("~/discstoch/Analyses/LineageData/Lawless_Bayes_parameters_BayesDetExp_1901:1946.txt",header=TRUE))
total_rates_bayes=total_params[,2]
rates=total_rates_bayes

lineagearea=read.table("~/discstoch/Analyses/LineageData/Lawless_area_shortTC.txt",header=FALSE)
lineagetime=read.table("~/discstoch/Analyses/LineageData/Lawless_time_shortTC.txt",header=FALSE)
lineagedata=read.table("~/discstoch/Analyses/LineageData/Lawless_data_shortTC.txt",header=FALSE)
names(lineagedata)=c("genotype","pin","parent","blob")

sortedfolders=as.vector(unique(data$clonalcolony))

# pdf(file="PinGrowthRates.pdf",title="Observed and Simulated Population Growth",width=20,height=28)
# par(mfrow=c(8,8))
# plot(NULL)
sortedfolders=sortedfolders[18]
obs=30
conv=92.5
it=100
Comp=PlotPinObsSim(area,times,sortedfolders,obs,conv,it,leg=FALSE)
print(Comp$tp)
print(mean(Comp$lagr))
print(mean(Comp$expr))
# dev.off()


