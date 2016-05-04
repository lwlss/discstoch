# Analysing non well-behaved growth curves 
library(data.table)

dataset<-function(x,ou.rm=FALSE){ #can choose to take area where outliers have been removed
  if (x == "Lawless"){
    # DataSet1: Lawless
    if(ou.rm==FALSE){area=fread("Lawless_area.txt",header=FALSE)}else{area=fread("Lawless_area_or.txt",header=FALSE)}
    times=fread("Lawless_time.txt",header=FALSE)
    data=fread("Lawless_data.txt",header=TRUE) #3rd column (Identifier) => strain_parentcolony 
    residuals=fread("Lawless_residuals.txt",header=FALSE)
    return(list("area"=area,"data"=data,"times"=times,"residuals"=residuals))
  }
  else if (x == "Levy"){
    # DataSet2: Levy
    if(ou.rm==FALSE){area=fread("Levy_area.txt",header=FALSE)}else{area=fread("Levy_area_or.txt",header=FALSE)}
    times=fread("Levy_time.txt",header=FALSE)
    data=fread("Levy_data.txt",header=TRUE) #3rd column (Identifier) => replicate
    residuals=fread("Levy_residuals.txt",header=FALSE)
    return(list("area"=area,"data"=data,"times"=times,"residuals"=residuals))
  }
  else if (x == "Ziv"){
    # DataSet3: Ziv
    if(ou.rm==FALSE){area=fread("Ziv_area.txt",header=FALSE)}else{area=fread("Ziv_area_or.txt",header=FALSE)}
    times=fread("Ziv_time.txt",header=FALSE)
    data=fread("Ziv_data.txt",header=TRUE) #3rd column (Identifier) => colony
    residuals=fread("Ziv_residuals.txt",header=FALSE)
    return(list("area"=area,"data"=data,"times"=times,"residuals"=residuals))
  }
  else {print("Not a valid dataset")}
}
plot_growth<-function(a,t,s,c){ 
  if (min(a,na.rm=TRUE)==0){
    y_range=c(1,max(a,na.rm=TRUE))
  }else{y_range=range(a,na.rm=TRUE)}
  plot(1,type='n', xlim=range(t,na.rm=TRUE), ylim=y_range,xlab="Time (h)", 
       ylab="Microcolony Area log(px)",main=paste(s),cex.lab=1.2, log='y')
  for (i in 1:dim(a)[1]){
    lines(as.numeric(t[i,]),as.numeric(a[i,]),col=c,lwd=2)
  }
}

datsetname="Lawless"
normtests=fread(paste(datsetname,"_residuals_normtests.txt",sep=""),header=TRUE)

#Durbing Watson Test: autocorrelation in the residuals 
#H0: there is no correlation among residuals (i.e. there are independent)
DW_indices=which(normtests$DurbinWatson<0.05)
#For Ziv there all data points are independent! (pre-filtering)

#Breusch Pagan Test: testing for homoskedasticity
#H0: the data is homoskedastic
BP_indices=which(normtests$BreuschPagan<0.05)

#A lot of the growth curves (~1000, Lawless) which are autocorrelated are also homoskedastic
#although one may very well be the result of the other... 
print(length(intersect(DW_indices,BP_indices)))
intersect1=intersect(DW_indices,BP_indices)

#Lilliefors Test: testing for normality in the residuals
#H0: the residuals are normally distributed 
LF_indices=which(normtests$Lilliefors.p<0.05) #less stringent; more likely to be normal 
#Shapiro-Wilk Test: testing for normality in the residuals 
#H0: the residuals are normally distributed 
SW_indices=which(normtests$ShapiroWilk.p<0.05)

#Almost all of the Lilliefors tested against normality also test against normality with the Shapiro-Wilk
print(length(intersect(LF_indices,SW_indices)))
intersect2=intersect(LF_indices,SW_indices)

#The following growth curves have failed all of the above tests:
indices=intersect(intersect1,intersect2)
print(length(indices))

# indices=intersect(intersect2,BP_indices) #Ziv
# # indices=BP_indices #Levy
# indices=sample(indices,100) #choosing to plot only one hundred of them 

#Extracting Growth Curves which have failed all of the above tests 
x=dataset(datsetname)
area=x$area[indices,]
area[area==0]=NA
times=x$times[indices,]
residuals=x$residuals[indices,]


colours=rainbow(length(indices))

pdf(height = 16, width = 16, file = paste(datsetname,"_NonWellBehaved.pdf",sep=""))
op=par(mfrow=c(3,2))
for(i in 1:length(indices)){
  plot_growth(area[i,],times[i,],paste("Growth Curve", toString(indices[i])),colours[i])
  plot(NULL,main=paste("Residuals of Growth Curve",toString(indices[i])),
       xlim=range(times[i,],na.rm=TRUE), ylim=range(residuals[i,],na.rm=TRUE))
  lines(as.numeric(times[i,]),as.numeric(residuals[i,]),type="p",col=colours[i])
}
par(op)
dev.off()


