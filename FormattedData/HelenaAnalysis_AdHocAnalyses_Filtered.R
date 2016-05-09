#ADD-HOC ANALYSES 

library(nortest)
library(moments)
library(data.table)
library(car)
library(lmtest)

############################################################################################################

# Analysing the lag phase at the population level
#x=fread("PopSimStats.txt",header=TRUE)

############################################################################################################

# Analysing the residuals at the single lineage level

####Functions####
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
subset_strain<-function(d,a,t,strain){ 
  s_name=toString(strain)
  indices=which(d$genotype == strain)
  s_area=a[indices,]; s_times=t[indices,]; s_data=d[indices,]
  return(list("area"=s_area,"times"=s_times,"data"=s_data, "name"=s_name, "indices"=indices))
}
subset_colony<-function(d,a,t,colony){
  c_name=toString(colony)
  indices=which(d$clonalcolony == colony)
  c_area=a[indices,]; c_times=t[indices,]; c_data=d[indices,]
  return(list("area"=c_area,"times"=c_times,"data"=c_data, "name"=c_name))
}
subset_identifier<-function(d,a,t,identifier){
  i_name=toString(identifier)
  indices=which(d$identifier == identifier)
  i_area=a[indices,]; i_times=t[indices,]; i_data=d[indices,]
  return(list("area"=i_area,"times"=i_times,"data"=i_data, "name"=i_name))
}
testing_normality<-function(r){
  tskew=c()
  tktosis=c()
  tpval1=c()
  tpval2=c()
  tTF1=c()
  tTF2=c()
  for(i in 1:dim(r)[1]){
    print(i)
    ri=as.numeric(r[i,])
    if(sum(ri,na.rm=TRUE)==0){
      print("Step1")
      skew=NA
      ktosis=NA
      pval1=NA
      TF1=TRUE #don't want to discard a perfect fit; although this will always be a straight line (all 1s or all 0s)
      pval2=NA
      TF2=TRUE 
      tskew=c(tskew,skew)
      tktosis=c(tktosis,ktosis)
      tpval1=c(tpval1,pval1)
      tpval2=c(tpval2,pval2)
      tTF1=c(tTF1,TF1)
      tTF2=c(tTF2,TF2)
    }
    else if(sum(!is.na(ri))<=4){
      print("Step2")
      skew=NA
      ktosis=NA
      pval1=NA
      TF1=FALSE #lineages of 4 or less data points will be discarded! 
      pval2=NA
      TF2=FALSE 
      tskew=c(tskew,skew)
      tktosis=c(tktosis,ktosis)
      tpval1=c(tpval1,pval1)
      tpval2=c(tpval2,pval2)
      tTF1=c(tTF1,TF1)
      tTF2=c(tTF2,TF2)
    }
    else{
      print("Step3")
      skew=skewness(ri,na.rm=TRUE)
      ktosis=kurtosis(ri,na.rm=TRUE)
      test1=lillie.test(ri) #H0: data are normally distributed; p>0.1: we do not reject H0
      test2=shapiro.test(ri)
      pval1=test1$p.value
      if (pval1<0.1){
        TF1=FALSE
      }else{TF1=TRUE}
      pval2=test2$p.value
      if (pval2<0.1){
        TF2=FALSE
      }else{TF2=TRUE}
      tskew=c(tskew,skew)
      tktosis=c(tktosis,ktosis)
      tpval1=c(tpval1,pval1)
      tpval2=c(tpval2,pval2)
      tTF1=c(tTF1,TF1) #p>0.1
      tTF2=c(tTF2,TF2)
    }
  }
  return(list("skewness"=tskew,"kurtosis"=tktosis,"Lilliefors.p"=tpval1, "Lilliefors.H0"=tTF1, "ShapiroWilk.p"=tpval2, "ShapiroWilk.H0"=tTF2))
}
LevTest<-function(pick){
  #Formatting the data
  area_list=c(t(pick$area[1:dim(pick$area)[1],]))
  categor=matrix(0,nrow=dim(pick$area)[1],ncol=dim(pick$area)[2])
  for (i in 1:dim(pick$area)[1]){
    categor[i,]=c(rep(i,dim(pick$area)[2]))
  }
  categor_list=factor(c(t(categor)))
  #Testing for equality of variance in "pick" population
  lev_data=data.frame("area_x"=log(area_list),"lineage"=categor_list)
  ltest=leveneTest(lm(area_x~lineage,data=lev_data),na.rm=TRUE)
  return(ltest)
}

######Data#####
datsetname="Lawless"
area=fread("Lawless_area_filtered1.txt",header=FALSE)
times=fread("Lawless_times_filtered1.txt",header=FALSE)
data=fread("Lawless_data_model3.txt",header=TRUE)
data=data.matrix(data)
residuals=data[,6:dim(data)[2]]
data=fread("Lawless_data_model3.txt",header=TRUE)

#Testing for normality 
normvals=testing_normality(residuals)
print(dim(residuals)[1])
print(length(which(normvals$Lilliefors.H0==TRUE)))
print(length(which(normvals$ShapiroWilk.H0==TRUE)))

#Storing the output in a file
op_residuals=cbind(data.frame(normvals),residuals)
writeout=paste(datsetname,"_residuals_normtests_model3.txt",sep="")
write.table(op_residuals,writeout,col.names=TRUE,row.names=FALSE)
x=fread(writeout,header=TRUE)

#Testing for normality in the residuals
#H0: the residuals are normally distributed 
normtests=fread("Lawless_residuals_normtests_model3.txt",header=TRUE) 
LF_indices=which(normtests$Lilliefors.p<0.05) #less stringent; more likely to be normal 
SW_indices=which(normtests$ShapiroWilk.p<0.05)
print(length(LF_indices))
print(length(SW_indices))
print(length(intersect(LF_indices,SW_indices)))

total_indices=unique(c(LF_indices,SW_indices))
area=area[-total_indices,]
residuals=residuals[-total_indices,]
times=times[-total_indices,]
data=data[-total_indices,]

x=read.table("ToBeRemoved_GrowthCurves.txt",header=TRUE)
print(length(intersect(x$model2,total_indices)))
indices_combined=unique(x$model2,total_indices)
print(length(indices_combined))
write.table(indices_combined,"ToBeRemoved_GrowthCurves.txt",col.names=FALSE,row.names=FALSE)

# # Filter Residuals
# # residuals must lie in the interval -50 to 50
# filter_indices=c()
# for (i in 1:dim(area)[1]){
#   if (max(area[i,],na.rm=TRUE)<1000){
#     filter_indices=c(filter_indices,i) 
#   }
# }
# 
# area=area[filter_indices,]
# residuals=residuals[filter_indices,]
# times=times[filter_indices,]
# data=data[filter_indices,]


############### Mean Residuals #################

# Analysing the mean residual fit for the exponential model
meanRes=rowSums(residuals, na.rm=TRUE)/dim(residuals)[2]
op=par(mfrow=c(1,1))
total_mean=mean(meanRes)
plot(meanRes, main=paste("Mean Residuals Values for Each Lineage, N = ",toString(dim(residuals)[1]),sep=""))
par(op)
print(mean(meanRes))


############################################################################################################

#Analysing  the intercept values of the exponential model

dat=data
int=dat$intercept

int=int[which(int<200)]

hist_cells<-function(dat,interval){
  lo=trunc(min(dat,na.rm=TRUE)*10)/10-interval #rounding down
  hi=trunc(max(dat,na.rm=TRUE)*10)/10+interval #rounding up
  cells=seq(lo,hi,interval)
  return(cells)
}

cells=hist_cells(int,1)
op=par(mfrow=c(1,1))
histo=hist(int, breaks=cells, xlab="log(Area)", main=paste(datsetname, "; Intercept Values for the Linear Regression Model", sep=""),cex.lab=1.4, cex.main=1.2,col="lightblue")
valguess=(histo$breaks[which.max(histo$counts)]+histo$breaks[which.max(histo$counts)+1])/2
abline(v=valguess,col="red",cex=1.5)
abline(v=median(int),col="green")
abline(v=mean(int),col="darkorange")
legend("topleft",legend=c("LME guess","Median","Mean"), col=c("red","green","darkorange"),pch=20)
par(op)

valguess=(histo$breaks[which.max(histo$counts)]+histo$breaks[which.max(histo$counts)+1])/2
print(paste("Estimated area of one cell is", valguess))
print(paste("Median starting area is",median(int)))
print(paste("Mean starting area is",mean(int)))


# Generating plots of well-behaved growth curves only 

# # Model 2 
# pdf(height = 16, width = 16, file = "Exponential_Model2_Fit_WellBehaved.pdf")
# par(mfrow=c(4,4))
# 
# area=data.matrix(area)
# times=data.matrix(times)
# AList=c()
# rList=c()
# residuals=matrix(NA,nrow=dim(area)[1],ncol=dim(area)[2]) #residuals maybe of different lengths; thus NA for when no data point is available
# for (i in 1:dim(area)[1]){
#   if(length(unique(area[i,]))>1){
#     inocguess=min(area[i,],na.rm=TRUE) #minimum data point;
#     # Getting rid of NA values
#     indices=as.numeric(which(area[i,]>0))
#     i_area=as.numeric(area[i,indices])
#     i_times=as.numeric(times[i,indices])
#     # Fit exponential model
#     result=nls(y~A+r*x,data=data.frame(x=i_times,y=log(i_area)),start=list(A=log(inocguess),r=0.1))
#     resid=residuals(result)
#     res=result$m$getPars()
#     A=exp(as.numeric(res[1])) #converting back from logs
#     r=(as.numeric(res[2])) #this still includes negative growth rates which have no biological meaning
#     plot(i_times,predict(result),type='l',main=paste("Growth Curve:", i),ylab="log(Area)",xlab="Time",cex.main=1.4,cex.lab=1.4,log='y')
#     lines(i_times,log(i_area),type='p')
#     plot(resid,type='p',main=paste("Residuals for Growth Curve",i))
#   }else{A=1;r=0;residuals=rep(NA,dim(area)[2])}
#   AList=c(AList,A)
#   rList=c(rList,r)
#   residuals[i,indices]=resid
# }
# 
# dev.off()



# Model 3
pdf(height = 16, width = 16, file = "Exponential_Model3_Fit_WellBehaved.pdf")
par(mfrow=c(4,4))

area=data.matrix(area)
times=data.matrix(times)
AList=c()
rList=c()
residuals=matrix(NA,nrow=dim(area)[1],ncol=dim(area)[2]) #residuals maybe of different lengths; thus NA for when no data point is available
for (i in 1:dim(area)[1]){
  print(i)
  if(length(unique(area[i,]))>1){
    inocguess=min(area[i,],na.rm=TRUE) #minimum data point;
    # Getting rid of NA values
    indices=as.numeric(which(area[i,]>0))
    i_area=as.numeric(area[i,indices])
    i_times=as.numeric(times[i,indices])
    # Fit exponential model
    result=nls(y~A*exp(r*x),data=data.frame(x=i_times,y=i_area),start=list(A=inocguess,r=0.1))
    resid=residuals(result)
    res=result$m$getPars()
    A=(as.numeric(res[1]))
    r=(as.numeric(res[2])) #this still includes negative growth rates which have no biological meaning
    # Producing the plot
    plot(i_times,predict(result),type='l',main=paste("Growth Curve:", i),ylab="Area",xlab="Time",cex.main=1.4,cex.lab=1.4)
    lines(i_times,i_area,type='p')
    plot(resid,type='p',main=paste("Residuals for Growth Curve",i))
  }else{A=1;r=0;residuals=rep(NA,dim(area)[2])}
  AList=c(AList,A)
  rList=c(rList,r)
  residuals[i,indices]=resid
}

dev.off()

