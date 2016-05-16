# Image Analysis 
library(gnumeric)
categories=read.gnumeric.sheet("HelenaAnalysis_ManualGrowthCurveCheck.ods",head=TRUE)

data=categories[2:7]

for (i in 1:length(data)){
  print(names(data)[i])
  dirt=which(data[,i]==0)
  print(paste("contains", length(dirt), "blobs of dirt"))
  side=which(data[,i]==2)
  print(paste("contains", length(side), "blobs which are cut-off at the side"))
  side0=which(data[,i]==20)
  print(paste("contains", length(side0), "blobs which are cut-off at the side and which exhibit 0 growth"))
  good=which(data[,i]==1)
  print(paste("contains", length(good), "blobs which look well-behaved!!!"))
  good0=which(data[,i]==10)
  print(paste("contains", length(good0), "blobs which look well-behaved with zero growth!!!"))
  multiple=which(data[,i]==4)
  print(paste("contains", length(multiple), "blobs which consist of multiple blobs"))
  print(" ")
}