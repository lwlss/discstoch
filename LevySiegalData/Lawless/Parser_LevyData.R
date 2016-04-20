#PARSER FOR LEVY DATA 
#To bring all three data sets into the same format 
#so that a single script can be used to extract growth curves 

load("d.110601.Rfile")
rep1=d
load("d.110612.Rfile")
rep3=d
load("d.110607.Rfile") 
rep2=d
load("d.110625.Rfile")
rep5=d
load("d.110706.Rfile") 
rep4=d

replicates=list(rep1,rep2,rep3,rep4,rep5)

for (r in 1:length(replicates)){
  rep=replicates[[r]]
  strain_names=rep$condition.names
  area=rep$areas
  times=rep$times
  identifier=rep(r,dim(area)[1])
  genotype=rep(0,dim(area)[1])
  clonalcolony=rep(0,dim(area)[1])
  for (i in 1:length(strain_names)){ #96 wells; find area indices for all of them
    indices=rep$well.list[[i]]
    genotype[indices]=rep$condition.names[[i]]
    clonalcolony[indices]=paste(toString(r),"_",toString(i), sep="")
  }
  SubData=data.frame(genotype,clonalcolony,identifier)
  write.table(area,"Levy_area.txt",append=TRUE,col.names=FALSE,row.names=FALSE)
  write.table(times,"Levy_time.txt",append=TRUE,col.names=FALSE,row.names=FALSE)
  write.table(SubData,"Levy_data.txt",append=TRUE,col.names=FALSE,row.names=FALSE)
}