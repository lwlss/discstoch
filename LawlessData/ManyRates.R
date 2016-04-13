

# Convert microQFA pinid to QFA pinid (plate is upside down)
pinid=function(pid){
	Col=as.numeric(substr(pid,5,6))
	NewCol=24-(Col-1)
	return(paste(substr(pid,1,4),sprintf("%02d",NewCol),sep=""))
}

# Read in library description file describing which strain at each spot
makeORFGetter=function(libraries,lib="LiqCurve",plate=1){
	# Open the library descriptions
	libs=read.delim(libraries,sep="\t",header=TRUE,stringsAsFactors=FALSE)
	libs$ORF=toupper(libs$ORF)
	# What library are we using? e.g. library="AndrewsOE_v2" library="SDL_v2"
	liblist=unique(libs$Library)
	NLIB=length(liblist)
	libdict=1:NLIB
	names(libdict)=liblist
	#libs=libs[libs$Library==library,]
	# Make an array for storing info about spots
	NROW=max(libs$Row)
	NCOL=max(libs$Column)
	NPLATE=max(libs$Plate)
	SPOTARRAY=array("missing",dim=c(NLIB,NPLATE,NROW,NCOL))
	# Fill the spot array object
	for (x in 1:length(libs$Plate)) SPOTARRAY[libdict[[libs[x,"Library"]]],libs[x,"Plate"],libs[x,"Row"],libs[x,"Column"]]=libs[x,"ORF"]
	getORF<-function(row,col) SPOTARRAY[libdict[[lib]],plate,row,col]
}


plotStuff=function(fname){
	strname=pinid(substr(fname,1,6))
	row=as.numeric(substr(strname,2,3))
	col=as.numeric(substr(strname,5,6))
	strgene=getORF(row,col)
	print(strgene)
	microdat=read.delim(fname,header=FALSE,sep="\t")
	numcolonies=dim(microdat)[2]

	cramp=colorRamp(c("grey","black","grey"),space="Lab")

	tim=(microdat[,1]-microdat[1,1])/(3600)
	#op=par(mfrow=c(1,2))
	plot(NULL,type="l",xlab="Time (h)",ylab="Microcolony Area (px)",ylim=c(0,3000),xlim=c(0,10),main=strgene,cex.lab=1.25)
	AList=c()
	rList=c()
	for(x in 2:numcolonies){
		dat=microdat[,x]
		dat=sapply(dat,max,1)
		if(length(unique(dat))>1){
			inocguess=max(min(dat),1)
			# Fit exponential model
			result=nls(y~A+r*x,data=data.frame(x=tim,y=log(dat)),start=list(A=log(inocguess),r=0.1))
			res=result$m$getPars()
			A=exp(as.numeric(res[1]))
			r=max(0,as.numeric(res[2]))
		}else{A=1;r=0}
		AList=c(AList,A)
		rList=c(rList,r)
		lines(tim,dat,type="l",lwd=0.2)#,col=rgb(cramp(min(1,r/0.3))/255.0))
	}
	brks=(0:500)/(500*0.3)
	hist(rList,breaks=brks,xlab="r (1/h)",main=paste(numcolonies,"microcolonies"),xlim=c(0,0.3),ylim=c(0,12.5),cex.lab=1.25)
	#par(op)

	sorter=order(rList,decreasing=TRUE)
	rsort=rList[sorter][1:(length(sorter)*0.10)]
	QFAsim=mean(rsort)
	
	return(data.frame(Pin=strname,A=AList,r=rList,QFA=rep(QFAsim,length(AList)),Strain=rep(strgene,length(AList))))
}


plotBoth=function(fnameA,fnameB){
	strnameA=pinid(substr(fnameA,1,6))
	row=as.numeric(substr(strnameA,2,3))
	col=as.numeric(substr(strnameA,5,6))
	strgeneA=getORF(row,col)
	print(strgeneA)
	microdatA=read.delim(fnameA,header=FALSE,sep="\t")
	numcoloniesA=dim(microdatA)[2]

	strnameB=pinid(substr(fnameB,1,6))
	row=as.numeric(substr(strnameB,2,3))
	col=as.numeric(substr(strnameB,5,6))
	strgeneB=getORF(row,col)
	print(strgeneB)
	microdatB=read.delim(fnameB,header=FALSE,sep="\t")
	numcoloniesB=dim(microdatB)[2]

	cramp=colorRamp(c("grey","black","grey"),space="Lab")
	tcol=c(rgb(1,0,0,0.3),rgb(0,1,0,0.3),rgb(1,1,0,0.3),rgb(0,0,1,0.3))
	ccol=c(rgb(1,0,0,1.0),rgb(0,1,0,1.0),rgb(1,1,0,1.0),rgb(0,0,1,1.0))

	timA=(microdatA[,1]-microdatA[1,1])/(3600)
	timB=(microdatB[,1]-microdatB[1,1])/(3600)
	op=par(mfrow=c(1,2))
	plot(NULL,type="l",xlab="Time (h)",ylab="Microcolony Area (px)",ylim=c(0,3000),xlim=c(0,10),cex.lab=1.25,cex.axis=0.85)
	rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "grey")
	rListA=c()
	for(x in 2:numcoloniesA){
		dat=microdatA[,x]
		dat=sapply(dat,max,1)
		if(length(unique(dat))>1){
			inocguess=max(min(dat),1)
			# Fit exponential model
			result=nls(y~A+r*x,data=data.frame(x=timA,y=log(dat)),start=list(A=log(inocguess),r=0.1))
			res=result$m$getPars()
			A=exp(as.numeric(res[1]))
			r=max(0,as.numeric(res[2]))
		}else{A=1;r=0}
		rListA=c(rListA,r)
		lines(timA,dat,type="l",lwd=1,col=tcol[1])
	}

	rListB=c()
	for(x in 2:numcoloniesB){
		dat=microdatB[,x]
		dat=sapply(dat,max,1)
		if(length(unique(dat))>1){
			inocguess=max(min(dat),1)
			# Fit exponential model
			result=nls(y~A+r*x,data=data.frame(x=timB,y=log(dat)),start=list(A=log(inocguess),r=0.1))
			res=result$m$getPars()
			A=exp(as.numeric(res[1]))
			r=max(0,as.numeric(res[2]))
		}else{A=1;r=0}
		rListB=c(rListB,r)
		lines(timB,dat,type="l",lwd=1,col=tcol[4])
	}



	brks=(0:500)/(500*0.3)
	plot(NULL,ylim=c(0,12),xlim=c(0,0.40),xlab="r (1/h)",ylab="Frequency",main="",cex.lab=1.5)
	rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "grey")
	hist(rListA,breaks=brks,xlab="r (1/h)",xlim=c(0,0.3),ylim=c(0,12.5),cex.lab=1.25,col=tcol[1],border=ccol[1],add=TRUE)
	hist(rListB,breaks=brks,xlab="r (1/h)",xlim=c(0,0.3),ylim=c(0,12.5),cex.lab=1.25,col=tcol[4],border=ccol[4],add=TRUE)

	legend("topright",c(paste(strgeneA,"(N =",numcoloniesA,")"),paste(strgeneB,"(N =",numcoloniesB,")")),col=c(tcol[1],tcol[4]),pch=15,pt.lwd=1,bg="white",cex=0.85)
	par(op)

	return(wilcox.test(rListA,rListB))
}



# Fit exponential model to all microcolonies
flist=list.files(pattern="*_OUT.txt")
dat=plotStuff(flist[1])
for (f in flist[2:length(flist)]){
	params=plotStuff(f)
	dat=rbind(dat,params)
}
dat$Gene=substr(dat$Strain,1,4)

pdf("ColourVersion.pdf",width=10,height=5)
plotBoth(flist[29],flist[31])
dev.off()

getORF=makeORFGetter("AUXILIARY\\LibraryDescription.txt")

svg("CompareAllMicrocolonies.pdf")

tcol=c(rgb(1,0,0,0.3),rgb(0,1,0,0.3),rgb(1,1,0,0.3),rgb(0,0,1,0.3))
ccol=c(rgb(1,0,0,1.0),rgb(0,1,0,1.0),rgb(1,1,0,1.0),rgb(0,0,1,1.0))
nbreaks=75
rmax=max(dat$r,na.rm=TRUE)
brks=rmax*(0:nbreaks)/nbreaks

NH3=length(dat$Gene[dat$Gene=="HIS3"])
NZ1=length(dat$Gene[dat$Gene=="HTZ1"])

dat$DT=log(2)/dat$r

plot(NULL,ylim=c(0,405),xlim=c(0,0.40),xlab="r (1/h)",ylab="Frequency",main="",cex.lab=1.5)
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "grey")
hist(dat$r[dat$Gene=="HIS3"],breaks=brks,add=TRUE,col=tcol[1],border=ccol[1],lwd=3)
hist(dat$r[dat$Gene=="HTZ1"],breaks=brks,add=TRUE,col=tcol[4],border=ccol[4],lwd=3)
legend("topright",c(paste("his3 (N =",NH3,")"),paste("htz1 (N =",NZ1,")")),col=c(tcol[1],tcol[4]),pch=15,pt.lwd=1,bg="white")

dev.off()

pdf("CompareTwoStrains.pdf")
op=par(mfrow=c(2,2),mai=c(0.65,0.65,0.4,0.2))
tt=plotStuff(flist[29])
tt=plotStuff(flist[31])
par(op)
dev.off()


op=par(ask=TRUE,mfrow=c(1,2))
for(x in 2:174){
	dat=microdat[,x]
	A=AList[x-1]
	r=rList[x-1]
	plot(tim,dat,ylim=c(0,3000),xlim=c(0,11))
	curve(A*exp(r*x),from=0,to=11,col="red",add=TRUE)
	plot(tim,log(dat),ylim=c(log(1),log(3000)),xlim=c(0,11))
	curve(log(A)+r*x,from=0,to=11,col="red",add=TRUE)
}

par(op)