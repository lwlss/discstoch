require(segmented)
strains=read.delim("Report.txt",header=TRUE,stringsAsFactors=FALSE)
snames=unique(strains$Culture)
snames=snames[order(snames)]
start=min(strains$Time) # Should replace this with actual spotting time...
strains$Time=(strains$Time-start)/(60*60)
#areamax=max(strains$TotalArea)
#areamin=min(strains$TotalArea)
areamin=2e3
areamax=6e5

tmax=20
tmin=10

brokenStick=function(x,init,s1,brk,s2){
	if(x<=brk){
		return(init+x*s1)
	}else{
		return(brokenStick(brk,init,s1,brk,s2)+(x-brk)*s2)
	}
}

bStick=function(x,init,s1,brk,s2) sapply(x,brokenStick,init,s1,brk,s2)

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

strains$OrigCulture=pinid(strains$Culture)
strains$Row=as.numeric(substr(strains$OrigCulture,2,3))
strains$Col=as.numeric(substr(strains$OrigCulture,5,6))

getORF=makeORFGetter("AUXILIARY\\LibraryDescription.txt")
strains$ORF=""
strains$Gene=""
for (i in 1:length(strains[,1])){
	strains$ORF[i]=getORF(strains$Row[i],strains$Col[i])
	strains$Gene[i]=strsplit(strains$ORF[i],"_")[[1]][1]
}

ORFs=c();Genes=c();Inoc=c();EarlyDT=c();LateDT=c();Break=c()

#pdf("microQFAPlots.pdf")

#svg(filename = "microQFAPlates%03d.svg",
#    width = 7, height = 7, pointsize = 12,
#    onefile = FALSE, family = "sans", bg = "white",
#    antialias = c("default", "none", "gray", "subpixel"))

# Fit exponential growth model to the data
#op<-par(mfrow=c(8,8),oma=c(4,4,4,4),mar=c(1.5,0,1.5,0),mgp=c(3,1,0),cex=0.7)
mkplot=c("R03C05","R03C06","R03C07","R03C08","R03C09","R03C10")
for (strain in snames){
	curr=strains[(strains$Culture==strain)&(strains$Time<tmax),]
	ORF=curr$ORF[1]
	Gene=curr$Gene[1]
	curr$LogArea=log(curr$TotalArea)
	try({
	res=lm(LogArea~Time,data=curr)
	seg=segmented.lm(res,seg.Z=~Time,psi=list(Time=tmax/2),control=seg.control(display=FALSE))
	init=as.numeric(seg$coefficients[1])
	s1=as.numeric(slope(seg)$Time[1,1])
	brk=as.numeric(seg$psi[2])
	s2=as.numeric(slope(seg)$Time[2,1])
	dt1=log(2)/s1
	dt2=log(2)/s2
	text(2,5e5,paste("Early DT:\n",signif(dt1,4),"(h)"),cex=1)
	text(17,3e3,paste("Late DT:\n",signif(dt2,4),"(h)"),cex=1)
	ORFs=c(ORFs,ORF)
	Genes=c(Genes,Gene)
	Inoc=c(Inoc,exp(init))
	EarlyDT=c(EarlyDT,dt1)
	LateDT=c(LateDT,dt2)
	Break=c(Break,brk)
	})
	if( strain%in%mkplot){
		plot(curr$Time,curr$TotalArea,type="n",main=ORF,xlab="Time since inoculation (h)",ylab="Log Culture Area (px)",xlim=c(0,tmax),
			ylim=c(areamin,areamax),log="y",axes=TRUE,cex.main=1.5,cex.lab=1.5)
		#axis(1,cex.lab=0.5,tck=-.1)
		abline(v=brk,col="blue",lwd=2)
		curve(exp(bStick(x,init,s1,brk,s2)),from=0,to=tmax,col="red",add=TRUE,lwd=3,log="y")
		points(curr$Time,curr$TotalArea,pch=16,cex=0.5)
	}
}
#title(main="microQFA Growth Curves",xlab="Time since inoculation (h)",line=1.5,ylab="Cell Density (AU)",cex.main=1,cex.lab=1,outer=TRUE)
#par(op)

data=data.frame(ORFs,Genes,Inoc,Break,EarlyDT,LateDT)
tst=wilcox.test(data$EarlyDT[data$Genes=="HIS3"],data$EarlyDT[data$Genes=="HTZ1"],conf.int=TRUE)
mlab=paste("Early (lag phase)  p-val:",signif(tst$p.value,3),"Diff:",signif(abs(tst$estimate),3))
plot(density(data$EarlyDT[data$Genes=="HIS3"]),xlim=c(2,8),ylim=c(0,0.7),xlab="Doubling Time (h)",main=mlab,lwd=3)
points(density(data$EarlyDT[data$Genes=="HTZ1"]),col="red",type="l",lwd=3)
legend("topleft",c("HIS3","HTZ1"),col=c("black","red"),lwd=3)

tst=wilcox.test(data$LateDT[data$Genes=="HIS3"],data$LateDT[data$Genes=="HTZ1"],conf.int=TRUE)
mlab=paste("Late (exponential phase)  p-val:",signif(tst$p.value,3),"Diff:",signif(abs(tst$estimate),3))
plot(density(data$LateDT[data$Genes=="HIS3"]),xlim=c(2.5,3.5),ylim=c(0,6),xlab="Doubling Time (h)",main=mlab,lwd=3)
points(density(data$LateDT[data$Genes=="HTZ1"]),col="red",type="l",lwd=3)
legend("topleft",c("HIS3","HTZ1"),col=c("black","red"),lwd=3)

tst=wilcox.test(data$Break[data$Genes=="HIS3"],data$Break[data$Genes=="HTZ1"],conf.int=TRUE)
mlab=paste("p-val:",signif(tst$p.value,3),"Diff:",signif(abs(tst$estimate),3))
plot(density(data$Break[data$Genes=="HIS3"]),xlim=c(0,15),ylim=c(0,0.9),xlab="Length of lag phase (h)",main=mlab,lwd=3)
points(density(data$Break[data$Genes=="HTZ1"]),col="red",type="l",lwd=3)
legend("topleft",c("HIS3","HTZ1"),col=c("black","red"),lwd=3)

tst=wilcox.test(data$Inoc[data$Genes=="HIS3"],data$Inoc[data$Genes=="HTZ1"],conf.int=TRUE)
mlab=paste("Inoculum density  p-val:",signif(tst$p.value,3),"Diff:",signif(abs(tst$estimate),3))
plot(density(data$Inoc[data$Genes=="HIS3"]),xlab="Inoculum area (px)",xlim=c(2000,25000),main=mlab,lwd=3)
points(density(data$Inoc[data$Genes=="HTZ1"]),col="red",type="l",lwd=3)
legend("topleft",c("HIS3","HTZ1"),col=c("black","red"),lwd=3)

dev.off()

write.table(data,"Phenotypes.txt",row.names=FALSE,quote=FALSE,sep="\t")



