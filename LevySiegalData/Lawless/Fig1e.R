## Figure 1e
source("plot.cum.histo2.R")
source("get.cdata.R")
#install.packages(c("psych","lawstat"))

##############################################################
load("d.110601.Rfile")
load("gr.110601.Rfile")
load("md.110601.Rfile")
gr$colparam[which(gr$colparam[,1] < 0),] <- 0
cf <- md$cenfield
md <- md$mindist
x <- md; x[] <- 0
x[which(md < 35)] <- 1
x[which(cf == 0)] <- 1
gr$colparam[which(x == 1),1] <- NA


d.110601 <- d
gr.110601 <- gr

if(!file.exists("cdata.110601.Rfile")){
  cdata <- get.cdata(d,gr)
  cdata.110601 <- cdata
  save(cdata.110601, file = "cdata.110601.Rfile")
}else{
  load("cdata.110601.Rfile")
}

pID <- pID.110601 <- d$pID
cn <- cn.110601 <- d$condition.names
cnu <- cnu.110601 <- unique(cn)

##############################################################
load("d.110612.Rfile")
load("gr.110612.Rfile")
load("md.110612.Rfile")
#load("cdata.110612.Rfile")
gr$colparam[which(gr$colparam[,1] < 0),] <- 0
cf <- md$cenfield
md <- md$mindist
x <- md; x[] <- 0
x[which(md < 35)] <- 1
x[which(cf == 0)] <- 1
gr$colparam[which(x == 1),1] <- NA


d.110612 <- d
gr.110612 <- gr

if(!file.exists("cdata.110612.Rfile")){
  cdata <- get.cdata(d,gr)
  cdata.110612 <- cdata
  save(cdata.110612, file = "cdata.110612.Rfile")
}else{
  load("cdata.110612.Rfile")
}

pID <- pID.110612 <- d$pID
cn <- cn.110612 <- d$condition.names
cnu <- cnu.110612 <- unique(cn)

##############################################################
load("d.110607.Rfile")
load("gr.110607.Rfile")
load("md.110607.Rfile")
gr$colparam[which(gr$colparam[,1] < 0),] <- 0
cf <- md$cenfield
md <- md$mindist
x <- md; x[] <- 0
x[which(md < 35)] <- 1
x[which(cf == 0)] <- 1
gr$colparam[which(x == 1),1] <- NA


d.110607 <- d
gr.110607 <- gr

if(!file.exists("cdata.110607.Rfile")){
  cdata <- get.cdata(d,gr)
  cdata.110607 <- cdata
  save(cdata.110607, file = "cdata.110607.Rfile")
}else{
  load("cdata.110607.Rfile")
}

  
pID <- pID.110607 <- d$pID
cn <- cn.110607 <- d$condition.names
cnu <- cnu.110607 <- unique(cn)



##############################################################
load("d.110625.Rfile")
load("gr.110625.Rfile")
load("md.110625.Rfile")
gr$colparam[which(gr$colparam[,1] < 0),] <- 0
cf <- md$cenfield
md <- md$mindist
x <- md; x[] <- 0
x[which(md < 35)] <- 1
x[which(cf == 0)] <- 1
gr$colparam[which(x == 1),1] <- NA


d.110625 <- d
gr.110625 <- gr

if(!file.exists("cdata.110625.Rfile")){
  cdata <- get.cdata(d,gr)
  cdata.110625 <- cdata
  save(cdata.110625, file = "cdata.110625.Rfile")
}else{
  load("cdata.110625.Rfile")
}

pID <- pID.110625 <- d$pID
cn <- cn.110625 <- d$condition.names
cnu <- cnu.110625 <- unique(cn)


##############################################################
load("d.110706.Rfile")
load("gr.110706.Rfile")
load("md.110706.Rfile")
gr$colparam[which(gr$colparam[,1] < 0),] <- 0
cf <- md$cenfield
md <- md$mindist
x <- md; x[] <- 0
x[which(md < 35)] <- 1
x[which(cf == 0)] <- 1
gr$colparam[which(x == 1),1] <- NA


d.110706 <- d
gr.110706 <- gr

if(!file.exists("cdata.110706.Rfile")){
  cdata <- get.cdata(d,gr)
  cdata.110706 <- cdata
  save(cdata.110706, file = "cdata.110706.Rfile")
}else{
  load("cdata.110706.Rfile")
}
  
pID <- pID.110706 <- d$pID
cn <- cn.110706 <- d$condition.names
cnu <- cnu.110706 <- unique(cn)


##############################################################
#make a combine well list (combining all replicates)

#rep1
cn <- cn.110601
cnu <- cnu.110601
cdata <- cdata.110601
wl <- list()
for(j in 1:length(cnu)){

  twl <- names(cn)[which(cn == cnu[j])]
  w <- cdata$test.norm[[ names(cdata$condition.names[twl[[1]]]) ]]
  for(i in 2:length(twl)){
    w <- c(w, cdata$test.norm[[ names(cdata$condition.names[twl[[i]]]) ]])
  }
  #remove outliers
   w <- w[which(w < 1.5)]
  wl[[j]] <- w
  names(wl)[j] <- cnu[j]
}
pwl <- pwl.110601 <- wl

#rep2
cn <- cn.110607
cnu <- cnu.110607
cdata <- cdata.110607
wl <- list()
for(j in 1:length(cnu)){

  twl <- names(cn)[which(cn == cnu[j])]
  w <- cdata$test.norm[[ names(cdata$condition.names[twl[[1]]]) ]]
  for(i in 2:length(twl)){
    w <- c(w, cdata$test.norm[[ names(cdata$condition.names[twl[[i]]]) ]])
  }
  #remove outliers
  w <- w[which(w < 1.5)]
  wl[[j]] <- w
  names(wl)[j] <- cnu[j]
}
pwl <- pwl.110607 <- wl


#rep3
cn <- cn.110612
cnu <- cnu.110612
cdata <- cdata.110612
wl <- list()
for(j in 1:length(cnu)){

  twl <- names(cn)[which(cn == cnu[j])]
  w <- cdata$test.norm[[ names(cdata$condition.names[twl[[1]]]) ]]
  for(i in 2:length(twl)){
    w <- c(w, cdata$test.norm[[ names(cdata$condition.names[twl[[i]]]) ]])
  }
  #remove outliers
  w <- w[which(w < 1.5)]
  wl[[j]] <- w
  names(wl)[j] <- cnu[j]
}
pwl <- pwl.110612 <- wl

#rep4
cn <- cn.110706
cnu <- cnu.110706
cdata <- cdata.110706
wl <- list()
for(j in 1:length(cnu)){

  twl <- names(cn)[which(cn == cnu[j])]
  wlt <- cdata$test.norm[twl]
  wl2 <- wlt
  for(i in 1:length(wlt)){
    if (length(wlt[[i]]) < 200){
      wl2[[i]] <- NULL
    }
  }
  twl <- names(wl2)
  w <- cdata$test.norm[[ names(cdata$condition.names[twl[[1]]]) ]]
  for(i in 2:length(twl)){
    w <- c(w, cdata$test.norm[[ names(cdata$condition.names[twl[[i]]]) ]])
  }
  #remove outliers
  w <- w[which(w < 1.5)]
  wl[[j]] <- w
  names(wl)[j] <- cnu[j]
}
pwl <- wl
names(pwl)[5] <- "PRS3"
pwl.110706 <- pwl




pdf(height = 16, width = 16, file = "ko1replicates.pdf")
par(mfrow = c(4,4))
for(j in 1:length(cnu)){
  twl <- names(cn)[which(cn == cnu[j])]
  wl <- cdata$test.norm[twl]
  wl2 <- wl
  for(i in 1:length(wl)){
    if (length(wl[[i]]) < 200){
      wl2[[i]] <- NULL
    }
  } 
  plot.cum.histo2(wl2, xlim = c(0.2, 1.6), text.x = 0.5,  main = cnu[j])
}
dev.off()


#nicely replicating
#KEM1 APQ12 YHL100C BEM1 RTT109 SCP160 SWA2 DIA2
#YHR095W

wl1 <- wl <- pwl.110601[c("YFR054C", "YHR095W", "HTZ1", "RAD50", "SNF6", "PET9", "YME1")]
wl2 <- wl<- pwl.110607[c("YFR054C", "YHR095W", "HTZ1", "RAD50", "SNF6", "PET9", "YME1")]
wl3 <- wl <- pwl.110612[c("YFR054C", "YHR095W", "HTZ1", "RAD50", "SNF6", "PET9", "YME1")]
wl4 <- pwl.110706[c("KEM1", "APQ12", "PRS3", "BEM1", "RTT109", "SCP160", "SWA2", "DIA2", "NOT5")]
wl4 <- pwl.110706[c("KEM1", "BEM1", "SCP160", "SWA2", "DIA2", "NOT5")]
bwl <- wl1
for (i in 1:length(bwl)){
  bwl[[i]] <- c(wl1[[i]], wl2[[i]], wl3[[i]])
}
bwl <- c(bwl, wl4)

#standardize by dubious orf
for( i in 1:length(bwl)){
  bwl[[i]] <- bwl[[i]]/((1.022 + 1.008)/2)
}

###
pdf(height = 5.25, width = 5.25, file = "fig1d.pdf")
mycol <- colorRampPalette(c(colors()[258], "gold", "blue3", "grey", colors()[520], "dark magenta"))
plot.cum.histo2(bwl[c("SWA2", "DIA2", "KEM1", "SNF6", "RAD50", "HTZ1", "SCP160", "YME1", "BEM1", "YHR095W", "YFR054C","PET9","NOT5")] , xlim = c(0, 1.6), text.x = 0.5, main = "", print.legend = F, colors = mycol(13))
legend(0,107, cex = .7, col = mycol(13), lwd=3, legend =  c(expression(paste(Delta , " SWA2")),expression(paste(Delta , " DIA2")), expression(paste(Delta , " KEM1")), expression(paste(Delta , " SNF6")), expression(paste(Delta , " RAD50")), expression(paste(Delta , " HTZ1")), expression(paste(Delta , " SCP160")), expression(paste(Delta , " YME1")), expression(paste(Delta , " BEM1")), expression(paste(Delta , " YHR095W")), expression(paste(Delta , " YFR054C")),expression(paste(Delta , " PET9")),expression(paste(Delta , " NOT5"))))
dev.off()


#Table S1
######################################################################
l.test3 <- function (a, b){
  library(lawstat)
  y <- c(a, b)
  z <- y
  z[1:length(a)] <- 1
  z[(length(a) + 1):length(z)] <- 2
  m <- levene.test(y, z, option="median", trim.alpha = 0)
  return(m)
}
library(psych)
######################################################################

x <- bwl[c("SWA2", "DIA2", "KEM1", "SNF6", "RAD50", "HTZ1", "SCP160", "YME1", "BEM1", "YHR095W", "YFR054C","PET9","NOT5")]
y <- matrix(NA, 13, 7)
rownames(y) <- c("SWA2", "DIA2", "KEM1", "SNF6", "RAD50", "HTZ1", "SCP160", "YME1", "BEM1", "YHR095W", "YFR054C","PET9","NOT5")
colnames(y) <- c("Mean", "Standard Deviation", "n", "Equality of Means*", "Equality of Variances**", "Percent Slow Growing", "Geometric Mean")
for(i in 1:length(x)){
  x[[i]] <- x[[i]][which(!x[[i]] == 0)]
  y[i,3] <- length(x[[i]])
  y[i,1] <- round(mean(x[[i]]), digits = 3)
  y[i,2] <- round(sd(x[[i]]), digits = 3)
  y[i,4] <- round(log(wilcox.test(x[[i]], x$YFR054C)[[3]], base = 10), digits = 1)#equality of means
  y[i,5] <- round(log(l.test3(x[[i]], x$YFR054C)[[2]], base = 10), digits = 1) #equality of variances
  y[i,6] <- round(length(which(x[[i]] < (median(x[[i]])/2)))/length(x[[i]])*100, digits = 3)
  y[i,7] <- round(geometric.mean(x[[i]]), digits = 3)
}
write.table(y, sep = "\t", file = "tableS1.txt", quote = F)


#levene test for equality of variances



#plotting multiple replicates on one plot
cn1 <- cn.110625
cnu1 <- cnu.110625
cdata1 <- cdata.110625
cn2 <- cn.110607
cnu2 <- cnu.110607
cdata2 <- cdata.110607
cn3 <- cn.110612
cnu3 <- cnu.110612
cdata3 <- cdata.110612
wl <- list()

pdf(height = 8, width = 12, file = "KO.Well.replicates.pdf")
par(mfrow = c(2,3))
u <- c("YFR054C", "PET9", "YME1", "HTZ1", "RAD50", "SNF6")
for(j in 1:length(u)){
  twl <- names(cn1)[which(cn1 == u[j])]
  w1 <- cdata1$test.norm[twl]
  twl <- names(cn2)[which(cn2 == u[j])]
  w2 <- cdata2$test.norm[twl]
  twl <- names(cn3)[which(cn3 == u[j])]
  w3 <- cdata3$test.norm[twl]
  w <- c(w1, w2, w3)
  

  for (i in 1:length(w)){
    w[[i]] <-  w[[i]][which(w[[i]] < 1.4)]
    w[[i]] <- w[[i]]/((1.038 + 1.046)/2) #standardize by dubious orf
  }
  wl[[j]] <- w
  names(wl)[j] <- u[j]
  mycol <- colors()[c(488,488,488,488,488,488,81,81,81,81,81,81,556,556,556,556,556,556)]
  plot.cum.histo2(wl[[j]], xlim = c(0, 1.5), text.x = 0.5, main = names(wl)[j], print.legend = F, colors = mycol, LWD = 1)
}
dev.off()



# Growth curves
x=d$areas[d$well.list$B3,]
t=d$times[d$well.list$B3,]

index=1751
plot(t[index,],x[index,],xlab="Time (h)",ylab="Area (px)",log="y")

