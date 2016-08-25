# Generates six synthetic lineage growth rate distribution for analysis

library(detstocgrowth)
#library(stats)
library(LambertW)
#library(bcp)
#library(segmented)

# Synthetic Distributions
set.seed(1)
# Sharp Peak
dist1=rnorm(1000000,0.3,0.015)
# Wide Peak
dist2=rnorm(1000000,0.3,0.05)
# Long right-hand tail
dist3=rLambertW(n=1000000, distname = "t", beta = c(0.2,0.035,100), delta = 0.08)
dist3[which(dist3<0)]=rnorm(1,0.2,0.05)
unimodals=list(dist3,dist1,dist2)
# Even Bimodal
dist4=c(rnorm(500000,0.2,0.03),rnorm(500000,0.4,0.03))
# Decreasing Bimodal
dist5=c(rnorm(800000,0.2,0.03),rnorm(200000,0.4,0.03))
# Increasing Bimodal
dist6=c(rnorm(200000,0.2,0.03),rnorm(800000,0.4,0.03))
bimodals=list(dist4,dist5,dist6)

# Plotting the Distributions
par(mfrow=c(1,2))
colours=c("red","blue","green")
colours=c("green","green3","darkgreen")
plot(NULL,xlim=c(0,0.8),ylim=c(0,30),
     xlab="Growth Rate (1/h)", ylab="Density", cex.lab=1.4,
     main="Unimodal Growth Rate Distributions")
for (i in 1:length(unimodals)){
  strain_rates=unimodals[[i]]
  lines(density(strain_rates,from=0),col=colours[i],lwd=4)
  polygon(density(strain_rates)$x,density(strain_rates)$y,col=adjustcolor(colours[i],0.2))
}
legend("topright",col=colours,legend=c("Sharp Peak","Moderate Peak","Long Tail"),lty=c(1,1,1),lwd=c(2,2,2),cex=0.85)

colours=c("red","blue","green")
plot(NULL,xlim=c(0,0.8),ylim=c(0,15),
     xlab="Growth Rate (1/h)", ylab="Density", cex.lab=1.4,
     main="Bimodal Growth Rate Distributions")
for (i in 1:length(bimodals)){
  strain_rates=bimodals[[i]]
  lines(density(strain_rates,from=0),col=colours[i],lwd=4)
  polygon(density(strain_rates)$x,density(strain_rates)$y,col=adjustcolor(colours[i],0.2))
}
legend("topright",col=colours,legend=c("Even Modality","Decreasing Modality","Increasing Modality"),lty=c(1,1,1),lwd=c(2,2,2),cex=0.85)

# # Plotting mini distributions to insert into the simulation plots
# allcolours=c("green","green3","darkgreen","cyan","cornflowerblue","navy")
# par(mfrow=c(2,3))
# for (i in 1:length(alldist)){
#   plot(density(alldist[[i]],from=0),col=allcolours[i],lwd=4,lty=1,xlab=" ",ylab=" ",main=" ")
#   polygon(density(alldist[[i]],from=0)$x,density(alldist[[i]],from=0)$y,col=adjustcolor(allcolours[i],0.2))
# }

# Generating and plotting the population simulations
alldist=list(dist1,dist2,dist3,dist4,dist5,dist6)
N=100
iter=100
time=seq(0,96,1)
par(mfrow=c(2,3))
for (i in 1:length(alldist)){
  r=alldist[[i]]
  syntheticplots=plot_pop_sim_dat(strain=0,r,time,iter,yl=c(N,10^25),
                   col=(colorRampPalette(c("lightyellow","yellow","orange","red"))(n = 181)),
                   colrange=seq(0,0.9,0.005))
}

# Getting Lag Phase Duration Estimates for different N
N=c(seq(100,10000,100))
#N=c(100,1000,10000)
time=seq(0,96,1)
iterations=100
a=lagduration(strain=0,dist1,time,N,iter)
b=lagduration(strain=0,dist2,time,N,iter)
b=lagduration(strain=0,dist3,time,N,iter)
d=lagduration(strain=0,dist4,time,N,iter)
e=lagduration(strain=0,dist5,time,N,iter)
f=lagduration(strain=0,dist6,time,N,iter)

# Plotting lag duration for the six different distributions
par(mfrow=c(1,1))
plot(c(N),c(d$m),ylim=c(19,48),xlab="Population Start Size (No. of Cells)",
     ylab="Apparent Lag Phase",type='l',cex.lab=1.4,col="cyan",lwd=2,
     main="Break Point Estimates over Time in \nPopulation Simulation of    and     ")
lines(c(N),d$uc,lty=2,col=adjustcolor("red",0.6))
lines(c(N),d$lc,lty=2,col=adjustcolor("red",0.6))
lines(c(N),e$m,lty=1,col="cornflowerblue",lwd=2)
lines(c(N),e$uc,lty=2,col=adjustcolor("red",0.6))
lines(c(N),e$lc,lty=2,col=adjustcolor("red",0.6))
lines(c(N),f$m,lty=1,col="navy",lwd=2)
lines(c(N),f$uc,lty=2,col=adjustcolor("red",0.6))
lines(c(N),f$lc,lty=2,col=adjustcolor("red",0.6))
lines(c(N),a$m,lty=1,col="green",lwd=2)
lines(c(N),a$uc,lty=2,col=adjustcolor("red",0.6))
lines(c(N),a$lc,lty=2,col=adjustcolor("red",0.6))
lines(c(N),b$m,lty=1,col="green3",lwd=2)
lines(c(N),b$uc,lty=2,col=adjustcolor("red",0.6))
lines(c(N),b$lc,lty=2,col=adjustcolor("red",0.6))
lines(c(N),c$m,lty=1,col="darkgreen",lwd=2)
lines(c(N),c$uc,lty=2,col=adjustcolor("red",0.6))
lines(c(N),c$lc,lty=2,col=adjustcolor("red",0.6))
