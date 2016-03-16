

directory = '/Volumes/X/out/GR/100604GR'
setwd(directory)
source("/Volumes/X/code/R/Functions.R")
load(paste(directory, "/d.Rfile", sep=""))
load(paste(directory, "/gr6.Rfile", sep="")) #growth rate based on 6 hours

load( paste(directory, "/md.Rfile", sep=""))
cf <- md$cenfield
md <- md$mindist
x <- md; x[] <- 0
x[which(md < 35)] <- 1
x[which(cf == 0)] <- 1
gr$colparam[which(x == 1),1] <- NA


a <- log(d$areas[,1:6]) #log areas matrix
b <- a
for(i in 1:(ncol(a)-1)){
  b[,i] <- a[,i+1] - a[,i]
}
b <- b[,1:5]
b[which(is.infinite(b))] <- NA
b[which(is.nan(b))] <- NA
n <- which(is.na(gr$colparam[,1]))


f <- d$fint[,1:6]
f <- f/d$areas[,1:6]
fm <- apply(f, 1, mean, na.rm=T)
fm[which(is.nan(fm))] <- NA
b[n,] <- NA; a[n,] <- NA; f[n,] <- NA

n <- unique(d$condition.names)
n <- n[which(! n == 0)]
w <- which(d$condition.names == n[1])
wl <- NULL
for(j in 1:length(w)){
  wl <- c(wl, d$well.list[[names(w)[j]]])
}
tsl <- wl
w <- which(d$condition.names == n[2])
wl <- NULL
for(j in 1:length(w)){
  wl <- c(wl, d$well.list[[names(w)[j]]])
}
tma <- wl


#binning TSL1
bins <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7)
cp <- gr$colparam[,1]
wl <- list()
wl[[1]] <- which(! is.na(cp) & 1:length(cp) %in% tsl)
for(i in 1:(length(bins) - 1)){
  wl[[i+1]] <- which(cp >= bins[i] & cp < bins[i +1] & 1:length(cp) %in% tsl)
  names(wl)[i+1] <- paste(bins[i],"-",bins[i+1])
}
names(wl)[1] <- "All"


tsl.colors <- colorRampPalette(c(colors()[258],"grey", "dark magenta"))
source("/Volumes/X/code/R/Barplot.R")
pdf(file = "~/Documents/Papers/Bethedging/figureparts/fig2c.pdf", height = 6, width = 5.6)
par(mai = c(2.3,1.5,.5,.5), lwd = 2)
i=6
Barplot(f[,i]*0.45, wl, asterix = T, cex.names = 1.3, ylim = c(200, 375), comparison = 1, space = 0.2, color = c("grey", tsl.colors(6)), density = NA, ylab = expression("TSL1-GFP " (counts / mu*m^2)), xlab = expression("Specific growth rate " (h^-1)), horizontal.labels = F, xlab.adjustment = 4, ylab.adjustment = -1, cex.lab = 1.3, asterix.bins = c(.01, 1e-5, 1e-10, 0))
dev.off()

