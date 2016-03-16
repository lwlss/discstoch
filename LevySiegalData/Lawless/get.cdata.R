get.cdata <- function(d, gr, control.color=2, test.color =1, well.variation = TRUE, contrast.stretching =TRUE, controlrate = 0.8507, conditions.file= NULL,  min.control.colonies = 200){

  #source("/Volumes/X/code/R/Functions.R")
  
  
#ID poorly imaged fields
  fieldscore <- 1:max(d$imnumber)
  for(j in 1:length(fieldscore)){
    x <- gr$areas.raw[which(d$imnumber == j),]
    y <- x
    if(is.matrix(x)){
      for(i in 2:ncol(x)){
        y[,i-1] <- x[,i] - x[,i-1]
      }
      y <- y[,1:(ncol(y) - 1)]
      z <- y
      z[] <- 0
      z[ which(y < -50)] <- 1
      a <- apply(z, 1, sum)
    }else{
      for(i in 2:length(x)){
        y[i-1] <- x[i] - x[i-1]
      }
      y <- y[1:(length(y) - 1)]
      z <- y
      z[] <- 0
      z[ which(y < -50)] <- 1
      a <- sum(z)
    }
    fieldscore[j] <- length(which(a > 1))/length(a)
  }
  
  fieldscore[which(is.nan(fieldscore))] <- 1
  for(i in 1:length(fieldscore)){
    if(fieldscore[i] > .1){
      gr$colparam[which(d$imnumber == i),] <- NA
    }
  }

  
  #Get concensus color
  #if(! file.exists(paste(directory, "/colony.color.Rfile", sep=""))){
  c <- d$color
  c[which(c==0)] <- NA
  cc <- 1:nrow(c) #concensus color
  
  
  for(i in 1:nrow(c)){
    x <- table(c[i,])
    if(max(x) > 3){
      cc[i] <- as.numeric(names(x)[which(x==max(x))])
    }else{
      cc[i] <- 0
    }
  }
  colony.color <- cc
  
  #Make color-defined well lists and colparams
  wl0 <- wl1 <- wl2 <- list()
  for(i in 1:length(d$well.list)){
    wl1[[i]] <- d$well.list[[i]][which(colony.color[d$well.list[[i]]] == test.color)]
    wl2[[i]] <- d$well.list[[i]][which(colony.color[d$well.list[[i]]] == control.color)]
    wl0[[i]] <- d$well.list[[i]][which(colony.color[d$well.list[[i]]] == 0)] 
  }
  names(wl0) <- names(wl1) <- names(wl2) <- names(d$well.list)

  cp1 <- cp2 <-  gr$colparam
  for(i in 1:length(d$well.list)){
    cp1[c(wl0[[i]], wl2[[i]]),] <- NA
    cp2[c(wl1[[i]], wl0[[i]]),] <- NA
    #cp3[c(wl2[[i]],wl1[[i]]),] <- NA
  }
  
  controlwl <- wl2
  testwl <- wl1
  controlcp <- cp2
  testcp <- cp1

  
  ########################################################################################
  #Get control info
  ########################################################################################
  y <- s <- ccount <- 1:96; y[ ] <- NA
  w <- d$wells.used
  for (i in 1:length(w)){
    z <- gr$colparam[controlwl[[i]], 1]
    ccount[w[i]] <- length(which(! is.na(z)))
    zsd <- sd(z, na.rm=T)
    s[w[i]] <- zsd
    y[w[i]] <- mean(z, na.rm=T)
  }
  controlmean <- y
  controlcount <- ccount
  controlsd <- s

  
  ########################################################################################
  #Get test info
  ########################################################################################
  y <- s <- ccount <- 1:96; y[ ] <- NA
  w <- d$wells.used
  for (i in 1:length(w)){
    z <- gr$colparam[testwl[[i]], 1]
    ccount[w[i]] <- length(which(! is.na(z)))
    zsd <- sd(z, na.rm=T)
    s[w[i]] <- zsd
    y[w[i]] <- mean(z, na.rm=T)
  }
  testmean <- y
  testcount <- ccount
  testsd <- s

   ########################################################################################
  
  #Make list of test values normalized by the control mean
  wm <- list()
  wcp <- gr$colparam
  wcp[,] <- NA

  if(contrast.stretching){
    e <- 1/(1+ (mean(controlmean, na.rm=T)/controlmean)^1.7) #contrast streching
  }else{
    e <- controlmean
  }
  for(i in 1:length(w)){
    if(controlcount[w[i]] > min.control.colonies && names(w[i]) %in% names(testwl)){ #need a minimum number of control colonies
      z <- gr$colparam[testwl[[names(w[i])]], 1]/e[w[i]]
      z <- z/controlrate
      wm[[names(w[i])]] <- z
      ccount[w[i]] <- length(which(! is.na(z)))
      zsd <- sd(z, na.rm=T)
      s[w[i]] <- zsd
      y[w[i]] <- mean(z, na.rm=T)
    }else{
      ccount[w[i]] <- s[w[i]] <- y[w[i]] <- NA
      wm[[names(w[i])]] <- NA
    }
  }
  testmean2 <- y
  testcount2 <- ccount
  testsd2 <- s
  testcp2 <- testcp
  names(wm) <- names(d$well.list) #the growth rates of the test in each well
  for(i in 1:length(d$well.list)){
    testcp2[ testwl[[i]] ,1] <- wm[[i]]
  }

  
 #colony data to use later
  cdata <- list()
  g <- d$well.list
  gt <- 1:length(d$well.list)
  for(i in 1:length(g)){
    g[[i]] <- testcp2[g[[i]], 1]
    g[[i]] <- g[[i]][which(! is.na(g[[i]]))]
    gt[i] <- length(g[[i]])
  }
  cdata[[1]] <- g
  gc <- 1:length(d$well.list)
  g <- d$well.list
  for(i in 1:length(g)){
    g[[i]] <- controlcp[g[[i]], 1]
    g[[i]] <- g[[i]][which(! is.na(g[[i]]))]
    gc[i] <- length(g[[i]])
  }
  cdata[[2]] <- g
  g <- d$well.list
  for(i in 1:length(g)){
    g[[i]] <- testcp[g[[i]], 1]
    g[[i]] <- g[[i]][which(! is.na(g[[i]]))]
  }
  cdata[[3]] <- g
  cdata[[4]] <- testwl
  cdata[[5]] <- controlwl
  cdata[[6]] <- testcp
  cdata[[7]] <- controlcp
  cdata[[8]] <- gt
  cdata[[9]] <- gc
  cdata[[10]] <- d$condition.names
  
  names(cdata) <- c("test.norm", "control", "test", "test.wl", "control.wl",  "test.cp", "control.cp", "test.count", "control.count", "condition.names")
  return(cdata)
}
