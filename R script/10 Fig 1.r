###################################################################################
#Results of simulated catch data with the catch-based method applied to 
#estimate status. 
#Created by Trevor A. Branch tbranch@uw.edu
#Starting 2 November 2011
#Latest modification: 4 December 2013 for Advanced R FISH553
###################################################################################

####stackpolyTB2#######################################
#Modified version of stackpoly
########################################################
set.seed(123)
stackpolyTB2 <- function (x, y = NULL, main = "", xlab = "", ylab = "", xat = NA, 
    xaxlab = NA, xlim = NA, ylim = NA, lty = 1, border = NA, 
    col = NA, staxx = FALSE, stack = FALSE, axis4 = TRUE, ...) 
{
    ydim <- dim(y)
    if (is.null(y[1])) {
        y <- x
        ydim <- dim(y)
        if (is.null(ydim)) 
            x <- 1:length(y)
        else x <- matrix(rep(1:ydim[1], ydim[2]), ncol = ydim[2])
    }
    if (stack) 
        y <- t(unlist(apply(as.matrix(y), 1, cumsum)))
    if (is.na(xlim[1])) 
        xlim <- range(x)
    if (is.na(ylim[1])) 
        ylim <- range(y)
    plot(0, main = main, xlab = xlab, ylab = ylab, xlim = xlim, 
        ylim = ylim, type = "n", xaxs = "i", yaxs = "i", axes = FALSE, 
        ...)
    box()
    if (is.matrix(y) || is.list(y)) {
        plotlim <- par("usr")
        if (is.na(xat[1])) {
            xat <- x[, 1]
            if (is.na(xaxlab[1])) 
                xaxlab <- xat
        }
        #if (staxx) 
        #    staxlab(at = xat, labels = xaxlab)
        #else axis(1, at = xat, labels = xaxlab)
        #axis(2,las=1)
        #if (axis4) 
        #    axis(4)
        if (is.na(col[1])) 
            col = rainbow(ydim[2])
        else if (length(col) < ydim[2]) 
            col <- rep(col, length.out = ydim[2])
        if (length(lty) < ydim[2]) 
            lty <- rep(lty, length.out = ydim[2])
        for (pline in seq(ydim[2], 1, by = -1)) {
            if (pline == 1) {
                polygon(c(x[1], x[, pline], x[ydim[1]]), c(plotlim[3], 
                  y[, pline], plotlim[3]), border = border, col = col[pline], 
                  lty = lty[pline])
            }
            else polygon(c(x[, pline], rev(x[, pline - 1])), 
                c(y[, pline], rev(y[, pline - 1])), border = border, 
                col = col[pline], lty = lty[pline])
        }
    }
    else {
        polygon(c(min(x), x, max(x), 0), c(0, y, 0, 0), border = border, 
            col = col, lty = lty)
        if (is.na(xat[1])) {
            xat <- x
            if (is.na(xaxlab[1])) 
                xaxlab <- xat
        }
        if (staxx) 
            staxlab(at = xat, labels = xaxlab)
        else axis(1, at = xat, labels = xaxlab)
        #axis(2,las=1)
        #if (axis4) 
        #    axis(4)
    }
}


#############################################################
#Stationary catches over time with mean mu but also with 
#autocorrelation with correlation rho. 
#note: simulates nyears+10 years and then deletes the first 10 years
#############################################################
st.autocorrel.catch <- function(nyears=100, mu=1000,rho=0.5,CV=0.2,plot=F) {
   sim.years <- nyears+10
   ln_SD <- sqrt(log(CV^2+1))  #SD parameter of lognormal distribution with desired variability
   ln_SD_e <- sqrt(ln_SD^2/(1-rho^2))  #SD for first value
   time_series_temp <- rnorm(mean=0,sd=1,n=sim.years)
   time_series <- vector(length=sim.years)
   time_series[1] <- ln_SD_e*time_series_temp[1]
   for (i in 1:(sim.years-1)) {
      time_series[i+1] <- rho*time_series[i]+time_series_temp[i]*ln_SD
      #time_series[i+1] <- rho*time_series[i]+time_series_temp[i]*ln_SD_e
   }
   time_series <- time_series + log(mu) - 0.5*ln_SD_e^2
   #time_series <- time_series + log(mu) - 0.5*ln_SD^2
   time_series <- exp(time_series)
   if (plot==T) {
      plot(time_series,type="l",ylim=c(0,1.05*max(time_series)),las=1,yaxs="i")
   }
   return(time_series[-c(1:10)])
}

##########################################################
#Plot some series of random catches.
##########################################################
lognorm.catch <- function(nyear=100, nstocks=1,meanlog=100, sdlog=0.2, noise="white",autocorrel=0.2, plot.cols,plot.yaxis=F) {
  nyear <- nyear
  nstocks <- nstocks
  catches <- matrix(nrow=nstocks,ncol=nyear)
  status <- matrix(nrow=nstocks,ncol=nyear)
  before.max <- vector(length=nyear)
  for (i in 1:nstocks) { 
     if (noise=="white") {  #catches will be lognormal and independent from one year to the next
        catches[i,] <- rlnorm(n=nyear,meanlog=log(meanlog),sdlog=sdlog)
     }
     if (noise=="Wilberg") {
        catches[i,] <- st.autocorrel.catch(nyears=nyear, mu=meanlog,rho=autocorrel,CV=sdlog)
     }
     max.c <- max(catches[i,])
     max.year <- (1:nyear)[catches[i,]==max.c]
     before.max[] <- F
     before.max[1:nyear < max.year] <- T

     status <- vector(length=nyear)
     
     #developing
     temp <- before.max & catches[i,] < 0.5*max.c
     status[temp] <- 1

     #fully exploited
     temp <- catches[i,] >= 0.5*max.c
     status[temp] <- 2

     #overexploited
     temp <- !before.max & catches[i,] >= 0.1*max.c & catches[i,] < 0.5*max.c
     status[temp] <- 3
     
     #collapsed
     temp <- !before.max & catches[i,] < 0.1*max.c
     status[temp] <- 4
     
     status.cols <- plot.cols[status]
     
     plot(catches[i,],ylim=c(0,1.05*max.c),type="l",ylab="",xlab="",yaxs="i",las=1,axes=F,col="black")  #set up plot area
     abline(h=0.5*max.c,lty=2,col="grey")
     abline(h=0.1*max.c,lty=2,col="grey")

     par(new=T)
     plot(catches[i,],ylim=c(0,1.05*max.c),type="p",ylab="",xlab="",yaxs="i",las=1,axes=F,pch=21,col="black",bg=status.cols)  #set up plot area
     box(col="black")
     if (plot.yaxis==T) {
        axis(2,at=c(0.1*max.c, 0.5*max.c, max.c), labels=c(0.1,0.5,"max"),las=1)
     }
  }
  invisible(catches)
}


##########################################################
#Repeat the plot in Pauly 2007: lognormally distributed random catches
##########################################################
random.status.pauly.2007 <- function(nyear=100, nstocks=1,meanlog=100, sdlog=0.2, noise="white", 
               autocorrel=0.2, plot.cols) {
  year.vec <- 0:nyear
  nyear <- length(year.vec)
  nstocks <- nstocks
  catches <- matrix(nrow=nstocks,ncol=nyear)
  status <- matrix(nrow=nstocks,ncol=nyear)
  before.max <- vector(length=nyear)
  for (i in 1:nstocks) { 
     if (noise=="white") {  #catches will be lognormal and independent from one year to the next
        catches[i,] <- rlnorm(n=nyear,meanlog=log(meanlog),sdlog=sdlog)
     }
     if (noise=="Wilberg") {
        catches[i,] <- st.autocorrel.catch(nyears=nyear, mu=meanlog,rho=autocorrel,CV=sdlog)
     }
     
     max.c <- max(catches[i,])
     max.year <- (1:nyear)[catches[i,]==max.c]
     before.max[] <- F
     before.max[1:nyear < max.year] <- T
     
     temp <- before.max & catches[i,] < 0.5*max.c
     status[i,temp] <- 1
     temp <- catches[i,] >= 0.5*max.c
     status[i,temp] <- 2
     temp <- !before.max & catches[i,] >= 0.1*max.c & catches[i,] < 0.5*max.c
     status[i,temp] <- 3
     temp <- !before.max & catches[i,] < 0.1*max.c
     status[i,temp] <- 4
     
  }
  grp.status <- matrix(nrow=4,ncol=nyear)
  for (j in 1:nyear) {
     for (s in 1:4) {
       grp.status[s,j] <- sum(status[,j]==s)
     }
  }
  
  for (j in 1:nyear) {
      grp.status[,j] <- grp.status[,j] / sum(grp.status[,j])*100
  }

  yvals <- apply(grp.status,MARGIN=2,cumsum)
  year.vecs <- rbind(year.vec,year.vec,year.vec,year.vec)
  stackpolyTB2(x=t(year.vecs),y=t(yvals),col=plot.cols, ylim=c(0,100),xlab="",ylab="", 
            xaxlab=seq(0,50,10),xat=seq(0,50,10), axis4=F)

  invisible(catches)
}

#####simulated.plot###############################
#Sets up a bunch of examples of simulated catch time series
##################################################
simulated.plot2 <- function(autocorrel=0.6, nsims=10000) {
    plot.order <- c(1,2,3,16,8,9,10,
                    4,5,6,16,11,12,13,
                    15,15,15,16,17,17,17,
                    7,7,7,16,14,14,14)

    mat<-matrix(nrow=4,ncol=7,byrow=T,plot.order)
    layout(mat=mat, widths=c(1,1,1,0.2,1,1,1), heights=c(0.8,0.8,0.18,2.5))
    par(mar=c(0.5,0,0,0),oma=c(4.5,5,3,1))
    plot.cols <- c("white","gray67","gray33","black")

    for (i in 1:14) {
       if (i %in% 1:6) {
           if (i %in% c(1,4)) { plot.yaxis<-T } else {plot.yaxis<-F}
           x <- lognorm.catch(nyear=50, nstocks=1,meanlog=100,sdlog=0.2,autocorrel=autocorrel,noise="Wilberg",plot.cols=plot.cols, plot.yaxis=plot.yaxis)
           if (i==1) {
              par(xpd=NA)
              text(x=-23,y=0,"Proportion of catch",srt=90,cex=1.45)
              text(x=0,y=1.2*max(x),pos=4,cex=1.3,"(a)")
              par(xpd=T)
           }
       }
       if (i==7) {
          x <- random.status.pauly.2007(nyear=50, nstocks=nsims,meanlog=100,sdlog=0.2, autocorrel=autocorrel,noise="Wilberg",
             plot.cols=plot.cols)
          par(xpd=NA)
          text(x=0,y=104,pos=4,cex=1.3,"(c)")
          par(xpd=T)
      }

       if (i %in% 8:13) {
           x <- lognorm.catch(nyear=50, nstocks=1,meanlog=100,sdlog=0.6,autocorrel=autocorrel,noise="Wilberg",plot.cols=plot.cols)
           if (i==8) {
               par(xpd=NA)
               text(x=0,y=1.2*max(x),pos=4,cex=1.3,"(b)")
               par(xpd=T)
           }
       }
       if (i==14) {
          random.status.pauly.2007(nyear=50, nstocks=nsims,meanlog=100,sdlog=0.6, autocorrel=autocorrel,noise="Wilberg",
             plot.cols=plot.cols)
          legend(x="bottomleft",legend=rev(c("Developing","Fully exploited","Overexploited","Collapsed")), 
             pch=22,pt.bg=rev(plot.cols), col="black",pt.cex=2.5, cex=1.18,bty="n")
          par(xpd=NA)
          text(x=0,y=104,pos=4,cex=1.3,"(d)")
          par(xpd=T)
       }

       box()
       if (i==2) {
          mtext(outer=F,line=1,"Low variability")
       }
       if (i==9) {
          mtext(outer=F,line=1,"High variability")
       }
       if (i==7) {
          mtext(side=2,outer=F, line=3,"Percentage of stocks")
          axis(2,las=1)
       }
       if (i %in% c(7,14)) {
          axis(1,las=1)
       }
    }

    mtext(side=1,outer=T, line=2.5,"Year")
}
pdf("Fig 1 v1.pdf",width=8,height=6)
simulated.plot2(autocorrel=0.5, nsims=2000)
dev.off()
#########source("10 Fig 1.r")