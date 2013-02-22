"mlms.calibration" <- function(d, prodate="2008-02-15", pres=120, temp=13, sensor="EMS3552", showPlots=FALSE) {
    
    ### mlms.calibration(d, prodate="2007-08-15", pres=120, temp=13, showPlots=TRUE)
    
# additional functions (subroutines)
    
    linfun <- function(x, x1, x2, m1, m2, p) {
        y1 <- as.vector(predict(m1, data.frame(MPRES=p)))
        y2 <- as.vector(predict(m2, data.frame(MPRES=p)))
        y <- ((((y2 - y1) / (x2 - x1))) * (x - x1)) + y1
        
        dat <- as.data.frame(cbind(p, y))
        names(dat) <- c("MPRES", "ERROR")
        
        m <- lm(ERROR ~ MPRES + I(MPRES^2)+ I(MPRES^3), dat)
        
        if(showPlots) {
            lines(p, y1, col="red")
            lines(p, y2, col="blue")
            lines(p, y,  col="black")
        }
        
        m
    }
    
    
    
# main program
    
    if(is.null(d$sen[[sensor]])) stop("Sensor not found")
    
    cal <- d$sen[[sensor]]$CAL
    
    caldate <- sapply(cal, function(i) as.Date(i$CALDATE))
    reftemp <- sapply(cal, function(i) i$REFTEMP)
    prepost <- sapply(cal, function(i) i$PREPOST)
    dat     <- lapply(cal, function(i) i$DAT)
    model   <- lapply(cal, function(i) i$MODEL)
    
  # simplify dates
    
    prodate <- as.numeric(as.Date(prodate))
    
    if(prodate < min(caldate) || prodate > max(caldate)) stop("Profile date (prodate) outside of calibration dates")
    
    dates <- c(max(caldate[caldate < prodate]), min(caldate[caldate > prodate]))
    
    cal1 <- caldate == dates[1] & prepost == "POST"
    cal2 <- caldate == dates[2] & prepost == "PRE"
    
    cals <- cal1 | cal2
    
  # simplify temperature
    
    temp <- as.numeric(temp)
    
    temps <- unique(reftemp[cals])
    
    if(temp < min(temps) || temp > max(temps)) stop("Temperature (temp) outside of calibration reference temperatures")
    
  # determine error range
    
    ylim <- NULL
    idxs <- (1:length(dat))[cals]
    for(idx in idxs) 
        ylim <- range(c(ylim, dat[[idx]]$ERROR), na.rm=TRUE)
    
    
  # temperature correction
    
    tempmodel <- list()
    j <- 0
    pps <- NULL
    
    for(i in dates) {
        
        idxs <- (1:length(dat))[cals & caldate %in% i]
        
        x <- temp
        x1 <- reftemp[idxs[1]]
        x2 <- reftemp[idxs[2]]
        m1 <- model[[idxs[1]]]
        m2 <- model[[idxs[2]]]
        
        mpres <- error <- NULL
        for(idx in idxs) {
            mpres <- c(mpres, dat[[idx]]$MPRES)
            error <- c(error, dat[[idx]]$ERROR)
        }
        
        xlim <- range(mpres)
        
        if(pres < xlim[1] || pres > xlim[2]) stop("Pressure (pres) outside of calibration range")
        
        p <- seq(xlim[1], xlim[2], length.out=100)
        
        pp <- unique(prepost[idxs])
        
        if(showPlots) {
            x11()
            plot(NA, xlim=xlim, ylim=ylim, type="n", 
                xlab="MEASURED PRESSURE IN PSIA", ylab="MEASUREMENT ERROR IN PSI", 
                main=paste("TEMPERATURE CORRECTION ON ", as.Date(i, origin="1970-01-01"), 
                     " ", pp, "-CALIBRATION DATA", sep=""), cex.main=0.95)
            abline(h=0, col="lightgray")
            points(dat[[idxs[1]]]$MPRES, dat[[idxs[1]]]$ERROR, col="red")
            points(dat[[idxs[2]]]$MPRES, dat[[idxs[2]]]$ERROR, col="blue")
        }
        
        pps <- tolower(c(pps, pp))
        
        j <- j + 1
        tempmodel[[j]] <- linfun(x, x1, x2, m1, m2, p)
        
        if(showPlots) {
            leg.txt <- c(paste("Temperature=", format(x1), "; R^2=", format(summary(m1)$r.squared), sep=""), 
                         paste("Temperature=", format(x2), "; R^2=", format(summary(m2)$r.squared), sep=""),
                         paste("Temperature=", x, sep=""))
            legend(x="bottomright", leg.txt, lty=1, pch=c(1, 1, -1), cex=0.8, pt.cex=1, col=c("red", "blue", "black"))
        }
    }
    
  # date correction using temperature corrected models
    
    if(showPlots) {
        x11()
        plot(NA, xlim=xlim, ylim=ylim, type="n", 
            xlab="MEASURED PRESSURE IN PSIA", ylab="MEASUREMENT ERROR IN PSI", 
            main="DATE CORRECTION ON TEMPERATURE CORRECTED DATA", cex.main=0.95)
        abline(h=0, col="lightgray")
    }
    
    if(length(tempmodel) > 1) {
        x <- prodate
        x1 <- dates[1]
        x2 <- dates[2]
        m1 <- tempmodel[[1]]
        m2 <- tempmodel[[2]]
        datetempmodel <- linfun(x, x1, x2, m1, m2, p)
    } else {
        datetempmodel <- tempmodel[[1]]
    }
    
 # predict calibration error
    
    corError <- as.numeric(predict(datetempmodel, data.frame(MPRES=pres)))
    truePres <- pres - corError
    
    if(showPlots) {
        leg.txt <- c(paste(as.Date(dates[1], origin="1970-01-01"), " ", pps[1], "-calibration; Temperature=", format(temp), sep=""), 
                     paste(as.Date(dates[2], origin="1970-01-01"), " ", pps[2], "-calibration; Temperature=", format(temp), sep=""),
                     paste(as.Date(prodate,  origin="1970-01-01")))
        legend(x="bottomright", leg.txt, lty=1, pch=-1, cex=0.8, pt.cex=1, col=c("red", "blue", "black"))
        
        points(pres, corError, pch=19, col="green")
        text(pres, corError, cex=0.8, paste(
            as.Date(prodate, origin="1970-01-01"), "; ", "Temperature=", temp, "\n",
            "Measured Pressure - Error = True Pressure\n", 
            format(pres), " - ", format(corError), " = ", format(truePres), " psia", sep=""), pos=1)
    }
    
    truePres
}
