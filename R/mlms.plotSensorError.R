"mlms.plotSensorError" <- function(d, sensor="EMS3552", xlim=NULL, ylim=NULL, CALDATE, TEMP) {


  # specific weight of water: temperature in degC

    specific.weight <- function(temp) {
        rho <- 1000 * (1 - (temp + 288.9414) / (508929.2 * (temp + 68.12963)) * (temp - 3.9863)^2) # kg/m^3
        g <- 9.80665  # m/s^2
        sw <- rho * g # N/m^3 or kg/(m^2*s^2)
        sw <- sw * 0.22480894387 / 35.314666721 # lb/ft^3
        sw
    }

  # pressure head: pressure in psi, temp in degC

    press.head <- function(press, temp) {
        psi <- press # psi
        psf <- psi * 144 # lb/ft^2
        p <- psf / specific.weight(temp) # ft
        p
    }



    mlab <- paste(sensor, "CALIBRATION")
    xlab <- "REFERENCE PRESSURE IN PSI"

    y1lab <- "ERROR IN PSI"
    y2lab <- "ERROR IN FEET"

    sen <- d$sen[[sensor]]

    if(missing(CALDATE)) CALDATE <- names(sen$CAL)

    txt <- NULL
    for(i in CALDATE)
        txt <- append(txt, paste(i, sen$CAL[[i]]$REFTEMP, sen$UNITTEMP))
    slab <- paste(txt, collapse=", ")

    n <- length(CALDATE)
    cols <- if(n == 1) "black" else c("red", "blue", "green", "brown", "black")[1:n]

  # axes limits

    x <- y1 <- y2 <- NULL
    for(i in CALDATE) {
        x <- c(x, sen$CAL[[i]]$DAT$TPRES)
        y1 <- c(y1, sen$CAL[[i]]$DAT[,"ERROR"])
        y2 <- c(y2, sen$CAL[[i]]$DAT[,"ERRORHEAD"])
    }
    if(is.null(xlim)) xlim <- extendrange(x,  f=0.02)
    if(is.null(ylim)) ylim <- extendrange(y1, f=0.02)




    cat(paste("\nRange of reference pressure: ", min(x), " to ", max(x),
        " ", sen$UNITPRESS, "\n", sep=""))
    cat("Error analysis:\n")
    cat(paste("Min =", min(y1),  "psi \n"))
    cat(paste("Max =", max(y1),  "psi \n"))
    cat(paste("Ave =", mean(y1), "psi \n"))
    cat(paste("SD =",  sd(y1),   "psi \n\n"))
    cat(paste("Min =", min(y2),  "ft \n"))
    cat(paste("Max =", max(y2),  "ft \n"))
    cat(paste("Ave =", mean(y2), "ft \n"))
    cat(paste("SD =",  sd(y2),   "ft \n\n"))


  # window setup

    x11(width=8, height=5)
    op <- par(mfrow=c(1, 1), bg="white", mar=c(4, 4, 3, 4) + 0.1)
    on.exit(par(op))

  # initialize plot

    plot(NA, xlim=xlim, ylim=ylim, type="n", xaxs="i", yaxs="i", ann=FALSE, cex.axis=0.7)
    axis(1:2)

    title(main=mlab, line=1.5)
    mtext(slab, side=3, line=0.5, cex=0.8, col="dark gray")
    title(ylab=y1lab, xlab=xlab, line=2.5)

    abline(h=0, col="darkgray", lty=1)

  # plot

    for(i in 1:n) {
        cal <- sen$CAL[[CALDATE[i]]]

        x <- cal$DAT$TPRES
        y <- cal$DAT[,"ERROR"]

        lines(x, y, col=cols[i])
        points(x, y, pch=20, cex=1, col=cols[i])
    }


  # add second y-axis

    par(new=TRUE)

    ylim <- press.head(ylim, TEMP)

    plot.window(xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")

    axis(4, at=pretty(range(y2)), padj=-0.7, cex.axis=0.7)

    mtext(y2lab, side=4, line=2.6)



  # legend

    if(n > 1)
        legend(x="topright", CALDATE, lty=1, pch=20, col=cols, cex=0.75,
            pt.cex=1, x.intersp=0.8, bty="n")
}
