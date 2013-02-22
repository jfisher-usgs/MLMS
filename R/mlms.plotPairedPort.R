"mlms.plotPairedPort" <- function(d, id="USGS 103", yvar="TOTHEADDIFF", stimes, ylim=NULL, showStats=FALSE) {

### require(ggplot2)

    vars <- list(
            POSALTDIFF  = "VERTICAL DISTANCE BETWEEN PAIRED PORTS IN FEET",
            TOTHEADDIFF = "TOTAL HEAD DIFFERENCE BETWEEN PAIRED PORTS IN FEET",
            PRESSDIFF   = "FLUID PRESSURE DIFFERENCE BETWEEN PAIRED PORTS IN PSI",
            BARODIFF    = "BARO. PRESSURE DIFFERENCE BETWEEN PAIRED PORTS IN PSI",
            TEMPDIFF    = "TEMP. DIFFERENCE BETWEEN PAIRED PORTS IN °C"
            )

    lab.x1 <- ""
    lab.y1 <- as.character(vars[yvar])

    id <- toupper(gsub(" ", "", id))

    mps <- d$mps[[id]]

    main <- paste(mps$STATION.C12, " (", mps$SITEID.C1, ")", sep="")
    mlab <- paste("Lat. ", mps$LAT.C9,
               ", Long. ", mps$LONG.C10,
               ", Datum ", mps$DATUM.C36, sep="")

    idxs <- grep(id, names(d$dat))

    if(!missing(stimes)) {
        plot.idxs <- NULL
        for(i in idxs)
            if(format(d$dat[[i]]$STIME, "%Y-%m-%d") %in% stimes)
                plot.idxs <- append(plot.idxs, i)
        idxs <- idxs[idxs %in% plot.idxs]
    }

    n <- length(idxs)

    x1 <- list()
    y1 <- list()

    lim.x1  <- lim.y1  <- periods <- intervals <- NULL

    for(i in 1:n) {
        dat <- d$dat[[idxs[i]]]

        x1[[i]] <- dat$ZONEDIFF[,"ZONEID"]
        y1[[i]] <- dat$ZONEDIFF[,yvar]

        lim.x1 <- range(append(lim.x1, x1[[i]]), na.rm=TRUE)
        lim.y1 <- range(append(lim.y1, y1[[i]]), na.rm=TRUE)

        periods <- append(periods, format(dat$STIME, "%Y-%m-%d"))

        intervals <- unique(append(intervals, x1[[i]]))
    }

  # initialize plot

    plot.new()

    op <- par(oma=c(0,0,0,2), mar=c(3,4,4,3) + 0.1)

    xlim <- extendrange(c(1, length(intervals)), f=0.1)
    if(is.null(ylim))
        ylim <- extendrange(lim.y1, f=0.02)

    plot.window(xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")

    x.at <- 1:length(intervals)
    y.at <- seq(par("yaxp")[1], par("yaxp")[2], length.out=par("yaxp")[3] + 1)

  # axis

    txt <- paste("Interval ", intervals, "\n(",
           round(mps$ZONE[intervals, "BOTALT"], 1), "-",
           round(mps$ZONE[intervals, "TOPALT"], 1), ")", sep="")

    axis(1, las=1, cex.axis=0.6, at=x.at, labels=txt, cex=0.8)
    axis(2, las=0, cex.axis=0.7)

    abline(h=y.at, v=x.at, col="lightgray", lty=3)
    abline(h=0, col="darkgray", lty=1)

  # titles

    title(main=main, xlab=lab.x1, ylab=lab.y1)
    mtext(mlab, side=3, line=0.5, cex=0.8, col="dark gray")

    box()

  # plot data

    cols <- colorRampPalette(c("red", "greenyellow", "skyblue", "darkgrey"))(n)

    x <- y <- col <- NULL
    for(i in 1:n) {
        hld <- x1[[i]]
        for(j in 1:length(intervals))
            hld[hld == intervals[j]] <- j
        x <- append(x, hld)
        y <- append(y, y1[[i]])
        col <- append(col, cols[rep(i, length(hld))])
    }
    points(x, y, type="p", lty=1, pch=21, col="black", bg=col, cex=1)

    statMin <- statMax <- statMean <- statSD <- NULL
    for(i in 1:length(intervals)) {
        val <- y[x == i]
        statMin[i]  <-  min(val, na.rm=TRUE)
        statMax[i]  <-  max(val, na.rm=TRUE)
        statMean[i] <- mean(val, na.rm=TRUE)
        statSD[i]   <-   sd(val, na.rm=TRUE)
    }

    statTbl <- cbind(intervals, statMin, statMax, statMean, statSD)
    write.table(statTbl, col.names=c("Interval", "Min", "Max", "Mean", "SD"),
        row.names=FALSE, quote=FALSE)


    if(showStats) {
        for(i in 1:length(intervals)) {
            points(i + 0.1, statMean[i], col="black", pch=15, cex=0.5)
            ucl <- statMean[i] + 2 * statSD[i]
            lcl <- statMean[i] - 2 * statSD[i]
            arrows(i + 0.1, ucl, i + 0.1, lcl, length=.05, angle=90, code=3, col="black")
        }
    }

  # add legend

    if(n > 1)
        periods <- c(paste(head(periods, -1), " ", sep=""), tail(periods, 1))

    legend(x=grconvertX(1, 'ndc'), y=grconvertY(0.5, 'nfc'), periods,
        xjust=1, yjust=0.5, pch=21, col="black", pt.bg=cols[1:n], cex=0.75, pt.cex=1,
        xpd=NA, x.intersp=0.8, bty="n")

    par(op)
}
