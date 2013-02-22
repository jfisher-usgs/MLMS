"mlms.plotProfile" <- function(d, id="USGS 103", xvar="TOTHEAD", stimes, brklns,
                 xlim=NULL, ylim=NULL, zone=FALSE, normalize=FALSE, corr=FALSE) {

    vars <- list(
            POS     = "DEPTH OF TRANSDUCER IN FEET BLS",
            POSALT  = "ALTITUDE OF TRANSDUCER IN FEET AMSL",
            BHEAD   = "BAROMETRIC HEAD IN FEET",
            PHEAD   = "ABSOLUTE PRESSURE HEAD",
            TOTHEAD = "HYDRAULIC HEAD IN FEET",
            PRSHEAD = "PRESSURE HEAD IN FEET",
            TEMP    = "TEMPERATURE IN °C",
            PRESS   = "FLUID PRESSURE IN PSIA",
            BARO    = "BAROMETRIC PRESSURE IN PSIA"
            )

    lab.x1 <- as.character(vars[xvar])
    lab.y1 <- as.character(vars["POSALT"])
    lab.y2 <- as.character(vars["POS"])

    id <- toupper(gsub(" ", "", id))

    mps <- d$mps[[id]]
    bcalt <- mps$BCALT.C16

    main <- paste(mps$STATION.C12, " (", mps$SITEID.C1, ")", sep="")
    mlab <- paste("Lat. ", mps$LAT.C9, ", Long. ", mps$LONG.C10,
            ", Datum ", mps$DATUM.C36, sep="")

    idxs <- grep(id, names(d$dat))
    if(length(idxs) == 0) stop(call.=FALSE, "no match on well ID")

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
    y2 <- list()

    lim.x1 <- lim.y1 <- lim.y2 <- NULL

    periods <- NULL

    lm.tbl <- NULL
    lm.tbl.names <- c("STime", "Bottom", "Top", "Intercept", "dx/dy", "R-squared", "NumPoints", "xrange", "xregrange", "xdiff")

    for(i in 1:n) {
        dat <- d$dat[[idxs[i]]]

        x1[[i]] <- dat$SAMPLE[,xvar]
        y1[[i]] <- dat$SAMPLE[,"POSALT"]
        y2[[i]] <- dat$SAMPLE[,"POS"]

        lim.x1 <- range(append(lim.x1, x1[[i]]), finite=TRUE)
        lim.y1 <- range(append(lim.y1, y1[[i]]), finite=TRUE)
        lim.y2 <- range(append(lim.y2, y2[[i]]), finite=TRUE)

        if(zone) {
            minlogic <- lim.y1[1] <= mps$ZONE[,"TOPALT"] & lim.y1[1] >= mps$ZONE[,"BOTALT"]
            maxlogic <- lim.y1[2] <= mps$ZONE[,"TOPALT"] & lim.y1[2] >= mps$ZONE[,"BOTALT"]
            if(any(minlogic)) {
                lim.y1[1] <- mps$ZONE[minlogic, "BOTALT"]
                lim.y2[1] <- mps$ZONE[minlogic, "BOT"]
            }
            if(any(maxlogic)) {
                lim.y1[2] <- mps$ZONE[maxlogic, "TOPALT"]
                lim.y2[2] <- mps$ZONE[maxlogic, "TOP"]
            }
        }

        st <- format(dat$STIME, "%Y-%m-%d")
        periods <- append(periods, st)

        if(!missing(brklns)) {
            if(length(brklns) < 2 | !is.numeric(brklns))
                stop(call.=FALSE, "Problem with brklns argument.")

            nsegs <- length(brklns) - 1

            for(j in 1:nsegs) {
                logic <- y1[[i]] >= brklns[j] & y1[[i]] <= brklns[j+1]

                y <- x1[[i]][logic]
                x <- y1[[i]][logic]

                if(length(y) < 2)
                    reg <- rep(NA, 3)
                else {
                    mdl <- lm(y ~ x)
                    reg  <- as.numeric(c(coef(mdl), summary(mdl)$r.squared))
                }

                numpts <- length(na.omit(y))

                if(numpts == 1) {
                    ran <- regran <- dif <- NA
                }
                else {
                    ran <- diff(range(y, na.rm=TRUE))
                    regran <- abs(reg[2] * diff(range(x, na.rm=TRUE)))
                    dif <- ran - regran
                }

                r <- data.frame(st, brklns[j], brklns[j+1], reg[1], reg[2], reg[3],
                     numpts, ran, regran, dif)
                names(r) <- lm.tbl.names
                lm.tbl <- rbind(lm.tbl, r)
            }
        }
    }

    if(corr) {
        dat.x1 <- NULL
        for(i in 1:n) {
            dat.x1 <- cbind(dat.x1, x1[[i]])
        }
        colnames(dat.x1) <- periods

        r <- cor(dat.x1, method="pearson", use="complete.obs")

        cat("Min. correlation coef. = ", min(r), "\n", sep="")
        print(r)
    }

    if(normalize) {
        slidex1 <- list()
        for(i in 1:n) {
            slidex1[[i]] <- mean(x1[[i]], na.rm=TRUE)
        }
        lim.x1 <- NULL
        for(i in 1:n) {
            x1[[i]] <- x1[[i]] - slidex1[[i]]
            lim.x1 <- range(c(lim.x1, x1[[i]]), na.rm=TRUE)
        }
        lab.x1 <- paste("ZERO-MEAN NORMALIZED", lab.x1)
    }


    cols <- c("Min", "Max", "Mean", "Var")
    tbl.summary <- matrix(NA, nrow=n, ncol=length(cols), dimnames=list(periods, cols))
    for(i in 1:n) {
        tbl.summary[i,] <- c(min(x1[[i]], na.rm=TRUE),
                             max(x1[[i]], na.rm=TRUE),
                            mean(x1[[i]], na.rm=TRUE),
                             var(x1[[i]], na.rm=TRUE))
    }
    print(tbl.summary)


    depth <- y2[[1]]
    xmin <- xmax <- x1[[1]]
    if(n > 1) {
        for(i in 2:n) {
            xmin <- pmin(xmin, x1[[i]], na.rm=TRUE)
            xmax <- pmax(xmax, x1[[i]], na.rm=TRUE)
        }
        maxdiff <- abs(xmin - xmax)
        write.table(cbind(depth, maxdiff), row.names=FALSE, quote=FALSE)
    }

    if(!missing(brklns)) {
        cat("Linear Regression Statistics Within Vertical Segments:\n")
        print(lm.tbl)
    }

  # set outer margin areas (only necessary in order to plot extra y-axis)

    op <- par(oma=c(2,0,0,0), mar=c(5,4,4,4) + 0.1)

  # initialize plot

    if(is.null(xlim)) xlim <- extendrange(lim.x1, f=0.02)
    if(is.null(ylim)) ylim <- extendrange(lim.y1, f=0.02)

    plot.new()
    plot.window(xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")

    axis(1, las=1, cex.axis=0.7)
    axis(2, las=0, cex.axis=0.7)

    title(main=main, xlab=lab.x1, ylab=lab.y1)
    mtext(mlab, side=3, line=0.5, cex=0.8, col="dark gray")

    if(!missing(brklns))
        abline(h=brklns, col="lightgray", lty=3)

    box()

  # plot data

    cols <- colorRampPalette(c("red", "skyblue", "greenyellow", "black"))(n)

    if(zone) {
        pchs <- 20
        for(i in 1:n) {
            for(j in mps$ZONE$ZONEID) {
                ucl <- mps$ZONE[j, "TOPALT"]
                lcl <- mps$ZONE[j, "BOTALT"]
                logic <- y1[[i]] >= lcl & y1[[i]] <= ucl
                if(!any(logic)) next
                x <- x1[[i]][logic]
                y <- y1[[i]][logic]
                xm <- mean(x1[[i]][logic], na.rm=TRUE)
                points(x, y, lty=1, pch=pchs, cex=0.8, col=cols[i])
                arrows(xm, ucl, xm, lcl, length=.05, angle=90, code=3, col=cols[i])

            }
        }
    }
    else {
        pchs <- as.integer(substr(periods, 1, 4))
        pchs <- 15 + (pchs - min(pchs))

        for(i in 1:n)
            points(x1[[i]], y1[[i]], type="o", lty=1, pch=pchs[i], col=cols[i], cex=0.8)
    }

  # add second y-axis

    par(new=TRUE)

    ylim <- bcalt - ylim

    plot.window(xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")

    axis(4, at=pretty(range(y2)), padj=-0.7, cex.axis=0.7)

    mtext(lab.y2, side=4, line=2.6)

  # add legend

    if(n > 1)
        periods <- c(paste(head(periods, -1), " ", sep=""), tail(periods, 1))

    n.col <- if(n < 8) n else 7

    legend(x=grconvertX(0.5, 'ndc'), y=grconvertY(0, 'nfc'), periods,
        xjust=0.5, yjust=1, pch=pchs, lty=1, col=cols, cex=0.7, pt.cex=0.8,
        xpd=NA, x.intersp=0.8, bty="n", ncol=n.col)

    par(op)
}
