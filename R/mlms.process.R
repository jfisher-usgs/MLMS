"mlms.process" <- function() {
    
# additional functions (subroutines)
    
  # altitude
    
    alt <- function(bcalt, pos) {
        bcalt - pos
    }
    
  # specific weight of water: temperature in degC
    
    specific.weight <- function(temp) {
        rho <- 1000 * (1 - (temp + 288.9414) / (508929.2 * (temp + 68.12963)) * (temp - 3.9863)^2) # kg/m^3
        g <- 9.80665  # m/s^2
        sw <- rho * g # N/m^3 or kg/(m^2*s^2)
        sw <- sw * 0.22480894387 / 35.314666721 # lb/ft^3
        sw
    }
    
  # pressure head: pressure in psia, temp in degC
    
    press.head <- function(press, temp) {
        psi <- press # psia
        psf <- psi * 144 # lb/ft^2
        p <- psf / specific.weight(temp) # ft
        p
    }
    
    
# main program
    
  # read files
    
    path <- tclvalue(tkchooseDirectory(initialdir=Data("default.dir"), 
            title="Choose Multilevel Monitoring Data Directory..."))
    
    mps.files <- list.files(path, full.names=TRUE, recursive=TRUE, pattern="[.][m][p][s]$")
    sen.files <- list.files(path, full.names=TRUE, recursive=TRUE, pattern="[.][s][e][n]$")
    dat.files <- list.files(path, full.names=TRUE, recursive=TRUE, pattern="[.][d][a][t]$")
    
    for(i in c(mps.files, sen.files, dat.files)) 
        cat(i, "\n")
    
  # process multiport system files (*.mps)
    
    mps <- list()
    for(f in mps.files) {
        d <- mlms.import(f)
        for(id in names(d)) {
            
            bcalt   <- d[[id]]$BCALT.C16
            top     <- d[[id]]$ZONE$TOP
            bot     <- d[[id]]$ZONE$BOT
            pos     <- d[[id]]$PORT$POS
            
            topalt <- alt(bcalt, top)
            botalt <- alt(bcalt, bot)
            posalt <- alt(bcalt, pos)
            len <- abs(topalt - botalt)
            
            d[[id]]$ZONE$TOPALT <- topalt
            d[[id]]$ZONE$BOTALT <- botalt
            d[[id]]$PORT$POSALT <- posalt
            d[[id]]$ZONE$LENGTH <- len
            d[[id]]$MPSFILE     <- f
        }
        mps <- append(mps, d)
    }
    
  # process sensor files (*.sen)
    
    sen <- list()
    for(f in sen.files) {
        d <- mlms.import(f)
        
        for(id in 1:length(d)) {
            
            for(idx in 1:length(d[[id]]$CAL)) {
                
                sw <- specific.weight(d[[id]]$CAL[[idx]]$REFTEMP)
                
                tpres <- d[[id]]$CAL[[idx]]$DAT[,"TPRES"]
                
                error <- d[[id]]$CAL[[idx]]$DAT[,"ERROR"]
                
                mpres <- tpres + error # measurement error = measured pressure - true pressure
                
                ehead <- error * 144 / sw # psia to ft
                
                d[[id]]$CAL[[idx]]$DAT[,"ERRORHEAD"] <- ehead
                
                d[[id]]$CAL[[idx]]$DAT[,"MPRES"] <- mpres
                
                dat <- d[[id]]$CAL[[idx]]$DAT
                
                model <- lm(ERROR ~ MPRES + I(MPRES^2)+ I(MPRES^3), dat)
                d[[id]]$CAL[[idx]]$MODEL <- model
                
                ### plot(ERROR ~ MPRES, dat)
                ### lines(seq(0,500,1), predict(model,data.frame(MPRES=seq(0,500,1))), col="blue")
                ### summary(model)$r.squared
                
            }
        }
        sen <- append(sen, d)
    }
    
  # process sampling data files (*.dat)
    
    dat <- list()
    for(f in dat.files) {
        d <- mlms.import(f)
        
        for(idx in 1:length(d)) {
        
            station <- d[[idx]]$STATION.C12
            portid  <- d[[idx]]$SAMPLE$PORTID
            press   <- d[[idx]]$SAMPLE$PRESS
            baro    <- d[[idx]]$SAMPLE$BARO
            temp    <- d[[idx]]$SAMPLE$TEMP
            
            id <- toupper(gsub(" ", "", station))
            
            bcalt  <- mps[[id]]$BCALT.C16
            zoneid <- mps[[id]]$ZONE[,"ZONEID"]
            top    <- mps[[id]]$ZONE[,"TOP"]
            bot    <- mps[[id]]$ZONE[,"BOT"]
            
            
            pos <- NULL
            for(i in portid) {
                pos <- c(pos, mps[[id]]$PORT[,"POS"][mps[[id]]$PORT[,"PORTID"] %in% i])
            }
            
            
          # total head for each sample
            
            z <- alt(bcalt, pos) # z in ft
            p <- press.head(press, temp) # ft
            b <- press.head(baro, temp) # ft
            H <- z + (p - b)
            
          # zone identification and pressure differences between samples within the same zone
            
            zonediff <- NULL
            zone <- rep(NA, length(pos))
            
            for(i in 1:length(zoneid)) {
                logic <- pos >= top[i] & pos <= bot[i]
                zone[logic] <- zoneid[i]
                
                if(sum(as.integer(logic)) > 1) {
                    sort.idxs <- sort(pos[logic], decreasing=TRUE, index.return=TRUE)$ix
                    
                    zdiff <- ave(diff(z[logic][sort.idxs]))
                    hdiff <- ave(diff(H[logic][sort.idxs]))
                    pdiff <- ave(diff(press[logic][sort.idxs]))
                    bdiff <- ave(diff(baro[logic][sort.idxs]))
                    tdiff <- ave(diff(temp[logic][sort.idxs]))
                    
                    zonediff <- rbind(zonediff, c(zoneid[i], zdiff, hdiff, pdiff, bdiff, tdiff))
                }
            }
            
            if(!is.null(zonediff)) 
                colnames(zonediff) <- c("ZONEID", "POSALTDIFF", "TOTHEADDIFF", 
                                        "PRESSDIFF", "BARODIFF", "TEMPDIFF")
            
            d[[idx]]$SAMPLE$POS     <- pos
            d[[idx]]$SAMPLE$POSALT  <- z
            d[[idx]]$SAMPLE$PHEAD   <- p
            d[[idx]]$SAMPLE$BHEAD   <- b
            d[[idx]]$SAMPLE$TOTHEAD <- H
            d[[idx]]$SAMPLE$ZONE    <- zone
            d[[idx]]$ZONEDIFF       <- zonediff
            d[[idx]]$DATFILE        <- f
        }
        dat <- append(dat, d)
    }
    
  # return processed information in a single list object
    
    list(mps=mps, sen=sen, dat=dat)
}
