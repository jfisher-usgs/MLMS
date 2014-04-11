"mlms.import" <- function(file=NULL) {
    
## function to import multiport source files
## example: dat <- mlms.import()
    
# additional functions (subroutines)
    
  # scan a single row in the input file
    
    rowScan <- function(con) {
        scan(con, comment.char="#", strip.white=TRUE, nlines=1, 
            what="character", sep="\t", quiet=TRUE)
    }
    
  # read multiport system file (*.mps)
    
    readMps <- function(con) {
        
        vars <- NULL
        while(length(vars) == 0) 
            vars <- rowScan(con)
        
        dat <- list()
        while(TRUE) {
            d <- list()
            
            d$STATION.C12 <- vars[1]
            d$SITEID.C1   <- vars[2]
            
            vars <- rowScan(con)
            d$FORMATDATE <- vars[1]
            d$UNITLEN    <- vars[2]
            d$UNITPRESS  <- vars[3]
            
            vars <- rowScan(con)
            d$DATUM.C36 <- vars[1]
            d$LAT.C9    <- as.numeric(vars[2])
            d$LONG.C10  <- as.numeric(vars[3])
            
            vars <- rowScan(con)
            d$DATUM.C22 <- vars[1]
            d$BCALT.C16 <- as.numeric(vars[2])
            d$ACC.C18   <- as.numeric(vars[3])
            
            vars <- rowScan(con)
            d$STICKUP       <- as.numeric(vars[1])
            d$BOREDEPTH.C27 <- as.numeric(vars[2])
            d$COMDEPTH.C28  <- as.numeric(vars[3])
            d$INSTDATE.C60  <- strptime(vars[4], d$FORMATDATE)
            d$DRILLCOMP.C21 <- strptime(vars[5], d$FORMATDATE)
            
            vars <- rowScan(con)
            d$NZONE <- as.integer(vars[1])
            
            d$ZONE <- read.table(con, flush=TRUE, sep="\t", strip.white=TRUE,
                      nrows=d$NZONE, col.names=c("ZONEID", "TOP", "BOT"),
                      colClasses=c("integer", rep("numeric", 2)))
            
            vars <- rowScan(con)
            d$NPORT <- as.integer(vars[1])
            
            d$PORT <- read.table(con, flush=TRUE, sep="\t", strip.white=TRUE, 
                      nrows=d$NPORT, col.names=c("PORTID", "POS", "ID"),
                      colClasses=c("integer", "numeric", "character"))
            
            id <- toupper(gsub(" ", "", d$STATION.C12))
            dat[[id]] <- d
            
            vars <- rowScan(con)
            if(length(vars) == 0) break
        }
        dat
    }
    
  # read sensor file (*.sen)
    
    readSen <- function(con) {
        
        vars <- NULL
        while(length(vars) == 0) 
            vars <- rowScan(con)
        
        dat <- list()
        while(TRUE) {
            d <- list()
            
            d$SENSOR <- vars[1]
            d$SERIAL <- as.integer(vars[2])
            d$BRAND  <- vars[3]
            d$MODEL  <- vars[4]
            
            vars <- rowScan(con)
            d$FORMATDATE <- vars[1]
            d$UNITPRESS  <- vars[2]
            d$UNITTEMP   <- vars[3]
            
            d$TYPE <- rowScan(con)
            
            vars <- rowScan(con)
            d$TEMPLOWER <- as.numeric(vars[1])
            d$TEMPUPPER <- as.numeric(vars[2])
            d$PRESLOWER <- as.numeric(vars[3])
            d$PRESUPPER <- as.numeric(vars[4])
            
            vars <- rowScan(con)
            d$PRESRES   <- as.numeric(vars[1])
            d$PRESACC   <- as.numeric(vars[2])
            d$REPEATACC <- as.numeric(vars[3])
            d$HYSTERACC <- as.numeric(vars[4])
            
            d$RESPTIME <- rowScan(con)
            
            d$NC <- as.integer(rowScan(con))
            
            if(d$NC > 0) {
                d$CAL <- list()
                for(i in 1:d$NC) {
                    vars <- rowScan(con)
                    d$CAL[[i]] <- list()
                    d$CAL[[i]]$CALDATE  <- strptime(vars[1], d$FORMATDATE)
                    d$CAL[[i]]$PREPOST  <- toupper(vars[2])
                    d$CAL[[i]]$REFTEMP  <- as.numeric(vars[3])
                    d$CAL[[i]]$TSDATE   <- strptime(vars[4], d$FORMATDATE)
                    d$CAL[[i]]$LABSTAND <- vars[5]
                    
                    d$CAL[[i]]$NR <- as.integer(rowScan(con))
                    
                    d$CAL[[i]]$DAT <- read.table(con, flush=TRUE, sep="\t", strip.white=TRUE, 
                        nrows=d$CAL[[i]]$NR, col.names=c("TPRES", "ERROR"),
                        colClasses=rep("numeric", 2))
                }
            }
            
            id <- toupper(gsub(" ", "", d$SENSOR))
            dat[[id]] <- d
            
            vars <- rowScan(con)
            if(length(vars) == 0) break
        }
        dat
    }
    
  # read multiport data file (*.dat)
    
    readDat <- function(con) {
        
        vars <- NULL
        while(length(vars) == 0) 
            vars <- rowScan(con)
        
        dat <- list()
        while(TRUE) {
            d <- list()
            
            d$STATION.C12  <- vars[1]
            
            vars <- rowScan(con)
            d$FORMATDT  <- vars[1]
            d$UNITLEN   <- vars[2]
            d$UNITPRESS <- vars[3]
            d$UNITTEMP  <- vars[4]
            
            vars <- rowScan(con)
            d$STIME <- strptime(vars[1], d$FORMATDT)
            d$ETIME <- strptime(vars[2], d$FORMATDT)
            
            vars <- rowScan(con)
            d$OPERATORS <- vars[1]
            d$COMMENT   <- vars[2]
            
            vars <- rowScan(con)
            d$BAROID    <- vars[1]
            d$BAROSTART <- as.numeric(vars[2])
            d$BAROEND   <- as.numeric(vars[3])
            
            vars <- rowScan(con)
            d$PRESID    <- vars[1]
            d$PRESSTART <- as.numeric(vars[2])
            d$PRESEND   <- as.numeric(vars[3])
            
            vars <- rowScan(con)
            d$TEMPID    <- vars[1]
            d$TEMPSTART <- as.numeric(vars[2])
            d$TEMPEND   <- as.numeric(vars[3])
            
            vars <- rowScan(con)
            d$NPS <- as.integer(vars[1])
            
            d$SAMPLE <- read.table(con, flush=TRUE, sep="\t", strip.white=TRUE,
                        nrows=d$NPS, col.names=c("PORTID", "BARO", "PRESS", "TEMP", "TIME", "COMMENT"),
                        colClasses=c("integer", rep("numeric", 3), rep("character", 2)))
            
            id <- paste(toupper(gsub(" ", "", d$STATION.C12)), 
                  format(d$STIME, "%Y%m%d%H%M"), sep="_")
            dat[[id]] <- d
            
            vars <- rowScan(con)
            if(length(vars) == 0) break
        }
        dat
    }
    
    
# main program
    
  # read file
    
    f <- RSurvey::GetFile(cmd="Open", initialdir=RSurvey::Data("default.dir"), 
         win.title=paste("Open Multiport File"), file=file)
    if(is.null(f)) return()
    
    con <- file(f, "r")
    
    dat <- NULL
    if(tolower(attr(f, "extension")) == "mps") dat <- readMps(con)
    if(tolower(attr(f, "extension")) == "sen") dat <- readSen(con)
    if(tolower(attr(f, "extension")) == "dat") dat <- readDat(con)
    
    close(con)
    
    dat
}
