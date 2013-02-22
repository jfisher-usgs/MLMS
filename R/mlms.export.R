"mlms.export" <- function(d, ids, stime, etime, outfile) {
    
### mlms.export(d, etime="2009-03-31")
### mlms.export(d, stime="2009-04-01")
    
    if(missing(d)) stop(call.=FALSE, "no data provided")
    
    if(missing(ids)) 
        ids <- names(d$mps)
    
    ids <- toupper(gsub(" ", "", ids))
    idxs <- NULL
    for(id in ids) 
        idxs <- append(idxs, grep(id, names(d$dat)))
    if(length(idxs) == 0) stop(call.=FALSE, "no match on well ID")
    
    sts <- NULL
    for(i in idxs) 
        sts <- append(sts, format(d$dat[[i]]$STIME, "%Y-%m-%d"))
    sts <- strptime(sts, "%Y-%m-%d", tz="")
    
    stime <- if(missing(stime)) min(sts) else strptime(stime, "%Y-%m-%d", tz="")
    etime <- if(missing(etime)) max(sts) else strptime(etime, "%Y-%m-%d", tz="")
    idxs <- idxs[sts >= stime & sts <= etime]
    
    
    if(missing(outfile)) {
        f <- GetFile(cmd="Save As", exts="txt", win.title="Save MLMS Data As", defaultextension="txt")
        if(is.null(f)) stop(call.=FALSE, "no output file")
        outfile <- f
    }
    
    dat <- NULL
    
    for(i in idxs) {
        
        id <- strsplit(names(d$dat)[i], "_")[[1]][1]
        siteid <- d$mps[[id]]$SITEID.C1
        bcalt <- d$mps[[id]]$BCALT.C16
        
        ports<- d$dat[[i]]$SAMPLE[,"PORTID"]
        n <- length(ports)
        if(id == "USGS103") ports <- ports + 1 # needs to be re-examined
        sitenums <- paste(substr(siteid, 1, nchar(siteid) - 2), sprintf("%02d", ports), sep="")
        
        dts <- strsplit(d$dat[[i]]$SAMPLE[,"TIME"], " ")
        dates <- gsub("-", "", sapply(dts, function(x) x[1]))
        times <- sapply(dts, function(x) x[2])
        
        wls <- sprintf("%.2f", bcalt - d$dat[[i]]$SAMPLE[,"TOTHEAD"])
        
        accrs <- rep(2, n)
        meths <- rep("F", n)
        srcs <- rep("S", n)
        agncs <- rep("USGS", n)
        
        wells <- rep(id, n)
        
        zones <- d$dat[[i]]$SAMPLE[,"ZONE"]
        ports <- d$dat[[i]]$SAMPLE[,"PORTID"]
        
        posalts <- sprintf("%.2f", d$dat[[i]]$SAMPLE[,"POSALT"])
        posdpts <- sprintf("%.2f", bcalt - d$dat[[i]]$SAMPLE[,"POSALT"])
        
        temps <- sprintf("%.2f", d$dat[[i]]$SAMPLE[,"TEMP"])
        
        baros <- sprintf("%.3f", d$dat[[i]]$SAMPLE[,"BARO"])
        
        press <- sprintf("%.2f", d$dat[[i]]$SAMPLE[,"PRESS"])
        
        heads <- sprintf("%.2f", d$dat[[i]]$SAMPLE[,"TOTHEAD"])
        
        cmmts <- d$dat[[i]]$SAMPLE[,"COMMENT"]
        cmmts[is.na(cmmts)] <- ""
        
        dat <- rbind(dat, cbind(sitenums, dates, times, wls, accrs, meths, srcs, agncs, 
               wells, zones, ports, posalts, posdpts, temps, baros, press, heads, cmmts))
    }
    
    head1 <- c("SITE NUMBER", "DATE", "TIME", "WL BLW LSD", "ACCURACY", "METHOD", "SOURCE", "AGENCY", 
             "WELL NAME", "ZONE NUMBER", "PORT NUMBER", "PORT ALT", "PORT DEPTH", 
             "FLUID TEMPERATURE", "BAROMETRIC PRESSURE", "FLUID PRESSURE", "HYDRAULIC HEAD", "COMMENT")
    head2 <- c("C1", "C235", "C709", "C237", "C276", "C239", "C244", "C247",
             "", "", "", "ft amsl", "ft bls", 
             "°C", "psi", "psi", "ft amsl", "")
    
    dat <- rbind(head1, head2, dat)
    
    con <- file(outfile, "w")
    write.table(dat, file=con, append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
    if(exists("con")) close(con)
    
}
