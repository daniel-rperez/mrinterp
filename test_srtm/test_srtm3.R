rm(list=ls())

library(parallel)
library(foreach)
library(doParallel)

library(sp)
library(raster)
library(rgdal)

#auxiliary functions
hypot <- function(a,b) sqrt( a*a+b*b )
mod <- function(x,a) (x-floor(x/a)*a)

#compute data weights with respect to cross-validation and set outlayers to NaN
#z0: depth data
#zm, sz: cross-validation depth and standar deviation estimation
#dz0: intrinsic error estimation (ex. dz0=0.5 m)
#rerr: relative error cutoff (ex. rerr = 0.1 = 10% of depth)
#fces: Tukey fences factor (ex. fces = 1.5 outlayer; fces = 3 far out)
compWgtsF <- function(z0, zm, sz, dz0=0.0, rerr=0.0, fces=Inf) {
  dz <- z0-zm
  
  #weight
  wgt <- 1.0 / ( 1.0 + (zm/sz)^2 )
  
  #outlayer 1: cross-validation error more than rerr of the depth
  m <- ( sz > sqrt(dz0^2 + (rerr*zm)^2) )
  
  #outlayer 2: Tukey fences
  r <- dz / zm
  f <- quantile(r, c(0.25,0.75), na.rm=T)
    f <- f + fces*c( -diff(f), diff(f) )
  m <- m | (( r < f[1]) | ( f[2] < r ))
  
  wgt[m] <- NaN
  
  return(wgt)
}


#Parallel setup
foldK <- 10
nCor <- 3 #5 cores por motivos de memoria; además es la longitud de seq(4,8)
pClust <- makeCluster(nCor, type = "PSOCK")
registerDoParallel(cl=pClust)
#getDoParRegistered()
#getDoParWorkers()

#Load SRTM file
if( ! file.exists( paste("data", "srtm2subset.tif", sep="/") ) ) {
  dxy <- 90
  s <- raster(paste("data", "srtm2.tif", sep="/"))
  b <- extent( c(-69.6,-65.7, -48.1,-44.2) )
  s <- crop(s, b)
  s <- projectRaster(s, crs=CRS("+proj=utm +zone=19 +south +datum=WGS84 +units=m +no_defs"), res=dxy)
  writeRaster(s, paste("data", "srtm2subset.tif", sep="/"))
} else {
  s <- raster(paste("data", "srtm2subset.tif", sep="/"))
}
bb <- extent( s )
dxy <- (bb@xmax-bb@xmin)/ncol(s)
#plot(s)

#sTmpFiles <- showTmpFiles()
#removeTmpFilesBut <- function(f) {
#  lfs <- showTmpFiles()
#  for( l in lfs ) {
#    if( l %in% f )
#      unlink( paste(tmpDir(), l, sep="") )
#  }
#}

#Define sampling strategies
samp_random <- function(r, p) {
  er <- extent(r)
  N <- r@ncols * r@nrows
  Np <- round(p*N)
  xyz <- data.frame(x=c(), y=c(), z=c())
  while( dim(xyz)[1] < Np ) {
    xy <- cbind( x=runif(Np-dim(xyz)[1], min=er@xmin, max=er@xmax),
                 y=runif(Np-dim(xyz)[1], min=er@ymin, max=er@ymax) )
    z <- extract(r, SpatialPoints(xy, proj4string=crs(r)))
    m <- !is.na(z)
#    cat("sum(m) = ", sum(m), "\n")
    if( any(m) ) xyz <- rbind( xyz, data.frame( x=xy[m, "x"], y=xy[m,"y"], z=z[m] ) )
  }
  
  xyz$ntr <- seq(1, dim(xyz)[1])
  return(xyz)
}

samp_randlow <- function(r, p, e=1) {
  er <- extent(r)
  N <- r@ncols * r@nrows
  Np <- round(p*N)
  xyz <- data.frame(x=c(), y=c(), z=c())
  h <- hist(r, breaks=100, plot=F)
  rejprob <- approxfun( h$mids, (1-cumsum(h$counts)/sum(h$counts))^e, rule=2 )
  while( dim(xyz)[1] < Np ) {
    xy <- cbind( x=runif(Np-dim(xyz)[1], min=er@xmin, max=er@xmax),
                 y=runif(Np-dim(xyz)[1], min=er@ymin, max=er@ymax) )
    z <- extract(r, SpatialPoints(xy, proj4string=crs(r)))
    m <- !is.na(z)
    q <- rejprob( z )
    m <- m & ( runif(length(q)) < q )
#    cat("sum(m) = ", sum(m), "\n")
    if( any(m) ) xyz <- rbind( xyz, data.frame( x=xy[m, "x"], y=xy[m,"y"], z=z[m] ) )
  }
  
  xyz$ntr <- seq(1, dim(xyz)[1])
  return(xyz)
}

minSide <- function(e) min( (e@xmax-e@xmin), (e@ymax-e@ymin) )
linspace <- function(a,b, by=1) {
 l <- hypot( b[1]-a[1], b[2]-a[2] )
 cbind(
   x=a[1]+(b[1]-a[1]) * seq(0,1,by=by/l),
   y=a[2]+(b[2]-a[2]) * seq(0,1,by=by/l) )
}
samp_transects <- function(r, p, l=minSide(extent(r))/10) {
  er <- extent(r)
  dxy <- ( (er@xmax-er@xmin)/ncol(s) + (er@ymax-er@ymin)/nrow(s) ) / 2
  N <- r@ncols * r@nrows
  Np <- round(p*N)
  ntr <- c()
  while( length(ntr) < Np) {
    t <- runif(1, min=0, max=pi)
    xyp <- c(x=runif(1, min=er@xmin, max=er@xmax), y=runif(1, min=er@ymin, max=er@ymax) )
    xyc <-linspace(xyp+c(-0.5*l*cos(t),-0.5*l*sin(t)), xyp+c(0.5*l*cos(t),0.5*l*sin(t)), by=dxy)
    cxy <- cellFromXY(r, xyc)
      cxy <- cxy[ !duplicated(cxy) & !is.na(cxy) ]
    
    if( n==1 ) {
      ntr <- rep(n, length(cxy))
      ccs <- cxy
    } else {
      ntr <- c(ntr, rep(n, length(cxy)))
      ccs <- c(ccs, cxy)
    }
  }

  xy <- xyFromCell(r,ccs)
  z  <- r[ccs]

  m <- !is.na(z)
  
  xyz <- data.frame(ntr=ntr[m], x=xy[m,"x"], y=xy[m,"y"], z=z[m] )
    
  return(xyz)
}

#
samp_partrans <- function(r, p, l=minSide(extent(r))/10, t=runif(1, min=0, max=pi)) {
  er <- extent(r)
  dxy <- ( (er@xmax-er@xmin)/ncol(s) + (er@ymax-er@ymin)/nrow(s) ) / 2
  N <- r@ncols * r@nrows
  Np <- round(p*N)
  ntr <- c()
  while( length(ntr) < Np ) {
    xyp <- c(x=runif(1, min=er@xmin, max=er@xmax), y=runif(1, min=er@ymin, max=er@ymax) )
    xyc <-linspace(xyp+c(-0.5*l*cos(t),-0.5*l*sin(t)), xyp+c(0.5*l*cos(t),0.5*l*sin(t)), by=dxy)
    cxy <- cellFromXY(r, xyc)
      cxy <- cxy[ !duplicated(cxy) & !is.na(cxy) ]
    
    if( n==1 ) {
      ntr <- rep(n, length(cxy))
      ccs <- cxy
    } else {
      ntr <- c(ntr, rep(n, length(cxy)))
      ccs <- c(ccs, cxy)
    }
  }

  xy <- xyFromCell(r,ccs)
  z  <- r[ccs]

  m <- !is.na(z)
  
  xyz <- data.frame(ntr=ntr[m], x=xy[m,"x"], y=xy[m,"y"], z=z[m] )
    
  return(xyz)
}

#Loop over different samplings and densities
for( t_ty in c("transpar") ) { # "random", "transpar", "randlow", "transects"
  
  #Test different point densities
  foreach( n_p=seq(8,4) ) %dopar% {
  
    #cat("Inside parloop n_p=", n_p, "\n")
    
    p=2^-n_p
    
    sId0 <- paste("sample", t_ty, n_p, sep="-")
    
    #Recall libraries in workers
    library(sp)
    library(raster)
    library(rgdal)

    #if all data has been simulated: do not recompute
    if( !file.exists(paste("data", paste(sId0,"-0",".csv", sep=""), sep="/")) ) {

      #Create random sample
      bdat <- switch(t_ty,
        random = samp_random(s, p),
        randlow = samp_randlow(s, p, e=0.5),
        transects = samp_transects(s, p),
        transpar = samp_partrans(s, p)
      )
      sbdat <- SpatialPointsDataFrame( bdat[,c("x","y")], dat=bdat, proj4string=CRS("+proj=utm +zone=20 +south +datum=WGS84 +units=m +no_defs") )

      #clean memory
      rm(list=c("bdat"))
    
      #transect statistics
      trStats <- data.frame( ntr=seq(1,max(sbdat$ntr,na.rm=T)) )
      trStats$npt <- NA
      trStats$lgt <- NA
      #compute number of points per transect and transect spanned length
      na <- 1
      nta <- sbdat$ntr[1]
      npt <- 0
      for(n in 1:dim(sbdat)[1]) {
        nt <- sbdat$ntr[n]
        if( nt != nta ) {
          #add the finished one
          trStats$npt[nta] <- npt
          if( npt > 1 ) {
            trStats$lgt[nta] <- hypot(
                 diff(range(sbdat@coords[sbdat$ntr[na:n-1],1])),
                 diff(range(sbdat@coords[sbdat$ntr[na:n-1],2])) )  
          } else {
            trStats$lgt[nta] <- 0
          }
        
          #next one
          na <- n
          npt <- 1
          nta <- nt
        } else {
          npt <- npt + 1 
        }
      }
      #add the last one
      trStats$npt[nta] <- npt
      if( npt > 1 ) {
        trStats$lgt[nta] <- hypot(
             diff(range(sbdat@coords[sbdat$ntr[na:n],1])),
             diff(range(sbdat@coords[sbdat$ntr[na:n],2])) )  
      } else {
        trStats$lgt[nta] <- 0
      }
    
      m <- !is.na(trStats$npt[ sbdat$ntr ])
      sbdat <- subset(sbdat, m )
      
      #save trStats
      write.csv( trStats, file=paste("data", paste(sId0, "_tranStats",".csv", sep=""), sep="/"), quote=F, row.names=F)
    
      #sbdat weights
      #sbdat$wgt.fr <- 1
      #sbdat$wgt.ffr <- 1
    
      #random labeling 
      sampLab <- floor( 1+foldK*runif(length(trStats$ntr)) )
      xyzw <- data.frame(x=sbdat@coords[,1], y=sbdat@coords[,2], z=sbdat$z, s=sampLab[sbdat$ntr]) #, w=sbdat$wgt.fr[m], ww=sbdat$wgt.ffr[m])
      write.csv( xyzw, file=paste("data", paste(sId0,"-0",".csv", sep=""), sep="/"), quote=F, row.names=F)
    } else {
      trStats <- read.csv( paste("data", paste(sId0, "_tranStats",".csv", sep=""), sep="/") )
      xyzw <- read.csv( paste("data", paste(sId0,"-0",".csv", sep=""), sep="/") )
      sbdat <- SpatialPointsDataFrame( xyzw[,c("x","y")], dat=xyzw, proj4string=CRS("+proj=utm +zone=20 +south +datum=WGS84 +units=m +no_defs") )
    }

    #garbage collect before parloop
    rm(list=c("xyzw"))
    gc(verbose=F)
  
    #cross validation interpolation
    for( n_s in seq(1,foldK) ) {
      
      #Sample identifier
      sId <- paste(sId0, n_s, sep="-")
      cat("sId =", sId, "\n")
      
      if( ! file.exists( paste("data", paste(sId, "-ffr", ".tif", sep=""), sep="/") ) ) {
        
        cat("K-folding data...\n")
        
        #Export data
        if( !file.exists(paste("data", paste(sId,"csv", sep="."), sep="/")) ) {
        
          if( "ntr" %in% names(sbdat) )
            m <- ( sbdat$ntr %in% trStats$ntr[ sampLab != n_s] ) #& !is.na(sbdat$wgt.fr)
          else
            m <- ( sbdat$s == n_s )

          xyzw <- data.frame(x=sbdat@coords[m,1], y=sbdat@coords[m,2], z=sbdat$z[m]) #, w=sbdat$wgt.fr[m], ww=sbdat$wgt.ffr[m])
          write.csv( xyzw, file=paste("data", paste(sId,"csv", sep="."), sep="/"), quote=F, row.names=F)
          rm(list=c("m", "xyzw"))
          gc(verbose=F)
        
        }
        
        #Run interpolation (call octave)
        cat(paste("octave --path='src' --eval 'sId=\"", sId, "\"; utmZ=19; utmH=-1; x_min=", bb@xmin,"; x_max=", bb@xmax,"; y_min=", bb@ymin,"; y_max=", bb@ymax,"; dxy=", dxy, "; test_srtm2'\n", sep=""))
        if( !file.exists(paste("data", paste(sId,"-ffr", ".asc", sep=""), sep="/")) )
          system(paste("octave --path='src' --eval 'sId=\"", sId, "\"; utmZ=19; utmH=-1; x_min=", bb@xmin,"; x_max=", bb@xmax,"; y_min=", bb@ymin,"; y_max=", bb@ymax,"; dxy=", dxy, "; test_srtm2'", sep=""))
        
        fr <- raster( paste("data", paste(sId, "-fr", ".asc", sep=""), sep="/"),
                      crs=readLines( paste("data", paste(sId, "-fr", ".crs", sep=""), sep="/") ) )
        fr <- crop(fr, bb)
        
        ss <- projectRaster( raster(paste("data", "srtm2.tif", sep="/")), fr)
        fr[ is.na(ss) ] <- NA
        
        writeRaster(fr, paste("data", paste(sId, "-fr", ".tif", sep=""), sep="/"), overwrite=T )
        
        rm(list="fr")
        
        ffr <- raster( paste("data", paste(sId, "-ffr", ".asc", sep=""), sep="/"),
                       crs=readLines( paste("data", paste(sId, "-ffr", ".crs", sep=""), sep="/") ) )
        ffr <- crop(ffr, bb)

        ss <- projectRaster( raster(paste("data", "srtm2.tif", sep="/")), ffr)
        ffr[ is.na(ss) ] <- NA
        
        writeRaster(ffr, paste("data", paste(sId, "-ffr", ".tif", sep=""), sep="/"), overwrite=T )
        
        rm(list=c("ffr", "ss"))
        
        #unlink( paste("data", paste(sId, "-fr", ".asc", sep=""), sep="/" ) )
        #unlink( paste("data", paste(sId, "-ffr", ".asc", sep=""), sep="/" ) )
        
        gc(verbose=F)
      }
    
    }
    
    #load stack of simple interpolations
    fr.ls <- paste("data", paste(paste(sId0, seq(1,foldK), sep="-"), "-fr", ".tif", sep=""), sep="/")
    frs <- stack(fr.ls)
    
    ss <- projectRaster( raster(paste("data", "srtm2.tif", sep="/")), frs)
      frs[ is.na(ss) ] <- NA
    
    fr.mrerr <- calc(frs-ss, fun=function (x) { sqrt(mean(x^2, na.rm=T)) } )
    writeRaster(fr.mrerr, paste("data", paste(sId0, "-fr_mrerr", ".tif", sep=""), sep="/"), overwrite=T )
    rm(list=c("fr.mrerr", "ss"))

    fr.cvmdn <- calc(frs, fun=function (x) { median(x, na.rm=T) } )
    writeRaster(fr.cvmdn, paste("data", paste(sId0, "-fr_cvmdn", ".tif", sep=""), sep="/"), overwrite=T)
    rm(list=c("fr.cvmdn"))
    
    fr.cverr <- calc(frs, fun=function (x) { sqrt(sum( (x-mean(x, na.rm=T))^2 )) } )
    writeRaster(fr.cverr, paste("data", paste(sId0, "-fr_cverr", ".tif", sep=""), sep="/"), overwrite=T )
    rm(list=c("fr.cverr"))
        
    #compute data weights with respect to cross-validation and set outlayers to NaN
    #sbdat$wgt.fr <- compWgtsF(z0=sbdat$depthc,
    #                          zm=extract(fr.cvmdn, sbdat), sz=extract(fr.cverr, sbdat),
    #                          dz0=0.5, rerr=0.15, fces=2.0)
    
    #garbage collect
    rm(list=c("fr.ls", "frs"))
    gc(verbose=F)
    
    #load stack of fractal interpolations
    ffr.ls <- paste("data", paste(paste(sId0, seq(1,foldK), sep="-"), "-ffr", ".tif", sep=""), sep="/")
    ffrs <- stack(ffr.ls)
    
    ss <- projectRaster( raster(paste("data", "srtm2.tif", sep="/")), ffrs)
      ffrs[ is.na(ss) ] <- NA
    
    ffr.mrerr <- calc(ffrs-ss, fun=function (x) { sqrt(mean(x^2, na.rm=T)) } )
    writeRaster(ffr.mrerr, paste("data", paste(sId0, "-ffr_mrerr", ".tif", sep=""), sep="/"), overwrite=T )
    rm(list=c("ffr.mrerr", "ss"))
    
    ffr.cvmdn <- calc(ffrs, fun=function (x) { median(x, na.rm=T) } )
    writeRaster(ffr.cvmdn, paste("data", paste(sId0, "-ffr_cvmdn", ".tif", sep=""), sep="/"), overwrite=T )
    rm(list=c("ffr.cvmdn"))
    
    ffr.cverr <- calc(ffrs, fun=function (x) { sqrt(sum( (x-mean(x, na.rm=T))^2 )) } )
    writeRaster(ffr.cverr, paste("data", paste(sId0, "-ffr_cverr", ".tif", sep=""), sep="/"), overwrite=T )
    rm(list=c("ffr.cverr"))
        
    #garbage collect
    rm(list=c("ffr.ls", "ffrs"))
    gc(verbose=F)
    
  }
  
}

#stop the parallel cluster
stopCluster(cl = pClust)

