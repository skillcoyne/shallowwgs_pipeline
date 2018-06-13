options(bitmapType = "cairo")

args = commandArgs(trailingOnly=TRUE)

datadir = args[1]
#datadir = '/Volumes/fh/fast/reid_b/collab/Killcoyne/Data/PerPatient/597'
print(datadir)

### Lifted directly from ASCAT
#Perform MAD winsorization:
madWins <- function(x,tau,k){
  xhat <- medianFilter(x,k)
  d <- x-xhat
  SD <- mad(d)
  z <- tau*SD
  xwin <- xhat + psi(d, z)
  outliers <- rep(0, length(x))
  outliers[x > xwin] <- 1
  outliers[x < xwin] <- -1
  return(list(ywin=xwin,sdev=SD,outliers=outliers))
}

psi <- function(x,z){
  xwin <- x
  xwin[x < -z] <- -z
  xwin[x > z] <- z
  return(xwin)
}

medianFilter <- function(x,k){
  n <- length(x)
  filtWidth <- 2*k + 1
  #Make sure filtWidth does not exceed n
  if(filtWidth > n){
    if(n==0){
      filtWidth <- 1
    }else if(n%%2 == 0){
      #runmed requires filtWidth to be odd, ensure this:
      filtWidth <- n - 1
    }else {
      filtWidth <- n
    }
  }
  runMedian <- runmed(x,k=filtWidth,endrule="median")
  return(runMedian)
}
### ---------------------- ###

suppressPackageStartupMessages( source('lib/data_func.R') )

chr.info = get.chr.lengths()
if (is.null(chr.info)) stop("Failed to get chr info")

tf = list.files(datadir, 'tiled.txt', full.names = T)
if (length(which(sapply(tf, file.size ) < 4025)) <= 0) {
  message("Tiled files already exist")
} else {
    
  
  load(list.files(datadir, 'Rdata', full.names = T), verbose=T)
  segraw = ascat.output$segments_raw
  rm(ascat.pcf, ascat.gg, ascat.output)
  
  ## Adjust the sign of the log ratio so that CN gains result in a positive LRR.  This is more similar to what we get from sWGS
  segraw = segraw %>% rowwise() %>% dplyr::mutate( adjustedLRR = ifelse(nMajor+nMinor > 2, abs(medLRR), medLRR))
  
  
  ## Winsorize, per sample, the adjusted log ratio values
  segraw = segraw %>% group_by(sample) %>% dplyr::mutate( winsLRR = madWins(adjustedLRR,2.5,25)$ywin )
  
  #qqnorm(subset(segraw, sample == segraw$sample[1])$adjustedLRR)
  #qqnorm(subset(segraw, sample == segraw$sample[1])$winsLRR)
  
  allsamples = NULL
  allarms = NULL
  for (sample in unique(segraw$sample)) {
    print(sample)
    df = segraw[which(segraw$sample == sample),]
    tiled = tile.segmented.data(df[c('chr','startpos','endpos','winsLRR')], chr.info=chr.info, verbose=T)
    if (is.null(allsamples)) allsamples = tiled[c(1:3)]
    allsamples[,sample] = tiled[,4]
    
    tiled.arms = tile.segmented.data(df[c('chr','startpos','endpos','winsLRR')], size='arms', chr.info=chr.info, verbose=T)
    if (is.null(allarms)) allarms = tiled.arms[c(1:3)]
    allarms[,sample] = tiled.arms[,4]
    head(allsamples)
  }
  
  write.table(allsamples, sep='\t', quote=F, row.names=F, file=paste(datadir,'/', basename(datadir), '_wins_tiled.txt', sep=''))
  write.table(allarms, sep='\t', quote=F, row.names=F, file=paste(datadir,'/', basename(datadir), '_wins_arms_tiled.txt', sep=''))
}
print("Finished")


# mtx = segment.matrix(allsamples)
# for (i in 1:ncol(mtx))
#   mtx[,i] = unit.var(mtx[,i], z.mean[i], z.sd[i])
# 
# armmtx = segment.matrix(allarms)
# for (i in 1:ncol(armmtx))
#   armmtx[,i] = unit.var(armmtx[,i], z.arms.mean[i], z.arms.sd[i])
# 
# cx = score.cx(mtx, 1)
# 
# 
# arrayDf = subtract.arms(mtx, armmtx)
# arrayDf = cbind(arrayDf, 'cx'=unit.var(cx, mn.cx, sd.cx))
# 
# predict(fitV, newx=arrayDf, s=l, type='response')
# predict(fitV, newx=arrayDf, s=l, type='link')



