options(bitmapType = "cairo")

library(ggplot2)
library(gridExtra)

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 2)
  stop("Missing data directory or summary file")

datadir = args[1]
summary.file = args[2]


if (!dir.exists(datadir))
  stop(paste(datadir, "doesn't exist or isn't readable"))

print(datadir)
print(summary.file)

smy = read.table(summary.file, sep = '\t', header = T)
smy = subset(smy, Copy.number <= 12)

a = as.data.frame(smy[,c('Copy.number','var.LRR')])
colnames(a) = c('x','y')
b = as.data.frame(smy[,c('Copy.number','median.LRR')])
colnames(b) = c('x','y')

fitSD = lm(y~x, data=a)
fitM = lm(y~x, data=b)

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

# lrr.rules<-function(medLRR, totalCN) {
#   LRR = jitter(medlrrs[as.character(totalCN)])
#   if (totalCN <= 1) jitter(medlrrs['1'], factor=50)
#   if (totalCN == 2) LRR = abs(medLRR)*2
#   if (totalCN > 2 & totalCN < 10) LRR = (abs(medLRR)+medlrrs[as.character(totalCN)])*totalCN^2
#   if (totalCN >= 10) LRR = jitter(medlrrs['10']*10, factor=20)
#   LRR
# }
#load('/fh/fast/reid_b/collab/Killcoyne/SNP_R/allpts_ascat.Rdata', verbose=T)

# all = do.call(rbind, segments.list)
# all = all %>% rowwise() %>% dplyr::mutate( total=nMajor+nMinor )
# y = all %>% group_by(total) %>% dplyr::summarise( median(medLRR), sd(medLRR), var(medLRR) )
# medlrrs = y$`median(medLRR)`
# names(medlrrs) = as.character(y$total)


chr.info = get.chr.lengths(file='hg19_info.txt')[1:22,]

if (is.null(chr.info)) stop("Failed to get chr info")
chr.info$chr = factor(sub('chr','',chr.info$chrom), levels=c(1:22), ordered = T)

#tf = list.files(datadir, 'tiled.txt', full.names = T)
#if (length(tf) ==2 && length(which(sapply(tf, file.size ) < 4025)) <= 0) {
#  message("Tiled files already exist")
#} else {
  load(list.files(datadir, 'Rdata', full.names = T), verbose=T)
  segraw = ascat.output$segments_raw
  rm(ascat.pcf, ascat.gg, ascat.output)
  
  segraw = subset(segraw, chr %in% c(1:22))
  segraw = segraw %>% rowwise() %>% dplyr::mutate( total = nMajor+nMinor )
  
  
  ## Adjust the sign of the log ratio so that CN gains result in a positive LRR.  This is more similar to what we get from sWGS
  # Also, anything that's a deletion peg to -1
  #segraw = segraw %>% rowwise() %>% dplyr::mutate( adjustedLRR = ifelse(nMajor+nMinor > 2, abs(medLRR), medLRR))
  #segraw = segraw %>% rowwise() %>% dplyr::mutate( adjustedLRR = ifelse(nMajor+nMinor > 2, abs(medLRR)*(nAraw+nBraw), medLRR))
  #segraw = segraw %>% rowwise() %>% dplyr::mutate( adjustedLRR = ifelse(nMajor+nMinor <= 1, jitter(-2, amount=1), medLRR))
  #segraw = subset(segraw, chr %in% c(1:22))
  #segraw = segraw %>% dplyr::mutate( totalRaw = (nAraw+nBraw)/2 )
  
  # segraw = segraw %>% rowwise() %>% dplyr::mutate(
  #   adjustedLRR = lrr.rules(medLRR, nMajor+nMinor)
  # )

  #head(segraw)
  ## Winsorize, per sample, the adjusted log ratio values
  #segraw = segraw %>% group_by(sample) %>% dplyr::mutate( winsLRR = madWins(adjustedLRR,2.5,25)$ywin )
  
  segraw$adjLRR = NA
  for (cn in unique(segraw$total)) {
    rows = which(segraw$total == cn)
    values = segraw[rows, 'medLRR', drop=T]
    
    newM = predict(fitM, newdata=cbind.data.frame('x'=cn))
    newSD = predict(fitSD, newdata=cbind.data.frame('x'=cn))
    
    segraw[rows,][['adjLRR']] = newM + (values - mean(values)) * (newSD/sd(values))
  }  

  #df = cbind.data.frame('chrom'=fit.data$chrom[good.bins], 'position'=fit.data$end[good.bins], 'seg.cov'=window.depths.standardised[good.bins,col])
  
  segraw$chr = factor(segraw$chr, levels=c(1:22), ordered=T)
  
  plotdir = paste(datadir,'/plots', sep='')
  dir.create(plotdir, showWarnings = F, recursive = T)
  
  m = (melt(segraw, measure.vars=c('medLRR','adjLRR')))
  ggsave(filename= paste(plotdir,'medLRR.png',sep='/'),
         plot=ggplot(m, aes(total, value, group=total)) + facet_grid(~variable) + geom_boxplot(), 
         width=9, height=7)

  allsamples = NULL
  allarms = NULL
  for (sample in unique(segraw$sample)) {
    print(sample)
    df = segraw[which(segraw$sample == sample),]
    
    df$winsLRR = madWins(df$adjLRR,2.5,25)$ywin

    tiled = tile.segmented.data(df[c('chr','startpos','endpos','winsLRR')], chr.info=chr.info, verbose=T)
    if (is.null(allsamples)) allsamples = tiled[c(1:3)]
    allsamples[,sample] = tiled[,4]
    
    tiled.arms = tile.segmented.data(df[c('chr','startpos','endpos','winsLRR')], size='arms', chr.info=chr.info, verbose=T)
    if (is.null(allarms)) allarms = tiled.arms[c(1:3)]
    allarms[,sample] = tiled.arms[,4]
  }

  plist = list()  
  for (smp in unique(segraw$sample)) {
    print(smp)
    df = as.data.frame(allsamples[,c('chr','start','end',smp)])
    colnames(df)[4] = 'CN'
    head(df)
    plist[[smp]] = ggplot(chr.info, aes(x=1:chr.length)) + facet_grid(~chr, space='free_x', scales='free_x') + ylim(-1,1) +
      geom_segment(data=df, aes(x=start, xend=end, y=CN, yend=CN), color='green4', lwd=3) +
      labs(title=smp, x='tiled') + theme_bw() + theme(axis.text.x=element_blank(), panel.spacing.x=unit(0,'lines'))
  }
  ggsave(filename = paste(plotdir,'tiled.png',sep='/'), plot = do.call(grid.arrange, c(plist, ncol=1)), width=10, height=25)
  
  write.table(allsamples, sep='\t', quote=F, row.names=F, file=paste(datadir,'/', basename(datadir), '_wins_tiled.txt', sep=''))
  write.table(allarms, sep='\t', quote=F, row.names=F, file=paste(datadir,'/', basename(datadir), '_wins_arms_tiled.txt', sep=''))
#}
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



