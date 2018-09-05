options(bitmapType = "cairo")

library(ggplot2)
library(gridExtra)
library(dplyr)

#suppressPackageStartupMessages( source('lib/data_func.R') )
suppressPackageStartupMessages( source('~/workspace/shallowwgs_pipeline/lib/data_func.R') )

args = commandArgs(trailingOnly=TRUE)

#if (length(args) < 2)
#  stop("Missing data directory or summary file")

if (length(args) < 1)
  stop("Missing data directory ")


datadir = args[1]
#summary.file = args[2]

# datadir = '~/Data/Reid_SNP/PerPatient-preLM/512'


byEndo = F
if (length(args) > 2)
  byEndo = as.logical(args[3])

if (!dir.exists(datadir))
  stop(paste(datadir, "doesn't exist or isn't readable"))

print(datadir)
#print(summary.file)

#smy = read.table(summary.file, sep = '\t', header = T)
#smy = subset(smy, Copy.number <= 12)

# a = as.data.frame(smy[,c('Copy.number','var.LRR')])
# colnames(a) = c('CN','y')
# b = as.data.frame(smy[,c('Copy.number','mean.LRR')])
# colnames(b) = c('CN','y')
# 
# fitSD = lm(y~CN, data=a)
# fitM = lm(y~CN, data=b)
# #fitCN = lm(mean~CN, data=b)

#autoplot(fitM)

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


chr.info = get.chr.lengths(file='~/tmp/hg19_info.txt')[1:22,]
print(chr.info)
if (is.null(chr.info)) stop("Failed to get chr info")
chr.info$chr = factor(sub('chr','',chr.info$chrom), levels=c(1:22), ordered = T)

file = list.files(datadir,'^ascat.Rdata',full.names=T)
print(file)
load(file, verbose=T)
sample.ploidy = ascat.output$ploidy
segraw = ascat.output$segments_raw

rm(ascat.pcf, ascat.gg, ascat.output)

segraw = subset(segraw, chr %in% c(1:22))

head(segraw)

#segraw = segraw %>% rowwise() %>% dplyr::mutate( total = nMajor+nMinor )

# Adjust based on the summary values per CN that we get in WGS (mostly) NP Barrett's cases
adjust.segraw<-function(segraw, ploidy) {
  #segraw$adjLRR = NA
  segraw$totalRaw = round(with(segraw, nAraw+nBraw), 3)

  segraw = segraw %>% rowwise() %>% mutate(totalRaw = round(nAraw+nBraw, 3), adjRaw = totalRaw-ploidy[sample], ploidy = ploidy[sample] )
  
  #segraw$total = with(segraw, nMajor+nMinor)
  # for (cn in unique(segraw$total)) {
  #   rows = which(segraw$total == cn)
  #   values = segraw[rows, 'medLRR', drop=T]
  #   newSD = predict(fitSD, newdata=cbind.data.frame('CN'=cn))
  #   
  #   # if (cn == 2) {
  #   #   newM = 0
  #   # } else if (cn > 2) {
  #   #   newM = abs(newM)
  #   # }
  #   if (cn == 0) {
  #     newM = predict(fitM, newdata=cbind.data.frame('CN'=cn))
  #     segraw[rows,][['adjLRR']] = newM + (values - mean(values)) * (newSD/sd(values))
  #   } else if (cn == 2) {
  #     newSD = sd(values)
  #     newM = 0
  #     segraw[rows,][['adjLRR']] = newM + (values - mean(values)) * (newSD/sd(values))
  #   } else {
  #     newM = mean(values)
  #     segraw[rows,][['adjLRR']] = newM + (values - mean(values)) * (newSD/sd(values))
  #   }
  # }
  return(segraw)
}
  
segraw = adjust.segraw(segraw, sample.ploidy)

segraw$chr = factor(segraw$chr, levels=c(1:22), ordered=T)

m = (melt(segraw, measure.vars=c('adjRaw')))
p = ggplot(m, aes(sample, value, fill=grepl('BLD|gastric',sample))) + geom_jitter() + geom_boxplot(alpha=0.8) + labs(y='Ploidy Adj. CN', x='', title=basename(datadir)) + theme(legend.position = 'none', axis.text.x = element_text(angle=45, hjust=1))

ggplot(chr.info, aes(x=1:chr.length)) + facet_grid(sample~chr, space='free_y', scales='free_x') + 
  geom_segment(data=segraw, aes(x=startpos,xend=endpos,y=adjRaw,yend=adjRaw), color='darkgreen') + 
  #geom_hline(yintercept=c(median(segraw$adjRaw)+sd(segraw$adjRaw), median(segraw$adjRaw), median(segraw$adjRaw)-sd(segraw$adjRaw)), color='red') + 
  theme_bw() + theme(axis.text.x=element_blank(), panel.spacing.x=unit(0, 'lines')) 

# m = (melt(segraw, measure.vars=c('medLRR','adjLRR')))
# m$totalRaw = round(round(m$totalRaw,1))
# p = ggplot(m, aes(variable, value, group=variable, fill=variable, color=variable)) + facet_grid(~total) + geom_jitter(alpha=0.5) + geom_boxplot(outlier.colour = NA) + labs(title=basename(datadir))
plotdir = paste(datadir,'/plots', sep='')
if (dir.exists(plotdir)) unlist(plotdir, recursive=T)

dir.create(plotdir, showWarnings = F, recursive = T)
ggsave(filename= paste(plotdir,'rawCN.png',sep='/'), plot=p, width=length(unique(segraw$sample))*2, height=6,scale=1.5, limitsize = F)
save(segraw, file=paste(datadir, 'segments_raw.Rdata',sep='/'))

allsamples = NULL
allarms = NULL

samples = unique(segraw$sample)

info = do.call(rbind.data.frame, c(strsplit(samples, '_|\\.'), stringsAsFactors=F))[,1:4]
colnames(info) = c('PatientID','Sample','EndoID','Level')
rownames(info) = samples

rowsPerSample = lapply(samples, grep, segraw$sample)
names(rowsPerSample) = samples
print(samples)
if (byEndo) {
  normal = grep( paste(info$PatientID[1],'.*(BLD|Gastric)', sep=''), segraw$sample, ignore.case=T)
  endoMatch = c( paste(apply(unique(info[,c(1,3)]), 1, paste, collapse='_.*_'), '_\\d+', sep=''  ) )
  print(endoMatch)
  rowsPerSample = lapply(endoMatch, grep, segraw$sample)
  names(rowsPerSample) = apply(unique(info[,c(1,3)]), 1, paste, collapse='_')
  
  if (length(normal) > 0) {
    nm = info[grep(paste(info$PatientID[1],'.*(BLD|Gastric)', sep=''), samples, ignore.case=T),]
		rowsPerSample[[ paste(nm$PatientID, toupper(nm$Level), sep='_') ]] = normal
  }
  rowsPerSample = rowsPerSample[ which(sapply(rowsPerSample, length) > 0) ]
}


chr.info$chr = factor(sub('chr', '', chr.info$chrom), levels=c(1:22))

allsamples = NULL; allarms = NULL
plist = list()  
for (i in 1:length(rowsPerSample)) {
  print( names(rowsPerSample)[i] )
  df = segraw[rowsPerSample[[i]],]
  
  #df$winsLRR = madWins(df$adjLRR,2.5,25)$ywin
  valueCol = 'adjRaw'
  
  tiled = tile.segmented.data(df[c('chr','startpos','endpos',valueCol)], chr.info=chr.info, verbose=F)
  if (is.null(allsamples)) allsamples = tiled[c(1:3)]
  allsamples[,names(rowsPerSample)[i]] = tiled[,4]
  
  tiled.arms = tile.segmented.data(df[c('chr','startpos','endpos',valueCol)], size='arms', chr.info=chr.info, verbose=F)
  if (is.null(allarms)) allarms = tiled.arms[c(1:3)]
  allarms[,names(rowsPerSample)[i]] = tiled.arms[,4]
  
  #lims = c(-1,1)
  #if (max(tiled[[valueCol]],na.rm=T) > 1) lims[2] = max(tiled$winsLRR,na.rm=T)
  #if (min(tiled$winsLRR,na.rm=T) < -1) lims[1] = min(tiled$winsLRR,na.rT)
  lims = range(tiled[[valueCol]], na.rm = T)
  if (valueCol == 'adjRaw') lims=c(0,4)
  
  tiled$chr = factor(tiled$chr, levels=chr.info$chr)

  p = ggplot(chr.info, aes(x=1:chr.length)) + ylim(lims) + facet_grid(~chr,space='free_x',scales='free_x') + geom_segment(data=tiled, aes(x=start,xend=end,y=tiled[,4],yend=tiled[,4]), size=3, color='darkgreen') + theme_bw() + theme(axis.text.x=element_blank(), panel.spacing.x=unit(0, 'lines')) + labs(x='', y=valueCol, title=names(rowsPerSample)[i])
  
  plist[[ names(rowsPerSample)[i] ]] = p
}

ggsave(filename = paste(plotdir,'tiled.png',sep='/'), plot = do.call(grid.arrange, c(plist, ncol=1)), width=10, height=5*length(plist), limitsize = F)

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



