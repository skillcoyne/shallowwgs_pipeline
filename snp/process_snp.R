options(bitmapType = "cairo")

library(tidyverse)
library(gridExtra)


suppressPackageStartupMessages( source('~/workspace/shallowwgs_pipeline/lib/data_func.R') )

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3)
  stop("Usage: <data directory> <all patient ASCAT Rdata> <Combine by endoscopy T/F> <exlucde low SCA T/F")

chr.info = get.chr.lengths(file='/tmp/hg19_info.txt')

datadir = args[1]
datadir = '~/Data/Reid_SNP/PerPatient/1001'
# datadir =  "~/Data/Reid_SNP/PerPatient-preLM/512"
allpts.file = args[2]
# allpts.file = '~/Data/Reid_SNP/PerPatient/allpts_ascat.Rdata'
byEndo = as.logical(args[3])
# byEndo = T
excludeLowSCA = as.logical(args[4])
# excludeLowSCA = T

if (!file.exists(allpts.file))
  stop(paste(allpts.file, "doesn't exist or isn't readable"))

load(allpts.file, verbose=T)
sampletype<-function(level) {
  type = 'BE'
  if (grepl('BLD',level,ignore.case=T)) {
    type = 'Blood Normal'
  } else if (grepl('gastric',level,ignore.case=T)) {
    type = 'Gastric Normal'
  }
  return(type)
}

qcdata = qcdata %>% rowwise %>% dplyr::mutate(
  PatientID = unlist(strsplit(Samplename, '_'))[1],
  SampleID = unlist(strsplit(Samplename, '_'))[2],
  EndoID = unlist(strsplit(Samplename, '_'))[3],
  Level = unlist(strsplit(Samplename, '_'))[4], 
  Samplename = sub('_$','',Samplename), 
  SampleType = sampletype(Level)
)
segments.list = lapply(segments.list, function(pt) {
  pt$sample = sub('\\.LogR','',pt$sample)
  pt
})

qcdata$`ASCAT SCA Ratio` = apply(qcdata,1,function(s) {
  smp = subset(segments.list[[ s[['PatientID']] ]], sample == s[['Samplename']] & chr %in% c(1:22))
  smp = smp %>% rowwise %>% dplyr::mutate(
    #'Total' = nAraw + nBraw,
    'Total' = nMajor + nMinor,
    'CNV' = round(Total) - round(as.numeric(s[['Ploidy']])) )
  
  x = subset(smp, CNV != 0 & chr %in% c(1:22))
  sum(x$endpos - x$startpos) / chr.info[22,'genome.length',drop=T]
})

lowSCACutoff = median((qcdata %>% filter(SampleType != 'BE') %>% select(`ASCAT SCA Ratio`))$`ASCAT SCA Ratio`)


if (!dir.exists(datadir))
  stop(paste(datadir, "doesn't exist or isn't readable"))
print(datadir)

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


chr.info = chr.info[1:22,]
print(chr.info)
if (is.null(chr.info)) stop("Failed to get chr info")
chr.info$chr = factor(sub('chr','',chr.info$chrom), levels=c(1:22), ordered = T)

file = list.files(datadir,'^ascat.Rdata',full.names=T)
print(file)
load(file, verbose=T)
sample.ploidy = round(ascat.output$ploidy,2)
segraw = as_tibble(ascat.output$segments_raw)
segraw = segraw %>% rowwise %>% dplyr::mutate( sample = sub('\\.LogR','',sample))
names(sample.ploidy) =  sub('\\.LogR','',names(sample.ploidy))
rm(ascat.pcf, ascat.gg, ascat.output)

segraw = segraw %>% filter(chr %in% c(1:22)) %>% rowwise() %>% 
  dplyr::mutate(
    ploidy = sample.ploidy[sample],
    # Adjust based on the summary values per CN that we get in WGS (mostly) NP Barrett's cases
    totalRaw = round(nAraw+nBraw,3),
    adjRaw = totalRaw-ploidy, 
    chr = factor(chr, levels=c(1:22), ordered=T)
) %>% arrange(sample, chr, startpos)

if (excludeLowSCA) 
  qcdata = qcdata %>% filter(SampleType == 'BE' & Samplename %in% intersect(Samplename, segraw$sample) & `ASCAT SCA Ratio` > lowSCACutoff  )

m = (melt(segraw, measure.vars=c('adjRaw')))
p = ggplot(m, aes(sample, value, fill=grepl('BLD|gastric',sample))) + geom_jitter() + geom_boxplot(alpha=0.8) + labs(y='Ploidy Adj. CN', x='', title=basename(datadir)) + theme(legend.position = 'none', axis.text.x = element_text(angle=45, hjust=1))

p2 = ggplot(chr.info, aes(x=1:chr.length)) + facet_grid(sample~chr, space='free_y', scales='free_x') + 
  geom_segment(data=segraw, aes(x=startpos,xend=endpos,y=adjRaw,yend=adjRaw), color='darkgreen') + 
  #geom_hline(yintercept=c(median(segraw$adjRaw)+sd(segraw$adjRaw), median(segraw$adjRaw), median(segraw$adjRaw)-sd(segraw$adjRaw)), color='red') + 
  theme_bw() + theme(axis.text.x=element_blank(), panel.spacing.x=unit(0, 'lines')) 

plotdir = paste(datadir,'/plots', sep='')
if (dir.exists(plotdir)) unlist(plotdir, recursive=T)

dir.create(plotdir, showWarnings = F, recursive = T)
ggsave(filename= paste(plotdir,'rawCN.png',sep='/'), plot=p, width=length(unique(segraw$sample))*2, height=6,scale=1.5, limitsize = F)
ggsave(filename= paste(plotdir,'adjCN.png',sep='/'), plot=p2, height=length(unique(segraw$sample))*2, width=8,scale=1.5, limitsize = F)

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

  p = ggplot(chr.info, aes(x=1:chr.length)) + ylim(lims) + facet_grid(~chr,space='free_x',scales='free_x') + 
    geom_segment(data=tiled, aes(x=start,xend=end,y=adjRaw,yend=adjRaw), size=3, color='darkgreen') + 
    theme_bw() + theme(axis.text.x=element_blank(), panel.spacing.x=unit(0, 'lines')) + labs(x='', y=valueCol, title=names(rowsPerSample)[i])

  plist[[ names(rowsPerSample)[i] ]] = p
}

ggsave(filename = paste(plotdir,'tiled.png',sep='/'), plot = do.call(grid.arrange, c(plist, ncol=1)), width=10, height=5*length(plist), limitsize = F)

filename1 = paste(basename(datadir), '_wins_tiled',sep='') 
filename2 = paste(basename(datadir), '_wins_arms_tiled',sep='') 
if (excludeLowSCA) {
  filename1 = paste(filename1, '_exLOW',sep='') 
  filename2 = paste(filename2, '_exLOW',sep='') 
}
filename1 = paste(filename1, '.txt',sep='') 
filename2 = paste(filename2, '.txt',sep='') 

write.table(allsamples, sep='\t', quote=F, row.names=F, file=filename1)
write.table(allarms, sep='\t', quote=F, row.names=F, file=filename2)

print("Finished")

