
args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3)
  stop("Missing required params: <data dir> <out dir> <blacklist file> <overwrite: default F>")

suppressPackageStartupMessages(source('~/workspace/shallowWGSpipeline//lib/fastPCF.R'))
suppressPackageStartupMessages(source('~/workspace/shallowWGSpipeline//lib/data_func.R'))
suppressPackageStartupMessages(library(ggplot2))

chrlen = get.chr.lengths(file='/tmp/hg19_info.txt')

dir = args[1]
outdir = args[2]
exclude.file = args[3]
overwrite = F
if (length(args) == 4) overwrite =  as.logical(args[4])

dir = '~/Data/Ellie/Cassandra/20180409_KosmidouC_RF_sWGS/qdnaseq_15kb/'
#dir = '~/Data/Ellie/Cassandra/20180604_KosmidouC_RF_sWGS/qdnaseq'
#outdir = '~/Data/Ellie/Analysis/VAL_Cohort/sWGS'
#exclude.file = '~/Data/Ellie/QDNAseq/qDNAseq_blacklistedRegions.txt'

dir.create(outdir, F)
#f = 'LP2000110-DNA_A01.binSize15.fittedReadCounts.txt'
#r = 'LP2000110-DNA_A01.binSize15.readCounts.txt'

if ( !file.exists(exclude.file) )
  stop(paste("Missing necessary exclusion file from 'qDNAseq_blacklistedRegions.txt' in",dir))

blacklisted.regions = read.table(exclude.file,sep="\t",header=T,stringsAsFactors=F)

samples = unique(sub('\\.binSize.*', '', list.files(dir)))

existing = sub('_raw.txt','',list.files(outdir))
existing = grep('\\.Rdata', existing, invert=T, value=T)

if (!overwrite) 
  samples = setdiff(samples, existing)

plotdir = paste(outdir, 'plots',sep='/')
logdir = paste(outdir,'logs',sep='/')

for (p in c(plotdir,logdir)) dir.create(p, recursive=T, showWarnings=F)

info = NULL
for (i in 1:length(samples)) {
  sampleFiles = grep(samples[i], list.files(dir), value=T)
  print(samples[i])
  
  fit = grep('\\.(fitted|corrected)ReadCounts.txt', sampleFiles, value=T)
  raw = grep('\\.readCounts.txt', sampleFiles, value=T)
  
  raw.data = read.table(paste(dir, raw, sep='/'), header=T,stringsAsFactors=F)
  fit.data = read.table(paste(dir, fit, sep='/'), header=T,stringsAsFactors=F)
  head(raw.data)
  head(fit.data)
  ncol(raw.data) == ncol(fit.data)
  
  binned = binSWGS(raw.data, fit.data, blacklisted.regions, min.probes=67, plot.dir=plotdir, metrics.file=paste(logdir, '/', samples[i],'.out',sep=''))
  res = binned$segvals

  head(res)
  p = ggplot(res, aes(res[,6])) + geom_histogram(bins=20) + 
    geom_vline(xintercept = median(res[,6])) +
    geom_vline(xintercept = c(median(res[,6])+sd(res[,6]), 
                              median(res[,6])-sd(res[,6])), color='red') +
    geom_vline(xintercept = c(median(res[,6])+sd(res[,6])*2, 
                            median(res[,6])-sd(res[,6])*2), color='blue') + labs(title=colnames(res)[6])
  ggsave(filename=paste(plotdir, '/', colnames(res)[6], '-hist.png', sep=''), plot=p, height=5,width=5,units='in')
  
  rownames(binned$info) = samples[i]
  info = rbind(info, binned$info)

  write.table(res, sep='\t', quote=F, row.names=F, file=paste(outdir, paste(samples[i], '_raw.txt',sep=''), sep='/'))
}

write.table(round(info,4), sep='\t', quote=F, row.names=T, col.names=NA, file=paste(outdir, 'sample_info.txt',sep='/'))
print("Finished")