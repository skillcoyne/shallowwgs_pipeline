
args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3)
  stop("Missing required params: <data dir> <out dir> <blacklist file> <overwrite: default F>")

suppressPackageStartupMessages(source('../lib/fastPCF.R'))
suppressPackageStartupMessages(source('../lib/data_func.R'))
suppressPackageStartupMessages(library(ggplot2))

chrlen = get.chr.lengths()

dir = args[1]
outdir = args[2]
exclude.file = args[3]
overwrite = F
if (length(args) == 4) overwrite =  as.logical(args[4])

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
  
  res = binSWGS(raw.data, fit.data, blacklisted.regions, min.probes=67, plot.dir=plotdir, metrics.file=paste(logdir, '/', samples[i],'.out',sep=''))

  write.table(res, sep='\t', quote=F, row.names=F, file=paste(outdir, paste(samples[i], '_raw.txt',sep=''), sep='/'))
}
print("Finished")