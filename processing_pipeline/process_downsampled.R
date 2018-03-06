
source('lib/fastPCF.R')
source('lib/data_func.R')

chrlen=get.chr.lengths()

dir = '~/Data/Ellie/QDNAseq/DownsampledBE/BEST2_NDBE'
outdir = '~/Data/Ellie/Analysis/downsampled5G_NDBE'
dir.create(outdir)
#f = 'LP2000110-DNA_A01.binSize15.fittedReadCounts.txt'
#r = 'LP2000110-DNA_A01.binSize15.readCounts.txt'

exclude.file = grep('qDNAseq_blacklistedRegions.txt', list.files('~/Data/Ellie',recursive=T,full.names=T), value=T)
if (length(exclude.file) <= 0)
  stop(paste("Missing necessary exclusion file from 'qDNAseq_blacklistedRegions.txt' in",data))

blacklisted.regions = read.table(exclude.file,sep="\t",header=T,stringsAsFactors=F)

samples = unique(sub('\\.binSize.*', '', list.files(dir)))

for (i in 1:length(samples)) {
  sampleFiles = grep(samples[i], list.files(dir), value=T)
  print(samples[i])
  
  fit = grep('\\.fittedReadCounts.txt', sampleFiles, value=T)
  raw = grep('\\.readCounts.txt', sampleFiles, value=T)
  
  raw.data = read.table(paste(dir, raw, sep='/'), header=T,stringsAsFactors=F)
  fit.data = read.table(paste(dir, fit, sep='/'), header=T,stringsAsFactors=F)
  head(raw.data)
  head(fit.data)
  ncol(raw.data) == ncol(fit.data)
  
  res = binSWGS(raw.data, fit.data, blacklisted.regions)
  write.table(res, sep='\t', quote=F, row.names=F, file=paste(outdir, paste(samples[i], '_raw.txt',sep=''), sep='/'))
}