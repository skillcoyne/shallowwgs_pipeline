
args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3)
  stop("Missing required params: <data dir> <out dir> <blacklist file>")

suppressPackageStartupMessages(source('lib/fastPCF.R'))
suppressPackageStartupMessages(source('lib/data_func.R'))


chrlen = get.chr.lengths()

dir = args[1]
outdir = args[2]
exclude.file = args[3]

# dir = '~/Data/Ellie/QDNAseq/DownsampledBE/20180206_KillcoyneS_RF_BarrettsCN/qdnaseq/'
# outdir = '~/Data/Ellie/Analysis/downsampled5G_BEAdjacent/'
#exclude.file = '~/Data/Ellie/QDNAseq/qDNAseq_blacklistedRegions.txt'

dir.create(outdir, F)
#f = 'LP2000110-DNA_A01.binSize15.fittedReadCounts.txt'
#r = 'LP2000110-DNA_A01.binSize15.readCounts.txt'

if ( !file.exists(exclude.file) )
  stop(paste("Missing necessary exclusion file from 'qDNAseq_blacklistedRegions.txt' in",dir))

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