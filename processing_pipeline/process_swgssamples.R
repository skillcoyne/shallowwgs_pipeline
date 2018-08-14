
args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3)
  stop("Missing required params: <data dir> <out dir> <blacklist file> <sample info file> <overwrite: default F>")

suppressPackageStartupMessages(source('~/workspace/shallowWGSpipeline//lib/fastPCF.R'))
suppressPackageStartupMessages(source('~/workspace/shallowWGSpipeline//lib/data_func.R'))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(readxl))

chrlen = get.chr.lengths(file='/tmp/hg19_info.txt')

dir = args[1]
outdir = args[2]
exclude.file = args[3]
info.file = args[4]

gamma2=250

overwrite = F
if (length(args) == 5) overwrite =  as.logical(args[5])

# dir = '~/Data/Ellie/QDNAseq/SLX-15638'
# outdir = '~/tmp/VAL'
# exclude.file = '~/Data/Ellie/QDNAseq/qDNAseq_blacklistedRegions.txt'
# info.file = '~/Data/Ellie/Cassandra/ValidationCohort.xlsx'

if ( !file.exists(exclude.file) )
  stop(paste("Missing necessary exclusion file 'qDNAseq_blacklistedRegions.txt'"))
if ( !file.exists(info.file) )
  stop(paste("Missing necessary sample info file"))

patient.info = read_xlsx(info.file, sheet='All')
colnames(patient.info) = gsub('\\s', '', colnames(patient.info))
patient.info$Samplename = with(patient.info, paste(sub('-', '_', IndexSequence), '_', `SLX-ID`, sep=''))
patient.info = arrange(patient.info, PatientID, `SLX-ID`, IndexSequence)

blacklisted.regions = read.table(exclude.file,sep="\t",header=T,stringsAsFactors=F)


files = list.files(dir, pattern='^(raw|fitted)ReadCounts.txt', full.names=T, recursive=T)
print(files)

if ( length(grep('rawReadCounts', files)) != 1 ) {
  message("Generating merged fitted/raw files")
  files = grep('SLX',list.files(dir, pattern='[raw|fitted]ReadCounts.txt', full.names=T, recursive=T),value=T)

  merged.raw = data.frame(matrix(ncol=4,nrow=0,dimnames=list(c(),c('location','chrom','start','end') ))); 
  merged.fitted = merged.raw
  
  mraw = grep('raw',files, value=T)
  mfitted = grep('fitted',files,value=T)

  fixSamplename<-function(nm) {
    nm = sub('\\.binSize.*','',nm)
    nm = sub('SLX[-|_]','',nm)
    nm2 = unlist(strsplit(nm,'\\.'))
    return(paste(nm2[2], nm2[1], sep='_'))
  }
    
  for (f in mraw) {
    print(f)
    rw = read.table(f,header=T,stringsAsFactors=F)
    colnames(rw)[grep('count',colnames(rw))] = fixSamplename(basename(f))
    
    ft = read.table(grep(sub('\\.binSize.*','',basename(f)), mfitted, value=T), header=T, stringsAsFactors=F)
    colnames(ft)[grep('count',colnames(ft))] = fixSamplename(basename(f))
    
    
    merged.raw = base::merge(merged.raw, rw, all=T)    
    merged.fitted = base::merge(merged.fitted, ft, all=T)    
    print(nrow(merged.raw))
  }
  merged.raw$chrom = factor(merged.raw$chrom, levels=c(1:22,'X','Y'))
  merged.fitted$chrom = factor(merged.fitted$chrom, levels=c(1:22,'X','Y'))
  
  merged.raw = arrange(merged.raw, chrom, start)
  merged.fitted = arrange(merged.fitted, chrom, start)
  
  write.table(merged.raw, sep='\t', quote=F, row.names=F, file=paste(dir, 'rawReadCounts.txt',sep='/') )
  write.table(merged.fitted, sep='\t', quote=F, row.names=F, file=paste(dir, 'fittedReadCounts.txt',sep='/') )
  files = c(paste(dir, 'rawReadCounts.txt',sep='/') , paste(dir, 'fittedReadCounts.txt',sep='/') )
}


fitted.file = grep('fitted', files, value=T)
raw.file = grep('raw',files,value=T)

message(paste("Loading raw file", raw.file))
raw.reads = read.table(raw.file, header=T, stringsAsFactors = F)
message(paste("Loading fitted file", fitted.file))
fitted.reads = read.table(fitted.file, header=T, stringsAsFactors = F)
# sort...
fitted.reads = fitted.reads[,colnames(raw.reads)]

cleanDataCols <- function(vec) {
  if (length(which(grepl('SLX',vec))) < 1) return(vec)
  vec = sub('SLX[\\.|-|_]', '', vec)  
  vec = gsub('tp', '', vec)
  
  vec2 = do.call(rbind, strsplit(vec,'\\.|-'))
  if ( length(unique(vec2[,1])) < length(vec)) {
    vec2 = paste(vec2[,2], vec2[,1], sep= '_')
  } else {
    vec2 = paste(vec2[,1], vec2[,2], sep= '_')
  }
  return(vec2)
}
colnames(fitted.reads)[-c(1:4)] = cleanDataCols(colnames(fitted.reads)[-c(1:4)])
colnames(raw.reads)[-c(1:4)] = cleanDataCols(colnames(raw.reads)[-c(1:4)])

head(fitted.reads)

for (pid in unique(patient.info$PatientID)) {
  plotdir = paste(outdir, pid,'plots',sep='/')
  dir.create(plotdir, recursive = T, showWarnings = F)
  logdir = paste(outdir,pid,'logs',sep='/')
  dir.create(logdir, recursive = T, showWarnings = F)
  
  message(paste('Patient',pid))
  info = subset(patient.info, PatientID == pid)
  noreads = grep('Insufficient Reads', info$QCIssues)
  if (length(noreads) > 0) {
    warning(paste('Patient',pid,':', length(noreads), 'samples have insufficient reads for segmentation.'))
    info = info[-noreads,]
    if (nrow(info) <= 0) next
  }
  
  samples = intersect(colnames(fitted.reads)[-c(1:4)], info$Samplename)
  if (length(samples) <= 0) next
  
  fit.data = fitted.reads[,c(colnames(fitted.reads)[1:4],samples)]
  raw.data = raw.reads[,c(colnames(raw.reads)[1:4],samples)]
  
  res.filename = paste(outdir, pid, paste('pid',info$PatientID[1], '_segmentedCoverage_fitted_gamma2_', gamma2, '.txt',sep=''), sep='/')
  segvals = NULL
  if (file.exists(res.filename) & file.size(res.filename)/1e3 > 1) { 
    message(paste("Reading segmented values from",res.filename))
    segvals = read.table(res.filename, sep='\t', header=T)
  }
  
  # metrics.file = paste(logdir, '/pid', info$PatientID[1],'.out',sep='')
  binned = binSWGS(raw.data, fit.data, blacklist=blacklisted.regions, min.probes=67, plot.dir=plotdir, gamma2=gamma2, metrics.file=paste(logdir, '/pid', info$PatientID[1],'.out',sep=''), seg.matrix = segvals)
  ggsave(filename=paste(plotdir, '/segmentedCoverage_pid', pid, '_gamma',gamma2,'.png',sep=''), plot=do.call(grid.arrange, c(binned$plots, ncol=1, top=paste('ID',pid))), height=4*length(samples), width=12)
  
  residuals = binned$resvar
  save(residuals, file=paste(outdir, pid, 'residuals.Rdata',sep='/'))
  
  write.table(binned$segvals, sep='\t', quote=F, row.names=F, file=res.filename)
}
#}

#write.table(seg.info, sep='\t', quote=F, row.names=F, file=paste(outdir, 'sample_info.txt',sep='/'))
#write.table(raw.summary, sep='\t', quote=F, row.names=F, file=paste(outdir, 'raw_data_summary.txt',sep='/'))
print("Finished")