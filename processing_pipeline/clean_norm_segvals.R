#!/usr/bin/env Rscript

# Step 3

args = commandArgs(trailingOnly=TRUE)


fileI = grep('file', args)
outI = grep('out', args)
probeI = grep('probes', args)
tileI = grep('tile=', args)

if (length(fileI) < 1 | length(outI) < 1)
  stop("Usage:\n
       \tfile=<segmented & fitted values REQ>\n
       \tout=<output dir REQ>\n 
       \tmin.probes=<DEF 67>\n
       \ttile=<DEF 5e6>\n
       ")

suppressPackageStartupMessages( source('~/workspace/shallowWGSpipeline/lib/data_func.R') )

args = sapply(args, strsplit, '=')

fitted_file = args[[fileI]][2] 
outdir = args[[outI]][2] 

#fitted_file='~/Data/Ellie/Analysis/multipcf_plots_fitted_perPatient/PR1_HIN_042/PR1_HIN_042_segmentedCoverage_fitted_gamma250.txt'
#outdir='~/Data/Ellie/Cleaned'


resid_file = grep('residuals.Rdata', list.files(dirname(fitted_file), full.names = T), value=T)
if (length(resid_file) <= 0)
  stop(paste("No residuals file to evaluate in ", dirname(fitted_file)))
load(resid_file, verbose=T)

residuals = sample.residual.variance(residuals)

samples = as.character(subset(residuals, varMAD_median <= 0.011)$sample)
message(paste(length(samples),'/',nrow(residuals), ' samples passed QC',sep=''))

outdir = paste(outdir, sub('\\..*', '', basename(fitted_file)), sep='/')

# metrics.file=paste(outdir, '/', dirname(fitted_file),'.out',sep='')
# if (file.exists(metrics.file)) {
#   mfile = file(metrics.file, 'a')
# } else {
#   mfile = file(metrics.file, 'w')
# }

print(paste("Writing to ", outdir, sep=''))

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive=T)
  message('Outdir created')
} else {
  message(paste(outdir, 'exists'))
}

min.probes=67 # equiv to 1MB
if (length(probeI) > 0)
  min.probes = as.numeric(args[[probeI]][2])

tile.w = 5e6
if (length(tileI) > 0)
  tile.w = as.numeric(args[[tileI]][2])

message(paste("Reading file:", fitted_file, '\nMin probes:', min.probes, ' Tile width:', tile.w, sep=''))

#write(paste("Reading file:", fitted_file, '\nMin probes:', min.probes, ' Tile width:', tile.w, sep=''), mfile, append=T)


gain.threshold = 1.1; loss.threshold = 0.9

if (!file.exists(fitted_file))
  stop(paste("File", fitted_file, "doesn't exist or is unreadable."))

segvals = as_tibble(read.table(fitted_file,sep="\t",stringsAsFactors=F,header=T))

# Select samples that passed QC
segvals = segvals[,c('chrom','arm','start.pos','end.pos','n.probes',samples)]

head(segvals)


probesCol = grep('probe', colnames(segvals), ignore.case=T)
probes = length(which(segvals[,probesCol] < min.probes))
message( paste(probes, ' probes (', round(probes/nrow(segvals),2)*100, '%)',' below the minimum probe count (',min.probes,')', sep=''))

#write( paste(probes, ' probes (', round(probes/nrow(segvals),2)*100, '%)',' below the minimum probe count (',min.probes,')', sep=''), mfile,append=T)

# get values
segvals = segvals[segvals[,probesCol] >= min.probes,]

res.summary<-function(df) {
  sm = signif(summary(df), 4)
  sm['Var.'] = signif(var(df), 4)
  sm
}

write.table(do.call(rbind, lapply(segvals[,-(1:5),drop=F], res.summary)), sep='\t', quote=F, col.names=NA, file=paste(outdir, 'probefiltered_summary.txt', sep='/'))

write.table(segvals, quote=F, sep='\t', row.names = F, file=paste(outdir, 'probefiltered_segvals.txt', sep='/'))

if (nrow(segvals) < 100) 
  warning(paste(fitted_file, 'has fewer than 100 genomic regions with segmented values:', nrow(segvals)))

sample.cols = intersect(grep('chr|arm|pos|probes', colnames(segvals), invert=T, ignore.case=T), which(sapply(segvals, is.numeric)))

chr.lengths = get.chr.lengths(file='/tmp/hg19_info.txt')

message(paste('genome coverage: ', round(sum(with(segvals, end.pos-start.pos))/chr.lengths[22,'genome.length']*100,1), '%', sep=''))

segs <- tile.segmented.data(data=segvals, size=tile.w, chr.info=chr.lengths)
arms <- tile.segmented.data(data=segvals, size='arms', chr.info=chr.lengths)

if (ncol(segs) != ncol(arms))
  stop("Segments and arms don't match.")

write.table(head(segs), row.names = F, col.names = T, quote=F)

write.table(segs, row.names=T, col.names=T, quote=F, sep='\t', file=paste(outdir, paste(tile.w,'_cleaned_tiled_segvals.txt', sep=''), sep='/'))
write.table(arms, row.names=T, col.names=T, quote=F, sep='\t', file=paste(outdir, 'arms_cleaned_tiled_segvals.txt', sep='/'))

print("Finished")
