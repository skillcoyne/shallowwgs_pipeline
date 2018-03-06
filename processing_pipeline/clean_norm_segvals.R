#!/usr/bin/env Rscript

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

suppressPackageStartupMessages( source('../lib/data_func.R') )

args = sapply(args, strsplit, '=')

fitted_file = args[[fileI]][2] 
outdir = args[[outI]][2] 

outdir = paste(outdir, sub('\\..*', '', basename(fitted_file)), sep='/')
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

print(paste("Reading file:", fitted_file, '\nMin probes:', min.probes, ' Tile width:', tile.w, sep=''))

gain.threshold = 1.1; loss.threshold = 0.9

if (!file.exists(fitted_file))
  stop(paste("File", fitted_file, "doesn't exist or is unreadable."))

segvals = as_tibble(read.table(fitted_file,sep="\t",stringsAsFactors=F,header=T))
head(segvals)

probesCol = grep('probe', colnames(segvals), ignore.case=T)
probes = length(which(segvals[,probesCol] < min.probes))
message( paste(probes, ' probes (', round(probes/nrow(segvals),2)*100, '%)',' below the minimum probe count (',min.probes,')', sep=''))

# get values
segvals = segvals[segvals[,probesCol] >= min.probes,]

write.table(segvals, quote=F, sep='\t', row.names = F, file=paste(outdir, 'raw_probefiltered_segvals.txt', sep='/'))


if (nrow(segvals) < 100) 
  warning(paste(fitted_file, 'has fewer than 100 genomic regions with segmented values:', nrow(segvals)))

sample.cols = intersect(grep('chr|arm|pos|probes', colnames(segvals), invert=T, ignore.case=T), which(sapply(segvals, is.numeric)))

if (length(sample.cols) > 0) {
  if (length(sample.cols) > 2) {  
    #HC = hclust(dist(t(segvals[,sample.cols])))
    #gg1 = ggdendrogram(HC, rotate=T) + labs(title = paste(basename(fitted_file),": fitted segment coverage"),x="")
    #print(gg1)
  }
  
  # Normalize  (value-mean(value))/sd(value)
  normalised.segvals = segvals[,sample.cols]
  if (length(sample.cols) > 1) {  
    for(c in 1:nrow(normalised.segvals)) 
      normalised.segvals[c,] = (normalised.segvals[c,]-mean(unlist(normalised.segvals[c,])))/sd(unlist(normalised.segvals[c,]))
    
    #HC = hclust(dist(t(normalised.segvals)))
    #gg2 = ggdendrogram(HC, rotate=T) + labs(title=paste(basename(fitted_file),": normalised fitted segment coverage", x=""))
    
    #grid.arrange(gg1, gg2, ncol=2)
  }
  
}

write.table(normalised.segvals, quote=F, sep='\t', row.names = F, file=paste(outdir, 'normalized_probefiltered_segvals.txt', sep='/'))


chr.lengths = get.chr.lengths()

segs <- tile.segmented.data(data=segvals, size=tile.w, chr.info=chr.lengths)
arms <- tile.segmented.data(data=segvals, size='arms', chr.info=chr.lengths)

if (ncol(segs) != ncol(arms))
  stop("Segments and arms don't match.")

head(segs)
write.table(head(segs), row.names = F, col.names = T, quote=F)

write.table(segs, row.names=T, col.names=T, quote=F, sep='\t', file=paste(outdir, paste(tile.w,'_cleaned_tiled_segvals.txt', sep=''), sep='/'))
write.table(arms, row.names=T, col.names=T, quote=F, sep='\t', file=paste(outdir, 'arms_cleaned_tiled_segvals.txt', sep='/'))

print("Finished")
