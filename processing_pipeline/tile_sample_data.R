#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#args = c()
#args[1] = '~/Data/Ellie/Analysis/All_patients.Rdata'
#args[2] = '~/Data/Ellie/Analysis/tile_patients_All'
#args[3] = 'arms'

if (length(args) < 3) {
  stop("Rdata file contining patient data, output filename, and tile size (def 1e7) required")
}

print(args)

suppressPackageStartupMessages( library(GenomicRanges) )
#suppressPackageStartupMessages( library(plyr) )
suppressPackageStartupMessages( library(dplyr) )

source('lib/load_patient_metadata.R')

load(args[1], verbose=T)

if (!exists('patient.data') | is.null(patient.data) | length(patient.data) <= 0)
  stop("Patient data file containing no information")

tileFile = args[2]
tile.w = as.numeric(args[3])
if (is.na(tile.w))
  tile.w = args[3]

tileFile = paste(tileFile, '_', as.character(tile.w), '.Rdata', sep='')
message(paste('Writing data to', tileFile))

apclust.data<-function(segdata, samples) {
  samples = intersect(colnames(segdata), samples)
  x1 = segdata[,samples]
  
  if (length(samples) <= 1) { 
    x1 = as.data.frame(x1)
    colnames(x1) = samples
  } else {
    x1 = x1[, samples]
  }
  rownames(x1) = (segdata[,c(1:4)] %>%
                    rowwise() %>%
                    mutate(location=paste(paste(chrom, arm, sep=''), '.', start.pos, '-', end.pos, sep='')))$location
  return(x1)  
}

tile.sample.segments<-function(sample.data) {
  message("Tiling patient data...")
  
  grList = lapply(sample.data, function(df) {
    print(unique(df$info$Hospital.Research.ID))
    df$info = plyr::arrange(df$info, Endoscopy.Year, Pathology)
    
    x1 = apclust.data(df$seg.vals, df$info$Samplename)
    
    segnames = as.data.frame(do.call(rbind, sapply(rownames(x1), strsplit, '[p|q].|-')))
    colnames(segnames) = c('chr','start','end')
    segnames[c('start','end')] = lapply(segnames[c('start','end')], function(x) as.numeric(as.character(x)))
    
    pt = unique(df$info$Hospital.Research.ID)
    return( makeGRangesFromDataFrame(cbind(segnames,pt,x1), keep.extra.columns=T) )
  })
  
  sampleNames = unlist(sapply(sample.data, function(df) {
    df$info = plyr::arrange(df$info, Endoscopy.Year, Pathology)
    x1 = apclust.data(df$seg.vals, df$info$Samplename)
    colnames(x1)
  }))
  #double check tiling in patient data
  tile.segments<-function(tiles, gr, mergedSegments) {
    meanSegs = c()
    ov = findOverlaps(tiles, gr) # Each tile is a chromosome, so the overlaps are only per chromosome here
    for (chr in unique(queryHits(ov))) { # for each chromosome get overlaps
      currentChr = tiles[[chr]]
      curov = findOverlaps(currentChr, gr)
      
      for (i in unique(queryHits(curov))) {
        bin = currentChr[i]

        segments = gr[subjectHits(curov[queryHits(curov) == i])]
        weights = lapply(segments,  function(r) {
          width(pintersect(bin, r))/width(bin)
        })
        
        rows = with(mergedSegments, which( chr==as.character(seqnames(bin)) & start == start(bin) & end == end(bin)))
        if (length(segments) > 1) 
          message(paste("chr", chr, "bin", bin, "has", length(segments), "matching segments"))
        
        meanSegs = c(meanSegs, length(segments))
        
        # weight means by the coverage of the bin  
        values = apply(as.data.frame(elementMetadata(segments)[-1]), 2, weighted.mean, w=weights)
        mergedSegments[rows, names(values)] = values
      }
    }
    message(paste("Mean number of CN segments per genome bin:", round(mean(meanSegs, na.rm=T), 2), "median:", round(median(meanSegs, na.rm=T), 2)))
    return(mergedSegments)  
  }
  
  mergedDf = do.call(rbind, lapply(tiles, function(tile) { 
    cbind('chr'=as.character(seqnames(tile)), as.data.frame(ranges(tile))[1:2]) 
  }) )
  mergedDf[,sampleNames] = NA
  
  for (gr in grList) {
    print(unique(gr$pt))
    mergedDf = tile.segments(tiles, gr, mergedDf)
  }
  
  mergedDf[is.na(mergedDf)] = 0
  rownames(mergedDf) = with(mergedDf, paste(chr, ':', start, '-', end, sep=''))
  return(mergedDf)
}


chr.lengths = get.chr.lengths()
chr.lengths$chrom = sub('chr','',chr.lengths$chrom)
chr.lengths$start = 1

genome = makeGRangesFromDataFrame(chr.lengths[1:22,], seqnames.field = 'chrom', end.field='chr.length')
if (is.numeric(tile.w)) {
  tiles = tile(genome, width=tile.w)
} else if (tile.w == 'arms') {
  parms = as.data.frame(chr.lengths %>% rowwise %>% summarise(
    'chr'=chrom,
    'start'=1, 'end'=chr.cent-cent.gap,
    'arm'='p'))

  qarms = as.data.frame(chr.lengths %>% rowwise %>% summarise(
    'chr'=chrom,
    'start'=chr.cent+cent.gap, 'end'=chr.length,
    'arm'='q'))

  genome = makeGRangesFromDataFrame(rbind(parms,qarms), seqnames.field = 'chr', keep.extra.columns = T)
  tiles = lapply(unique(seqnames(genome)), function(seq) genome[seqnames(genome)==seq] )
  tiles = GRangesList(tiles)
  if (length(tiles) != length(unique(seqnames(genome)))) stop("Tile and genome length doesn't match")
} else if (tile.w == 'chr') {
  tiles = lapply(unique(seqnames(genome)), function(seq) genome[seqnames(genome)==seq] )
  tiles = GRangesList(tiles)
  if (length(tiles) != length(unique(seqnames(genome)))) stop("Tile and genome length doesn't match")
}

mergedDf = tile.sample.segments(patient.data)
dim(mergedDf)

save(mergedDf, file=tileFile)
