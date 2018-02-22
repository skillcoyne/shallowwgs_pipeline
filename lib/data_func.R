require(GenomicRanges) 
require(dplyr)
require(data.table)
require(tibble)


prep.matrix<-function(dt) {
  output = list()
  
  means = apply(dt,2, mean)
  sd = apply(dt,2, sd)
  
  output[['z.mean']] = means
  output[['z.sd']] = sd
  
  dt = apply(dt, 2, unit.var)
  
  output[['matrix']] = dt
  return(output)
}

load.segment.matrix<-function(segFile) {
  dt = read.table(segFile, sep='\t', header=T)
  
  if (length(which(grepl('chr|start|end', colnames(dt)))) < 3)
    stop("Cannot load without chr, start, end columns")
  
  chrCol = grep('chr',colnames(dt))
  startCol = grep('start',colnames(dt))
  endCol = grep('end',colnames(dt))
  
  rows = paste(dt[,chrCol], ':', dt[,startCol], '-', dt[,endCol], sep='')
  cols = colnames(dt)[-c(chrCol, startCol, endCol)]
  
  dt = as.matrix(dt[,-c(chrCol, startCol, endCol)])
  rownames(dt) = rows
  colnames(dt) = cols
  return( t(dt) )
}


unit.var <- function(x, mean=NULL, sd=NULL) {
  if (is.null(mean) | is.null(sd)) {
    warning("No mean or sd provided.")
    if (length(x) == 1 | sd(x, na.rm=T) == 0) {
      warning("Unit normalization can't be performed with less than 2 samples or SD was 0")
      return(x)
    } else {
      return((x-mean(x,na.rm=T))/sd(x,na.rm=T) )
    }
  }
  
  uv = (x-mean)/sd 
  if (sd == 0) uv = 0
  
  return(uv)
}

score.cx <- function(df, MARGIN) {  
  cx = apply(df, MARGIN, function(x) length(which(x >= sd(x)*2 | x <= -sd(x)*2)))
  return(cx)
}

.tile.genome<-function(tile.w=5e6, chr.info=NULL, incGender=F) {
  
  if (is.null(chr.info))
    chr.info = read.table('data/chr_hg19.txt', sep='\t', header=T, stringsAsFactors=F)
  
  chr.info$chrom = sub('^chr', '', chr.info$chrom)
  
  chrs = c(1:22)
  if (incGender) chrs = c(1:22, 'X','Y')
  chr.info = subset(chr.info, chrom %in% chrs)
  chr.info$chrom = factor(chr.info$chrom, levels=chrs, ordered=T)
  
  chr.info$start = 1
  
  genome = makeGRangesFromDataFrame(chr.info, seqnames.field = 'chrom', end.field='chr.length')
  if (is.numeric(tile.w)) {
    tiles = tile(genome, width=tile.w)
  } else if (tile.w == 'arms') {
    parms = as.data.frame(chr.info %>% rowwise %>% dplyr::summarise(
      'chr'=chrom,
      'start'=1, 'end'=chr.cent-cent.gap,
      'arm'='p'))
    qarms = as.data.frame(chr.info %>% rowwise %>% dplyr::summarise(
      'chr'=chrom,
      'start'=chr.cent+cent.gap, 'end'=chr.length,
      'arm'='q'))
    genome = makeGRangesFromDataFrame(rbind(parms,qarms), seqnames.field = 'chr', keep.extra.columns = T)
    tiles = lapply( levels(seqnames(genome)), function(seq) genome[seqnames(genome)==seq] )
    tiles = GRangesList(tiles)
    if (length(tiles) != length(unique(seqnames(genome)))) stop("Tile and genome length doesn't match")
  } else if (tile.w == 'chr') {
    tiles = lapply(levels(seqnames(genome)), function(seq) genome[seqnames(genome)==seq] )
    tiles = GRangesList(tiles)
    if (length(tiles) != length(unique(seqnames(genome)))) stop("Tile and genome length doesn't match")
  }
  return(tiles)
}

# Presumes that the segmented data has already been filtered for num of probes, sd etc
tile.segmented.data<-function(data, size=5e6, chr.info=NULL) {
  
  if (!is.tibble(data))
    data = as_tibble(data)
  
  descCols = sort(union(grep('chr|arm|start|end|probes', colnames(data)), which(!sapply(data, is.numeric))))
  dataCols = c((descCols[length(descCols)]+1):ncol(data))
  
  chrCol = grep('chr',colnames(data),value=T)
  armCol = grep('arm',colnames(data),value=T)
  startPos = grep('start',colnames(data),value=T)
  endPos = grep('end',colnames(data),value=T)
  
  data = data[which(!data[[chrCol]] %in% c('X','Y')),]
  x1 = data[,c(chrCol, startPos, endPos, colnames(data)[dataCols])]
  
  tiles = .tile.genome(size, chr.info, incGender=length(which(grepl('X|Y', unique(data[[chrCol]])))) > 0)
  gr = makeGRangesFromDataFrame(x1, keep.extra.columns=T, start.field = startPos, end.field = endPos  ) 
  
  mergedDf = (do.call(rbind, lapply(tiles, function(tile) { 
    cbind('chr'=as.character(seqnames(tile)), as.data.frame(ranges(tile))[1:2]) 
  }) ))
  mergedDf[colnames(data)[dataCols]] = NA
  
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
      
      rows = with(mergedDf, which( chr==as.character(seqnames(bin)) & start == start(bin) & end == end(bin)))
      if (length(segments) > 1) 
        message(paste("chr", chr, "bin", bin, "has", length(segments), "matching segments"))
      
      meanSegs = c(meanSegs, length(segments))
      
      # weight means by the coverage of the bin  
      values = apply(as.data.frame(elementMetadata(segments)), 2, weighted.mean, w=weights)
      mergedDf[rows, names(values)] = values
    }
  }
  message(paste("Mean number of CN segments per genome bin:", round(mean(meanSegs, na.rm=T), 2), "median:", round(median(meanSegs, na.rm=T), 2)))
  
  mergedDf[is.na(mergedDf)] = 0
  
  return(as_tibble(mergedDf))
}

## Subtract arms from segs
subtract.arms<-function(segments, arms) {
  get.loc<-function(df) {
    locs = do.call(rbind.data.frame, lapply(colnames(df), function(x) unlist(strsplit( x, ':|-'))))
    colnames(locs) = c('chr','start','end')
    locs[c('start','end')] = lapply(locs[c('start','end')], function(x) as.numeric(as.character(x)))
    locs$chr = factor(locs$chr, levels=c(1:22), ordered=T)
    locs
  }
  
  if (is.null(segments) | is.null(arms))
    stop("Two matrices required")
  
  if (nrow(segments) != nrow(arms))
    stop(paste("Segment matrix cannot be adjusted by an arm matrix with different numbers of samples:", nrow(segments), ",", nrow(arms)))
  
  seg.loc = get.loc(segments)
  arm.loc = get.loc(arms)
  
  armsDF = makeGRangesFromDataFrame(arm.loc)
  segDF = makeGRangesFromDataFrame(seg.loc)
  
  tmp = segments
  # subtract arms from 5e6 and merge both
  ov = findOverlaps(armsDF, segDF)
  for (hit in unique(queryHits(ov))) {
    cols = subjectHits(ov)[which(queryHits(ov) == hit)]
    for (i in 1:nrow(tmp)) {
      tmp[i,cols] = tmp[i,cols] - arms[i,hit]
    }
  }
  mergedDf = cbind(tmp, arms)  
  #head(mergedDf)
  #summary(mergedDf[,1])
  
  #plot(apply(segs,2,mean))  
  #plot(apply(mergedDf, 2, mean))
  return(mergedDf)  
}



get.chr.lengths<-function(chrs = paste('chr', c(1:22, 'X','Y'), sep=''), build='hg19') {
  require(plyr)
  
  chr.lengths = read.table(paste('http://genome.ucsc.edu/goldenpath/help/', build, '.chrom.sizes',sep='') , sep='\t', header=F)
  colnames(chr.lengths) = c('chrom','chr.length')
  chr.lengths = subset(chr.lengths, chrom %in% chrs)
  chr.lengths$chrom = factor(chr.lengths$chrom, levels=chrs)
  chr.lengths = arrange(chr.lengths, chrom)
  
  cytoband.url = 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz'
  cytoband.file = paste('/tmp', basename(cytoband.url), sep='/')
  
  if (!file.exists(cytoband.file)) {
    success = download.file(cytoband.url, cytoband.file)
    if (!success)
      stop(paste("Could not download", cytoband.url))
  }
  
  cytobands = read.table(cytoband.file, sep='\t', header=F)
  colnames(cytobands) = c('chrom','start','end','band','attr')
  cytobands$chrom = factor(cytobands$chrom, levels=chrs)
  
  centromeres = subset(cytobands, attr == 'acen')
  
  chr.lengths$chr.cent = ddply(centromeres, .(chrom), summarise, cent=mean(range(start, end)) )$cent
  chr.lengths$cent.gap = ddply(centromeres, .(chrom), summarise, gap = (max(end)-min(start))/2 )$gap
  chr.lengths$genome.length = cumsum(as.numeric(chr.lengths$chr.length))
  
  return(chr.lengths)
}