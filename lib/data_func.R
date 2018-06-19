require(GenomicRanges) 
require(dplyr)
require(data.table)
require(tibble)


#source("fastPCF.R")
logTV<-function(x) {
  x = log2(x)
  x[x==-Inf] = 0
  x
}

logTransform<-function(x,inf=F) {
  if (is.matrix(x) | is.data.frame(x)) {
    x[] = apply(x, 2, logTV)
  } else {
    x = logTV(x)
  }
  return(x)
}
  

binSWGS<-function(raw.data, fit.data, blacklist, logTransform=F, min.probes=67, plot.dir='.', metrics.file='binSWGS.out') {
  require(copynumber)
  message(paste('Logging to',metrics.file))
  
  if (file.exists(metrics.file)) {
    mfile = file(metrics.file, 'a')
  } else {
    mfile = file(metrics.file, 'w')
  }
  
  if (is.null(blacklist)) 
    stop('QDNASeq blacklisted regions missing')
    
  countCols = grep('loc|feat|chr|start|end', colnames(fit.data), invert=T)
  
  write(paste('fit data ',paste(dim(fit.data), collapse=':')), mfile, append=T)
  write(paste('raw data ',paste(dim(raw.data), collapse=':')), mfile, append=T)

  rows = which(fit.data[,1] %in% raw.data[,1])
  #raw.data = raw.data[rows,]
  fit.data = fit.data[rows,]

  window.depths = as.vector(as.matrix(raw.data[,countCols]))/as.vector(as.matrix(fit.data[,countCols]))
  
  #window.depths = raw.data[,countCols]/fit.data[,countCols]
  # QDNAseq does this in 'correctBins' but it's too late to use this now, so will adjust ours the same way (we weren't 14-May-2018)
  #corrected <- counts / fit
  #corrected[fit <= 0] <- 0
  negs = which(as.vector(as.matrix(fit.data[,countCols])) <= 0)
  if (length(negs) > 0) window.depths[negs] = 0
  
  window.depths = matrix(window.depths, ncol=length(countCols))  
  
  
  message(paste(nrow(blacklist), "genomic regions in the exclusion list."))
  write(paste(nrow(blacklist), "genomic regions in the exclusion list."),mfile,append=T)

  fit.data$in.blacklist = F
  for(r in 1:nrow(blacklist)) {
    fit.data$in.blacklist[ fit.data$chrom == blacklist$chromosome[r] & 
                           fit.data$start >= blacklist$start[r] & 
                           fit.data$end <= blacklist$end[r] ] = T
  }
  message(paste("# excluded probes = ",sum(fit.data$in.blacklist),sep=""))
  write(paste("# excluded probes = ",sum(fit.data$in.blacklist),sep=""),mfile,append=T)
  
  if (sum(fit.data$in.blacklist) <= 0)
    stop("There's an issue with the excluded list, no probes excluded.")
  
  if (length(countCols) == 1)
    window.depths = as.data.frame(window.depths)
  
  window.depths.standardised = as.data.frame(window.depths[which(!fit.data$in.blacklist),])
  if (logTransform)
    window.depths.standardised = log2( window.depths.standardised+abs(min(window.depths.standardised, na.rm=T))+1 )

  fit.data = fit.data[!fit.data$in.blacklist,-ncol(fit.data)]
  sdevs = sapply(c(1:length(countCols)), function(s) {
    getMad( window.depths.standardised[!is.na(window.depths.standardised[,s]),s], k=25 )
  })
  sdevs[sdevs==0] = NA
  sdev = exp(mean(log(sdevs[!is.na(sdevs)]))) 
  
  write(paste('sdev=',signif(sdev,3),sep=''),mfile,append=T)
  
  good.bins = which(!is.na(rowSums(as.data.frame(window.depths.standardised[,!is.na(sdevs)]))))
  gamma2=250 

  data = cbind(fit.data[good.bins,c('chrom','start')],window.depths.standardised[good.bins,!is.na(sdevs)])
  
  if (ncol(data) < 4) { # Single sample
    write("Segmenting single sample",mfile,append=T)
    res = pcf( data=data, gamma=gamma2*sdev, fast=F, verbose=T, return.est=F)
    colnames(res)[grep('mean', colnames(res))] = colnames(raw.data)[countCols]
    res$sampleID = NULL
  } else { # for most we have multiple samples
    write(paste("Segmenting", (ncol(data)-4), "samples"),mfile,append=T)
    res = multipcf( data=data, gamma=gamma2*sdev, fast=F, verbose=T, return.est=F)
  }
  
  if (!is.null(min.probes) & !is.na(min.probes)) {
    probes = which(res$n.probes < min.probes)
    message(paste(round(length(probes)/nrow(res), 3), 'segments with fewer than', min.probes, '"probes"'))
    write(paste(round(length(probes)/nrow(res), 3), 'segments with fewer than', min.probes, '"probes"'), mfile, append=T)
  }
  res = res[-probes,]
  
  chr.info = get.chr.lengths()
  ylims = c(floor(range(res[,-c(1:5)])[1]), ceiling(range(res[,-c(1:5)])[2]))
  for(col in which(!is.na(sdevs))) {
    sample = colnames(res)[col+5]
  
    df = cbind.data.frame('chrom'=fit.data$chrom[good.bins], 'position'=fit.data$end[good.bins], 'seg.cov'=window.depths.standardised[good.bins,col])
    df2 = cbind.data.frame('chrom'=res$chrom, 'start'=res$start.pos, 'end'=res$end.pos, 'seg.val'=res[,col+5])
    
    df$chrom = factor(df$chrom, levels=c(1:22), ordered=T)
    df2$chrom = factor(df2$chrom, levels=c(1:22), ordered=T)

    p = ggplot(chr.info, aes(x=1:chr.length)) + ylim(ylims) + facet_grid(~chrom, space='free_x', scales='free_x') +
      geom_point(data=df, aes(x=position, y=seg.cov), color='darkred', alpha=.4) +
      geom_segment(data=df2, aes(x=start, xend=end, y=seg.val, yend=seg.val), color='green4', lwd=5) +
      labs(y='segmented coverage', title=sample) + theme_bw() + theme(axis.text.x=element_blank(), panel.spacing.x=unit(0,'lines'))

    res$chrom = factor(res$chrom, levels=c(1:22), ordered=T)
    ggsave(filename=paste(plotdir, paste(sample,'_segmented.png',sep=''), sep='/'), plot=p, height=8, width=10)
  }
  close(mfile)
  return(res)
}

## Should NAs be replace by per-sample mean value instead??
prep.matrix<-function(dt, na.replace=0) {
  output = list()
  
  means = apply(dt,2, mean, na.rm=T)
  sd = apply(dt,2, sd, na.rm=T)
  
  output[['z.mean']] = means
  output[['z.sd']] = sd
  
  dt = apply(dt, 2, unit.var)
  
  if (!is.null(na.replace))
    dt[is.na(dt)] = na.replace

  output[['matrix']] = dt
  return(output)
}

segment.matrix<-function(dt) {
  dt = as.data.frame(dt)
  
  chrCol = grep('chr',colnames(dt),ignore.case=T)
  startCol = grep('start',colnames(dt),ignore.case=T)
  endCol = grep('end',colnames(dt),ignore.case=T)
  
  if (length(chrCol)<=0 | length(startCol)<=0 | length(endCol)<=0)
    stop("Missing columns chr|start|end")
  
  rows = paste(dt[,chrCol], ':', dt[,startCol], '-', dt[,endCol], sep='')
  cols = colnames(dt)[-c(chrCol, startCol, endCol)]
  
  dt = as.matrix(dt[,-c(chrCol, startCol, endCol)])
  rownames(dt) = rows
  colnames(dt) = cols
  return( t(dt) )
}

load.segment.matrix<-function(segFile) {
  dt = read.table(segFile, sep='\t', header=T)
  
  if (length(which(grepl('chr|start|end', colnames(dt)))) < 3)
    stop("Cannot load without chr, start, end columns")
  
 return(segment.matrix(dt))
}

unit.var <- function(x, mean=NULL, sd=NULL, warn=F) {
  if ((is.null(mean) | is.null(sd)) || (is.na(mean) | is.na(sd))) {
    warning("No mean or sd provided.")
    if (length(x) == 1 | length(which(is.na(x))) == length(x) | sd(x, na.rm=T) == 0) {
      if (warn) warning("Unit normalization can't be performed with less than 2 samples or SD was 0")
      return(x)
    } else {
      return( (x-mean(x,na.rm=T))/sd(x,na.rm=T) )
    }
  } else {
    uv = (x-mean)/sd 
    if (sd == 0) uv = 0
    return(uv)
  }
}

score.cx <- function(df, MARGIN) {  
  cx = apply(df, MARGIN, function(x) {
    length(which(x >= sd(x,na.rm=T)*2 | x <= -sd(x,na.rm=T)*2))
  })
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
tile.segmented.data<-function(data, size=5e6, chr.info=NULL, verbose=F) {
  if (!is.tibble(data)) data = as_tibble(data)
  
  descCols = sort(union(grep('chr|arm|start|end|probes', colnames(data), ignore.case=T), which(!sapply(data, is.numeric))))
  dataCols = c((descCols[length(descCols)]+1):ncol(data))
  
  chrCol = grep('chr',colnames(data),ignore.case=T, value=T)
  armCol = grep('arm',colnames(data),ignore.case=T,value=T)
  startPos = grep('start',colnames(data),ignore.case=T,value=T)
  endPos = grep('end',colnames(data),ignore.case=T,value=T)
  
  data = data[which(!data[[chrCol]] %in% c('X','Y')),]
  x1 = data[,c(chrCol, startPos, endPos, colnames(data)[dataCols])]
  
  tiles = .tile.genome(size, chr.info, incGender=length(which(grepl('X|Y', unique(data[[chrCol]])))) > 0)
  gr = makeGRangesFromDataFrame(x1, keep.extra.columns=T, start.field=startPos, end.field=endPos  ) 
  
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
      weights = sapply(as(segments,'GRangesList'),function(r) {
        width(pintersect(bin, r))/width(bin)
      })
      
      rows = with(mergedDf, which( chr==as.character(seqnames(bin)) & start == start(bin) & end == end(bin)))
      if (length(segments) > 1 & verbose) 
        message(paste("chr", chr, "bin", bin, "has", length(segments), "matching segments"))
      
      meanSegs = c(meanSegs, length(segments))
      
      # weight means by the coverage of the bin  
      values = apply(as.data.frame(elementMetadata(segments)), 2, weighted.mean, w=weights, na.rm=T)
      mergedDf[rows, names(values)] = values
    }
  }
  message(paste("Mean number of CN segments per genome bin:", round(mean(meanSegs, na.rm=T), 2), "median:", round(median(meanSegs, na.rm=T), 2)))
  
  # Not sure if this should be NA or 0
  #mergedDf[is.na(mergedDf)] = 0
  
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
  # subtract arms from 5e6 and merge both (per sample)
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

get.chr.lengths<-function(chrs = paste('chr', c(1:22, 'X','Y'), sep=''), build='hg19', file=NULL) {
  require(plyr)
  
  if (!is.null(file) && file.exists(file)) {
    message("Reading from local file")
    chr.lengths = read.table(file, header=T, sep='\t')
  } else {
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
  }
  return(chr.lengths)
}


modelMatrixSetup <- function(rs, chr.info, segMeans=NULL, segSDs=NULL, armMeans=NULL, armSDs=NULL) {
  segs <- tile.segmented.data(rs, size=5e6, chr.info=chrlen)
  segM = as.matrix(segs[,-c(1:3)])
  rownames(segM) = paste(segs$chr, ':', segs$start, '-', segs$end, sep='')
  segM = t(segM)
  
  if (!is.null(segMeans) && !is.null(segSDs)) {
    for (i in 1:ncol(segM)) 
      segM[,i] = unit.var(segM[,i], segMeans[i], segSDs[i])
  } else {
    segM = unit.var(segM)
  }
  
  arms <- tile.segmented.data(rs, size='arms', chr.info=chrlen)
  armsM = as.matrix(arms[,-c(1:3)])
  rownames(armsM) = paste(arms$chr, ':', arms$start, '-', arms$end, sep='')
  armsM = t(armsM)
  
  if (!is.null(armMeans) && !is.null(armSDs)) {
    for (i in 1:ncol(armsM)) 
      armsM[,i] = unit.var(armsM[,i], z.arms.mean[i], z.arms.sd[i])
  } else {
    armsM = unit.var(armsM)
  }
  
  nrow(armsM) == nrow(segM)
  
  cx.score = score.cx(segM,1)
  
  mergedDf = subtract.arms(segM, armsM)
  mergedDf = cbind(mergedDf, 'cx' = unit.var(cx.score, mn.cx, sd.cx))
  return(mergedDf)
}