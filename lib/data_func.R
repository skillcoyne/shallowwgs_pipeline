
## Should NAs be replace by per-sample mean value instead??
prep.matrix<-function(dt, na.replace=0) {
  output = list()
  
  means = apply(dt,2, mean, na.rm=T)
  sd = apply(dt,2, sd, na.rm=T)
  
  output[['z.mean']] = means
  output[['z.sd']] = sd
  
  dt = apply(dt, 2,  BarrettsProgressionRisk:::unit.var)
  
  if (!is.null(na.replace))
    dt[is.na(dt)] = na.replace

  output[['matrix']] = dt
  return(output)
}

segment.matrix<-function(dt, na.replace='mean',samplenames=NULL) {
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
  if (!is.null(samplenames) & length(samplenames) == ncol(dt))
    colnames(dt) = samplenames
  
  dt = t(dt)
  
  if (na.replace == 'mean') {
    dt[is.na(dt)] = mean(dt,na.rm=T)
  } else if (na.replace == 'median') {
    dt[is.na(dt)] = median(dt,na.rm=T)
  } else {
    dt[is.na(dt)] = 0
  }

  return( dt )
}

load.segment.matrix<-function(segFile) {
  dt = read.table(segFile, sep='\t', header=T)
  
  if (length(which(grepl('chr|start|end', colnames(dt)))) < 3)
    stop("Cannot load without chr, start, end columns")
  
 return(segment.matrix(dt))
}



get.chr.lengths<-function(chrs = paste('chr', c(1:22, 'X','Y'), sep=''), build='hg19', file=NULL) {
  require(plyr)
  
  if (!is.null(file) && file.exists(file)) {
    message("Reading from local file")
    chr.lengths = read.table(file, header=T, sep='\t', stringsAsFactors=F)
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
  write.table(chr.lengths, sep='\t', row.names=F, file=paste('/tmp/', build, '_info.txt', sep=''))

  return(chr.lengths)
}


