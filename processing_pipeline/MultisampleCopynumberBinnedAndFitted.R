## Step 1

library(plyr)
library(ggplot2)
library(ggdendro)
# this version removes qDNAseq 'blacklisted' regions and uses 'fitted' read counts

#args = commandArgs(trailingOnly=TRUE)

#if (length(args) < 3)
#  stop("Missing required arguments: <qdna data dir> <patient spreadsheet> <output dir>")

source('lib/load_patient_metadata.R')

#data = args[1]
#patient.file = args[2]
#outdir = args[3]
data = '~/Data/Ellie/QDNAseq'
patient.file = paste(data, 'All_patient_info.xlsx', sep='/')
outdir = '~/Data/Ellie/Analysis'


all.patient.info = read.patient.info(patient.file, set='All')
if (is.null(all.patient.info))
  stop(paste("Failed to read patient file", patient.file))
all.patient.info$Patient = gsub('/', '_',all.patient.info$Patient)


outdir = paste(outdir, "plots_fitted", sep='/')
if ( !dir.exists(outdir) ) 
  dir.create(outdir)

data.dirs = list.dirs(data, full.names=T, recursive=F)
if (length(data.dirs) <= 0)
  stop(paste("No directories in", data))

# Load the data, binding each new data file into single data table
fit.data = NULL
raw.data = NULL
samplenames=NULL
for(dir in data.dirs) {
  print(paste("Reading directory", dir))
  files = list.files(dir,pattern="[fitted|corrected]ReadCounts.txt", full.names=T)
  slx = sub('SLX-','', basename(dir))
  
  # Files that do not contain individual indecies
  if (length(files) == 1 || !grepl('D\\d+_D\\d+', files)) {
    for (file in files) {
      print(paste("Reading file", file))
      spl = strsplit(basename(file),"\\.")
      
      fitted.file = file
      raw.file = list.files(dir, pattern='\\.(raw)?readCounts', full.names=T)
      if (length(raw.file) > 1)
        stop(paste("Error finding raw read counts in", dir))
      
      message( paste('FIT:', basename(fitted.file), '  RAW:', basename(raw.file)) )
      
      fit = read.table(fitted.file,sep="\t",header=T,stringsAsFactors=F)  # fitted reads
      f.cols = grep('D\\d', colnames(fit))
      colnames(fit)[f.cols] = paste(colnames(fit)[f.cols], slx, sep="_") 
      
      raw = read.table(raw.file,sep="\t",header=T,stringsAsFactors=F) # raw reads
      r.cols = grep('D\\d', colnames(raw))
      colnames(raw)[r.cols] = paste(colnames(raw)[r.cols], slx, sep="_") 
      raw = raw[,colnames(fit)]	    
      
      if (is.null(fit.data)) {
        fit.data = fit
        raw.data = raw
      } else {
        fit.data = cbind(fit.data,fit[,f.cols])
        raw.data = cbind(raw.data,raw[,f.cols])
      }
    }
  } else { # individual index
    for(file in files){
      print(paste("Reading file", file))
      spl = strsplit(basename(file),"\\.")
      
      fitted.file = file
      raw.file = sub('fitted', 'raw', fitted.file)
      
      message( paste('FIT:', basename(fitted.file), '  RAW:', basename(raw.file)) )
      
      if (!file.exists(fitted.file) || !file.exists(raw.file)) 
        stop("Missing fitted/raw file ")
      
      fit = read.table(fitted.file,sep="\t",header=T,stringsAsFactors=F)  # fitted reads
      colnames(fit)[5] = paste(spl[[1]][2],slx,sep="_")
      
      raw = read.table(raw.file,sep="\t",header=T,stringsAsFactors=F) # raw reads
      colnames(raw)[5] = paste(spl[[1]][2],slx,sep="_")
      
      if (is.null(fit.data)) {
        fit.data = fit
        raw.data = raw
      } else { # add to table
        fit.data = cbind(fit.data,fit[,5])
        raw.data = cbind(raw.data, raw[,5])
        
        colnames(fit.data)[ncol(fit.data)] = colnames(fit)[5]
        colnames(raw.data)[ncol(raw.data)] = colnames(raw)[5]
      }	  
    }
  }
  print(dim(fit.data))
}

#save.image(file=paste(data,"MultisampleCopynumberBinnedAndFitted1.RData", sep='/'))
save(fit.data, raw.data, file=paste(data,"MultisampleCopynumberBinnedAndFitted1.RData", sep='/'))

chrs = ddply(fit.data, .(chrom), summarise, length=max(end))
chrs$chrom = factor(chrs$chrom, c((1:22),"X","Y","M"), ordered=T)
chrs = chrs[order(chrs$chrom),]
chrs$cum.lengths = cumsum(as.numeric(chrs$length))
for (c in 1:nrow(chrs)) {
  chrs[c,'l.pos'] = ifelse(c==1, as.numeric(chrs$cum.lengths[c]/2), as.numeric((chrs$cum.lengths[c-1]+chrs$cum.lengths[c])/2))
}

# moving the position information from chromosomes to the genome position (cumulutive sum of the lengths of the chrs)
fit.data$genome.pos = as.numeric(fit.data$start)
for (c in 2:nrow(chrs)) {
  fit.data$genome.pos[fit.data$chrom == chrs[c, 'chrom']] = fit.data$genome.pos[fit.data$chrom == chrs[c, 'chrom']] + chrs[c-1, 'cum.lengths']   #cum.lengths[c-1]
}

# "count" columns
depth.cols = grep('^D\\d',colnames(fit.data))
no.samples = length(depth.cols)
# adjusting the raw counts by the fitted counts (need to find out what "fitted" is in this context)
window.depths = raw.data[,depth.cols]/fit.data[,depth.cols] 
# QDNAseq does this in 'correctBins' but it's too late to use this now, so will adjust ours the same way (we weren't 14-May-2018)
#corrected <- counts / fit
#corrected[fit <= 0] <- 0
negs = apply(fit.data[,depth.cols], 2, function(x) which(x<=0))
window.depths[] = lapply(names(negs), function(n) {
  window.depths[negs[[n]],n] = 0
  window.depths[,n]
})
window.depths = cbind(fit.data[,c('chrom','start','end','genome.pos')], window.depths)

samplenames = colnames(fit.data)[depth.cols]
chr.info = colnames(window.depths)[1:depth.cols[1]-1]
#save.image(file=paste(data,"MultisampleCopynumberBinnedAndFitted2.RData", sep='/'))
save(fit.data, raw.data, window.depths, samplenames, chrs, depth.cols, chr.info, file=paste(data,"MultisampleCopynumberBinnedAndFitted2.RData", sep='/'))


plot.dir = paste(outdir, 'coverage_binned_fitted', sep='/')
dir.create(plot.dir, showWarnings=F)

for (s in samplenames) {
  print(s)
  pt = paste('SLX-', subset(all.patient.info, Samplename == s)$SLX.ID, sep='')
  
  #if (length(pt) > 0) {
  pt.plot.dir = strip.whitespace(paste(plot.dir, pt, sep='/'))
  dir.create(pt.plot.dir, showWarnings=F)
  png(paste(pt.plot.dir,paste(s,"coverage_binned_fitted.png",sep="_"),sep="/"),width=1600,height=800)
  #} else {
  #  png(paste(plot.dir,paste(s,"coverage_binned_fitted.png",sep="_"),sep="/"),width=1600,height=800)
  #}
  print(paste("plotting",pt, s))
  
  gg = ggplot(window.depths[,c(chr.info, s)],aes_string(x="genome.pos", y=s)) + 
    lims(y=c(round(min(window.depths[[s]],na.rm=T), 2),round(max(window.depths[[s]],na.rm=T), 2)), x =c(0,max(window.depths$genome.pos))) + 
    geom_point(color='darkred', alpha=0.5) +
    geom_vline(data=chrs, xintercept=c(0,chrs$cum.lengths), color='darkgrey') + 
    geom_text(data=chrs, aes(x=l.pos, y=rep(3, length(chrom)), label=chrom)) + 
    ggtitle(s) + labs(x='genome pos', y="corrected depth/15KB")
  print(gg)
  dev.off()
}

hc = hclust(dist(t( window.depths[depth.cols] )))
ggd = ggdendrogram(hc)

png(paste(outdir,"/hierarchical_clustering_plot_CNA_binned_fitted.png",sep=""),width=2400, height=1600)
print(ggd)
dev.off()

q(save="no")
