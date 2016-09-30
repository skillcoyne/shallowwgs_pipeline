library(plyr)
library(ggplot2)
library(ggdendro)
# this version removes qDNAseq 'blacklisted' regions and uses 'fitted' read counts

#setwd("/lustre/scratch112/sanger/cgppipe/PanCancerDownloads/workspace/dw9/Ellie")

#setwd("~/Documents/MRC-OwnCloud/skillcoyne/Documents/Ellie BE Project/")

data = "~/Data/Ellie/QDNAseq"

data.dirs = list.dirs(data, full.names=T, recursive=F)
#data.dirs = paste(data, c("SLX-10725_10729","SLX-9242","SLX-9246"), sep="/")

outdir = "~/tmp/plots_fitted"
if ( !dir.exists(outdir) ) 
  dir.create(outdir)

# Load the data, binding each new data file into single data table
fit.data = NULL
raw.data = NULL
samplenames=NULL
for(dir in data.dirs) {
  print(dir)
	files = list.files(dir,pattern="[fitted|corrected]ReadCounts.txt")

	if (length(files) == 1) {
	  fit = read.table(paste(dir,"/",files[1],sep=""),sep="\t",header=T,stringsAsFactors=F)  # fitted reads
	  cols = grep('D\\d', colnames(fit))
	  colnames(fit)[cols]= paste(colnames(fit)[cols], '10722_11823', sep="_") # hardcoded...fix this from spreadsheets later
	  
	  raw = read.table(paste(dir,"/",gsub("fitted|correctedR","r",files[1]),sep=""),sep="\t",header=T,stringsAsFactors=F) # raw reads
	  colnames(raw)[cols]= paste(colnames(raw)[cols], '10722_11823', sep="_") # hardcoded...fix this from spreadsheets later
	  
	  if (is.null(fit.data)) {
	    fit.data = fit
	    raw.data = raw
	  } else {
	    fit.data = cbind(fit.data,fit[,cols])
	    raw.data = cbind(raw.data,raw[,cols])
	  }

	} else {
  	short.name = gsub("_fitted","",gsub("^.*/SLX-","",dir))
	  for(file in files){
		  print(file)
		  spl = strsplit(file,"\\.")
		
  		fitted.file = paste(dir,"/",file,sep="")
  		raw.file = paste(dir,"/",gsub("fitted","raw",file),sep="")
  		if (!file.exists(fitted.file) | !file.exists(raw.file)) 
  		  stop("Missing fitted/raw file ")

  		fit = read.table(fitted.file,sep="\t",header=T,stringsAsFactors=F)  # fitted reads
  		colnames(fit)[5] = paste(spl[[1]][2],short.name,sep="_")
  		
  		raw = read.table(raw.file,sep="\t",header=T,stringsAsFactors=F) # raw reads
  		colnames(raw)[5] = paste(spl[[1]][2],short.name,sep="_")

	  	if(is.null(fit.data)){
	  		fit.data = fit
	  		raw.data = raw
	  	}else{ # add to table
	  		fit.data = cbind(fit.data,fit[,5])
	  		raw.data = cbind(raw.data, raw[,5])
	  		
	  		colnames(fit.data)[ncol(fit.data)] = colnames(fit)[5]
	  		colnames(raw.data)[ncol(raw.data)] = colnames(raw)[5]
	  	}	  
	  }
  }
}
save.image(file=paste(data,"MultisampleCopynumberBinnedAndFitted1.RData", sep='/'))


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
save.image(file=paste(data,"MultisampleCopynumberBinnedAndFitted2.RData", sep='/'))

# "count" columns
depth.cols = grep('^D\\d',colnames(fit.data))
no.samples = length(depth.cols)
# adjusting the raw counts by the fitted counts (need to find out what "fitted" is in this context)
window.depths = raw.data[,depth.cols]/fit.data[,depth.cols] 
window.depths = cbind(fit.data[,c('chrom','start','end','genome.pos')], window.depths)

samplenames = colnames(fit.data)[depth.cols]
chr.info = colnames(window.depths)[1:depth.cols[1]-1]
for (s in samplenames) {
  print(paste("plotting ",s,sep=""))
  png(paste(outdir,paste(s,"coverage_binned_fitted.png",sep="_"),sep="/"),width=1600,height=800)
  gg = ggplot(window.depths[,c(chr.info, s)],aes_string(x="genome.pos", y=s)) + lims(y=c(-0.25,3), x =c(0,max(window.depths$genome.pos))) + 
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
