library(plyr)
# this version removes qDNAseq 'blacklisted' regions and uses 'fitted' read counts

setwd("/lustre/scratch112/sanger/cgppipe/PanCancerDownloads/workspace/dw9/Ellie")

setwd("~/Documents/MRC-OwnCloud/skillcoyne/Documents/Ellie BE Project/")

data = "~/Data/Ellie/QDNAseq"
data.dirs = paste(data, c("SLX-10725_10729","SLX-9242","SLX-9246"), sep="/")

outdir = "plots_fitted"
if ( !dir.exists(outdir) ) 
  dir.create(outdir)

# Load the data, binding each new data file into single data table
fit.data = NULL
raw.data = NULL
samplenames=NULL
for(dir in data.dirs){
  print(dir)
	files = list.files(dir,pattern="fittedReadCounts.txt")
	short.name = gsub("_fitted","",gsub("^.*/SLX-","",dir))
	for(file in files){
		print(file)
		spl = strsplit(file,"\\.")
		
		fitted = paste(dir,"/",file,sep="")
		raw = paste(dir,"/",gsub("fitted","raw",file),sep="")
		if (!file.exists(fitted)) {
		  stop("Missing fitted file ")
		}
		if (!file.exists(raw))
		  stop("Missing raw file")
		  
		next
		
		if(is.null(fit.data)){
			fit.data = read.table(paste(dir,"/",file,sep=""),sep="\t",header=T,stringsAsFactors=F)  # fitted reads
			samplenames = paste(spl[[1]][2],short.name,sep="_")
			raw.data = read.table(paste(dir,"/",gsub("fitted","raw",file),sep=""),sep="\t",header=T,stringsAsFactors=F) # raw reads
		}else{ # add to table
			fit.data = cbind(fit.data,(read.table(paste(dir,"/",file,sep=""),sep="\t",header=T,stringsAsFactors=F))[,5])
			samplenames = c(samplenames,paste(spl[[1]][2],short.name,sep="_"))
			raw.data = cbind(raw.data,(read.table(paste(dir,"/",gsub("fitted","raw",file),sep=""),sep="\t",header=T,stringsAsFactors=F))[,5])
		}
	}
}
names(fit.data)[-(1:4)] = samplenames
names(raw.data)[-(1:4)] = samplenames
save.image(file="MultisampleCopynumberBinnedAndFitted1.RData")

depth.cols = (1:ncol(fit.data))[-(1:4)] # "count" columns
no.samples = length(samplenames)
# adjusting the raw counts by the fitted counts (need to find out what "fitted" is in this context)
window.depths = raw.data[,depth.cols]/fit.data[,depth.cols] 


chrs = ddply(fit.data, .(chrom), summarise, length=max(end))
chrs$chrom = factor(chrs$chrom, c((1:22),"X","Y","M"), ordered=T)
chrs = chrs[order(chrs$chrom),]
chrs$cum.lengths = cumsum(as.numeric(chrs$length))

# chrs = unique(fit.data$chrom)
# chr.lengths = vector(length = length(chrs),mode="integer")
# for(c in 1:length(chrs)){
# 	chr = chrs[c]
# 	chr.lengths[c] = max(fit.data$end[fit.data$chrom==chr])
# }

#cum.lengths = cumsum(as.numeric(chr.lengths))
fit.data$genome.pos = as.numeric(fit.data$start)
# moving the position information from chromosomes to the genome position (cumulutive sum of the lengths of the chrs)
#for(c in 2:length(chr.lengths)){
for (c in 2:nrow(chrs)) {
	fit.data$genome.pos[fit.data$chrom == chrs[c, 'chrom']] = fit.data$genome.pos[fit.data$chrom == chrs[c, 'chrom']] + chrs[c-1, 'cum.lengths']   #cum.lengths[c-1]
}

save.image(file="MultisampleCopynumberBinnedAndFitted2.RData")

png(paste(outdir,"/coverage_binned_fitted.png",sep=""),width=24000,height=4000)
par(mfrow=c(12,10))
for(s in 1:length(samplenames)){
	print(paste("plotting ",samplenames[s],sep=""))
	plot(fit.data$genome.pos,window.depths[,s],pch=20,col="red",cex=0.5,xlab="",ylab="corrected depth/15KB",xaxt="n",main=samplenames[s],ylim=c(0,3))
	abline(v=0,lty=1,col="lightgrey",lwd=2)
	for(c in chrs){
	  abline(v=chrs[c, 'cum.lengths'],lty=1,col="lightgrey",lwd=2)
		#abline(v=cum.lengths[c],lty=1,col="lightgrey",lwd=2)
		if(c==1){
			text(cum.lengths[c]/2,3,chrs[c],pos=1,cex=2)
		}else{
			text((cum.lengths[c-1]+cum.lengths[c])/2,3,chrs[c],pos=1,cex=2)
		}
	}
}
dev.off()

hclust.data = window.depths
names(hclust.data)=samplenames
HC = hclust(dist(t(hclust.data)))
png(paste(outdir,"/hierarchical_clustering_plot_CNA_binned_fitted.png",sep=""),width=2400)
plot(HC,xlab="",main="",sub="")
dev.off()

q(save="no")
