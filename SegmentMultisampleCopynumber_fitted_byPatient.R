#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
args = c('PR1/HIN/044')

if (length(args) == 0) {
  stop("Patient name required as argument")
}

source("lib/fastPCF.R")

#source("/lustre/scratch110/sanger/dw9/haplotype_pipeline/fastPCF.R")
library(copynumber)

#load(file="MultisampleCopynumberBinnedFittedAndStandardised.RData")
load(file="MultisampleCopynumberBinnedAndFitted2.RData")

#200216 filter dodgy regions
blacklisted.regions = read.table("data/qDNAseq_blacklistedRegions.txt",sep="\t",header=T,stringsAsFactors=F)
fit.data$in.blacklist = F

for(r in 1:nrow(blacklisted.regions)) {
	print(r)
	fit.data$in.blacklist[ fit.data$chrom==blacklisted.regions$chromosome[r] & 
	                        fit.data$start>=blacklisted.regions$start[r] & 
	                        fit.data$end<=blacklisted.regions$end[r] ] = T
}

print(paste("# blacklisted probes = ",sum(fit.data$in.blacklist),sep=""))
window.depths.standardised = window.depths[!fit.data$in.blacklist,]
names(window.depths.standardised) = samplenames
fit.data = fit.data[!fit.data$in.blacklist,-ncol(fit.data)]

#just use samples from one patient
patient.info = read.table("data/All_patient_info.txt",sep="\t",header=T,stringsAsFactors=F)
head(patient.info)
patient.info$SLX_ID[patient.info$SLX_ID=="SLX-10729"] = "SLX-10725_10729"
patient.info$SLX_ID[patient.info$SLX_ID=="SLX-10725"] = "SLX-10725_10729"
patient.info$samplename = paste(patient.info$Index,gsub("SLX-","",patient.info$SLX_ID),sep="_")

# I think this is collating the necessary information for a single patient
#patient.index = 2 #as.numeric(commandArgs(T)[1])
patient.name = as.character(args[1])   #unique(patient.info$Patient)[patient.index]
print(patient.name)
patient.info = subset(patient.info, Patient == patient.name)
if (nrow(patient.info) <= 0)
  stop(paste("No patient information on patient", patient.name))

sample.indices = match(patient.info$samplename,samplenames)
samplenames = samplenames[sample.indices]
window.depths.standardised = window.depths.standardised[,sample.indices]

# there are mutiple samples per patient too
no.samples=length(samplenames)

# get the median absolution deviation for the observed window depths in each sample
sdevs = vector(mode="numeric",length=no.samples)
for(s in 1:no.samples) {
	sdevs[s] <- getMad( window.depths.standardised[!is.na(window.depths.standardised[,s]),s], k=25)
}
sdevs[sdevs==0] = NA
#samples 46 and 47 have no good bins, sdev=NA for these samples
sdev = exp(mean(log(sdevs[!is.na(sdevs)]))) # geometric mean??

print("sdevs")
print(sdevs)
print(sdev)

good.bins = which(!is.na(rowSums(window.depths.standardised[,!is.na(sdevs)])))
#good.bins = which(!is.na(rowSums(window.depths.standardised)))

patient.name = gsub("/","_",patient.name)

if (!dir.exists('multipcf_plots_fitted_perPatient/'))
  dir.create('multipcf_plots_fitted_perPatient/')

for(gamma2 in c(5,10,25,50,100,250,500,1000)){ # this could be parallelized but not really necessary right now I suppose
  message(paste("gamma2:",gamma2))
	# columns 2 and 3 contain chr and pos, which multipcf requires
  # ?? gamma is the important parameter here, is the loop to help find an appropriate gamma? 
  ## Note that the pre processing of the data was to provide Windsorized values
	res = multipcf( cbind(fit.data[good.bins,c('chrom','start')],window.depths.standardised[good.bins,!is.na(sdevs)]), gamma=gamma2*sdev, fast=F)
	write.table(res,paste("multipcf_plots_fitted_perPatient/",patient.name,"_segmentedCoverage_fitted_gamma",gamma2,".txt",sep=""),sep="\t",quote=F)
	print(dim(res))
	chrs=unique(res$chrom)
	chr.lengths = vector(mode="integer",length=length(chrs))
	for(c in 1:length(chrs)){
		chr.lengths[c] = max(res$end.pos[res$chrom==chrs[c]])+999
	}
	cum.lengths = c(0,cumsum(chr.lengths))
# what is the rest of this loop for?
	res$chrpos.start = NA
	res$chrpos.end = NA
	res$chrpos.start = res$start.pos + cum.lengths[match(res$chrom,chrs)]
	res$chrpos.end = res$end.pos + cum.lengths[match(res$chrom,chrs)]
}

#reload and plot segmented vals
#do 250 first
for(gamma2 in c(250,5,10,25,50,100,500,1000)){
	segvals = read.table(paste("multipcf_plots_fitted_perPatient/",patient.name,"_segmentedCoverage_fitted_gamma",gamma2,".txt",sep=""),sep="\t",stringsAsFactors=F,header=T)
	plotdir = paste("multipcf_plots_fitted_perPatient",patient.name, paste("gamma2", gamma2, sep="_"), sep="/")
	if (!dir.exists(plotdir))
	  dir.create(plotdir, recursive=T)
	
	for(chr in 1:22){
		print(chr)
		indices = intersect(which(fit.data$chrom==chr),good.bins)
		indices2 = which(segvals$chrom==chr)
		#take out whitespace
		wid=12000
		if(no.samples<10){
			wid=12000*no.samples/10
		}
		ht=5000
		n.rows = ceiling(no.samples / 10)
		if(n.rows < 12){
			ht = 5000 * n.rows / 12
		}
		
		paste(plotdir, "/", patient.name,"_segmentedCoverage_chr",chr,"_gamma",gamma2,".png",sep="")
		png(paste(plotdir, "/", patient.name,"_segmentedCoverage_chr",chr,"_gamma",gamma2,".png",sep=""),width=wid,height=ht)
		par(mfrow=c(n.rows,min(no.samples,10)))
		for(col in which(!is.na(sdevs))){
			plot(fit.data$end[indices],window.depths.standardised[indices,col],col="red",pch=20,cex=1,cex.axis=2,cex.lab=2,cex.main=2,xlab="pos",ylab="segmented coverage",main=paste("chr",chr,", ",samplenames[col],sep=""),ylim=c(0,3))
			for(ind in indices2){
				lines(c(segvals$start.pos[ind],segvals$end.pos[ind]),rep(segvals[ind,col+5],2),col="green",lwd=4)
			}
		}
		dev.off()
	}
	
	sds = apply(segvals[,-(1:5)],1,sd)
	means = apply(segvals[,-(1:5)],1,mean)
	CofV = sds/means

	shapiro.wilk = sapply(1:nrow(segvals),function(i) {
	  row.data=unlist(segvals[i,-(1:5)])
	  ifelse(all(row.data==row.data[1]), NA, shapiro.test(row.data)$statistic)
	})
	
	#out = cbind(segvals[,1:5],CofV,shapiro.wilk)
	#names(out)[6:7] = c("coeff_of_var","ShapiroWilk")
	#use SDs rather than CofV
	out = cbind(segvals[,1:5],sds,shapiro.wilk)
	names(out)[6:7] = c("st_dev","ShapiroWilk")
	
	write.table(out,paste("multipcf_plots_fitted_perPatient/",patient.name,"_allSamples_fitted_multipcf_variability_gamma",gamma2,".txt",sep=""),sep="\t",row.names=F,quote=F)
	
	png(paste(plotdir, "/", patient.name,"_multipcf_variability_clean_gamma",gamma2,".png",sep=""))
	plot(out[,6],out[,7],cex=log(segvals$n.probes)+1,bg=rgb(1,0,0,0.5),pch=21,xlab="standard deviation",ylab="Shapiro-Wilk statistic",xlim=c(0,1))
	dev.off()
	#DENS=density(CofV[!is.na(CofV) & !is.infinite(CofV)])
	#use SDs rather than CofV
	DENS=density(sds[!is.na(sds) & !is.infinite(sds)],from=0,to=1)
	#png(paste("multipcf_plots_fitted_perPatient/",patient.name,"_multipcf_CoeffOfVar_density_gamma",gamma2,".png",sep=""))
	png(paste(plotdir, "/",patient.name,"_multipcf_StDev_density_gamma",gamma2,".png",sep=""))
	plot(DENS,xlab="standard deviation",main="")
	dev.off()
	DENS=density(shapiro.wilk[!is.na(shapiro.wilk) & !is.infinite(shapiro.wilk)])
	png(paste(plotdir, "/",patient.name,"_multipcf_ShapiroWilk_density_gamma",gamma2,".png",sep=""))
	plot(DENS,xlab="Shapiro-Wilk statistic",main="")
	dev.off()

	sd.threshold=0.08
	n.threshold=50
	
	#plot density just for segments with at least 50 probes
	if(length(sds[!is.na(sds) & !is.infinite(sds) & segvals$n.probes>=n.threshold])>=2) {
		DENS=density(sds[!is.na(sds) & !is.infinite(sds) & segvals$n.probes>=n.threshold],from=0,to=1)
		png(paste(plotdir, "/",patient.name,"_multipcf_StDev_density_gamma",gamma2,"_50probes.png",sep=""))
		plot(DENS,xlab="standard deviation",main="")
		dev.off()
	}
	
	##thresholds selected from density plots
	#190216 This is done, but we may want to try different thresholds
	variable.region.indices = which(sds>=sd.threshold & segvals$n.probes>=n.threshold)
	print("variable regions:")
	print(segvals[variable.region.indices,1:5])
	write.table(segvals[variable.region.indices,1:5],paste("multipcf_plots_fitted_perPatient/",patient.name,"_variableRegions_gamma",gamma2,"_sd",sd.threshold,"_n",n.threshold,".txt",sep=""),sep="\t",row.names=F,quote=F)
	
	no.variable.regions = length(variable.region.indices)
	if(no.variable.regions>0){
		hclust.data = segvals[variable.region.indices,-(1:5)]	
		colnames(hclust.data) = samplenames
		HC = hclust(dist(t(hclust.data)))
		png(paste(plotdir, "/",patient.name,"_hierarchical_clustering_plot_multipcf_gamma",gamma2,"_sd",sd.threshold,"_n",n.threshold,".png",sep=""))
		plot(HC,xlab="",sub="")
		dev.off()
	}
}

q(save="no")
