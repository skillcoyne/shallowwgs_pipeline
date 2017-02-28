#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3) {
  stop("Data directory, patient name, and plot output directory required as arguments")
}

suppressPackageStartupMessages( library(ggplot2) )
suppressPackageStartupMessages( library(ggdendro) )
suppressPackageStartupMessages( library(copynumber) )

source("lib/fastPCF.R")
source('lib/load_patient_metadata.R')


data = args[1]
patient.name = args[2]
outdir = args[3]

#data = '~/Data/Ellie/QDNAseq'
#patient.name = 'AD0098' # 'PR1/HIN/042' #'PR1/HIN/044'
#outdir = '~/Data/Ellie/Analysis'
data.dirs = list.dirs(data, full.names=T, recursive=F)
data.files = list.files(data, full.names=T, recursive=T)

plot.dir = paste(outdir, 'multipcf_plots_fitted_perPatient', sep='/')
print(paste("Plot directory:", plot.dir))
if (!dir.exists(plot.dir))  
  dir.create(plot.dir, recursive=T)

## Patient info file
patient.file = grep('All_patient_info.xlsx', list.files(data, full.names=T), value=T)
if (length(patient.file) != 1)
  stop(paste("Missing patient info file in", data))

all.patient.info = read.patient.info(patient.file)
all.patient.info$Patient = gsub('\\/', '_', all.patient.info$Patient)
patient.name = gsub('\\/', '_', patient.name)

head(all.patient.info)
#print(unique(all.patient.info$Patient))

patient.plot.dir = paste(plot.dir,patient.name, sep='/')
if (!dir.exists(patient.plot.dir))  
  dir.create(patient.plot.dir, recursive=T)


print(patient.name)
patient.info = subset(all.patient.info, Patient == patient.name)
if (nrow(patient.info) <= 0)
  stop(paste("No patient information on patient", patient.name))

## Data file
data.file = grep('MultisampleCopynumberBinnedAndFitted2.RData', data.files, value=T)
if (length(data.file) <= 0)
  stop("Missing necessary Rdata file from Multisample script")
load(file=data.file, verbose=T)

#200216 filter dodgy regions
exclude.file = grep('qDNAseq_blacklistedRegions.txt', data.files, value=T)
if (length(data.file) <= 0)
  stop(paste("Missing necessary exclusion file from 'qDNAseq_blacklistedRegions.txt' in",data))

blacklisted.regions = read.table(exclude.file,sep="\t",header=T,stringsAsFactors=F)
print(paste(nrow(blacklisted.regions), "genomic regions in the exclusion list."))
fit.data$in.blacklist = F

for(r in 1:nrow(blacklisted.regions)) {
  #if (r %% 50 ==0) print(r)
	#print(blacklisted.regions[r,])
	fit.data$in.blacklist[ fit.data$chrom==blacklisted.regions$chromosome[r] & 
	                        fit.data$start>=blacklisted.regions$start[r] & 
	                        fit.data$end<=blacklisted.regions$end[r] ] = T
}

print(paste("# excluded probes = ",sum(fit.data$in.blacklist),sep=""))
if (sum(fit.data$in.blacklist) <= 0)
  stop("There's an issue with the excluded list, no probes excluded.")
window.depths.standardised = window.depths[!fit.data$in.blacklist,]
fit.data = fit.data[!fit.data$in.blacklist,-ncol(fit.data)]

# Some of the indicies were not run for various reasons. At the moment the patient file does not include that information
missingIndicies = setdiff(patient.info$samplename, colnames(window.depths.standardised))
if (length(missingIndicies > 0))
  warning(paste("Missing several samples: ", paste(missingIndicies, collapse=', '), sep=''))
# Grab just those specific samples
window.depths.standardised = window.depths.standardised[,na.omit( match(patient.info$Samplename, colnames(window.depths.standardised)) )]

# there are mutiple samples per patient too
no.samples=ncol(window.depths.standardised)
if (no.samples <= 0)
  stop(paste("No samples represented in the data for", patient.name))

# get the median absolution deviation for the observed window depths in each sample
sdevs = sapply(c(1:no.samples), function(s) getMad( window.depths.standardised[!is.na(window.depths.standardised[,s]),s], k=25))
sdevs[sdevs==0] = NA
sdev = exp(mean(log(sdevs[!is.na(sdevs)]))) # geometric mean??

print("sdevs")
print(sdevs)
print(sdev)

good.bins = which(!is.na(rowSums(window.depths.standardised[,!is.na(sdevs)])))

for(gamma2 in c(250,5,10,25,50,100,500,1000)) { 
  #gamma2=250
  message(paste("gamma2:",gamma2))
	# columns 2 and 3 contain chr and pos, which multipcf requires
  # ?? gamma is the important parameter here, is the loop to help find an appropriate gamma? 
  ## Note that the pre processing of the data was to provide Windsorized values
	res = multipcf( data=cbind(fit.data[good.bins,c('chrom','start')],window.depths.standardised[good.bins,!is.na(sdevs)]), 
	                gamma=gamma2*sdev, fast=F, verbose=T)
	
	write.table(res,paste(patient.plot.dir, '/', patient.name,"_segmentedCoverage_fitted_gamma",gamma2,".txt",sep=""),sep="\t",quote=F)
#read.table(paste(patient.plot.dir, '/', patient.name,"_segmentedCoverage_fitted_gamma",gamma2,".txt",sep=""), sep="\t", header=T)->res 
	
	print(dim(res))
	chrs=unique(res$chrom)
	chr.lengths = vector(mode="integer",length=length(chrs))
	for(c in 1:length(chrs))
		chr.lengths[c] = max(res$end.pos[res$chrom==chrs[c]])+999
	
	cum.lengths = c(0,cumsum(chr.lengths))
  # what is the rest of this loop for?
	res$chrpos.start = NA
	res$chrpos.end = NA
	res$chrpos.start = res$start.pos + cum.lengths[match(res$chrom,chrs)]
	res$chrpos.end = res$end.pos + cum.lengths[match(res$chrom,chrs)]
}

#reload and plot segmented vals
#do 250 first
for(gamma2 in c(250,5,10,25,50,100,500,1000)) {
  print(gamma2)
	segvals = read.table(paste(patient.plot.dir, "/",patient.name,"_segmentedCoverage_fitted_gamma",gamma2,".txt",sep=""),sep="\t",stringsAsFactors=F,header=T)
	gamma.plot = paste(patient.plot.dir, paste("gamma2", gamma2, sep="_"), sep="/")
	print(paste("Plots in", gamma.plot))
	if (!dir.exists(gamma.plot))
	  dir.create(gamma.plot, recursive=T, showWarnings=F)
	
	for(chr in 1:22){
		print(paste('chr',chr,sep=''))
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
		
		print(paste(gamma.plot, "/", "segmentedCoverage_chr",chr,"_gamma",gamma2,".png",sep=""))
		png(paste(gamma.plot, "/", "segmentedCoverage_chr",chr,"_gamma",gamma2,".png",sep=""),width=wid,height=ht)
		par(mfrow=c(n.rows,min(no.samples,10)))
		for(col in which(!is.na(sdevs))){
		  #print(colnames(window.depths.standardised)[col])
			plot(fit.data$end[indices],window.depths.standardised[indices,col],col="red",pch=20,cex=1,cex.axis=2,cex.lab=2,cex.main=2,xlab="pos",ylab="segmented coverage",main=paste("chr",chr,", ",colnames(window.depths.standardised)[col],sep=""),ylim=c(0,3))
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
	  ifelse(all(row.data==row.data[1]) | length(row.data) > 5000 , NA, shapiro.test(row.data)$statistic)
	})
	if (length(which(is.na(shapiro.wilk))) < 1) {
		#out = cbind(segvals[,1:5],CofV,shapiro.wilk)
		#names(out)[6:7] = c("coeff_of_var","ShapiroWilk")
		#use SDs rather than CofV
		out = cbind(segvals[,1:5],sds,shapiro.wilk)
		names(out)[6:7] = c("st_dev","ShapiroWilk")
	
		print(paste(gamma.plot, paste(patient.name,"_allSamples_fitted_multipcf_variability_gamma.txt",sep=""),sep="/"))
		write.table(out,paste(gamma.plot, paste(patient.name,"_allSamples_fitted_multipcf_variability_gamma.txt",sep=""),sep="/"),sep="\t",row.names=F,quote=F)
	
		png(paste(gamma.plot, "/multipcf_variability_clean_gamma",gamma2,".png",sep=""))
		plot(out[,6],out[,7],cex=log(segvals$n.probes)+1,bg=rgb(1,0,0,0.5),pch=21,xlab="standard deviation",ylab="Shapiro-Wilk statistic",xlim=c(0,1))
		dev.off()

		DENS=density(shapiro.wilk[!is.na(shapiro.wilk) & !is.infinite(shapiro.wilk)])
    png(paste(gamma.plot, "/multipcf_ShapiroWilk_density_gamma",gamma2,".png",sep=""))
    plot(DENS,xlab="Shapiro-Wilk statistic",main="")
  	dev.off()
}
	#DENS=density(CofV[!is.na(CofV) & !is.infinite(CofV)])
	#use SDs rather than CofV
	DENS=density(sds[!is.na(sds) & !is.infinite(sds)],from=0,to=1)
	png(paste(gamma.plot, "/multipcf_StDev_density_gamma",gamma2,".png",sep=""))
	plot(DENS,xlab="standard deviation",main="")
	dev.off()
	
	sd.threshold = 0.08
	min.probes = 67 #50
	
	#plot density just for segments with at least 50 probes
	if(length(sds[!is.na(sds) & !is.infinite(sds) & segvals$n.probes >= min.probes]) >= 2) {
		DENS = density(sds[!is.na(sds) & !is.infinite(sds) & segvals$n.probes >= min.probes],from=0,to=1)
		png(paste(gamma.plot, "/multipcf_StDev_density_gamma",gamma2,".png",sep=""))
		plot(DENS,xlab="standard deviation",main="")
		dev.off()
	}
	
	##thresholds selected from density plots
	#190216 This is done, but we may want to try different thresholds
	variable.region.indices = which(sds >= sd.threshold & segvals$n.probes >= min.probes)
	
	print("variable regions:")
	print(segvals[variable.region.indices,1:5])

	write.table(segvals[variable.region.indices,1:5],
	            paste(gamma.plot, "/", patient.name, "_variableRegions_gamma",gamma2,"_sd", sd.threshold, "_n", min.probes, ".txt",sep="")
	            ,sep="\t",row.names=F,quote=F)
	
	
	if(length(variable.region.indices) > 0) {
		hclust.data = segvals[variable.region.indices,-(1:5)]	
		#colnames(hclust.data) = samplenames
		HC = hclust(dist(t(hclust.data)))
		png(paste(gamma.plot, "/hierarchical_clustering_plot_multipcf_gamma",gamma2,"_sd",sd.threshold,"_n",min.probes,".png",sep=""))
		print(ggdendrogram(HC, rotate=T) + ggtitle(paste(patient.name, 'gamma2:', gamma2))  )
		dev.off()
	}
}
print("Finished")

#q(save="no")
