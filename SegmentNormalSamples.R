#!/usr/bin/env Rscript

# Step 2
args = commandArgs(trailingOnly=TRUE)

if (length(args) < 2) {
  stop("Data directory and plot output directory required as arguments")
}

suppressPackageStartupMessages( library(ggplot2) )
suppressPackageStartupMessages( library(gridExtra) )
suppressPackageStartupMessages( library(ggdendro) )
suppressPackageStartupMessages( library(copynumber) )
suppressPackageStartupMessages( library(dplyr) )

suppressPackageStartupMessages( source("~/workspace/shallowWGSpipeline/lib/fastPCF.R") )
suppressPackageStartupMessages( source('~/workspace/shallowWGSpipeline/lib/load_patient_metadata.R') )
suppressPackageStartupMessages( source('~/workspace/shallowWGSpipeline/lib/data_func.R') )


data = args[1]
#patient.name = args[2]
outdir = args[2]

# data = '~/Data/Ellie/QDNAseq'
# outdir = '~/Data/Ellie/Analysis/Normals'
data.dirs = list.dirs(data, full.names=T, recursive=F)
data.files = list.files(data, full.names=T, recursive=T)

plot.dir = paste(outdir, 'multipcf_plots_normals', sep='/')
print(paste("Plot directory:", plot.dir))
if (!dir.exists(plot.dir))  
  dir.create(plot.dir, recursive=T)

## Patient info file
filename = 'All_patient_info.xlsx'
patient.file = grep(filename, list.files(data, full.names=T), value=T)
if (length(patient.file) != 1)
  stop(paste("Missing patient info file in", data))
message(paste("Reading file", patient.file))

normals = read.patient.info(patient.file, set='all')$removed

head(normals)

## Data file
data.file = grep('MultisampleCopynumberBinnedAndFitted2.RData', data.files, value=T)
if (length(data.file) <= 0)
  stop("Missing necessary Rdata file from Multisample script")
load(file=data.file, verbose=T)

fit.data = fit.data[,c('location','chrom','start','end',normals$Samplename)]
raw.data = raw.data[,c('location','chrom','start','end',normals$Samplename)]
window.depths = window.depths[,c(normals$Samplename)]

rm(samplenames,chrs,depth.cols,chr.info)

dim(fit.data)
dim(raw.data)
dim(window.depths)


patient.plot.dir = plot.dir
write.table(fit.data[,c('location','chrom','start','end',normals$Samplename)], sep='\t', quote=F, file=paste(patient.plot.dir, 'raw_fitted.txt', sep='/'))
write.table(raw.data[,c('location','chrom','start','end',normals$Samplename)], sep='\t', quote=F, file=paste(patient.plot.dir, 'raw_counts.txt', sep='/'))


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

write.table(do.call(rbind, lapply(window.depths.standardised, summary)), sep='\t', quote=F, col.names=NA, file=paste(patient.plot.dir, '/', "normals_std_window_depths_Summary.txt",sep=""))

# there are mutiple samples per patient too - usually
no.samples = ncol(window.depths.standardised)
if (no.samples <= 0)
  stop(paste("No samples represented in the data for", patient.name))

# get the median absolution deviation for the observed window depths in each sample
sdevs = sapply(c(1:no.samples), function(s) getMad( window.depths.standardised[!is.na(window.depths.standardised[,s]),s], k=25))
sdevs[sdevs==0] = NA
#sdev = exp(mean(log(sdevs[!is.na(sdevs)]))) # geometric mean across all sample sds, doesn't work when there's only one...

print("sdevs")
print(sdevs)

good.bins = which(!is.na(rowSums(as.data.frame(window.depths.standardised[,!is.na(sdevs)]))))
patient.name = 'normals'
  gamma2=250 
  message(paste("gamma2:",gamma2))
  # columns 2 and 3 contain chr and pos, which multipcf requires
  # ?? gamma is the important parameter here, is the loop to help find an appropriate gamma? 
  ## Note that the pre processing of the data was to provide Windsorized values
  
  data = cbind(fit.data[good.bins,c('chrom','start')],window.depths.standardised[good.bins,!is.na(sdevs)])
  
  head(data)
  samples = grep('^D', colnames(data), value=T)
  sample.summaries = tibble()
  raw.sum = tibble()
  for (i in 1:length(samples)) {
    sdev = sdevs[i]
    sample = samples[i]
    print(paste(i,sample))
    
    patient.plot.dir = paste(plot.dir, sample, sep='/')
    dir.create(patient.plot.dir, showWarnings = F, recursive = T)

    filename = paste(patient.plot.dir, '/', "segmentedCoverage_fitted_gamma",gamma2,".txt",sep="")
    
    res = pcf( data=data[,c('chrom','start',sample)], gamma=gamma2*sdev, fast=F, verbose=T, return.est=F)
    colnames(res)[grep('mean', colnames(res))] = sample
    res$sampleID = NULL
    
    res.summary<-function(df) {
      sm = summary(df)
      sm['SD'] = sd(df)
      sm['High.Ratio'] = length(which(df > 1.5)) /length(df)
      sm
    }

    
    rs = as_tibble(do.call(rbind, lapply(data[,sample,drop=F], res.summary)))
    rs = cbind('Samplename'=sample, rs, 'n.segs'=nrow(res))

    raw.sum = rbind(raw.sum, rs)
    
    rs = as_tibble(do.call(rbind, lapply(res[,-(1:5),drop=F], res.summary)))
    sample.summaries = rbind(sample.summaries, rs)
    
    #plotChrom(data=data,segments=res,chrom=1, layout=c(ncol(data)-2,1))
    message(paste("Writing pcf output to",filename))
    write.table(res, file=filename, sep="\t", quote=F)
    
    # This used to be a separate loop but might as well plot at the same time
    segvals = res
    gamma.plot = paste(patient.plot.dir, paste("gamma2", gamma2, sep="_"), sep="/")
    dir.create(gamma.plot, showWarnings = F, recursive = T)
    
    print(paste("Plots in", gamma.plot))
    if (!dir.exists(gamma.plot))
      dir.create(gamma.plot, recursive=T, showWarnings=F)
    
    
    chr.info = get.chr.lengths(file='/tmp/hg19_info.txt')[1:22,]
    chr.info$chrom = sub('chr','', chr.info$chrom)
    chr.info$chrom = factor(chr.info$chrom, levels=c(1:22), ordered=T)
    
    p = plot.segmented.genome(fit.data[good.bins,c('location','chrom','start','end',sample)],
                              segvals, window.depths.standardised[good.bins,sample,drop=F]) + labs(title=sample)
    
      ggsave(filename=paste(gamma.plot, '/',"segmentedCoverage_gamma",gamma2,".png",sep=""),plot=p, width=12, height=4, units='in', scale=2)
    
}
  sample.summaries = cbind('Samplename'=samples, sample.summaries)
  
  write.table(raw.sum,  sep='\t', quote=F, row.names=F, file=paste(plot.dir, 'normals_raw_summary.txt',sep='/'))
  
  write.table(sample.summaries, sep='\t', quote=F, col.names=T, row.names=F, file=paste(plot.dir, '/', "normals_segmentedCoverage_gamma",gamma2,"_Summary.txt",sep=""))
  
  
print("Finished")

#q(save="no")
