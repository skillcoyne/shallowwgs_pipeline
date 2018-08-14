#!/usr/bin/env Rscript

# Step 2
args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3) {
  stop("Data directory, patient name, and plot output directory required as arguments")
}

suppressPackageStartupMessages( library(ggplot2) )
suppressPackageStartupMessages( library(gridExtra) )
suppressPackageStartupMessages( library(ggdendro) )
suppressPackageStartupMessages( library(copynumber) )
suppressPackageStartupMessages( library(dplyr) )

suppressPackageStartupMessages( source("~/workspace/shallowwgs_pipeline/lib/fastPCF.R") )
suppressPackageStartupMessages( source('~/workspace/shallowwgs_pipeline/lib/load_patient_metadata.R') )
suppressPackageStartupMessages( source('~/workspace/shallowwgs_pipeline/lib/data_func.R') )


data = args[1]
patient.name = args[2]
outdir = args[3]

# data = '~/Data/Ellie/QDNAseq'
# patient.name = 'PR1/HIN/042' #'AD0098' # 'PR1/HIN/042' #'PR1/HIN/044'
# outdir = '~/Data/Ellie/Analysis'
data.dirs = list.dirs(data, full.names=T, recursive=F)
data.files = list.files(data, full.names=T, recursive=T)

patient.name = gsub('\\/', '_', patient.name )

plot.dir = paste(outdir, 'multipcf_plots_fitted_perPatient', sep='/')
print(paste("Plot directory:", plot.dir))
if (!dir.exists(plot.dir))  
  dir.create(plot.dir, recursive=T)

## Patient info file
filename = 'All_patient_info.xlsx'
patient.file = grep(filename, list.files(data, full.names=T), value=T)
if (length(patient.file) != 1)
  stop(paste("Missing patient info file in", data))
message(paste("Reading file", patient.file))

all.patient.info = read.patient.info(patient.file, set='all')$info

head(subset(all.patient.info, Hospital.Research.ID == patient.name))
#print(unique(all.patient.info$Patient))

patient.plot.dir = paste(plot.dir,patient.name, sep='/')
if (!dir.exists(patient.plot.dir))  
  dir.create(patient.plot.dir, recursive=T)

print(patient.name)
patient.info = subset(all.patient.info, Hospital.Research.ID == patient.name)
if (nrow(patient.info) <= 0)
  stop(paste("No patient information on patient", patient.name))

## Data file
data.file = grep('MultisampleCopynumberBinnedAndFitted2.RData', data.files, value=T)
if (length(data.file) <= 0)
  stop("Missing necessary Rdata file from Multisample script")
load(file=data.file, verbose=T)

write.table(fit.data[,c('location','chrom','start','end',patient.info$Samplename)], sep='\t', quote=F, file=paste(patient.plot.dir, 'raw_fitted.txt', sep='/'))
write.table(raw.data[,c('location','chrom','start','end',patient.info$Samplename)], sep='\t', quote=F, file=paste(patient.plot.dir, 'raw_counts.txt', sep='/'))

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
missingIndicies = setdiff(patient.info$Samplename, colnames(window.depths.standardised))
if (length(missingIndicies > 0))
  warning(paste("Missing several samples: ", paste(missingIndicies, collapse=', '), sep=''))

# Grab just those specific samples
window.depths.standardised = as.data.frame(window.depths.standardised[,na.omit( match(patient.info$Samplename, colnames(window.depths.standardised)) )])
colnames(window.depths.standardised) = patient.info$Samplename

raw.summary = as.data.frame(do.call(rbind, lapply(window.depths.standardised, function(x) {
  x = x[complete.cases(x)]
  sm = summary(x)
  sm['Var.'] = var( x[which(x <= 2 & x >= 0)] )
  sm['High.Ratio'] = length(which(x > 1.5 | x < 0.5))/length(x) 
  sm
})))
#raw.summary$Excluded = raw.summary$High.Ratio >= 0.05
write.table(raw.summary, sep='\t', quote=F, col.names=NA,file=paste(patient.plot.dir, '/', 'raw_sample_summary.txt', sep=''),)
write.table(do.call(rbind, lapply(window.depths.standardised, summary)), sep='\t', quote=F, col.names=NA, file=paste(patient.plot.dir, '/', patient.name,"_std_window_depths_Summary.txt",sep=""))

#exclude = which(raw.summary$High.Ratio >= 0.05)
#samplesExcluded = colnames(window.depths.standardised)[exclude]

#window.depths.standardised[,exclude] = NA
head(window.depths.standardised)
ncol(window.depths.standardised)

# there are mutiple samples per patient too - usually
no.samples = ncol(window.depths.standardised)
if (no.samples <= 0)
  stop(paste("No samples represented in the data for", patient.name))

# get the median absolution deviation for the observed window depths in each sample
sdevs = sapply(c(1:no.samples), function(s) getMad( window.depths.standardised[!is.na(window.depths.standardised[,s]),s], k=25))
sdevs[sdevs==0] = NA
sdev = exp(mean(log(sdevs[!is.na(sdevs)]))) # geometric mean across all sample sds, doesn't work when there's only one...

print("sdevs")
print(sdevs)
print(sdev)

good.bins = which(!is.na(rowSums(as.data.frame(window.depths.standardised[,!is.na(sdevs)]))))
#for(gamma2 in c(250,5,10,25,50,100,500,1000)) { 
  {
  gamma2=250 
  message(paste("gamma2:",gamma2))
  # columns 2 and 3 contain chr and pos, which multipcf requires
  # ?? gamma is the important parameter here, is the loop to help find an appropriate gamma? 
  ## Note that the pre processing of the data was to provide Windsorized values
  filename = paste(patient.plot.dir, '/', patient.name,"_segmentedCoverage_fitted_gamma",gamma2,".txt",sep="")
  
  data = cbind(fit.data[good.bins,c('chrom','start')],window.depths.standardised[good.bins,!is.na(sdevs)])
  
  #data = cbind(fit.data[good.bins,c('chrom','start')],window.depths.standardised[good.bins,sample,drop=F])
  head(data)
  
#head(data)
  if (ncol(data) < 4) { # Single sample
    #colnames(data)[3] = patient.info$Samplename
    res = pcf( data=data, gamma=gamma2*sdev, fast=F, verbose=T, return.est=F)
    colnames(res)[grep('mean', colnames(res))] = patient.info$Samplename
    res$sampleID = NULL
  } else { # for most we have multiple samples
    res = multipcf( data=data, gamma=gamma2*sdev, fast=F, verbose=T, return.est=F)
  }
  
  residuals = calculateSegmentResiduals(res[res$n.probes >= 67,], data)
  residuals = residuals[which(!is.na(sdevs))]

  res.summary<-function(df) {
    sm = summary(df)
    sm['Var.'] = var(df)
    sm
  }

  seg.summary = as.data.frame(do.call(rbind, lapply(res[,-(1:5),drop=F], res.summary)))
  seg.summary$n.segs = nrow(res)
  seg.summary$n.filtered.segs = length(which(res$n.probes >= 67))
  write.table(do.call(rbind, lapply(res[,-(1:5),drop=F], res.summary)), sep='\t', quote=F, col.names=NA, file=paste(patient.plot.dir, '/', patient.name,"_segmentedCoverage_gamma",gamma2,"_Summary.txt",sep=""))

  save(residuals, file=paste(patient.plot.dir, '/', 'residuals.Rdata',sep=''))
    
  #plotChrom(data=data,segments=res,chrom=1, layout=c(ncol(data)-2,1))
  message(paste("Writing pcf output to",filename))
  write.table(res, file=filename, sep="\t", quote=F)
  write.table(data, file=paste(patient.plot.dir, 'filtered_fitted.txt', sep='/'), sep='\t', quote=F)
  
  
  #read.table(filename, sep="\t", header=T)->res 
  print(dim(res))

  # This used to be a separate loop but might as well plot at the same time
  segvals = res
  gamma.plot = paste(patient.plot.dir, paste("gamma2", gamma2, sep="_"), sep="/")
  print(paste("Plots in", gamma.plot))
  if (!dir.exists(gamma.plot))
    dir.create(gamma.plot, recursive=T, showWarnings=F)
  

  chr.info = get.chr.lengths(file='/tmp/hg19_info.txt')[1:22,]
  chr.info$chrom = sub('chr','', chr.info$chrom)
  chr.info$chrom = factor(chr.info$chrom, levels=c(1:22), ordered=T)
  
  plots = list()
  #for(col in which(!is.na(sdevs))) {
  samples = grep('^D',colnames(segvals))
  for (col in samples) {
    sample = colnames(segvals)[col] #[col+5]
    
    p = plot.segmented.genome(fit.data[good.bins,], segvals[,c(1:5,col)], window.depths.standardised[good.bins,sample,drop=F]) + labs(title=sample)
    
    #df = cbind.data.frame('chrom'=fit.data$chrom[good.bins], 'position'=fit.data$end[good.bins], 'seg.cov'=window.depths.standardised[good.bins,col])
    #df2 = cbind.data.frame('chrom'=segvals$chrom, 'start'=segvals$start.pos, 'end'=segvals$end.pos, 'seg.val'=segvals[,col+5])
    
    rows = which(segvals$n.probes < 67)
    df3 = cbind.data.frame('chrom'=segvals$chrom[rows], 'start'=segvals$start.pos[rows], 'end'=segvals$end.pos[rows], 'seg.val'=segvals[rows,col])
    p = p + geom_segment(data=df3, aes(x=start, xend=end, y=seg.val, yend=seg.val), color='blue', lwd=5) 

   ggsave(filename=paste(gamma.plot, '/',"segmentedCoverage_", sample, "_gamma",gamma2,".png",sep=""),plot=p, width=12, height=4, units='in', scale=2)

   plots[[sample]] = p
  }  
  #ggsave(filename=paste(patient.name,".png",sep=""),plot=do.call(grid.arrange, c(plots, ncol=1)), width=15, height=4*no.samples, units='in', scale=1.5, limitsize = F)
  
  ggsave(filename=paste(gamma.plot, '/',"segmentedCoverage_chr_gamma",gamma2,".png",sep=""),plot=do.call(grid.arrange, c(plots, ncol=1)), width=15, height=4*no.samples, units='in', scale=1.5, limitsize = F)

  if ( length(grep('^D\\d+', colnames(segvals))) > 1  ) {
    sds = apply(segvals[,-(1:5)],1,sd)
    means = apply(segvals[,-(1:5)],1,mean)
    CofV = sds/means
  
    if ((ncol(segvals)-5) > 3) {
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
        head(out)

        p = ggplot(out, aes(st_dev, ShapiroWilk, size=log(n.probes)+1),color='darkred') + geom_point(alpha=0.5, fill='indianred1',shape=21) + xlim(0,1) + scale_size(range=c(0,10)) + theme(legend.position = 'none') + labs(x='std.dev.', y='Shapiro-Wilk statistic', title=patient.name)
      
        ggsave(file=paste(gamma.plot, "/multipcf_variability_clean_gamma",gamma2,".png",sep=""), plot=p, width=5,height=5,units = 'in')

        p = ggplot(out) + geom_density(aes(ShapiroWilk),fill='indianred3', alpha=0.5) + labs(x='Shapiro-Wilk statistic', title=patient.name)
        ggsave(filename = paste(gamma.plot, "/multipcf_ShapiroWilk_density_gamma",gamma2,".png",sep=""), plot=p, height=5,width=5,units='in')
      }
    } else {
      warning(paste("Too few samples to run shapiro wilk", (ncol(segvals)-5)))
    }
    
    sd.threshold = 0.08
    min.probes = 67 #50
    
    #plot density just for segments with at least 50 probes
    # if(length(sds[!is.na(sds) & !is.infinite(sds) & segvals$n.probes >= min.probes]) >= 2) {
    #   png(paste(gamma.plot, "/multipcf_StDev_density_gamma",gamma2,".png",sep=""))
    #   DENS = density(sds[!is.na(sds) & !is.infinite(sds) & segvals$n.probes >= min.probes],from=0,to=1)
    #   plot(DENS,xlab="standard deviation",main=patient.name)
    #   dev.off()
    # }
    
    ##thresholds selected from density plots
    #190216 This is done, but we may want to try different thresholds
    variable.region.indices = which(sds >= sd.threshold & segvals$n.probes >= min.probes)
    
    print("variable regions:")
    #print(segvals[variable.region.indices,1:5])
    
    write.table(segvals[variable.region.indices,1:5],
                paste(gamma.plot, "/", patient.name, "_variableRegions_gamma",gamma2,"_sd", sd.threshold, "_n", min.probes, ".txt",sep="")
                ,sep="\t",row.names=F,quote=F)
    
    
    if (length(variable.region.indices) > 0) {
      hclust.data = segvals[variable.region.indices,-(1:5)]	
      colnames(hclust.data) = paste(colnames(hclust.data), apply(subset(patient.info, Samplename %in% colnames(hclust.data), select=c('Endoscopy.Year','Pathology')),1,paste,collapse=':'), sep='-')
      #colnames(hclust.data) = samplenames
      png(paste(gamma.plot, "/hierarchical_clustering_plot_multipcf_gamma",gamma2,"_sd",sd.threshold,"_n",min.probes,".png",sep=""))
      HC = hclust(dist(t(hclust.data)))
      print(ggdendro::ggdendrogram(HC, rotate=T) + ggtitle(paste(patient.name, 'gamma2:', gamma2))  )
      dev.off()
    }
  }
}

print("Finished")

#q(save="no")
