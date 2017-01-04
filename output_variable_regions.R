#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3) {
  stop("Missing required arguments: <data directory> <output directory> <cosmic tsv gene file>")
}

library(xlsx) 
library(GenomicRanges)

source("lib/fastPCF.R")
source('lib/load_patient_metadata.R')

Negate('%in%') -> `%nin%`

data = args[1]
outdir = args[2]
cosmic = args[3]

gamma2 = 250
if (length(args) == 4)
  gamma2 = args[4]


#data = '~/Data/Ellie/Analysis'
#gamma2 = 250
#outdir = "~/Documents/MRC-OwnCloud/Progressor Study/Variable Regions"
# cosmic = '~/Data/CosmicCensusGenes.tsv'

plot.dir = paste(data, 'multipcf_plots_fitted_perPatient', sep='/')
print(paste("Plot directory:", plot.dir))
if (!dir.exists(plot.dir))  
  stop("multipcf needs to be run first")

ccgenes = read.table(cosmic, sep='\t', header=T, stringsAsFactors=F)

noLoc = grep(':-',ccgenes$Genome.Location)
noLocGenes = ccgenes[noLoc,]
ccgenes = ccgenes[-noLoc,]

locs = as.data.frame(do.call(rbind, strsplit(ccgenes$Genome.Location, ':|-')))
colnames(locs) = c('chr','start','end')
locs[c('start','end')] = lapply(locs[c('start','end')], function(x) as.numeric(as.character(x)))

ccgenes = cbind(locs, ccgenes)
ccgenes = makeGRangesFromDataFrame(ccgenes, keep.extra.columns = T, start.field = 'start', end.field = 'end')


## Patient info file
patient.file = '~/Data/Ellie/QDNAseq/All_patient_info.xlsx'
if (length(patient.file) != 1)
  stop(paste("Missing patient info file in", data))

all.patient.info = read.patient.info(patient.file)
all.patient.info = arrange(all.patient.info, Status, Patient, Endoscopy.Year)
all.patient.info$Patient = gsub("/", "_", all.patient.info$Patient)
head(all.patient.info)

## Thresholds
sd.threshold = 0.08
# this is what's used in the plotting script... = 1mb
min.probes = 67  

gain.threshold = 1.1
loss.threshold = 0.9
## ------- ##

#filename = path.expand( paste(outdir, "variable_regions.xlsx", sep='/') )

cancer.consensus.genes<-function(regions, ccgenes) {
  reg = as_tibble(regions)
  gr = GenomicRanges::makeGRangesFromDataFrame(regions, start.field = 'start.pos', end.field = 'end.pos')

  ov = findOverlaps(ccgenes, gr)
  for (hit in unique(subjectHits(ov))) {
    reg[hit, 'Gene.Symbols'] =  paste(ccgenes[queryHits(ov)[which(subjectHits(ov) == hit)],]$Gene.Symbol, collapse=',')
  }
  
  regions = cbind(reg, regions[,-(1:5)])
  return(regions)
}

appendFile = F
for (patient.name in unique(all.patient.info$Patient)) {
  print(patient.name)
  
  patient.info = subset(all.patient.info, Patient == patient.name)
  if (nrow(patient.info) <= 0) {
    print(paste("No patient information on patient", patient.name))
    next
  }
  patient.plot.dir = paste(plot.dir,patient.name, sep='/')
  segvals = read.table(paste(patient.plot.dir, "/",patient.name,"_segmentedCoverage_fitted_gamma",gamma2,".txt",sep=""),sep="\t",stringsAsFactors=F,header=T)
  #head(segvals)
  
  sds = apply(segvals[,-(1:5)],1,sd)
  means = apply(segvals[,-(1:5)],1,mean)
  CofV = sds/means
  
  cg = cancer.consensus.genes(segvals[,(1:5)], ccgenes)
  segvals = cbind(cg$Gene.Symbols, segvals)
  colnames(segvals) = sub('cg\\$', '', colnames(segvals) )
  segvals = segvals[,c(2:6,1,7:ncol(segvals))] 
  
  variable.region.indices = which(sds >= sd.threshold & segvals$n.probes >= min.probes)

  variable.regions = segvals[variable.region.indices, patient.info$Samplename]
  colnames(variable.regions) = paste(paste(patient.info$'Endoscopy.Year', patient.info$'Pathology', sep='-'), " (",patient.info$Plate.Index, ")", sep='')
  
  head(segvals)
  
  filename = paste(outdir, "/", patient.name, "_variable_regions_gamma", gamma2, '.txt', sep='')  
  print(paste("variable regions:", nrow(variable.regions)))

  if (nrow(variable.regions) > 0) {
    gain_loss = variable.regions
    for (i in 1:nrow(gain_loss)) {
      gain_loss[i,which(variable.regions[i,] > gain.threshold)] = 'GAIN'
      gain_loss[i,which(variable.regions[i,] < loss.threshold)] = 'LOSS'
      gain_loss[i,which(variable.regions[i,] > loss.threshold & variable.regions[i,] < gain.threshold)] = 'NEUTRAL'
    }
    
    # write.xlsx2(cbind(segvals[variable.region.indices,c(1:5)], rbind( variable.regions, gain_loss ) ), file=filename, sheetName=patient.name, row.names=F, 
    #             append=appendFile)
    # 
    
    write.table(cbind(segvals[variable.region.indices,c(1:6)], rbind( variable.regions ) ), sep='\t', quote=T, row.names=F, file=filename, append = F)
  } else {
    #write.xlsx2(cbind(segvals[variable.region.indices,c(1:5)],  variable.regions ), file=filename, sheetName=patient.name, row.names=F, append=appendFile)
    
    write.table(cbind(segvals[variable.region.indices,c(1:6)],  variable.regions ), sep='\t', quote=F, row.names=F, file=filename, append = F)
  }
  
  # sink(filename, append=T)
  # cat("\n\n## Thresholds applied:")
  # cat("## gamma2:", gamma2, "sd:", sd.threshold, "  min.probes:", min.probes, "\n")
  # cat("## Gain > 1.1, Loss < 0.9")
  # sink()  
  appendFile = T
}

#write.xlsx2(cbind( 'gamma2'=gamma2, 'sd'=sd.threshold, 'min.probes'=min.probes, 'gain'=gain.threshold, 'loss'=loss.threshold), filename, sheetName='thresholds',row.names=F, append=T)



