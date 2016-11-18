#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#if (length(args) < 3) {
#  stop("Data directory, patient name, and plot output directory required as arguments")
#}

library(xlsx) 

source("lib/fastPCF.R")
source('lib/load_patient_metadata.R')

Negate('%in%') -> `%nin%`

data = args[1]
gamma2 = args[2]
outdir = args[3]

data = '~/Data/Ellie/Analysis'
gamma2 = 250
outdir = "~/Documents/MRC-OwnCloud/Progressor Study/Variable Regions"

plot.dir = paste(data, 'multipcf_plots_fitted_perPatient', sep='/')
print(paste("Plot directory:", plot.dir))
if (!dir.exists(plot.dir))  
  stop("multipcf needs to be run first")

## Patient info file
patient.file = grep('patient_info.xlsx', data.files, value=T)
if (length(patient.file) != 1)
  stop(paste("Missing patient info file in", data))

all.patient.info = read.patient.info(patient.file)
all.patient.info = arrange(all.patient.info, Status, Patient, Endoscopy.Year)
all.patient.info$Patient = gsub("/", "_", all.patient.info$Patient)
head(all.patient.info)
## TODO Mistake in earlier version of the patient file caused this will be fixed for the next run
snames = grep('_1072(5|9)$', all.patient.info$Samplename)
if ( length(snames > 0) ) {
  all.patient.info$Samplename[snames] =  paste( sub('-','_', all.patient.info$Plate.Index[snames]), '10725_10729', sep='_' )
}

## Thresholds
sd.threshold = 0.08 
min.probes = 67  # this is what's used in the plotting script...

gain.threshold = 1.1
loss.threshold = 0.9
## ------- ##

filename = path.expand( paste(outdir, "variable_regions.xlsx", sep='/') )

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
  
  variable.region.indices = which(sds >= sd.threshold & segvals$n.probes >= min.probes)

  variable.regions = segvals[variable.region.indices, patient.info$Samplename]
  colnames(variable.regions) = paste(paste(patient.info$'Endoscopy.Year', patient.info$'Pathology', sep='-'), " (",patient.info$Plate.Index, ")", sep='')
  #filename = paste(outdir, "/", patient.name, "_variable_regions_gamma", gamma2, '.txt', sep='')  
  print(paste("variable regions:", nrow(variable.regions)))

  if (nrow(variable.regions) > 0) {
    gain_loss = variable.regions
    for (i in 1:nrow(gain_loss)) {
      gain_loss[i,which(variable.regions[i,] > gain.threshold)] = 'GAIN'
      gain_loss[i,which(variable.regions[i,] < loss.threshold)] = 'LOSS'
      gain_loss[i,which(variable.regions[i,] > loss.threshold & variable.regions[i,] < gain.threshold)] = 'NEUTRAL'
    }
    
    write.xlsx2(cbind(segvals[variable.region.indices,c(1:5)], rbind( variable.regions, gain_loss ) ), file=filename, sheetName=patient.name, row.names=F, 
                append=appendFile)
    
    #write.table(cbind(segvals[variable.region.indices,c(1:5)], rbind( variable.regions, gain_loss ) ), sep='\t', quote=F, row.names=F, file=filename, append = F)
  } else {
    write.xlsx2(cbind(segvals[variable.region.indices,c(1:5)],  variable.regions ), file=filename, sheetName=patient.name, row.names=F, append=appendFile)
    #write.table(cbind(segvals[variable.region.indices,c(1:5)],  variable.regions ), sep='\t', quote=F, row.names=F, file=filename, append = F)
  }
  
  # sink(filename, append=T)
  # cat("\n\n## Thresholds applied:")
  # cat("## gamma2:", gamma2, "sd:", sd.threshold, "  min.probes:", min.probes, "\n")
  # cat("## Gain > 1.1, Loss < 0.9")
  # sink()  
  appendFile = T
}

write.xlsx2(cbind( 'gamma2'=gamma2, 'sd'=sd.threshold, 'min.probes'=min.probes, 'gain'=gain.threshold, 'loss'=loss.threshold), filename, sheetName='thresholds',row.names=F, append=T)



