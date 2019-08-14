# predict downsampled
args = commandArgs(trailingOnly=TRUE)

if (length(args) < 4)
  stop("Missing required params: <data dir> <model dir> <outdir> <info file path> <same name>")

library(tidyverse)
library(BarrettsProgressionRisk)
#source('~/workspace/shallowwgs_pipeline/lib/load_patient_metadata.R')

dir = args[1]
model.dir = args[2]
outdir = args[3]
info = args[4]
sample = args[5]

if (!file.exists(info)) stop(paste0("File does not exist or is not readable: ", info))
all.ds.info = readxl::read_xlsx(info)

x = list.files(model.dir, 'all.pt.alpha.Rdata', recursive = T, full.names = T)
load(x, verbose=F)
fit = models$`0.9`
s = performance.at.1se$`0.9`$lambda  

#dir = '~/Data/Ellie/QDNAseq/all_downsampled/20180206_KillcoyneS_RF_BarrettsCN/qdnaseq'
#outdir = '~/Data/Ellie/Analysis/downsampled'

qndaseq.files = list.files(dir, pattern=sample, recursive = T, full.names = T)

#files = list.files(dir, '.txt', full.names=T)
#sampleNames = unique(sub('\\.binsize.*', '', basename(files), ignore.case = T))

all.preds = list(); dspred = tibble()
#for (name in sampleNames) {
#print(name)
  
info = loadSampleInformation( all.ds.info %>% filter(`Illumina ID` == sample) %>% 
                                dplyr::mutate(Endoscopy = "01/01/2001") %>% dplyr::rename(Sample = `Illumina ID`) )

rawFile = grep('*.raw', qndaseq.files, value=T)
fittedFile = grep('.*fitted|corr', qndaseq.files, value=T)
  
if (length(rawFile) < 1 | length(fittedFile) < 1) stop(paste0("Missing raw or fitted counts file for ", sample))
  
pred.dir = paste(outdir, sample, sep='/')
dir.create(pred.dir, showWarnings = F, recursive = T)
  
segObj = BarrettsProgressionRisk::segmentRawData(info, rawFile, fittedFile, cutoff=0.011)
  
p1 = BarrettsProgressionRisk::plotSegmentData(segObj, 'plot') 

ggsave(paste0(pred.dir, '/', sample, '_segmentedCoverage.png'), plot=grid.arrange(p1),  width=20, height=6, units='in', limitsize=F)

prr = BarrettsProgressionRisk::predictRiskFromSegments(segObj, model = fit, s = s)

save(prr, file=paste0(pred.dir, '/predictions.Rdata'))

print("Finished'")

