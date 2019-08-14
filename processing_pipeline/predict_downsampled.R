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
  
prr = BarrettsProgressionRisk::predictRiskFromSegments(segObj, model = fit, s = s)
  
save(prr, file=paste0(pred.dir, '/predictions.Rdata'))

#save(dspred, all.preds, file=paste(model.dir, 'downsampled-predictions.Rdata', sep='/'))

# dspred$Sample = rownames(dspred)
# 
# y = c(0.3,0.5,0.8)
# names(y) = c('Low','Moderate','High')
# dspred$y = y[dspred$Risk]
# 
# pcal = showPredictionCalibration() + 
#   geom_point(data=dspred, aes(Prob,y), color='grey39', shape=18, size=5, show.legend = F) + 
#   geom_point(data=dspred, aes(Prob,y,color=Risk), shape=18, size=4, show.legend=F) +
#   geom_text_repel(data=dspred, aes(Prob,y,label=Sample,angle=45), show.legend=F )
# 
# ggsave(filename=paste(outdir, 'predictions_plot.png',sep='/'), plot=pcal, width=10, height=10, units='in')

print("Finished'")

