# predict downsampled
args = commandArgs(trailingOnly=TRUE)

if (length(args) < 4)
  stop("Missing required params: <data dir> <model dir> <outdir> <info file path> <same name> <alpha=0.9 DEF")

library(tidyverse)
library(gridExtra)
library(BarrettsProgressionRisk)
#source('~/workspace/shallowwgs_pipeline/lib/load_patient_metadata.R')

dir = args[1]
model.dir = args[2]
outdir = args[3]
info = args[4]

#sample = args[5]


# dir = '~/Data/BarrettsProgressionRisk/QDNAseq/all_downsampled/100kb'
# info = '~/Data/BarrettsProgressionRisk/QDNAseq/all_downsampled/downsampled_ids.xlsx'
# model.dir = '~/Data/BarrettsProgressionRisk/Analysis/models_5e6/100kb'
# outdir = '~/tmp/100kb'
# sample = 'LP6008280-DNA_A06'

select.alpha = '0.9'
if (length(args) == 5) {
  select.alpha = args[5]
  
  if (!select.alpha %in% c(0.0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0))
    stop("Alpha values available: 0.0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0")
}
outdir = paste0(outdir, '/', select.alpha)
dir.create(outdir, showWarnings = F, recursive = T)

if (!file.exists(info)) stop(paste0("File does not exist or is not readable: ", info))
all.ds.info = readxl::read_xlsx(info)

x = list.files(model.dir, paste0('loo_',select.alpha,'.Rdata'), recursive = T, full.names = T)
message(paste0('LOO: ', x))
load(x, verbose=F)
nz = nzcoefs
rm(plots,performance.at.1se,fits,pg.samp,coefs)

x = list.files(model.dir, 'all.pt.alpha.Rdata', recursive = T, full.names = T)
message(paste0('MODEL INFO: ',x))
load(x, verbose=T)
fit = models[[select.alpha]]
lambda = performance.at.1se[[select.alpha]]$lambda  
cvRR = BarrettsProgressionRisk:::cvRR(dysplasia.df, coefs[[select.alpha]])

x = list.files(model.dir, 'model_data.Rdata', recursive = T, full.names = T)
message(paste0('MODEL DATA: ',x))
load(x, verbose=T)
rm(dysplasia.df, coefs, labels)

be.model = BarrettsProgressionRisk:::be.model.fit(fit, lambda, 5e6, z.mean, z.arms.mean, z.sd, z.arms.sd, mn.cx, sd.cx, nz, cvRR, NULL)

qdnaseq.files = list.files(dir, 'raw|fitted', recursive = T, full.names = T)
#qndaseq.files = list.files(dir, pattern=sample, recursive = T, full.names = T)

all.preds = list(); dspred = tibble()
#for (name in sampleNames) {
#print(name)

rawFile = grep('raw', qdnaseq.files, value=T)
fittedFile = grep('.*fitted|corr', qdnaseq.files, value=T)
  
if (length(rawFile) < 1 | length(fittedFile) < 1) stop(paste0("Missing raw or fitted counts file for ", sample))
  
pred.dir = paste(outdir, sample, sep='/')
dir.create(pred.dir, showWarnings = F, recursive = T)


raw.data = read_tsv(rawFile, col_types = cols(chrom = col_character(), location = col_character(), .default = col_double()))
fit.data = read_tsv(fittedFile, col_types = cols(chrom = col_character(), location = col_character(), .default = col_double()))

for (sample in all.ds.info$`Illumina ID`) {
  message(paste0('Processing sample ', sample))
  info = loadSampleInformation( all.ds.info %>% filter(`Illumina ID` == sample) %>% dplyr::mutate(Endoscopy = "01/01/2001") %>% dplyr::rename(Sample = `Illumina ID`) )
  
  kb = as.integer(sub('kb','',basename(dir)))

  raw = raw.data %>% dplyr::select(location,chrom,start,end,!!sample)
  fit = fit.data %>% dplyr::select(location,chrom,start,end,!!sample)
    
  segmented = BarrettsProgressionRisk::segmentRawData(info, raw, fit, verbose=T, kb = kb, multipcf = F )

  p1 = BarrettsProgressionRisk::plotSegmentData(segmented, 'plot') 
  ggsave(paste0(pred.dir, '/', sample, '_segmentedCoverage.png'), plot=grid.arrange(p1),  width=20, height=6, units='in', limitsize=F)
  
  prr = BarrettsProgressionRisk::predictRiskFromSegments(segmented, be.model)
  p2 = copyNumberMountainPlot(prr,T,F)
  ggsave(paste0(pred.dir, '/', sample, '_copyNumberMtn.png'), plot=grid.arrange(p2),  width=20, height=6, units='in', limitsize=F)

  dspred = bind_rows(dspred, predictions(prr))
}
  
save(prr, file=paste0(pred.dir, '/predictions.Rdata'))

readr::write_tsv(predictions(prr), path=paste0(outdir,'/predictions.tsv'))

print("Finished'")

