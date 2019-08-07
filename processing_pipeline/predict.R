args = commandArgs(trailingOnly=TRUE)

if (length(args) < 1)
  stop("Missing required params: <data dir> <sample info file> <model dir> <outdir> ")

library(BarrettsProgressionRisk)
source('~/workspace/shallowwgs_pipeline/lib/load_patient_metadata.R')

datadir = args[1]
info = args[2]
modeldir = args[3]
outdir = args[4]

datadir = '~/Data/BarrettsProgressionRisk/Analysis/multipcf_perPatient/PR1_WSH_030/'
info = '~/Data/BarrettsProgressionRisk/QDNAseq/training/All_patient_info.xlsx'

pt = basename(datadir)

info = BarrettsProgressionRisk::loadSampleInformation(
  read.patient.info(info)$info %>% 
  dplyr::filter(Hospital.Research.ID == pt | Patient == pt) %>% 
    dplyr::select(matches('ID|Endoscopy|Path|p53|Samplename|Block'), -matches('SLX')) %>% 
    dplyr::rename(Endoscopy = Endoscopy.Year, Sample = Samplename, GEJ.Distance = Block, `P53 IHC` = p53.Status)
)


x = list.files(modeldir, 'all.pt.alpha.Rdata', recursive = T, full.names = T)
load(x, verbose=F)
fit = models$`0.9`
s = performance.at.1se$`0.9`$lambda  

loadRData<-function(filename) {
  load(filename, verbose=F)
  get(ls()[ls() != "filename"])
}

  
segFile = list.files(datadir, 'Rdata', recursive = T, full.names = T)


if (length(segFile) <= 0) {
  raw = list.files(datadir, 'raw', full.names = T, recursive = T)
  fit = list.files(datadir, 'fit', full.names = T, recursive = T)
  
  segs = BarrettsProgressionRisk::segmentRawData(info, raw, fit)
} else if (length(segFile) > 0) {
  message(paste0("Loading segment data: ", segFile))
  segs = loadRData(segFile)  
  
  segs$sample.info = info
}

prr = predictRiskFromSegments(segs, model = fit, s = s)
predictions(prr, 'sample')

#prr2 = predictRiskFromSegments(segs, verbose = T)
#predictions(prr2, 'sample')

write_tsv(predictions(prr, 'sample'), paste0(outdir, '/', pt, '_preds.tsv'))


print('Finished')


