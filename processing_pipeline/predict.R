args = commandArgs(trailingOnly=TRUE)

if (length(args) < 1)
  stop("Missing required params: <data dir> <sample info file> <model dir> <outdir> <alpha=0.9 DEF>")

library(BarrettsProgressionRisk)
source('~/workspace/shallowwgs_pipeline/lib/load_patient_metadata.R')

datadir = args[1]
info = args[2]
modeldir = args[3]
outdir = args[4]

select.alpha = '0.9'
if (length(args) == 5) {
  select.alpha = args[5]
 
  if (!select.alpha %in% c(0.0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0))
    stop("Alpha values available: 0.0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0")
}

outdir = paste0(outdir, '/', select.alpha)
dir.create(outdir, showWarnings = F, recursive = T)

print(paste0("Output path:", outdir))

#datadir = '~/Data/BarrettsProgressionRisk/Analysis/multipcf_perPatient/PR1_WSH_030/'
#info = '~/Data/BarrettsProgressionRisk/QDNAseq/training/All_patient_info.xlsx'

pt = basename(datadir)

info = BarrettsProgressionRisk::loadSampleInformation(
  read.patient.info(info)$info %>% 
  dplyr::filter(Hospital.Research.ID == pt | Patient == pt) %>% 
    dplyr::select(matches('ID|Endoscopy|Path|p53|Samplename|Block'), -matches('SLX')) %>% 
    dplyr::rename(Endoscopy = Endoscopy.Year, Sample = Samplename, GEJ.Distance = Block, `P53 IHC` = p53.Status)
)


x = list.files(modeldir, 'model_data.Rdata', recursive = T, full.names = T)
load(x, verbose=T)
rm(dysplasia.df, labels)

x = list.files(modeldir, 'all.pt.alpha.Rdata', recursive = T, full.names = T)
load(x, verbose=F)
fit = models[[select.alpha]]
s = performance.at.1se[[select.alpha]]$lambda  

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
prr = BarrettsProgressionRisk::predictRiskFromSegments(segs, model = fit, s = s, tile.mean = z.mean, tile.sd = z.sd, arms.mean = z.arms.mean, arms.sd = z.arms.sd, cx.mean = mn.cx, cx.sd = sd.cx, verbose = F)

pred.dir = paste0(outdir, '/predictions')
print(pred.dir)
if (!dir.exists(pred.dir)) dir.create(pred.dir, showWarnings=F, recursive = T)
save(prr, file=paste0(pred.dir, '/', pt, '.Rdata'))

predictions(prr, 'sample')

#prr2 = predictRiskFromSegments(segs, verbose = T)
#predictions(prr2, 'sample')

write_tsv(predictions(prr, 'sample'), paste0(outdir, '/', pt, '_preds.tsv'))


print('Finished')


