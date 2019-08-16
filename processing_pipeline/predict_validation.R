args = commandArgs(trailingOnly=TRUE)

if (length(args) < 1)
  stop("Missing required params: <data dir> <sample info file> <model dir> <outdir> <alpha=0.9 DEF>")

library(BarrettsProgressionRisk)
source('~/workspace/shallowwgs_pipeline/lib/load_patient_metadata.R')

datadir = args[1]
info.file = args[2]
model.dir = args[3]
outdir = args[4]

select.alpha = '0.9'
if (length(args) == 5) {
  select.alpha = args[5]
  
  if (!select.alpha %in% c(0.0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0))
    stop("Alpha values available: 0.0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0")
}

outdir = paste0(outdir, '/', select.alpha)
dir.create(outdir, showWarnings = F, recursive = T)


x = list.files(model.dir, 'all.pt.alpha.Rdata', recursive = T, full.names = T)
load(x, verbose=F)
fit = models[[select.alpha]]
s = performance.at.1se[[select.alpha]]$lambda  


sheets = readxl::excel_sheets(info.file)[9:13]
info = do.call(bind_rows, lapply(sheets, function(s) {
  print(s)
  readxl::read_xlsx(info.file, s, trim_ws = T) %>% 
    dplyr::select(`Hospital Research ID`, `Block ID`, Endoscopy, Pathology, `SLX-ID`, `Index Sequence`, `Path Notes`) %>% dplyr::mutate(Sample = paste0(`SLX-ID`,'.',`Index Sequence`))
}))


segFiles = grep(paste(unique(info$`Hospital Research ID`),collapse='|'), list.files(datadir, '[2|3|4]_segObj',  full.names = T, recursive = T), value=T)

preds = do.call(bind_rows, lapply(segFiles, function(f) {
  load(f)  
  segmented$sample.info = BarrettsProgressionRisk::loadSampleInformation(info %>% filter(Sample %in% segmented$sample.info$Sample) )
  
  prr = BarrettsProgressionRisk::predictRiskFromSegments(segmented, model = fit, s = lambda, verbose = F)
  predictions(prr)
}))

outdir = paste0(outdir, '/', select.alpha, '/val_predictions.tsv')

write_tsv(preds, path=file)

print('Finished')
