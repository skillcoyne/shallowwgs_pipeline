args = commandArgs(trailingOnly=TRUE)

if (length(args) < 1)
  stop("Missing required params: <data dir> <sample info file> <model dir> <outdir> <alpha=0.9 DEF>")

suppressPackageStartupMessages( library(BarrettsProgressionRisk) )
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

pt = basename(datadir)
print(paste0('Patient ', pt))


outdir = paste0(outdir, '/', select.alpha, '/', pt)
dir.create(outdir, showWarnings = F, recursive = T)

print(paste0("Output path:", outdir))


x = list.files(model.dir, 'loo.Rdata', recursive = T, full.names = T)
load(x, verbose=F)
rm(plots,performance.at.1se,fits,pg.samp,coefs)

x = list.files(model.dir, 'all.pt.alpha.Rdata', recursive = T, full.names = T)
load(x, verbose=T)
fit = models[[select.alpha]]
s = performance.at.1se[[select.alpha]]$lambda  
cvRR = BarrettsProgressionRisk:::cvRR(dysplasia.df, coefs[[select.alpha]])

x = list.files(model.dir, 'model_data.Rdata', recursive = T, full.names = T)
load(x, verbose=F)
rm(dysplasia.df, coefs, labels)


be.fit.obj = BarrettsProgressionRisk:::be.model.fit(model=fit, s=s, tile.size=5e6, 
      tile.mean=z.mean, arms.mean=z.arms.mean, tile.sd=z.sd, arms.sd=z.arms.sd, 
      cx.mean=mn.cx, cx.sd=sd.cx, per.pt.nzcoefs = nzcoefs, cvRR = cvRR)


sheets = readxl::excel_sheets(info.file)[8:13]
info = do.call(bind_rows, lapply(sheets, function(s) {
  readxl::read_xlsx(info.file, s, trim_ws = T) %>% 
    dplyr::select(`Hospital Research ID`, `Block ID`, Endoscopy, Pathology, `SLX-ID`, `Index Sequence`, `Path Notes`) %>% dplyr::mutate(Sample = paste0(`SLX-ID`,'.',`Index Sequence`))
}))


segFiles = list.files(datadir, '[2|3|4]_segObj',  full.names = T, recursive = T)
print(segFiles)
if (length(segFiles) <= 0) stop(paste0('No segmentation file found in ', datadir))

preds = do.call(bind_rows, lapply(segFiles, function(f) {
  load(f)  
  segmented$sample.info = BarrettsProgressionRisk::loadSampleInformation(info %>% filter(Sample %in% segmented$sample.info$Sample) )
  
  prr = BarrettsProgressionRisk::predictRiskFromSegments(segmented,be.model=be.fit.obj,verbose=F)
  save(prr, file = paste0(outdir, '/predictions.Rdata'))
  predictions(prr)
cnPlot = BarrettsProgressionRisk::copyNumberMountainPlot(prr,T,F) 
ht = nrow(segmented$sample.info)
ggsave(file=paste0(outdir,'/copyNumberMtn.png'), plot=gridExtra::grid.arrange(cnPlot, top=pt), width=25, height=3*ht, dpi=300, units='in',limitsize=F)
}))



write_tsv(preds, path=paste0(outdir, '/', pt, '_predictions.tsv'))

print('Finished')
