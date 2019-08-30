args = commandArgs(trailingOnly=TRUE)

if (length(args) < 1)
  stop("Missing required params: <data dir> <sample info file> <model dir> <outdir> <alpha=0.9 DEF>")

suppressPackageStartupMessages( library(BarrettsProgressionRisk) )
suppressPackageStartupMessages( library(mclust) )
source('~/workspace/shallowwgs_pipeline/lib/load_patient_metadata.R')

smooth.EM.loess<-function(x, span=0.1, plot.dir = '.') {

  BIC = mclustBIC(x) 
  ds = densityMclust(x, x=BIC) 
  
  message(paste0('Clusters ', ds$G))
  png(filename=paste0(plot.dir,'/density.png'), width=800, height=600)
  plot(ds, what = "density", data = x, breaks = 50)
  dev.off()
  
  clusters = which(round(ds$parameters$mean,1) <= 1)
  #sG = which.max(table(ds$classification)[clusters])
  
  x[ds$classification %in% clusters] = loess( x[ds$classification %in% clusters]~c(1:sum(table(ds$classification)[clusters])), span=span )$fitted
  
  # for (i in 1:G) {
  #   if (table(ds$classification)[i] < 10) next
  #   x[ds$classification == i] = 
  #     loess( x[ds$classification == i]~c(1:table(ds$classification)[i]), span=span )$fitted
  # }
  
  return(x)  
}


fourier<-function(x, path='.') {
  fx = fft(x)
  #plot(Re(fx), type='l', xlim=c(0,20))
  fx[20:length(fx)] = 0+0i # Other than looking at the plot manually, how woudl I determine this?
  fxx = fft(fx, inverse = TRUE)/length(fx) 

  png(filename=paste0(path, '/fourier_transform.png'), width=600, height=800, units='px')
  layout(matrix(1:4,4,1)) 
  plot(x, type="l", main="Original Data")
  hist(x, breaks=50, main='')
  plot(Re(fxx),type="l", main="Fourier Transform Filtering") 
  hist(Re(fxx), breaks=50, main='')
  dev.off()
  
  return(Re(fxx))
}

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
nz = nzcoefs
rm(plots,performance.at.1se,fits,pg.samp,coefs)

x = list.files(model.dir, 'all.pt.alpha.Rdata', recursive = T, full.names = T)
load(x, verbose=T)
fit = models[[select.alpha]]
lambda = performance.at.1se[[select.alpha]]$lambda  
cvRR = BarrettsProgressionRisk:::cvRR(dysplasia.df, coefs[[select.alpha]])

x = list.files(model.dir, 'model_data.Rdata', recursive = T, full.names = T)
load(x, verbose=F)
rm(dysplasia.df, coefs, labels)

be.model = BarrettsProgressionRisk:::be.model.fit(fit, lambda, 5e6, z.mean, z.arms.mean, z.sd, z.arms.sd, mn.cx, sd.cx, nz, cvRR, NULL)


sheets = readxl::excel_sheets(info.file)[8:13]
info = do.call(bind_rows, lapply(sheets, function(s) {
  readxl::read_xlsx(info.file, s, trim_ws = T) %>% 
    dplyr::select(`Hospital Research ID`, `Block ID`, Endoscopy, Pathology, `SLX-ID`, `Index Sequence`, `Path Notes`) %>% dplyr::mutate(Sample = paste0(`SLX-ID`,'.',`Index Sequence`))
}))


segFiles = list.files(datadir, '[2|3|4]_segObj',  full.names = T, recursive = T)
print(segFiles)
if (length(segFiles) <= 0) stop(paste0('No segmentation file found in ', datadir))

preds = tibble()
for (f in segFiles) { 
  load(f)
  index = sub('_.*', '', basename(f))
  segmented$sample.info = BarrettsProgressionRisk::loadSampleInformation(info %>% filter(Sample %in% segmented$sample.info$Sample) )
  
  segmented$seg.vals = segmented$seg.vals %>% mutate_at(vars(matches(paste(segmented$sample.info$Sample, collapse='|'))), list(~fourier(.,path=outdir))  )

  prr = BarrettsProgressionRisk::predictRiskFromSegments(segmented,be.model,verbose=T)
print(predictions(prr))  
save(prr, file = paste0(outdir, '/', index, '_predictions.Rdata'))
  preds = bind_rows(preds, predictions(prr))
print("generating mtn plot")
cnPlot = BarrettsProgressionRisk::copyNumberMountainPlot(prr,T,F) 
print('saving plot')
ht = nrow(segmented$sample.info)
ggsave(file=paste0(outdir,'/', index, '_copyNumberMtn.png'), plot=gridExtra::grid.arrange(cnPlot, top=pt), width=25, height=3*ht, dpi=300, units='in',limitsize=F)
}

write_tsv(preds, path=paste0(outdir, '/', pt, '_predictions.tsv'))

print('Finished')
