args = commandArgs(trailingOnly=TRUE)

if (length(args) < 1)
  stop("Missing required params: <data dir> <sample info file> <model dir> <outdir> <alpha=0.9 DEF>")

suppressPackageStartupMessages( library(BarrettsProgressionRisk) )
#suppressPackageStartupMessages( library(mclust) )
#source('~/workspace/shallowwgs_pipeline/lib/load_patient_metadata.R')

datadir = args[1]
# datadir = '~/Data/BarrettsProgressionRisk/Analysis/validation/multipcf/AHM1400'
info.file = args[2]
# info.file = '~/Data/BarrettsProgressionRisk/QDNAseq/validation/sWGS_validation_batches.xlsx'
model.dir = args[3]
# model.dir = '~/Data/BarrettsProgressionRisk/Analysis/5e6_arms'
outdir = args[4]
# outdir = '/tmp'
training.dir = args[5]
# training.dir = '~/Data/BarrettsProgressionRisk/Analysis/multipcf_perPatient'

select.alpha = '0.9'
if (length(args) == 6) {
  select.alpha = args[6]
  
  if (!select.alpha %in% c(0.0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0))
    stop("Alpha values available: 0.0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0")
}

# ----- Load ALL data

files = list.files(training.dir, '5e06_cleaned_tiled', recursive=T, full.names = T)
if (length(files) <= 0) stop(paste0("No tiled files in ", training.dir))
training.tiles.5mb = do.call(bind_rows, purrr::map(files, function(f) {
  read_tsv(f,col_types=c(.default=col_double()))
})) %>% dplyr::rename('Sample' = 'X1')

files = list.files(training.dir, 'arms_cleaned_tiled', recursive=T, full.names = T)
if (length(files) <= 0) stop(paste0("No tiled arm files in ", training.dir))
training.tiles.arms = do.call(bind_rows, purrr::map(files, function(f) {
  read_tsv(f,col_types=c(.default=col_double()))
})) %>% dplyr::rename('Sample' = 'X1')

if (nrow(training.tiles.5mb) != nrow(training.tiles.arms)) stop("5MB files don't match arm files in training directory.")

## Val
files = list.files(dirname(datadir), '5e06_tiles', recursive=T, full.names = T)
val.tiles.5e6 = do.call(bind_rows, purrr::map(files, function(f) {
  read_tsv(f,col_types=c(.default=col_double()))
})) %>% dplyr::rename('Sample' = 'X1')

files = list.files(dirname(datadir), 'arm_tiles', recursive=T, full.names = T)
val.tiles.arms = do.call(bind_rows, purrr::map(files, function(f) {
  read_tsv(f,col_types=c(.default=col_double())) %>% dplyr::filter(X1 != '')
})) %>% dplyr::rename('Sample' = 'X1')

if (nrow(val.tiles.5e6) != nrow(val.tiles.arms)) stop("5MB files don't match arm files in validation directory.")

# Rank adjust validation
find.rank<-function(r,x) {
  if ( nrow(x %>% filter(rank == r)) > 0 ) return(x %>% filter(rank == r))
  find.rank(r-1,x)
}

rank.adjust<-function(train,val) {
  for (n in 2:ncol(train)) {
    col = colnames(train)[n]
    
    rk = sapply(seq(1, length(train[[col]]), 4), function(s) mean(train[[col]][s:(s+4)], na.rm=T))
    tt.rank = tibble(value=rk, rank=rank(rk)) %>% arrange(rank)
    vt.rank = floor(rank(val[[col]]))
    
    for (i in 1:length(vt.rank)) {
      r = vt.rank[i]
      if (is.na(mean(find.rank(r,tt.rank)$value,na.rm=T))) stop(paste(col, r))
      val[i, col] = mean(find.rank(r,tt.rank)$value,na.rm=T)
    }
  }
return(val)
}

#print("Rank adjusting validation cohort")
#val.tiles.5e6 = rank.adjust(training.tiles.5mb, val.tiles.5e6)
#val.tiles.arms = rank.adjust(training.tiles.arms, val.tiles.arms)


# ------

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

sheets = readxl::excel_sheets(info.file)[8:13]
info = do.call(bind_rows, lapply(sheets, function(s) {
  readxl::read_xlsx(info.file, s, trim_ws = T) %>% 
    dplyr::select(`Hospital Research ID`, `Block ID`, Endoscopy, Pathology, `SLX-ID`, `Index Sequence`, `Path Notes`) %>% dplyr::mutate(Sample = paste0(`SLX-ID`,'.',`Index Sequence`))
}))

## Means!
val.mean = apply(as.matrix(val.tiles.5e6[,-1]),2,mean,na.rm=T)
val.sd = apply(as.matrix(val.tiles.5e6[,-1]),2,sd,na.rm=T)

arm.mean = apply(as.matrix(val.tiles.arms[,-1]),2,mean,na.rm=T)
arm.sd = apply(as.matrix(val.tiles.arms[,-1]),2,sd,na.rm=T)

info = info %>% dplyr::filter(`Hospital Research ID` == pt) %>% arrange(Sample)

val.tiles = val.tiles.5e6 %>% filter(Sample %in% info$Sample) %>% arrange(Sample)
val.tiles = as.matrix(val.tiles[,-1])
rownames(val.tiles) = info$Sample

for (i in 1:ncol(val.tiles))
  val.tiles[,i] = BarrettsProgressionRisk:::unit.var(val.tiles[,i], val.mean[i], val.sd[i])
print('val unit.var done')

val.arms = val.tiles.arms %>% filter(Sample %in% info$Sample) %>% arrange(Sample)
val.arms = as.matrix(val.arms[,-1])
rownames(val.arms) = info$Sample

for (i in 1:ncol(val.arms))
  val.arms[,i] = BarrettsProgressionRisk:::unit.var(val.arms[,i], arm.mean[i], arm.sd[i])
print('val arm unit.var done')

cx = BarrettsProgressionRisk:::scoreCX(as.matrix(val.tiles[,-1]),1)
print(paste('cx', cx))

val.df = cbind(BarrettsProgressionRisk::subtractArms(val.tiles, val.arms), 'cx'= BarrettsProgressionRisk:::unit.var(cx, mn.cx, sd.cx))
be.model = BarrettsProgressionRisk:::be.model.fit(fit, lambda, 5e6, val.mean, arm.mean, val.sd, arm.sd, mn.cx, sd.cx, nz, cvRR, NULL)

print("Loading segment file")
segFiles = list.files(datadir, '[2|3|4]_segObj',  full.names = T, recursive = T)
print(segFiles)
if (length(segFiles) <= 0) stop(paste0('No segmentation file found in ', datadir))

preds = tibble()
for (f in segFiles) { 
  load(f)
  index = sub('_.*', '', basename(f))
  segmented$sample.info = BarrettsProgressionRisk::loadSampleInformation(info %>% filter(Sample %in% segmented$sample.info$Sample) )
  
  tiles = BarrettsProgressionRisk:::tileSamples(segmented, verbose=F)
  tiles$tiles = val.df
  
  prr = BarrettsProgressionRisk:::predictRisk(segmented, tiles, be.model)
  predictions(prr)

  #segmented$seg.vals = segmented$seg.vals %>% mutate_at(vars(matches(paste(segmented$sample.info$Sample, collapse='|'))), list(~fourier(.,path=outdir))  )
  #prr = BarrettsProgressionRisk::predictRiskFromSegments(segmented,be.model,verbose=T)
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
