# predict downsampled
args = commandArgs(trailingOnly=TRUE)

if (length(args) < 1)
  stop("Missing required params: <data dir> <model dir> <outdir>")

library(BarrettsProgressionRisk)
source('~/workspace/shallowwgs_pipeline/lib/load_patient_metadata.R')

dir = args[1]
model.dir = args[2]
outdir = args[3]

logT = F

#dir = '~/Data/Ellie/QDNAseq/all_downsampled/20180206_KillcoyneS_RF_BarrettsCN/qdnaseq'
#outdir = '~/Data/Ellie/Analysis/downsampled'

files = list.files(dir, '.txt', full.names=T)

sampleNames = unique(sub('\\.binsize.*', '', basename(files), ignore.case = T))

all.ds.info = readxl::read_xlsx('~/Data/BarrettsProgressionRisk/QDNAseq/all_downsampled/downsampled_ids.xlsx')

all.preds = list(); dspred = tibble()
for (name in sampleNames) {
  print(name)
  
  info = loadSampleInformation(
    all.ds.info %>% filter(`Illumina ID` == name) %>% dplyr::mutate(Endoscopy = "01/01/2001") %>% dplyr::rename(Sample = `Illumina ID`) )

  rawFile = grep(paste(name,'.*raw', sep=''), files, value=T)
  fittedFile = grep(paste(name,'.*fitted|corr', sep=''), files, value=T)
  
  pred.dir = paste(outdir, name, sep='/')
  dir.create(pred.dir, showWarnings = F, recursive = T)
  
  if (length(rawFile) < 1 | length(fittedFile) < 1) {
    message(paste("Missing raw or fitted file for ", name, " skipping.",sep=''))
    next
  }
  segObj = BarrettsProgressionRisk::segmentRawData(info, rawFile, fittedFile,cutoff=0.011)
  
  prr = BarrettsProgressionRisk::predictRiskFromSegments(segObj)
  
  save(prr, file=paste0(pred.dir, '/predictions.Rdata'))

  #all.preds[[name]] = prr
    
  #dspred = bind_rows(dspred, predictions(prr))
}
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

