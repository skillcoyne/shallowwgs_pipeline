# predict downsampled
args = commandArgs(trailingOnly=TRUE)

if (length(args) < 1)
  stop("Missing required params: <data dir> <model dir> ")

library(BarrettsProgressionRisk)
source('~/workspace/shallowwgs_pipeline/lib/load_patient_metadata.R')

dir = args[1]
outdir = args[2]

logT = F
if (length(args) > 2) logT = as.logical(args[3])

#dir = '~/Data/Ellie/QDNAseq/all_downsampled/20180206_KillcoyneS_RF_BarrettsCN/qdnaseq'
#outdir = '~/Data/Ellie/Analysis/downsampled'

files = list.files(dir, '.txt', full.names=T)

sampleNames = unique(sub('\\.binsize.*', '', basename(files), ignore.case = T))
dspred = data.frame(matrix(ncol=6, nrow=length(sampleNames), dimnames=list(sampleNames,c('Prob','RR', 'Risk', 'Adj.Prob','Adj.RR', 'Adj.Risk'))))


grep('raw',files, value=T)
grep('fitted',files, value=T)


for (name in sampleNames) {
  print(name)

  rawFile = grep(paste(name,'.*raw', sep=''), files, value=T)
  fittedFile = grep(paste(name,'.*fitted|corr', sep=''), files, value=T)
  
  if (length(rawFile) < 1 | length(fittedFile) < 1) {
    message(paste("Missing raw or fitted file for ", name, " skipping.",sep=''))
    next
  }
  pr = BarrettsProgressionRisk::predictRisk(path=dir, raw.file.grep=paste(name,'.*raw', sep=''), corrected.file.grep=paste(name,'.*fitted|corr', sep='') )

  dspred[name, ] = cbind( predictions(pr)[,c(2:4)], predictions(BarrettsProgressionRisk::adjustRisk(pr, 'mean'))[,c(2:4)])
}

save(dspred, file=paste(outdir, 'predictions.Rdata', sep='/'))

dspred$Sample = rownames(dspred)

y = c(0.3,0.5,0.8)
names(y) = c('Low','Moderate','High')
dspred$y = y[dspred$Risk]

pcal = showPredictionCalibration() + 
  geom_point(data=dspred, aes(Prob,y), color='grey39', shape=18, size=5, show.legend = F) + 
  geom_point(data=dspred, aes(Prob,y,color=Risk), shape=18, size=4, show.legend=F) +
  geom_text_repel(data=dspred, aes(Prob,y,label=Sample,angle=45), show.legend=F )

ggsave(filename=paste(outdir, 'predictions_plot.png',sep='/'), plot=pcal, width=10, height=10, units='in')

print("Finished'")

