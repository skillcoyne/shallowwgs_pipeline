
args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3)
  stop("Missing required arguments: <qdna data dir> <patient spreadsheet> <output dir> <patient name OPT>")


library(tidyverse)
library(ggrepel)
library(BarrettsProgressionRisk)
source('~/workspace/shallowwgs_pipeline/lib/load_patient_metadata.R')

data = args[1]
patient.file = args[2]
outdir = args[3]

data = '~/Data/BarrettsProgressionRisk//QDNAseq/kit_test/'
patient.file = '~/Data/BarrettsProgressionRisk/QDNAseq/kit_test/kit_pilot.xlsx'
outdir = '~/Data/BarrettsProgressionRisk/Analysis/kit_test'

patient.name = NULL
if (length(args) == 4)
  patient.name = args[4]

patient.info = readxl::read_xlsx(patient.file)
patient.info$Samplename = gsub('-', '_', paste(patient.info$`Index Sequence`,strip.whitespace(patient.info$`SLX-ID`),sep="_"))

colnames(patient.info)[1] = 'Patient ID'

#patient.info

data.dirs = list.dirs(data, full.names=T, recursive=F)
if (length(data.dirs) <= 0)
  stop(paste("No directories in", data))

pts_slx = arrange(unique(patient.info[c('SLX-ID','Patient ID')]), `Patient ID`)

preds = tibble()
failedQC = tibble()

for (slx in unique(patient.info$`SLX-ID`)) {
  print(slx)
  dir = grep(slx,data.dirs,value=T)
  
  rawFile = grep('raw',list.files(dir,'txt',full.names=T,recursive=T), value=T)
  fittedFile = grep('fitted',list.files(dir,'txt',full.names=T,recursive=T), value=T)
  
  if (length(rawFile) < 1 | length(fittedFile) < 1) {
    message(paste("Missing raw or fitted file for SLX ID",slx,"in",dir))
    next
  } else if (length(rawFile) == 1) {
    raw.data = data.table::fread(rawFile)
    fitted.data = data.table::fread(fittedFile)
  } else {
    raw.data = NULL; fitted.data = NULL
    for (file in rawFile) {
      print(file)
      raw = data.table::fread(file)
      fit = data.table::fread(grep(sub('\\.binSize.*','',basename(file)), fittedFile, value=T))
      colnames(fit)[5] = sub('\\.binSize.*','',basename(file))
      colnames(raw)[5] = sub('\\.binSize.*','',basename(file))
      
      if (is.null(fitted.data)) {
        fitted.data = fit
        raw.data = raw
      } else {
        fitted.data = merge(fitted.data, fit, by=c('location','chrom','start','end'), all=T) 
        raw.data = merge(raw.data, raw, by=c('location','chrom','start','end'), all=T) 
      }
    }
  }

  if ( length(which(grepl('tp',colnames(raw.data)))) > 0) {
    colnames(raw.data) = gsub('tp','', colnames(raw.data))
    colnames(fitted.data) = gsub('tp','', colnames(fitted.data))
  }

  for (s in (patient.info %>% filter(`SLX-ID` == slx))$Sample) {
    #message(paste('Patient',id))
    message(paste('Sample',s))
    
    pid = patient.info %>% filter(Sample == s) %>% select(`Patient ID`) %>% pull
    
    samples = patient.info %>% filter(Sample == s & `SLX-ID` == slx) %>% mutate(`Index Sequence` = sub('-', '_', gsub('tp', '', `Index Sequence`)), Samplename = paste(`SLX-ID`, `Index Sequence`, sep='.') ) %>% select(Samplename) %>% pull
    
    rcols = grep(paste(samples,collapse='|'), colnames(raw.data))
    fcols = grep(paste(samples,collapse='|'), colnames(fitted.data))
      
    rd = raw.data[,c(1:4,rcols),with=F]
    fd = fitted.data[,c(1:4,fcols),with=F]
    
    plot.dir = paste(outdir, pid, 'plots',sep='/')
    dir.create(plot.dir, showWarnings = F, recursive = T)
    
    tryCatch({
      
      si = patient.info %>% filter(Sample == s) %>% mutate(Endoscopy = '2019/01/01')
      
      segmented = BarrettsProgressionRisk::segmentRawData(loadSampleInformation(si),rd,fd,verbose=F)
      residuals = BarrettsProgressionRisk::sampleResiduals(segmented)
      
      failedQC = bind_rows(failedQC, residuals)
      
      prr = BarrettsProgressionRisk::predictRiskFromSegments(segmented)
      preds = bind_rows(preds, predictions(prr))
      
      #plot = BarrettsProgressionRisk::plotCorrectedCoverage(segmented, 'plot')
      ggsave(paste(plot.dir, paste(s, 'segmentedCoverage.png',sep='_'), sep='/'),  plot=grid.arrange(BarrettsProgressionRisk::plotSegmentData(segmented,'plot')), height=4*length(rcols), width=20, units='in')
      
      readr::write_tsv(residuals, path=paste(dirname(plot.dir), paste(s,'residuals.txt',sep='_'), sep='/'))
      save(segmented, prr, paste(dirname(plot.dir), 'riskObj.Rdata',sep='/'))
      
      if (nrow(subset(residuals, !Pass)) > 0)
        failedQC = bind_rows(failedQC,cbind.data.frame('Patient ID'=id,subset(residuals, !Pass), stringsAsFactors=F))
      
      if (nrow(subset(residuals,Pass)) > 0) {
        pr = BarrettsProgressionRisk::predictRiskFromSegments(segmented, verbose=F)
        preds = predictions(pr)
        } else {
        message(paste('All samples failed QC for patient ID',id))
        }
    }, error = function(e) {
      message(paste("Error in risk prediction for patient",id,',skipping:\n\t',e))
    })
  }
}


y = c(0.3,0.5,0.8)
names(y) = c('Low','Moderate','High')
preds$y = y[preds$Risk]

pcal = showPredictionCalibration() + 
  geom_point(data=preds, aes(Probability,y), color='grey39', shape=18, size=5, show.legend = F) + 
  geom_point(data=preds, aes(Probability,y,color=Risk), shape=18, size=4, show.legend=F) +
  geom_text_repel(data=preds, aes(Probability,y,label=`Patient ID`,angle=45), show.legend=F )

ggsave(filename=paste(outdir, 'predictions_plot.png',sep='/'), plot=pcal, width=10, height=10, units='in')


write.table(preds, sep='\t', quote=F, row.names=F, file=paste(outdir, 'predictions.txt', sep='/'))
write.table(failedQC, sep='\t', quote=F, row.names=F, file=paste(outdir, 'failed_qc_samples.txt', sep='/'))
  
print("Finished")
