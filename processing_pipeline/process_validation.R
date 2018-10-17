
args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3)
  stop("Missing required arguments: <qdna data dir> <patient spreadsheet> <output dir> <patient name OPT>")


library(plyr)
library(ggplot2)
library(ggrepel)
library(BarrettsProgressionRisk)
source('~/workspace/shallowwgs_pipeline/lib/load_patient_metadata.R')

data = args[1]
patient.file = args[2]
outdir = args[3]

# data = '~/Data/Ellie/QDNAseq/validation'
# patient.file = '~/Data/Ellie/Cassandra/ValidationCohort.xlsx'
# outdir = '~/Data/Ellie/Analysis/VAL_Cohort'

patient.name = NULL
if (length(args) == 4)
  patient.name = args[4]

patient.info = readxl::read_xlsx(patient.file)
patient.info$Samplename = gsub('-', '_', paste(patient.info$`Index Sequence`,strip.whitespace(patient.info$`SLX-ID`),sep="_"))

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
  
  rawFile = grep('raw',list.files(dir,'txt',full.names=T), value=T)
  fittedFile = grep('fitted',list.files(dir,'txt',full.names=T), value=T)
  
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
  
  ids = unique(subset(patient.info, `SLX-ID` == slx)$`Patient ID`)
  for (id in ids) {
    message(paste('Patient',id))
    samples = sub('-','_',subset(patient.info, `Patient ID` == id & `SLX-ID` == slx)$`Index Sequence`)

    rcols = grep(paste(samples,collapse='|'), colnames(raw.data))
    fcols = grep(paste(samples,collapse='|'), colnames(fitted.data))
      
    rd = raw.data[,c(1:4,rcols),with=F]
    fd = fitted.data[,c(1:4,fcols),with=F]
    
    tryCatch({
      segmented = BarrettsProgressionRisk::segmentRawData(rd,fd,verbose=F)
      residuals = BarrettsProgressionRisk::sampleResiduals(segmented)
      
      if (nrow(subset(residuals, !Pass)) > 0)
        failedQC = bind_rows(failedQC,cbind.data.frame('Patient ID'=id,subset(residuals, !Pass), stringsAsFactors=F))
      
      if (nrow(subset(residuals,Pass)) > 0) {
        pr = BarrettsProgressionRisk::predictRiskFromSegments(segmented, verbose=F)
        adj = predictions(adjustRisk(pr, 'mean'))[2:4]
        colnames(adj) = paste('Adjusted', colnames(adj))
        preds = as_tibble(rbind(as.data.frame(preds), cbind('Patient ID'=id, predictions(pr), adj)))
    
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
