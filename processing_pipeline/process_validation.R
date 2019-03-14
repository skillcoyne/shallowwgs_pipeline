
args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3)
  stop("Missing required arguments: <qdna data dir> <patient spreadsheet> <output dir> <patient name OPT>")


library(tidyverse)
library(ggrepel)
library(gridExtra)
library(BarrettsProgressionRisk)
source('~/workspace/shallowwgs_pipeline/lib/load_patient_metadata.R')

data = args[1]
patient.file = args[2]
outdir = args[3]

#data = '~/Data/BarrettsProgressionRisk/QDNAseq/validation/'
#patient.file = '~/Data/BarrettsProgressionRisk/QDNAseq/validation/validation_samples.xlsx'
#outdir = '~/Data/BarrettsProgressionRisk/Analysis/val_test/'

patient.name = NULL
if (length(args) == 4)
  patient.name = args[4]

patient.info = readxl::read_xlsx(patient.file)
patient.info$Samplename = gsub('-', '_', paste(patient.info$`Index Sequence`,strip.whitespace(sub('SLX-','',patient.info$`SLX-ID`)),sep="_"))

#patient.info

data.dirs = list.dirs(data, full.names=T, recursive=F)
if (length(data.dirs) <= 0)
  stop(paste("No directories in", data))

pts_slx = arrange(unique(patient.info[c('SLX-ID','Patient ID')]), `Patient ID`)

preds = tibble(); failedQC = tibble()

merged.raw = NULL; merged.fit = NULL
for (slx in unique(pts_slx$`SLX-ID`)) {
#for (slx in unique(patient.info$`SLX-ID`)) {
  print(slx)
  dir = grep(slx,data.dirs,value=T)
  
  rawFile = grep('raw',list.files(dir,'txt',full.names=T,recursive=T), value=T)
  fittedFile = grep('fitted',list.files(dir,'txt',full.names=T,recursive=T), value=T)
  
  if (length(rawFile) < 1 | length(fittedFile) < 1) {
    message(paste("Missing raw or fitted file for SLX ID",slx,"in",dir))
    next
  } else if (length(rawFile) == 1) {
    raw.data = read_tsv(rawFile, col_types = cols('chrom'=col_character()))
    fitted.data = read_tsv(fittedFile, col_types = cols('chrom'=col_character()))
  } else {
    raw.data = NULL; fitted.data = NULL
    for (file in rawFile) {
      print(file)
      raw = read_tsv(file, col_types = cols('chrom'=col_character()))
      fit = read_tsv(grep(sub('\\.binSize.*','',basename(file)), fittedFile, value=T), col_types = cols('chrom'=col_character()))

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
  
  if (is.null(merged.fit)) {
    merged.fit = fitted.data
    merged.raw = raw.data
  } else {
    merged.fit = full_join(merged.fit, fitted.data, by=c('location','chrom','start','end'))
    merged.raw = full_join(merged.raw, raw.data, by=c('location','chrom','start','end'))
  }
  
  print(dim(merged.fit))
}

#patients = (pts_slx %>% filter(`SLX-ID` == slx)) %>% select(`Patient ID`) %>% pull
patients = pts_slx$`Patient ID`
for (pid in  patients) {
  message(paste('Patient',pid))
    
  si = patient.info %>% filter(`Patient ID` == pid) 
    
  plot.dir = paste(outdir, pid, 'plots',sep='/')
  if (length(list.files(plot.dir)) >= nrow(si)) next  # skip patients I've already done

  dir.create(plot.dir, showWarnings = F, recursive = T)
  si = si %>% mutate(Sample = paste(`SLX-ID`, gsub('-','_',si$`Index Sequence`), sep='.'))

  samples = si$Sample

  rcols = grep(paste(samples,collapse='|'), colnames(merged.raw))
  fcols = grep(paste(samples,collapse='|'), colnames(merged.fit))
    
  rd = merged.raw %>% select(location,chrom,start,end,!!rcols)
  fd = merged.fit %>% select(location,chrom,start,end,!!fcols)

  rd = rd %>% dplyr::rename_at(vars(matches(paste(samples,collapse='|'))), funs( str_match(.,paste(samples,collapse='|'))))
  fd = fd %>% dplyr::rename_at(vars(matches(paste(samples,collapse='|'))), funs( str_match(.,paste(samples,collapse='|'))))
    
  tryCatch({
      
    if (is.null(si[['Endoscopy']])) si = si %>% mutate(Endoscopy = '2019/01/01')

    segmented = BarrettsProgressionRisk::segmentRawData(loadSampleInformation(si),rd,fd,verbose=F)
    residuals = BarrettsProgressionRisk::sampleResiduals(segmented)
      
    failedQC = bind_rows(failedQC, residuals)
      
    #prr = BarrettsProgressionRisk::predictRiskFromSegments(segmented)
    #preds = bind_rows(preds, predictions(prr))
      
    #plots = BarrettsProgressionRisk::plotCorrectedCoverage(segmented, 'list')
    plots = BarrettsProgressionRisk::plotSegmentData(segmented, 'list')
    for (s in names(plots))
      ggsave(paste(plot.dir, paste(s, 'segmentedCoverage.png',sep='_'), sep='/'),  plot=plots[[s]], height=4, width=20, units='in')
      
    ggsave(paste(plot.dir, paste(pid, 'segmentedCoverage.png',sep='_'), sep='/'),  plot=do.call(grid.arrange, c(plots,ncol=1)), height=4*length(rcols), width=20, units='in')
      
    readr::write_tsv(residuals, path=paste(dirname(plot.dir), 'residuals.txt',sep='/'))
    save(segmented, file=paste(dirname(plot.dir), 'segObj.Rdata',sep='/'))
      
    if (nrow(subset(residuals, !Pass)) > 0)
      failedQC = bind_rows(failedQC, residuals %>% filter(!Pass) %>% mutate(`Patient ID` = pid) %>% select(`Patient ID`, everything()) )
      
    # if (nrow(subset(residuals,Pass)) > 0) {
    #   pr = BarrettsProgressionRisk::predictRiskFromSegments(segmented, verbose=F)
    #   preds = predictions(pr)
    #   } else {
    #   message(paste('All samples failed QC for patient ID',pid))
    #   }
  }, error = function(e) {
    message(paste("Error in risk prediction for patient",id,',skipping:\n\t',e))
  })
}



# y = c(0.3,0.5,0.8)
# names(y) = c('Low','Moderate','High')
# preds$y = y[preds$Risk]
# 
# pcal = showPredictionCalibration() + 
#   geom_point(data=preds, aes(Probability,y), color='grey39', shape=18, size=5, show.legend = F) + 
#   geom_point(data=preds, aes(Probability,y,color=Risk), shape=18, size=4, show.legend=F) +
#   geom_text_repel(data=preds, aes(Probability,y,label=`Patient ID`,angle=45), show.legend=F )
# 
# ggsave(filename=paste(outdir, 'predictions_plot.png',sep='/'), plot=pcal, width=10, height=10, units='in')
# 
# 
# write.table(preds, sep='\t', quote=F, row.names=F, file=paste(outdir, 'predictions.txt', sep='/'))
# write.table(failedQC, sep='\t', quote=F, row.names=F, file=paste(outdir, 'failed_qc_samples.txt', sep='/'))
#   
# print("Finished")
