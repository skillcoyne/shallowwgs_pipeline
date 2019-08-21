
args = commandArgs(trailingOnly=TRUE)

if (length(args) < 5)
  stop("Missing required arguments: <qdna data dir> <patient spreadsheet> <model.dir> <output dir> <patient>")


suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(BarrettsProgressionRisk))
#source('~/workspace/shallowwgs_pipeline/lib/load_patient_metadata.R')

data = args[1]
val.file = args[2]
model.dir = args[3]
outdir = args[4]
patient = args[5]


patient = str_replace( unlist(str_split(patient, ',| ')), ' ', '')
patient = str_replace_all( patient, '/', '_')

if (dir.exists(outdir))
  stop(paste0("Directory already exists: ", outdir))

dir.create(outdir, recursive = T, showWarnings = F)



select.alpha = '0.9'
if (length(args) == 6) {
  select.alpha = args[6]
  
  if (!select.alpha %in% c(0.0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0))
    stop("Alpha values available: 0.0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0")
}

x = list.files(model.dir, 'model_data.Rdata', recursive = T, full.names = T)
load(x, verbose=F)
rm(dysplasia.df, labels)

x = list.files(model.dir, 'all.pt.alpha.Rdata', recursive = T, full.names = T)
load(x, verbose=F)
fit = models[[select.alpha]]
lambda = performance.at.1se[[select.alpha]]$lambda  


# data = '~/Data/BarrettsProgressionRisk/QDNAseq/validation/'
# val.file = '~/Data/BarrettsProgressionRisk/QDNAseq/validation/sWGS_validation_batches.xlsx'
# outdir = '~/Data/BarrettsProgressionRisk/Analysis/validation'

#patients = c('OCCAMS_AH_100','AHM0185','AHM0555','AHM0348','AHM0281','AHM0500','AHM0546','AHM0584','AHM1432','AHM1471','AHM0567','AHM0570')


sheets = readxl::excel_sheets(val.file)[1:13]


all.val = do.call(bind_rows, lapply(sheets, function(s) {
  readxl::read_xlsx(val.file, s) %>% select(`Hospital Research ID`, matches('Status'), `Sample Type`, `SLX-ID`, `Index Sequence`, Cohort, Batch, RA) %>% mutate_at(vars(`SLX-ID`), list(as.character)) %>% dplyr::filter(!is.na(`SLX-ID`))
}))

pastefun<-function(x) {
  if ( !grepl('SLX-', x) ) x = paste0('SLX-',x)
  return(x)
}
all.val = all.val %>% rowwise %>% mutate_at(vars(`SLX-ID`), list(pastefun) ) %>% ungroup
#all.val = all.val %>% arrange(Batch, `Hospital Research ID`) %>% group_by(`Hospital Research ID`) %>% mutate(AID = group_indices()) %>% ungroup
all.val = all.val %>% mutate(`Hospital Research ID` = str_replace_all( str_remove_all(`Hospital Research ID`, " "), '/', '_'), `Index Sequence` = str_replace_all(`Index Sequence`, 'tp', ''))
all.val = all.val %>% mutate(Samplename = paste(`SLX-ID`, `Index Sequence`, sep='.'), RA = factor(RA))

if (!file.exists(  paste0(data, '/merged_raw_fit.Rdata')))
  stop("Missing merged data files, run merge_qdnaseq_data.R' first.")

load(file=paste0(data, '/merged_raw_fit.Rdata'))

si = si %>% filter(`Hospital Research ID` == patient)

if (nrow(si) <= 0)
  stop(paste0('No patient ', patient, ' found in spreadsheet.'))

message(paste('Patient',patient))

si$Sample = si$Samplename 
if (is.null(si[['Endoscopy']])) si = si %>% mutate(Endoscopy = '2019/01/01')

plot.dir = paste(outdir, patient, 'plots',sep='/')

dir.create(plot.dir, showWarnings = F, recursive = T)

residuals = tibble(); predictions = tibble()
for (sample in si$Samplename) {
  rcols = grep(paste(sample,collapse='|'), colnames(merged.raw))
  fcols = grep(paste(sample,collapse='|'), colnames(merged.fit))
  
  if (length(rcols) != length(samples) | length(fcols) != length(samples)) {
    warning(paste0(pid, ' from RA ', ra, ' samples do not match. Skipping'))
    break
  }
  
  rd = merged.raw %>% dplyr::select(location,chrom,start,end,!!rcols)
  fd = merged.fit %>% dplyr::select(location,chrom,start,end,!!fcols)
  
  tryCatch({
    segmented = BarrettsProgressionRisk::segmentRawData(loadSampleInformation(si),rd,fd,cutoff=0.011, verbose=T)
    residuals = bind_rows(residuals, BarrettsProgressionRisk::sampleResiduals(segmented))
    
    plots = BarrettsProgressionRisk::plotSegmentData(segmented, 'list')
    for (s in names(plots))
      ggsave(paste(plot.dir, paste(s, 'segmentedCoverage.png',sep='_'), sep='/'),  plot=plots[[s]], height=4, width=20, units='in')
    readr::write_tsv(residuals, path=paste(dirname(plot.dir), paste0(which(levels(all.val$RA) == ra),'_residuals.txt'),sep='/'), col_names = F, append=F)
    
    save(segmented, file=paste(dirname(plot.dir), paste0(sample, '_segmented.Rdata'),sep='/'))
    
    if (BarrettsProgressionRisk::sampleResiduals(segmented)$Pass) {
      prr = BarrettsProgressionRisk::predictRiskFromSegments(segmented, model=fit, s=lambda, tile.mean = z.mean, tile.sd = z.sd, arms.mean = z.arms.mean, arms.sd = z.arms.sd, cx.mean = mn.cx, cx.sd = sd.cx)
      
      predictions = bind_rows(predictions, predictions(prr))
      
      save(prr, file=paste(dirname(plot.dir), paste0(sample, '_predictions.Rdata'),sep='/'))
    }
    
  }, error = function(e) {
    message(paste("Error in segmentation for patient",pid,'from RA:', ra, ', skipping:\n\t',e))
  })
  
  write_tsv(residuals, path=paste0(dirname(plot.dir), '/residuals.tsv')) 
  write_tsv(predictions, path=paste0(dirname(plot.dir), '/predictions.tsv')) 
}

print("Finished")
