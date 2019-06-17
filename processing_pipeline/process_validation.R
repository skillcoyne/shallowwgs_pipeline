
args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3)
  stop("Missing required arguments: <qdna data dir> <patient spreadsheet> <output dir> <patient list OPT>")


suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(BarrettsProgressionRisk))
#source('~/workspace/shallowwgs_pipeline/lib/load_patient_metadata.R')

data = args[1]
val.file = args[2]
outdir = args[3]

# data = '~/Data/BarrettsProgressionRisk/QDNAseq/validation/'
# val.file = '~/Data/BarrettsProgressionRisk/QDNAseq/validation/sWGS_validation_batches.xlsx'
# outdir = '~/Data/BarrettsProgressionRisk/Analysis/validation'

patients = NULL
if (length(args) == 4) patients = args[4]

#patients = c('OCCAMS_AH_100','AHM0185','AHM0555','AHM0348','AHM0281','AHM0500','AHM0546','AHM0584','AHM1432','AHM1471','AHM0567','AHM0570')

if (!is.null(patients)) {
  patients = str_replace( unlist(str_split(patients, ',| ')), ' ', '')
  patients = str_replace_all( patients, '/', '_')
}


sheets = readxl::excel_sheets(val.file)[1:10]


all.val = do.call(bind_rows, lapply(sheets, function(s) {
  readxl::read_xlsx(val.file, s) %>% select(`Hospital Research ID`, matches('Status'), `Sample Type`, `SLX-ID`, `Index Sequence`, Cohort, Batch, RA) %>% mutate_at(vars(`SLX-ID`), list(as.character))
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

if (!is.null(patients))
  all.val = all.val %>% filter(`Hospital Research ID` %in% patients)

#patients = (pts_slx %>% filter(`SLX-ID` == slx)) %>% select(`Patient ID`) %>% pull
failedQC = tibble()

# Need to process the batches from CK and SA separately even if the patients overlap


for (ra in levels(all.val$RA) ) {
  print(ra)
  pts = all.val %>% filter(RA == ra) %>% select(`Hospital Research ID`) %>% unique() %>% pull  
  if (length(pts) <= 0) next
  
  for (pid in pts) {
    message(paste('Patient',pid))
    
    si = all.val %>% filter(`Hospital Research ID` == pid & RA == ra) 
    si$Sample = si$Samplename 
    if (is.null(si[['Endoscopy']])) si = si %>% mutate(Endoscopy = '2019/01/01')

    plot.dir = paste(outdir, pid, 'plots',sep='/')
    #if (length(list.files(plot.dir)) >= nrow( all.val %>% filter(`Hospital Research ID` == pid) )) next  
    if (file.exists(paste(dirname(plot.dir), paste0(which(levels(all.val$RA) == ra), '_segObj.Rdata'),sep='/'))) next # skip patients I've already done

    dir.create(plot.dir, showWarnings = F, recursive = T)
    #si = si %>% mutate(Sample = paste(`SLX-ID`, gsub('-','_',si$`Index Sequence`), sep='.'))

    samples = si %>% filter(RA == ra & `Hospital Research ID` == pid) %>% select(Samplename) %>% pull

    rcols = grep(paste(samples,collapse='|'), colnames(merged.raw))
    fcols = grep(paste(samples,collapse='|'), colnames(merged.fit))
    
    if (length(rcols) != length(samples) | length(fcols) != length(samples)) {
      warning(paste0(pid, ' from RA ', ra, ' samples do not match. Skipping'))
      break
    }
    
    rd = merged.raw %>% dplyr::select(location,chrom,start,end,!!rcols)
    fd = merged.fit %>% dplyr::select(location,chrom,start,end,!!fcols)

    tryCatch({
      
      segmented = BarrettsProgressionRisk::segmentRawData(loadSampleInformation(si),rd,fd,cutoff=0.011, verbose=F)
      residuals = BarrettsProgressionRisk::sampleResiduals(segmented)

      failedQC = bind_rows(failedQC, residuals)
      
      #prr = BarrettsProgressionRisk::predictRiskFromSegments(segmented)
      #preds = bind_rows(preds, predictions(prr))
      
      #plots = BarrettsProgressionRisk::plotCorrectedCoverage(segmented, 'list')
      plots = BarrettsProgressionRisk::plotSegmentData(segmented, 'list')
      for (s in names(plots))
        ggsave(paste(plot.dir, paste(s, 'segmentedCoverage.png',sep='_'), sep='/'),  plot=plots[[s]], height=4, width=20, units='in')

      ggsave(paste(plot.dir, paste(pid, ra, 'segmentedCoverage.png',sep='_'), sep='/'),  plot=do.call(grid.arrange, c(plots,ncol=1)), height=4*length(rcols), width=20, units='in')
      
      file.remove( paste(dirname(plot.dir), 'residuals.txt',sep='/' ) )
      readr::write_tsv(residuals, path=paste(dirname(plot.dir), paste0(which(levels(all.val$RA) == ra),'_residuals.txt'),sep='/'), col_names = F, append=F)
      save(segmented, file=paste(dirname(plot.dir), paste0(which(levels(all.val$RA) == ra), '_segObj.Rdata'),sep='/'))
      
      #if (nrow(subset(residuals, !Pass)) > 0)
      #  failedQC = bind_rows(failedQC, residuals %>% filter(!Pass) %>% mutate(`Hospital Research ID` = pid) %>% select(`Patient ID`, everything()) )
      
      # if (nrow(subset(residuals,Pass)) > 0) {
      #   pr = BarrettsProgressionRisk::predictRiskFromSegments(segmented, verbose=F)
      #   preds = predictions(pr)
      #   } else {
      #   message(paste('All samples failed QC for patient ID',pid))
      #   }
  }, error = function(e) {
    message(paste("Error in segmentation for patient",pid,'from RA:', ra, ', skipping:\n\t',e))
  })
}
}

print("Finished")
