args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3)
  stop("Missing required arguments: <qdna data dir> <patient spreadsheet> <output dir> <patient name OPT>")


library(tidyverse)
#library(ggrepel)
#library(gridExtra)
#library(BarrettsProgressionRisk)
#source('~/workspace/shallowwgs_pipeline/lib/load_patient_metadata.R')

data = args[1]
val.file = args[2]
outdir = args[3]

if (length(args) <= 0) {
  data = '~/Data/BarrettsProgressionRisk/QDNAseq/validation/'
  val.file = '~/Data/BarrettsProgressionRisk/QDNAseq/validation/sWGS validation batches.xlsx'
  outdir = '~/Data/BarrettsProgressionRisk/Analysis/validation'
}

sheets = readxl::excel_sheets(val.file)

all.val = do.call(bind_rows, lapply(sheets, function(s) {
  print(s)
  readxl::read_xlsx(val.file, s) %>% select(`Hospital Research ID`, matches('Status'), `Sample Type`, `SLX-ID`, `Index Sequence`, Cohort, Batch, RA) %>% mutate_at(vars(`SLX-ID`), list(as.character))
}))

pastefun<-function(x) {
  if ( !grepl('SLX-', x) ) x = paste0('SLX-',x)
  return(x)
}
all.val = all.val %>% rowwise %>% mutate_at(vars(`SLX-ID`), list(pastefun) ) %>% ungroup
all.val = all.val %>% arrange(Batch, `Hospital Research ID`) %>% group_by(`Hospital Research ID`) %>% mutate(AID = group_indices()) %>% ungroup
all.val = all.val %>% mutate(`Hospital Research ID` = str_replace_all( str_remove_all(`Hospital Research ID`, " "), '/', '_'), `Index Sequence` = str_replace_all(`Index Sequence`, 'tp', ''))

data.dirs = grep( paste(sub('SLX-', '', unique(all.val$`SLX-ID`)), collapse='|'), list.dirs(data, full.names=T, recursive=F), value = T)
if (length(data.dirs) <= 0)
  stop(paste("No directories in", data))

merged.raw = NULL; merged.fit = NULL
for (dir in data.dirs) {
  print(dir)
  
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

fix.index <-function(x) { 
  x = str_replace_all(x, '_', '-') 
  str_replace_all(x, '.HWFF5BBXX.s-4', '')
}
merged.fit = merged.fit %>% rename_at(vars(matches('^SLX-')), fix.index)
merged.raw = merged.raw %>% rename_at(vars(matches('^SLX-')), fix.index)


print(paste0("Writing merged raw and fit data to ", data, '/merged_raw_fit.Rdata'))

save(merged.fit, merged.raw, file=paste0(data, '/merged_raw_fit.Rdata'))

print("Finished.")
