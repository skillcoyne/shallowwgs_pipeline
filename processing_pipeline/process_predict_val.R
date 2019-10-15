args = commandArgs(trailingOnly=TRUE)

if (length(args) < 5)
  stop("Missing required arguments: <qdna data dir> <patient spreadsheet> <model dir> <output dir> <patient>")


suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(BarrettsProgressionRisk))

data = args[1]
val.file = args[2]
model.dir = args[3]
outdir = args[4]
pid = args[5]

# data = '~/Data/BarrettsProgressionRisk/QDNAseq/validation/'
# val.file = '~/Data/BarrettsProgressionRisk/QDNAseq/validation/sWGS_validation_batches.xlsx'
# model.dir = '~/Data/BarrettsProgressionRisk/Analysis/5e6_arms/'
# outdir = '~/Data/BarrettsProgressionRisk/Analysis/validation'
# pid = 'AHM1807'

pid = str_replace_all( pid, '/', '_')

#training.raw = read_tsv('~/Data/BarrettsProgressionRisk/Analysis/multipcf_perPatient/pre_seg_medians.tsv', col_types = 'cdd') %>% dplyr::rename(sd='SD',sample='Sample') %>% mutate(cohort = 'train')
var.cutoff<-function(x) { x > 0.143093 } #{ x > median(training.raw$sd) }


select.alpha = '0.9'

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

sheets = readxl::excel_sheets(val.file)[8:14]
all.val = do.call(bind_rows, lapply(sheets, function(s) {
  readxl::read_xlsx(val.file, s) %>% dplyr::select(`Hospital Research ID`, matches('Status'), `Sample Type`, `SLX-ID`, `Index Sequence`, Cohort, Batch, RA) %>% 
    dplyr::mutate_at(vars(`SLX-ID`), list(as.character)) %>% dplyr::filter(!is.na(`SLX-ID`))
}))

pastefun<-function(x) {
  if ( !grepl('SLX-', x) ) x = paste0('SLX-',x)
  return(x)
}
all.val = all.val %>% rowwise %>% mutate_at(vars(`SLX-ID`), list(pastefun) ) %>% ungroup
#all.val = all.val %>% arrange(Batch, `Hospital Research ID`) %>% group_by(`Hospital Research ID`) %>% mutate(AID = group_indices()) %>% ungroup
all.val = all.val %>% mutate(`Hospital Research ID` = str_replace_all( str_remove_all(`Hospital Research ID`, " "), '/', '_'), `Index Sequence` = str_replace_all(`Index Sequence`, 'tp', ''))
all.val = all.val %>% mutate(Samplename = paste(`SLX-ID`, `Index Sequence`, sep='.'), RA = factor(RA))

if ( length(grep('15|50|100', list.files(data, 'kb'))) < 3) 
  stop('Segmentation files for 15kb, 50kb, and 100kb are required.')

# First time use the 15kb file
data.15kb = paste0(data, '/15kb')

if (!file.exists(  paste0(data.15kb, '/merged_raw_fit.Rdata')))
  stop("Missing merged data files, run merge_qdnaseq_data.R' first.")

load(file=paste0(data.15kb, '/merged_raw_fit.Rdata'))

if (!is.null(pid))
  all.val = all.val %>% filter(`Hospital Research ID` %in% pid)

get.qdnaseg<-function(samples, dir) {
  if (!file.exists(paste0(dir, '/merged_raw_fit.Rdata'))) stop(paste0('Missing merged rdata file in ', dir))
  load(file=paste0(dir, '/merged_raw_fit.Rdata'))

  rcols = grep(paste(samples,collapse='|'), colnames(merged.raw))
  fcols = grep(paste(samples,collapse='|'), colnames(merged.fit))

  if (length(rcols) != length(samples) | length(fcols) != length(samples)) 
    stop(paste0(pid, ' samples do not match. Skipping'))
    
  rd = merged.raw %>% dplyr::select(location,chrom,start,end,!!rcols)
  fd = merged.fit %>% dplyr::select(location,chrom,start,end,!!fcols)

  return(list('rd'=rd,'fd'=fd))
}

  message(paste('Patient',pid))
  
  si = all.val %>% filter(`Hospital Research ID` == pid) 
  si$Sample = si$Samplename 
  if (is.null(si[['Endoscopy']])) si = si %>% mutate(Endoscopy = '2019/01/01')
  
  plot.dir = paste(outdir, pid, 'plots',sep='/')
  #if (length(list.files(plot.dir)) >= nrow( all.val %>% filter(`Hospital Research ID` == pid) )) next  
  #if (file.exists(paste(dirname(plot.dir), paste0(which(levels(all.val$RA) == ra), '_segObj.Rdata'),sep='/'))) next # skip patients I've already done
  
  dir.create(plot.dir, showWarnings = F, recursive = T)
  
  #for (sample in si$Samplename) {
  samples = si %>% filter(`Hospital Research ID` == pid) %>% dplyr::select(Samplename) %>% pull

  variance = tibble( 'Samples' = samples, '15kb'=NA, '50kb'=NA, '100kb'=NA)
  tryCatch({
    info = loadSampleInformation(si %>% filter(Samplename == samples))

    residuals = tibble(); predictions = tibble()
    
    objList = list()  
    # This is ridiculously kludgy but...  
    for (sample in info$Sample) {
      message('15kb...')
      qdna = get.qdnaseg(sample, paste0(data,'/15kb'))
      prepped = BarrettsProgressionRisk:::prepRawSWGS(qdna$rd,qdna$fd)
      raw.variance = prepped$data %>% dplyr::summarise_at(vars(-chrom,-start), list(~sd(.,na.rm=T)))
      kb = 15
      variance[variance$Samples == sample,] = variance %>% filter(Samples == sample) %>% mutate('15kb'=t(raw.variance)[,1])

      if (var.cutoff(raw.variance[[1]])) {
        message('50kb...')
        qdna = get.qdnaseg(sample, paste0(data,'/50kb'))
        prepped = BarrettsProgressionRisk:::prepRawSWGS(qdna$rd,qdna$fd)
        raw.variance = prepped$data %>% rename_at(vars(matches(paste(info$Samplename,collapse='|'))), funs(sub('.H\\dG.*','',.))) %>%
          dplyr::summarise_at(vars(-chrom,-start), list(~sd(.,na.rm=T)))
        variance[variance$Samples == sample,] = variance %>% filter(Samples == sample) %>% mutate('50kb'=t(raw.variance)[,1])
        kb = 50
      }
        
      if (var.cutoff(raw.variance[[1]])) {
        message('100kb...')
        qdna = get.qdnaseg(sample, paste0(data,'/100kb'))
        prepped = BarrettsProgressionRisk:::prepRawSWGS(qdna$rd,qdna$fd)
        raw.variance = prepped$data %>% rename_at(vars(matches(paste(info$Samplename,collapse='|'))), funs(sub('.H\\dG.*','',.))) %>%
          dplyr::summarise_at(vars(-chrom,-start), list(~sd(.,na.rm=T)))
        variance[variance$Samples == sample,] = variance %>% filter(Samples == sample) %>% mutate('100kb'=t(raw.variance)[,1])
        kb = 100
      } 
        
      qdna$rd = qdna$rd %>% select(matches('loc|chr|start|end'), matches(paste(info$Samplename,collapse='|'))) %>% 
        rename_at(vars(matches(paste(info$Samplename,collapse='|'))), funs(sub('.H\\dG.*','',.)))
      qdna$fd = qdna$fd %>% select(matches('loc|chr|start|end'), matches(paste(info$Samplename,collapse='|'))) %>% 
        rename_at(vars(matches(paste(info$Samplename,collapse='|'))), funs(sub('.H\\dG.*','',.)))

      segmented = BarrettsProgressionRisk::segmentRawData(info %>% filter(Sample == sample),qdna$rd, qdna$fd, kb=kb, cutoff = 0.011, multipcf=F,verbose=T)   
      residuals = bind_rows(residuals, BarrettsProgressionRisk::sampleResiduals(segmented) %>% mutate('kb' = kb) )
      objList[[sample]] = segmented
      
      plots = BarrettsProgressionRisk::plotSegmentData(segmented, 'list')
      for (s in names(plots))
        ggsave(paste(plot.dir, paste(s, 'segmentedCoverage.png',sep='_'), sep='/'),  plot=plots[[s]], height=4, width=20, units='in')
      
      pred.dir = paste0(dirname(plot.dir), '/predictions')
      dir.create(pred.dir,showWarnings = F)

      prr = predictRiskFromSegments(segmented, be.model = be.model)
      predictions = bind_rows(predictions, predictions(prr))
      
      plots = BarrettsProgressionRisk::copyNumberMountainPlot(prr,annotate = T,legend = F, as='list')
      for (s in names(plots))
        ggsave(paste(pred.dir, paste0(s,'_cnMtn.png'),sep='/'), plot=plots[[s]], height=4, width=20, units='in')
    }
    
    variance %>% write_tsv(paste0(plot.dir,'/variance.tsv'))
    predictions %>% write_tsv(path=paste0(pred.dir,'/preds.tsv'))
    residuals %>% write_tsv(paste0(dirname(plot.dir), '/residuals.tsv'))
  }, error = function(e) {
    message(paste("Error in segmentation/predictions for patient",pid,': ',e))
  })


message('Finished')







