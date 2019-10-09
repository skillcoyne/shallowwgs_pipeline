args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3)
  stop("Missing required arguments: <qdna data dir> <patient spreadsheet> <output dir> <patient list OPT>")


suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(BarrettsProgressionRisk))

data = args[1]
val.file = args[2]
model.dir = args[3]
outdir = args[4]
patients = args[5]


patients = str_replace( unlist(str_split(patients, ',| ')), ' ', '')
patients = str_replace_all( patients, '/', '_')

#training.raw = read_tsv('~/Data/BarrettsProgressionRisk/Analysis/multipcf_perPatient/pre_seg_medians.tsv', col_types = 'cdd') %>% dplyr::rename(sd='SD',sample='Sample') %>% mutate(cohort = 'train')

var.cutoff<-function(x) { x > 0.143093 } #{ x > median(training.raw$sd) }


# data = '~/Data/BarrettsProgressionRisk/QDNAseq/validation/'
# val.file = '~/Data/BarrettsProgressionRisk/QDNAseq/validation/sWGS_validation_batches.xlsx'
# model.dir = '~/Data/BarrettsProgressionRisk/Analysis/5e6_arms/'
# outdir = '~/Data/BarrettsProgressionRisk/Analysis/validation'

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

sheets = readxl::excel_sheets(val.file)[8:13]
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

if (!is.null(patients))
  all.val = all.val %>% filter(`Hospital Research ID` %in% patients)



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




for (pid in patients) {
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
    
    ## There's an issue using multipcf with the 100kb.  It works for single sample pcf.  So...the right answer may be that each sample is selected separately and the variance tested.  That means they also get segmented separately.  Need to update the BE code for this.
    
    
    qdna = get.qdnaseg(samples, paste0(data,'/15kb'))
    prepped = BarrettsProgressionRisk:::prepRawSWGS(qdna$rd,qdna$fd)
    raw.variance = prepped$data %>% dplyr::summarise_at(vars(-chrom,-start), list(~sd(.,na.rm=T)))
    kb = 15
    variance = variance %>% mutate('15kb'=t(raw.variance[variance$Samples])[,1])

    
    if (length(which( raw.variance %>% mutate_all(var.cutoff) %>% unlist )) > 0) {
      message('Using 50kb segments')
      qdna = get.qdnaseg(samples, paste0(data,'/50kb'))
      prepped = BarrettsProgressionRisk:::prepRawSWGS(qdna$rd,qdna$fd)
      raw.variance = prepped$data %>% dplyr::summarise_at(vars(-chrom,-start), list(~sd(.,na.rm=T)))
      kb = 50
      variance = variance %>% mutate('50kb'=t(raw.variance[variance$Samples])[,1])
    }

    if (length(which( raw.variance %>% mutate_all(var.cutoff) %>% unlist )) > 0) {
      message('Using 100kb segments')
      qdna = get.qdnaseg(samples, paste0(data,'/100kb'))
      prepped = BarrettsProgressionRisk:::prepRawSWGS(qdna$rd,qdna$fd)
      raw.variance = prepped$data %>% dplyr::summarise_at(vars(-chrom,-start), list(~sd(.,na.rm=T)))
      kb = 100
      variance = variance %>% mutate('100kb'=t(raw.variance[variance$Samples])[,1])
    }
    variance %>% write_tsv(paste0(plot.dir,'/variance.tsv'))

    segmented = BarrettsProgressionRisk::segmentRawData(info,qdna$rd,qdna$fd,kb=kb, multipcf=F,verbose=T)
    residuals = BarrettsProgressionRisk::sampleResiduals(segmented)

    plots = BarrettsProgressionRisk::plotSegmentData(segmented, 'list')
    for (s in names(plots))
      ggsave(paste(plot.dir, paste(s, 'segmentedCoverage.png',sep='_'), sep='/'),  plot=plots[[s]], height=4, width=20, units='in')
    
    ggsave(paste(plot.dir, paste(pid, ra, 'segmentedCoverage.png',sep='_'), sep='/'),  plot=do.call(grid.arrange, c(plots,ncol=1)), height=4*length(rcols), width=20, units='in')
    
    file.remove( paste(dirname(plot.dir), 'residuals.txt',sep='/' ) )
    readr::write_tsv(residuals, path=paste(dirname(plot.dir), paste0(which(levels(all.val$RA) == ra),'_residuals.txt'),sep='/'), col_names = F, append=F)
    save(segmented, file=paste(dirname(plot.dir), paste0(which(levels(all.val$RA) == ra), '_segObj.Rdata'),sep='/'))
    
    failed = sampleResiduals(segmented) %>% dplyr::filter(!Pass)
    # if (nrow(failed) < nrow(sampleResiduals(segmented))) {
    #   tiles = BarrettsProgressionRisk::tileSegments(segmented)
    #   arms = BarrettsProgressionRisk::tileSegments(segmented)
    #   
    #   write.table(tiles$tiles, sep='\t', quote=F, col.names=NA, row.names=T, file=paste0(dirname(plot.dir), '/5e06_cleaned_tiled.tsv'))
    #   write.table(arms$tiles, sep='\t', quote=F, col.names=NA, row.names=T, file=paste0(dirname(plot.dir), '/arms_cleaned_tiled.tsv'))
    # } 
    
    pred.dir = paste0(dirname(plot.dir), '/predictions')
    dir.create(pred.dir,showWarnings = F)
    
    prr = predictRiskFromSegments(segmented, be.model = be.model)
    plots = BarrettsProgressionRisk::copyNumberMountainPlot(prr,annotate = T,legend = F, as='list')
    for (s in names(plots))
      ggsave(paste(pred.dir, paste0(s,'_cnMtn.png'),sep='/'), plot=plots[[s]], height=4, width=20, units='in')

    predictions(prr) %>% write_tsv(path=paste0(pred.dir,'/preds.tsv'))
    
  }, error = function(e) {
    message(paste("Error in segmentation/predictions for patient",pid,': ',e))
  })
}


message('Finished')







