## load all of Ellie's data, merge, then break into per-patient raw/fitted files. Run segmentation and save results.

options(bitmapType = "cairo")


args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3)
  stop("Missing required arguments: <qdna data dir> <patient spreadsheet> <output dir> <patient name OPT>")

suppressPackageStartupMessages( library(BarrettsProgressionRisk) )
source('~/workspace/shallowwgs_pipeline/lib/load_patient_metadata.R')

data = args[1]
patient.file = args[2]
outdir = args[3]
# data = '~/Data/Ellie/QDNAseq/training/merged/50kb'
# patient.file = '~/Data/Ellie/QDNAseq/training/All_patient_info.xlsx'
# outdir = '~/Data/Ellie/Analysis/'


patient.name = NULL
if (length(args) == 4)
  patient.name = gsub('\\/', '_', args[4])

all.patient.info = read.patient.info(patient.file, set='All')$info
if (is.null(all.patient.info))
  stop(paste("Failed to read patient file", patient.file))

pts_slx = all.patient.info %>% dplyr::select(SLX.ID, Hospital.Research.ID) %>% distinct %>% arrange(Hospital.Research.ID)

#datafile = paste(data,'excluded_qdnaseq_output.Rdata',sep='/')
datafile = paste(data,"merged_raw_fit.Rdata", sep='/')
if (!file.exists(datafile)) stop("Run merge_qdnaseq_data.R first")

load(datafile, verbose=T)


multipcfdir = paste(outdir,"pcf_perPatient", basename(data), sep='/')
if ( !dir.exists(multipcfdir) ) 
  dir.create(multipcfdir, recursive = T)

if (!is.null(patient.name))
  pts_slx = pts_slx %>% filter(Hospital.Research.ID == patient.name)
  
kb = as.integer(sub('kb', '',basename(data)))

tiled = NULL; tile.MSE = NULL
arms.tiled = NULL; arm.MSE = NULL

for (pt in unique(pts_slx$Hospital.Research.ID)) {
  print(pt)
  pd = paste(multipcfdir, pt,sep='/')
  
  if (!file.exists(paste0(pd,'/segment.Rdata'))) {
    dir.create(pd,showWarnings = F)
    
    info = all.patient.info %>% filter(Hospital.Research.ID == pt) %>% 
      dplyr::select('Patient','Hospital.Research.ID', 'Samplename','Endoscopy.Year','Block','Pathology','p53.Status') %>% 
      dplyr::rename('Sample' = 'Samplename', 'GEJ.Distance' = 'Block', 'p53 IHC' = 'p53.Status', 'Endoscopy' = 'Endoscopy.Year') %>% rowwise %>% 
      dplyr::mutate( Endoscopy = paste(as.character(Endoscopy), '01', '01', sep='-')   )
    info = BarrettsProgressionRisk::loadSampleInformation(info, path=c('BE','ID','LGD','HGD','IMC'))
  
    cols = which(colnames(merged.fit) %in% info$Sample)
    
    segmented = BarrettsProgressionRisk::segmentRawData(info, merged.raw[,c(1:4,cols)],merged.fit[,c(1:4,cols)], verbose=T, kb = kb, multipcf = kb<100 )
    
    raw_dist = tibble(
      patient = segmented$sample.info$Hospital.Research.ID,
      sample = segmented$sample.info$Sample,

      min = (segmented$prepped.data %>% dplyr::summarise_at(vars(-chrom, -start), list(~min(.,na.rm=T))) %>% t())[,1],
      max = (segmented$prepped.data %>% dplyr::summarise_at(vars(-chrom, -start), list(~max(.,na.rm=T))) %>% t())[,1],

      mean = (segmented$prepped.data %>% dplyr::summarise_at(vars(-chrom, -start), list(~mean(.,na.rm=T))) %>% t())[,1],
      median = (segmented$prepped.data %>% dplyr::summarise_at(vars(-chrom, -start), list(~median(.,na.rm=T))) %>% t())[,1],
      stdev = (segmented$prepped.data %>% dplyr::summarise_at(vars(-chrom, -start), list(~sd(.,na.rm=T))) %>% t())[,1],
      var = (segmented$prepped.data %>% dplyr::summarise_at(vars(-chrom, -start), list(~var(.,na.rm=T))) %>% t())[,1],
      
      Q1 = (segmented$prepped.data %>% dplyr::summarise_at(vars(-chrom, -start), list(~quantile(.,probs=0.25,na.rm=T))) %>% t())[,1],
      Q3 = (segmented$prepped.data %>% dplyr::summarise_at(vars(-chrom, -start), list(~quantile(.,probs=0.75,na.rm=T))) %>% t())[,1]
    )
    raw_dist %>% write_tsv(path=paste0(pd,'/raw_dist.tsv'))
    
    residuals = BarrettsProgressionRisk::sampleResiduals(segmented) %>% add_column('patient'=pt, .before=1)
  
    save(segmented, file=paste0(pd,'/segment.Rdata'))
    readr::write_tsv(residuals, path=paste0(pd,'/residuals.tsv'))
  
    message("Saving plots")
  
    plotdir = paste0(pd,'/segmented_plots')
    dir.create(plotdir, showWarnings=F )
    plots = BarrettsProgressionRisk::plotSegmentData(segmented, 'list')
    for (sample in names(plots)) 
      ggsave(filename=paste(plotdir, '/', sample, '_segmented.png',sep=''), plot=plots[[sample]] + labs(title=paste(pt, sample)), width=20, height=6, units='in', limitsize=F)
  
    if (length(plots) <= 10)
      ggsave(filename=paste0(plotdir, '/', pt, '_segmented.png'), plot=do.call(gridExtra::grid.arrange, c(plots, ncol=1, top=pt)), width=20, height=6*length(plots), units='in', limitsize=F) 
  
  
    plotdir = paste0(pd,'/cvg_binned_fitted')
    dir.create(plotdir, showWarnings=F )
    plots = BarrettsProgressionRisk::plotCorrectedCoverage(segmented, 'list') 
    for (sample in names(plots))    
      ggsave(filename=paste0(plotdir, '/', sample, '_cvg_binned.png'), plot=plots[[sample]] + labs(title=paste(pt, sample)), width=20, height=6, units='in', limitsize = F)

    if (length(plots) <= 10)
      ggsave(filename=paste0(plotdir, '/', pt, '_cvg_binned.png'), plot=do.call(gridExtra::grid.arrange, c(plots, ncol=1, top=pt)), width=20, height=6*length(plots), units='in', limitsize=F) 
    
    
    
  } else {
    load(file = paste0(pd,'/segment.Rdata'))
  }
  
  tiles = BarrettsProgressionRisk::tileSegments(segmented, verbose=T)
  arms = BarrettsProgressionRisk::tileSegments(segmented, 'arms', verbose=T)
  
  write.table(arms$tiles, file=paste0(pd, '/arms_tiled_segvals.txt'), quote=F, sep='\t', row.names=T, col.names=NA)
  write.table(tiles$tiles, file=paste0(pd, '/5e06_tiled_segvals.txt'), quote=F, sep='\t', row.names=T, col.names=NA)

}

message("Finished")
q(save="no")
