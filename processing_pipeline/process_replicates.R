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
# outdir = '~/Data/Ellie/Analysis/replicates'


patient.name = NULL
if (length(args) == 4)
  patient.name = gsub('\\/', '_', args[4])

all.patient.info = read.patient.info(patient.file, set='All')$replicates
if (is.null(all.patient.info))
  stop(paste("Failed to read patient file", patient.file))

#all.patient.info = all.patient.info %>% mutate(Samplename = sub('_','-',Samplename))

pts_slx = all.patient.info %>% dplyr::select(SLX.ID, Hospital.Research.ID,Samplename) %>% distinct %>% arrange(Hospital.Research.ID)

datafile = paste(data,'merged_raw_fit.Rdata',sep='/')
if (!file.exists(datafile)) stop("Run merge_qdnaseq_data.R first")

load(datafile, verbose=T)

multipcfdir = paste(outdir,"pcf_perPatient", basename(data), sep='/')
if ( !dir.exists(multipcfdir) ) 
  dir.create(multipcfdir, recursive = T)

if (!is.null(patient.name))
  pts_slx = pts_slx %>% filter(Hospital.Research.ID == patient.name)
  
kb = as.integer(sub('kb', '',basename(data)))

#tiled = NULL; tile.MSE = NULL
#arms.tiled = NULL; arm.MSE = NULL

residuals = tibble()
for (i in 1:nrow(pts_slx)) {
  pt = dplyr::select(dplyr::slice(pts_slx,i), Hospital.Research.ID) %>% pull
  
  pd = paste(multipcfdir, pt,sep='/')
  
#  if (!file.exists(paste0(pd,'/segment.Rdata'))) { 
  dir.create(pd,showWarnings = F)
    
  sm = dplyr::select(dplyr::slice(pts_slx,i), Samplename) %>% pull
  print(sm)
  info = all.patient.info %>% filter(Samplename == sm) %>% 
      dplyr::select('Patient','Hospital.Research.ID', 'Samplename','Endoscopy.Year','Block','Pathology','p53.Status') %>% 
      dplyr::rename('Sample' = 'Samplename', 'GEJ.Distance' = 'Block', 'p53 IHC' = 'p53.Status', 'Endoscopy' = 'Endoscopy.Year') %>% rowwise %>% 
      dplyr::mutate( Endoscopy = paste(as.character(Endoscopy), '01', '01', sep='-'))
    
    info = BarrettsProgressionRisk::loadSampleInformation(info, path=c('Normal', 'NDBE','ID','LGD','HGD','IMC'))

    cols = which(colnames(merged.fit) %in% info$Sample)
#print(colnames(merged.fit)[cols])    
    segmented = BarrettsProgressionRisk::segmentRawData(info, merged.raw[,c(1:4,cols)],merged.fit[,c(1:4,cols)], verbose=T, kb = kb, multipcf = F )
#print(colnames(segmented$seg.vals))

    residuals = bind_rows(residuals, BarrettsProgressionRisk::sampleResiduals(segmented) %>% add_column('patient'=pt, .before=1))
  
    save(segmented, file=paste0(pd,'/',sm,'-segment.Rdata'))
  
    message("Saving plots")
  
    plotdir = paste0(pd,'/segmented_plots')
    dir.create(plotdir, showWarnings=F )
    plots = BarrettsProgressionRisk::plotSegmentData(segmented, 'list')
    for (sample in names(plots)) 
      ggsave(filename=paste(plotdir, '/', sample, '_segmented.png',sep=''), plot=plots[[sample]] + labs(title=paste(pt, sample)), width=20, height=6, units='in', limitsize=F)

}
readr::write_tsv(residuals, path=paste0(pd,'/residuals.tsv'))

message("Finished")
q(save="no")
