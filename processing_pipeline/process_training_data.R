## load all of Ellie's data, merge, then break into per-patient raw/fitted files. Run segmentation and save results.


args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3)
  stop("Missing required arguments: <qdna data dir> <patient spreadsheet> <output dir> <patient name OPT>")

library(BarrettsProgressionRisk)
source('~/workspace/shallowwgs_pipeline/lib/load_patient_metadata.R')

data = args[1]
patient.file = args[2]
outdir = args[3]
#data = '~/Data/Ellie/QDNAseq/training/'
#patient.file = paste(data, 'All_patient_info.xlsx', sep='/')
#outdir = '~/Data/Ellie/Analysis'

patient.name = NULL
if (length(args) == 4)
  patient.name = args[4]

all.patient.info = read.patient.info(patient.file, set='All')
if (is.null(all.patient.info))
  stop(paste("Failed to read patient file", patient.file))

data.dirs = list.dirs(data, full.names=T, recursive=F)
if (length(data.dirs) <= 0)
  stop(paste("No directories in", data))

pts_slx = arrange(unique(all.patient.info$info[c('SLX.ID','Hospital.Research.ID')]), Hospital.Research.ID)

data.dirs = grep(paste(unique(pts_slx$SLX.ID), collapse='|'), data.dirs, value=T)

datafile = paste(data,"merged_qdnaseq_output.Rdata", sep='/')
if (!file.exists(datafile)) {
  # Load the data, binding each new data file into single data table
  fit.data = NULL
  raw.data = NULL
  samplenames=NULL
  for(dir in data.dirs) {
    print(paste("Reading directory", dir))
    files = list.files(dir,pattern="[fitted|corrected]ReadCounts.txt", full.names=T)
    slx = sub('SLX-','', basename(dir))
    
    # Files that do not contain individual indecies
    if (length(files) == 1 || !grepl('D\\d+_D\\d+', files)) {
      for (file in files) {
        print(paste("Reading file", file))
        spl = strsplit(basename(file),"\\.")
        
        fitted.file = file
        raw.file = list.files(dir, pattern='\\.(raw)?readCounts', full.names=T)
        if (length(raw.file) > 1)
          stop(paste("Error finding raw read counts in", dir))
        
        message( paste('FIT:', basename(fitted.file), '  RAW:', basename(raw.file)) )
        
        #fit = read.table(fitted.file,sep="\t",header=T,stringsAsFactors=F)  # fitted reads
        fit = data.table::fread(fitted.file)
        f.cols = grep('D\\d', colnames(fit))
        colnames(fit)[f.cols] = paste(colnames(fit)[f.cols], slx, sep="_") 
        
        #raw = read.table(raw.file,sep="\t",header=T,stringsAsFactors=F) # raw reads
        raw = data.table::fread(raw.file)
        r.cols = grep('D\\d', colnames(raw))
        colnames(raw)[r.cols] = paste(colnames(raw)[r.cols], slx, sep="_") 
        raw = raw[,colnames(fit),with=F]	    
        
        if (is.null(fit.data)) {
          fit.data = fit
          raw.data = raw
        } else {
          fit.data = merge(fit.data, fit, by=c('location','chrom','start','end'), all=T) 
          raw.data = merge(raw.data, raw, by=c('location','chrom','start','end'), all=T) 
          
          #fit.data = cbind(fit.data,fit[,f.cols])
          #raw.data = cbind(raw.data,raw[,f.cols])
        }
      }
    } else { # individual index
      for(file in files){
        print(paste("Reading file", file))
        spl = strsplit(basename(file),"\\.")
        
        fitted.file = file
        raw.file = sub('fitted', 'raw', fitted.file)
        
        message( paste('FIT:', basename(fitted.file), '  RAW:', basename(raw.file)) )
        
        if (!file.exists(fitted.file) || !file.exists(raw.file)) 
          stop("Missing fitted/raw file ")
  
        fit = data.table::fread(fitted.file)      
        #fit = read.table(fitted.file,sep="\t",header=T,stringsAsFactors=F)  # fitted reads
        colnames(fit)[5] = paste(spl[[1]][2],slx,sep="_")
        
        raw = data.table::fread(raw.file)
        #raw = read.table(raw.file,sep="\t",header=T,stringsAsFactors=F) # raw reads
        colnames(raw)[5] = paste(spl[[1]][2],slx,sep="_")
        
        if (is.null(fit.data)) {
          fit.data = fit
          raw.data = raw
        } else { # add to table
          fit.data = merge(fit.data, fit, by=c('location','chrom','start','end'), all=T) 
          raw.data = merge(raw.data, raw, by=c('location','chrom','start','end'), all=T) 
          #fit.data = cbind(fit.data,fit[,5])
          #raw.data = cbind(raw.data, raw[,5])
          
          #colnames(fit.data)[ncol(fit.data)] = colnames(fit)[5]
          #colnames(raw.data)[ncol(raw.data)] = colnames(raw)[5]
        }	  
      }
    }
    print(dim(fit.data))
  }
 
  fit.data = as.data.frame(fit.data)
  raw.data = as.data.frame(raw.data)
   
  save(fit.data, raw.data, file=datafile)
} else {
  load(datafile, verbose=T)
}


plot.dir = paste(outdir, "coverage_binned_fitted", sep='/')
if ( !dir.exists(plot.dir) ) 
  dir.create(plot.dir)

multipcfdir = paste(outdir, "multipcf_perPatient", sep='/')
if ( !dir.exists(multipcfdir) ) 
  dir.create(multipcfdir)


if (!is.null(patient.name))
  pts_slx = subset(pts_slx, Hospital.Research.ID == patient.name)

tiled = NULL; tile.MSE = NULL
arms.tiled = NULL; arm.MSE = NULL

for (pt in unique(pts_slx$Hospital.Research.ID)) {
  print(pt)
  pd = paste(multipcfdir, pt,sep='/')
  dir.create(pd,showWarnings = F)
  
  cols = which(colnames(fit.data) %in% subset(all.patient.info$info, Hospital.Research.ID == pt)$Samplename)

  #write.table(fit.data[,c(1:4,cols)], sep='\t', quote=F, row.names=F, col.names=T, file=paste(pd,'fittedReadCounts.txt',sep='/')) 
  #write.table(raw.data[,c(1:4,cols)], sep='\t', quote=F, row.names=F, col.names=T, file=paste(pd,'rawReadCounts.txt',sep='/')) 
  
  segmented = BarrettsProgressionRisk::segmentRawData(raw.data[,c(1:4,cols)],fit.data[,c(1:4,cols)])
  tile = BarrettsProgressionRisk::tileSegments(segmented, size=5e6)
  arms = BarrettsProgressionRisk::tileSegments(segmented, size='arms')

  if (is.null(tiled)) {
    tiled = tile$tiles
    arms.tiled = arms$tiles
    tile.MSE = tile$error
    arm.MSE = arms$error
  } else {
    tiled = rbind(tiled, tile$tiles)
    arms.tiled = rbind(arms.tiled, arms$tiles)
    tile.MSE = rbind(tile.MSE, tile$error)
    arm.MSE = rbind(arm.MSE, arms$error)
  }
    
  #filename = paste(pd, '/', pt,"_probefiltered_segvals_gamma250.txt",sep="")
  #write.table(segmented$seg.vals, file=filename, sep="\t", quote=F)
  save(segmented, tile, arms, file=paste(pd, 'segment.Rdata', sep = '/'))

  message("Saving plots")
  ggsave(filename=paste(plot.dir, '/', pt,'_cvg_binned_fitted.png',sep=''), plot=BarrettsProgressionRisk::plotCorrectedCoverage(segmented) + labs(title=pt), width=20, height=4*length(cols), units='in', limitsize = F)
  
  ggsave(filename=paste(pd, '/segmented.png',sep=''), plot=BarrettsProgressionRisk::plotSegmentData(segmented), width=20, height=4*length(cols), units='in', limitsize=F)
}

write.table(tiled, sep='\t', quote=F, col.names=NA, row.names=T, file=paste(multipcfdir, '5e6_tiled.txt', sep='/'))
write.table(tile.MSE, sep='\t', quote=F, col.names=NA, row.names=T, file=paste(multipcfdir, '5e6_tiled_MSE.txt', sep='/'))

write.table(arms.tiled, sep='\t', quote=F, col.names=NA, row.names=T, file=paste(multipcfdir, 'arms_tiled.txt', sep='/'))
write.table(arm.MSE, sep='\t', quote=F, col.names=NA, row.names=T, file=paste(multipcfdir, 'arms_tiled_MSE.txt', sep='/'))


message("Finished")
q(save="no")
