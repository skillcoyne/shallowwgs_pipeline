
suppressPackageStartupMessages( source('~/workspace/shallowwgs_pipeline/lib/data_func.R') )

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 1)
  stop('Missing data dir')

data.dir = args[1]

ptdirs = list.dirs(data.dir, full.names=T, recursive=F)
patients = basename(ptdirs)
 
 
  mergedSegs = NULL
  mergedArms = NULL
  length(ptdirs)
  
  for (pt in ptdirs) {
    print(pt)
    if (length(list.files(pt, '*wins_tiled.txt', full.names=T)) <= 0) {
      message(paste("No tiled files for",pt))
      next
    }
    
    segvals = as.data.frame(data.table::fread(list.files(pt, '*wins_tiled.txt', full.names=T)))
    armvals = as.data.frame(data.table::fread(list.files(pt, '*arms_tiled.txt', full.names=T)))
    
    segvals = segment.matrix(segvals)
    segvals[is.na(segvals)] = mean(segvals,na.rm=T)
    
    armvals = segment.matrix(armvals)
    armvals[is.na(armvals)] = mean(armvals,na.rm=T)
    
    if (is.null(segvals) | is.null(armvals))
      stop(paste("Missing values in ", pt))
    
    if (is.null(mergedSegs)) {
      mergedSegs = segvals
      mergedArms = armvals
    } else {
      mergedSegs = rbind(mergedSegs, segvals)    
      mergedArms = rbind(mergedArms, armvals)    
    }
  }
  nrow(mergedSegs) == nrow(mergedArms)
  #setdiff(rownames(mergedSegs), rownames(mergedArms))
  dim(mergedSegs)
  
fileout = paste(data.dir, 'tmp_seg_pt.Rdata', sep='/')
  save(mergedSegs, mergedArms, file=fileout)

message(paste("Finished, saved to:",fileout))

