
args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3)
  stop("Missing required params: <file in> <model dir> <output dir>")

library(copynumber)
library(tibble)
library(dplyr)
#library(ggplot2)
#library(gridExtra)
library(glmnet)

source('lib/data_func.R')

fileIn = args[1]
modelDir = args[2]
outdir = args[3]

dir.create(outdir, showWarnings=F, recursive=T)

load(paste(modelDir, '/model_data.Rdata', sep=''), verbose=T)
rm(dysplasia.df, labels)
load(paste(modelDir,'/all.pt.alpha.Rdata', sep=''), verbose = T)
fitV = models$`0.9`
l = performance.at.1se$`0.9`$lambda
rm(plots,performance.at.1se,dysplasia.df,models,cvs,labels)


#filelist = list.files('~/Data/Ellie/arrays/SNP', pattern='txt', full.names=T)
samplelist = colnames(as_tibble(data.table::fread(fileIn, nrows=2)))

whole_epi = grep('*whole.*Log', samplelist, value=T)  
normals = grep('*normal.*Log', samplelist, value=T)  


fix.chr<- function(df) {
  df = subset(df, Chr %in% c(1:22))
  df$Chr = as.numeric(df$Chr)
  dplyr::arrange(df, Chr, Position)
}

process.patients<-function(df, outfile=NULL) {
  predictions = data.frame(matrix(nrow=0,ncol=2,dimnames=list(c(), c('Prob','RR'))))
  allptsdata = list()
  
  infocols = grep('Chr|Pos', colnames(df), value=T)
  pts = grep('SNP|Chr|Pos', colnames(df), invert=T, value=T)
  
  chr.len = get.chr.lengths()
  
  for (pt in pts)  {
    print(pt)    
    data = df[,c(infocols,pt)]
    segdata = copynumber::pcf(as.data.frame(data), gamma=28, fast=T, verbose=F)  
    #head(segdata)
    
    #ggplot(segdata, aes(x=1:nrow(segdata), y=mean)) + geom_point() + labs(title='LogR whole epi', x='', y='')
    
    # this is the LogR, doesn't need to be logged again
    filename = paste(dirname(outfile), '/', unique(segdata$sampleID), '.segmented_tiled.txt', sep='')
    if (file.exists(filename)) {
      tiled = load.segment.matrix(filename)
    } else {
      tiles = tile.segmented.data(segdata[,-1], chr.info=chr.len)
      write.table(tiles, sep='\t', quote=F, row.names=F, file=filename)
      tiled = segment.matrix(tiles)
      rownames(tiled) = unique(segdata$sampleID)
    }
    tiled.noUV = tiled
    
    for (i in 1:ncol(tiled)) 
      tiled[,i] = unit.var(tiled[,i], z.mean[i], z.sd[i]) 
    tiled[is.na(tiled)] = mean(tiled, na.rm=T) #0
    
    cx = score.cx(tiled, 1)
    
    filename = paste(dirname(outfile), '/', unique(segdata$sampleID), '.segmented_tiled_arms.txt', sep='')
    if (file.exists(filename)) {
      tiled.arms = load.segment.matrix(filename)
    } else {
      tiles = tile.segmented.data(segdata[,-1], size='arms', chr.info=chr.len)
      write.table(tiles, sep='\t', quote=F, row.names=F, file=filename)
      tiled.arms = segment.matrix(tiles)
      rownames(tiled.arms) = unique(segdata$sampleID)
    }
    tiled.armsNoUV = tiled.arms
    
    for (i in 1:ncol(tiled.arms)) 
      tiled.arms[,i] = unit.var(tiled.arms[,i], z.arms.mean[i], z.arms.sd[i]) 
    tiled.arms[is.na(tiled.arms)] = mean(tiled, na.rm=T) #0
    
    arrayDf = subtract.arms(tiled, tiled.arms)
    arrayDf = cbind(arrayDf, 'cx'=unit.var(cx, mn.cx ,sd.cx))
    
    # Whole epi
    predictions[pt,] = c(predict(fitV, newx=arrayDf, s=l, type='response')[,1], 'RR'=predict(fitV, newx=arrayDf, s=l, type='link')[,1] )
    
    allptsdata[[pt]] = list('seg'=segdata, 'tile'=tiled, 'tileNoUV'=tiled.noUV, 'tile.arms'=tiled.arms, 'tile.arms.noUV'=tiled.armsNoUV, 'arrayDf'=arrayDf, 'cx'=cx)
  }  
  if(!is.null(outfile))   save(allptsdata, predictions, file=outfile)
  
  return(predictions)
}

dt = as_tibble(data.table::fread(fileIn, showProgress=T))

epi = dt[,c(grep('SNP|Chr|Pos', colnames(dt), value=T),whole_epi)]
nm = dt[,c(grep('SNP|Chr|Pos', colnames(dt), value=T),normals)]

epi = fix.chr(epi)
nm = fix.chr(nm)
rm(dt)

tmpE = paste(outdir, sub('\\.txt','',basename(fileIn)), '_epi.Rdata', sep='')
tmpN = paste(outdir, sub('\\.txt','',basename(fileIn)), '_norm.Rdata', sep='')


all.preds = rbind(process.patients(epi, outfile=tmpE),process.patients(nm, outfile=tmpN))

print(all.preds)

write.table(all.preds, sep='\t', quote=F, row.names=F, file=paste(sub('\\.txt', '', fileIn), '-predictions.txt', sep=''))
