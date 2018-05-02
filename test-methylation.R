

library(ggplot2)
#library(Matrix)
library(glmnet)

source('lib/data_func.R')
#source('lib/load_patient_metadata.R')

file = '~/Data/Ellie/Analysis/5e6_arms_all_exAHM0320/model_data.Rdata'
if (!file.exists(file))
  stop(paste("Missing data file", file))
load(file, verbose=T)
file = '~/Data/Ellie/Analysis/5e6_arms_all_exAHM0320/all.pt.alpha.Rdata'
if (!file.exists(file))
  stop(paste("Missing data file", file))
load(file, verbose=T)
rm(dysplasia.df,coefs,plots,labels)

select.alpha = '0.9'
fitV = models[[select.alpha]]
lambda.opt = performance.at.1se[[select.alpha]][, 'lambda']


file = '~/Data/Ellie/Analysis/5e6_arms_all_exAHM0320/loo.Rdata'
if (!file.exists(file))
  stop(paste("Missing data file", file))
load(file, verbose=T)
rm(plots, coefs, nzcoefs, fits)

pg.samp = do.call(rbind, pg.samp)
pg.samp$Hospital.Research.ID = NULL

status = unique(pg.samp[,c('Patient','Status')])$Status

cases = table(status)['P']
m = round((table(status)['P']/(0.0225*100))*100)
offsetMean = log(cases/m)


adjustRisk <- function(RR, offset, type='risk') {
  if (type == 'prob') {
    x = 1/(1+exp(-RR+abs(offset)))
  } else {
    x = RR+offset
  }
  return(x)
}

# Predict function giving me difficulty when I have only a single sample, this ensures the dimensions are the same
predict.progression<-function(mergedDf) {
  sparsed_test_data <- Matrix(data=0, nrow=nrow(mergedDf),  ncol=ncol(mergedDf),
                              dimnames=list(rownames(mergedDf),colnames(mergedDf)), sparse=T)
  
  for(col in colnames(mergedDf)) sparsed_test_data[,col] = mergedDf[,col]
  
  preds = predict(fitV, newx=sparsed_test_data, s=lambda.opt, type='response')
  RR = predict(fitV, newx=sparsed_test_data, s=lambda.opt, type='link')
  return(list('pred'=preds, 'RR'=RR))
}

chrlen = get.chr.lengths()


files = list.files('~/Data/Ellie/arrays/methylation', full.names = T)
dspred = data.frame(matrix(ncol=5, nrow=length(files), dimnames=list( sub('\\.txt','', basename(files)),c('segRatio','Prob','RR', 'Adj.Prob','Adj.RR'))))

for (file in files) {
  methy = read.table(file, header=T, sep='\t', stringsAsFactors = F)
  id = methy[1,1]
  methy = methy[,c('chrom','loc.start','loc.end','seg.mean')]
  methy$chrom = sub('^chr','', methy$chrom)

  r = length( which(methy[3]-methy[2] >= 5e6) )/nrow(methy)
  segs <- tile.segmented.data(methy, size=5e6, chr.info=chrlen)

  segMethy = as.matrix(segs[,-c(1:3)])
  rownames(segMethy) = paste(segs$chr, ':', segs$start, '-', segs$end, sep='')
  segMethy = t(segMethy)

  segMethy = segMethy * 1/sd(segMethy)
  #for (i in 1:ncol(segMethy)) 
  #  segMethy[,i] = unit.var(segMethy[,i], z.mean[i], z.sd[i])
  
  arms <- tile.segmented.data(methy, size='arms', chr.info=chrlen)
  armsMethy = as.matrix(arms[,-c(1:3)])
  rownames(armsMethy) = paste(arms$chr, ':', arms$start, '-', arms$end, sep='')
  armsMethy = t(armsMethy)
  
  armsMethy = armsMethy * 1/sd(armsMethy)
  
  cx.score = score.cx(segMethy,1)
  
  mergedMethy = subtract.arms(segMethy, armsMethy)
  mergedMethy = cbind(mergedMethy, 'cx' = unit.var(cx.score, mn.cx, sd.cx))
  
  pp = predict.progression(mergedMethy)
  dspred[id,] = c(r, pp$pred[,1],pp$RR[,1],adjustRisk(pp$RR, offsetMean, 'prob'),adjustRisk(pp$RR, offsetMean))
}

