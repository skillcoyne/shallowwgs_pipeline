

library(ggplot2)
#library(Matrix)
library(glmnet)

source('lib/data_func.R')

file = '~/Data/Ellie/Analysis/5e6_arms_all_exAHM0320/model_data.Rdata'
if (!file.exists(file))
  stop(paste("Missing data file", file))
load(file, verbose=T)
file = '~/Data/Ellie/Analysis/5e6_arms_all_exAHM0320/all.pt.alpha.Rdata'
if (!file.exists(file))
  stop(paste("Missing data file", file))
load(file, verbose=T)
rm(dysplasia.df,coefs,plots)

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

quant.norm<-function(vecA, dfA) {
  ranked = apply(dfA, 2, rank, ties.method='min')
  sorted = apply(dfA, 2, sort)
  means <- apply(sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  df_final <- apply(ranked, 2, index_to_mean, my_mean=means)
  
  rankVecA = rank(vecA, ties.method='min')
  
  return(sapply(1:length(rankVecA), function(i) df_final[rankVecA[i],i]  ))
}


norm.norm<-function(vecA, dfA) {
  
  
  pnorm(dfA)
  
  
}


files = list.files('~/Data/Ellie/arrays/methylation', full.names = T)
dspred = data.frame(matrix(ncol=5, nrow=length(files), dimnames=list( sub('\\.txt','', basename(files)),c('segRatio','Prob','RR', 'Adj.Prob','Adj.RR'))))

#r = length( which(methy[3]-methy[2] >= 5e6) )/nrow(methy)

r = do.call(rbind, lapply(files, function(file) {
  methy = read.table(file, header=T, sep='\t', stringsAsFactors = F)
  c(length( which(methy[,'loc.end']-methy[,'loc.start'] >= 5e6) )/nrow(methy), nrow(methy))
}))
rownames(r) = sub('.txt','',basename(files))
colnames(r) = c('Ratio.above.5e6', 'Total.Segs')
r

segMethy = do.call(rbind, lapply(files, function(file) {
  methy = read.table(file, header=T, sep='\t', stringsAsFactors = F)
  id = sub('.txt','',basename(file))
  methy = methy[,c('chrom','loc.start','loc.end','seg.mean')]
  methy$chrom = sub('^chr','', methy$chrom)
  segs <- tile.segmented.data(methy, size=5e6, chr.info=chrlen)
  
  segM = as.matrix(segs[,-c(1:3)])
  rownames(segM) = paste(segs$chr, ':', segs$start, '-', segs$end, sep='')
  segM = t(segM)
  rownames(segM) = id
  segM
}))


armsMethy = do.call(rbind, lapply(files, function(file) {
  methy = read.table(file, header=T, sep='\t', stringsAsFactors = F)
  id = sub('.txt','',basename(file))
  methy = methy[,c('chrom','loc.start','loc.end','seg.mean')]
  methy$chrom = sub('^chr','', methy$chrom)
  segs <- tile.segmented.data(methy, size='arms', chr.info=chrlen)
  
  segM = as.matrix(segs[,-c(1:3)])
  rownames(segM) = paste(segs$chr, ':', segs$start, '-', segs$end, sep='')
  segM = t(segM)
  rownames(segM) = id
  segM
}))


#r = length( which(methy[3]-methy[2] >= 5e6) )/nrow(methy)

#segMethy = t(apply(segMethy, 1, quant.norm, dfA=allDf[,1:ncol(segMethy)]))
by = 1.5
hist(segMethy)
mx = mean(segMethy)+sd(segMethy)*by
mn = mean(segMethy)-sd(segMethy)*by
abline(v=c(mx, mn), col='red')

## Get the distribution as below across ALL methylation samples then sub the mean and sd for ALL instead of for one
  
newsd = sd(segMethy[segMethy <= mx & segMethy >= mn])
newmean = mean(segMethy[segMethy <= mx & segMethy >= mn])   
  
#xx = qnorm(pnorm(segMethy, mean(segMethy), sd(segMethy)), mean(allDf[,1:ncol(segMethy)]), sd(allDf[,1:ncol(segMethy)]))
#hist(xx)

xx = qnorm(pnorm(segMethy, newmean, newsd, lower.tail=F), mean(allDf[,1:ncol(segMethy)]), sd(allDf[,1:ncol(segMethy)]), lower.tail=F)
hist(xx)

apply(xx, 2, function(r) which(r==Inf | r== -Inf))


hist(armsMethy)
mxA = mean(armsMethy)+sd(armsMethy)*by
mnA = mean(armsMethy)-sd(armsMethy)*by
abline(v=c(mxA, mnA), col='red')


newsdA = sd(armsMethy[armsMethy <= mxA & armsMethy >= mnA])
newmeanA = mean(armsMethy[armsMethy <= mxA & armsMethy >= mnA])   


xxA = qnorm(pnorm(armsMethy, newmeanA, newsdA, lower.tail=F),mean(allDf[,(ncol(segMethy)+1):(ncol(allDf)-1)]), sd(allDf[,(ncol(segMethy)+1):(ncol(allDf)-1)]), lower.tail=F)
hist(xxA)
  
#armsMethy = armsMethy * 1/sd(armsMethy)
#for (i in 1:ncol(armsMethy)) 
#  armsMethy[,i] = unit.var(armsMethy[,i], z.arms.mean[i], z.arms.sd[i])
  
#cx.score = score.cx(segMethy,1)
cx.score = score.cx(xx,1)
  
#mergedMethy = subtract.arms(segMethy, armsMethy)
mergedMethy = subtract.arms(xx, xxA)
mergedMethy = cbind(mergedMethy, 'cx' = unit.var(cx.score, mn.cx, sd.cx))
range(mergedMethy)
  
pp = predict.progression(mergedMethy)
dspred = cbind( pp$pred[,1],pp$RR[,1],adjustRisk(pp$RR, offsetMean, 'prob'),adjustRisk(pp$RR, offsetMean))
colnames(dspred) = c('Prob','RR', 'Adj.Prob','Adj.RR')

dspred = merge(dspred, r, by='row.names')
dspred