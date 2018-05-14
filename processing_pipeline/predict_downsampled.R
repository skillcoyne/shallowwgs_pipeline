# predict downsampled

library(ggplot2)
library(Matrix)
library(glmnet)

source('lib/data_func.R')
source('lib/load_patient_metadata.R')

chrlen = get.chr.lengths()

adjustRisk <- function(RR, offset, type='risk') {
  if (type == 'prob') {
    x = 1/(1+exp(-RR+abs(offset)))
  } else {
    x = RR+offset
  }
  return(x)
}


dir = '~/Data/Ellie/Cleaned/Downsampled'
#dir = '~/Data/Ellie/Cleaned/AH0329_segmentedCoverage_fitted_gamma250/'
#dir = '~/Data/Ellie/Analysis/VAL_Cohort/'
files = list.files(dir, full.names=T)

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

sampleNames = sub('_raw.txt', '', basename(files))
dspred = data.frame(matrix(ncol=4, nrow=length(sampleNames), dimnames=list(sampleNames,c('Prob','RR', 'Adj.Prob','Adj.RR'))))

for (n in 1:length(files)) {
  print(n)
  print(sampleNames[n])
  
  rs = read.table(files[n], sep='\t',header=T)
  if (ncol(rs) < 4)
    stop("File needs to contain 4 columns: chr, start, end, sample_value")
  
  segs <- tile.segmented.data(rs, size=5e6, chr.info=chrlen)
  segM = as.matrix(segs[,-c(1:3)])
  rownames(segM) = paste(segs$chr, ':', segs$start, '-', segs$end, sep='')
  segM = t(segM)
  for (i in 1:ncol(segM)) 
    segM[,i] = unit.var(segM[,i], z.mean[i], z.sd[i])
  
  arms <- tile.segmented.data(rs, size='arms', chr.info=chrlen)
  armsM = as.matrix(arms[,-c(1:3)])
  rownames(armsM) = paste(arms$chr, ':', arms$start, '-', arms$end, sep='')
  armsM = t(armsM)
  for (i in 1:ncol(armsM)) 
    armsM[,i] = unit.var(armsM[,i], z.arms.mean[i], z.arms.sd[i])
  
  nrow(armsM) == nrow(segM)
    
  cx.score = score.cx(segM,1)
  
  mergedDf = subtract.arms(segM, armsM)
  mergedDf = cbind(mergedDf, 'cx' = unit.var(cx.score, mn.cx, sd.cx))
  dim(mergedDf)
  
  # Predict function giving me difficulty when I have only a single sample, this ensures the dimensions are the same
  sparsed_test_data <- Matrix(data=0, nrow=nrow(mergedDf),  ncol=ncol(mergedDf),
                              dimnames=list(rownames(mergedDf),colnames(mergedDf)), sparse=T)
  for(col in colnames(mergedDf)) sparsed_test_data[,col] = mergedDf[,col]
  
  preds = predict(fitV, newx=sparsed_test_data, s=lambda.opt, type='response')
  RR = predict(fitV, newx=sparsed_test_data, s=lambda.opt, type='link')

  dspred[n,] = c(preds[,1],RR[,1],adjustRisk(RR, offsetMean, 'prob'),adjustRisk(RR, offsetMean))
  print(dspred[1:n,])
}

save(dspred, file='dowsampled_predictions.Rdata')
save(dspred, file=paste(dir, 'predictions.Rdata', sep='/'))

load(paste(dir, 'predictions.Rdata', sep='/'), verbose=T)

riskPal = rev(RColorBrewer::brewer.pal(11, 'RdYlBu'))
unadjRisk = ggplot(pg.samp, aes(OR)) + geom_histogram(aes(fill=..x..), bins=10, show.legend = T) +
  scale_fill_gradientn(colors = riskPal,  name='RR') + 
  labs(y='n Samples', x='Relative Risk', title='Unadjusted relative risk') + theme_light(base_size = 14)
unadjRisk

m = melt(dspred, measure.vars=c('RR','Adj.RR'))
unadjRisk + geom_point(data=m, aes(x=value,y=Adj.Prob*30, color=variable))
  
ids = readxl::read_excel('~/Data/Ellie/Analysis/downsampled_ids.xlsx', sheet=1)


preds = base::merge(dspred, ids, by.x='row.names', by.y='Illumina ID', all=T) 
# ignore the three LGD cases, they're not clear NP or P in any case (downsampled just out of interest)
preds = subset(preds, Grade_biopsy != 'LGD')


preds[ which(subset(preds, `Expected Risk` == 'Low')$Adj.Prob > 0.3) ,]

preds[ which(subset(preds, `Expected Risk` == 'High')$Adj.Prob < 0.6) ,]






