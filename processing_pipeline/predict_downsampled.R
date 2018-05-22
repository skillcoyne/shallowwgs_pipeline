# predict downsampled
args = commandArgs(trailingOnly=TRUE)

if (length(args) < 1)
  stop("Missing required params: <data dir> <model dir> <out dir> ")


library(ggplot2)
library(Matrix)
library(glmnet)

source('../lib/data_func.R')
source('../lib/load_patient_metadata.R')

chrlen = get.chr.lengths()

adjustRisk <- function(RR, offset, type='risk') {
  if (type == 'prob') {
    x = 1/(1+exp(-RR+abs(offset)))
  } else {
    x = RR+offset
  }
  return(x)
}


dir = args[1]
modeldir = args[2]
#outdir = args[2]

logT = F
if (length(args) > 2) logT = as.logical(args[3])

#dir = '~/Data/Ellie/Analysis/VAL_Cohort'
#modeldir = '~/Data/Ellie/Analysis/5e6_arms_all_logR'

#dir = '~/Data/Ellie/Analysis/VAL_Cohort/'
files = list.files(dir, full.names=T)

file = paste(modeldir, 'model_data.Rdata', sep='/')
if (!file.exists(file))
  stop(paste("Missing data file", file))
load(file, verbose=T)
file = paste(modeldir, 'all.pt.alpha.Rdata', sep='/')
if (!file.exists(file))
  stop(paste("Missing data file", file))
load(file, verbose=T)
rm(dysplasia.df,coefs,plots,labels)

select.alpha = '0.9'
fitV = models[[select.alpha]]
lambda.opt = performance.at.1se[[select.alpha]][, 'lambda']


file = paste(modeldir, 'loo.Rdata', sep='/')
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
  
  tiled.segs <- tile.segmented.data(rs, size=5e6, chr.info=chrlen)
  tiled.segs = segment.matrix(tiled.segs)
  if (logT) tiled.segs = t(apply(tiled.segs, 1, logTransform))
  
  for (i in 1:ncol(tiled.segs)) 
    tiled.segs[,i] = unit.var(tiled.segs[,i], z.mean[i], z.sd[i])
  
  tiled.arms <- tile.segmented.data(rs, size='arms', chr.info=chrlen)
  tiled.arms = segment.matrix(tiled.arms)
  if (logT) tiled.arms = t(apply(tiled.arms, 1, logTransform))
  
  for (i in 1:ncol(tiled.arms)) 
    tiled.arms[,i] = unit.var(tiled.arms[,i], z.arms.mean[i], z.arms.sd[i])
  
  nrow(tiled.segs) == nrow(tiled.arms)
    
  cx.score = score.cx(tiled.segs,1)
  
  mergedDf = subtract.arms(tiled.segs, tiled.arms)
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

#load(paste(dir, 'predictions.Rdata', sep='/'), verbose=T)

riskPal = rev(RColorBrewer::brewer.pal(11, 'RdYlBu'))
unadjRisk = ggplot(pg.samp, aes(OR)) + geom_histogram(aes(fill=..x..), bins=10, show.legend = T) +
  scale_fill_gradientn(colors = riskPal,  name='RR') + 
  labs(y='n Samples', x='Relative Risk', title='Unadjusted relative risk') + theme_light(base_size = 14)
unadjRisk

m = melt(dspred, measure.vars=c('RR','Adj.RR'))
unadjRisk + geom_point(data=m, aes(x=value,y=Adj.Prob*30, color=variable))
  
ids = readxl::read_excel('~/Data/Ellie/Analysis/downsampled_ids.xlsx', sheet=1)

x = base::merge(ids, dspred, by.x='Illumina ID', by.y='row.names')

table(subset(x, `Expected Risk` == 'High')$Prob >= 0.7)
table(subset(x, `Expected Risk` == 'Low')$Prob <= 0.4)

table(subset(x, Tissue == 'Normal')$Prob <= 0.4)
subset(x, Tissue == 'Normal' & Prob > 0.4)

table(subset(x, Tissue == 'Barretts' & Study == 'OCCAMS')$Prob >= 0.7)

table(subset(x, `Expected Risk` == 'High' & Study == 'Annalise')$Prob >= 0.7)
table(subset(x, `Expected Risk` == 'Low' & Study == 'Annalise')$Prob <= 0.4)








