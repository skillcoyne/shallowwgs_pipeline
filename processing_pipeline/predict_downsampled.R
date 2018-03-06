# predict downsampled

library(ggplot2)
library(Matrix)
library(glmnet)

source('lib/data_func.R')
source('lib/load_patient_metadata.R')

chrlen = get.chr.lengths()

dir = '~/Data/Ellie/Analysis/downsampled5G'
files = list.files(dir, full.names=T)


file = '~/Data/Ellie/Analysis/5e6_arms_disc/model_data.Rdata'
if (!file.exists(file))
  stop(paste("Missing data file", file))
load(file, verbose=T)
file = '~/Data/Ellie/Analysis/5e6_arms_disc/all.pt.alpha.Rdata'
if (!file.exists(file))
  stop(paste("Missing data file", file))
load(file, verbose=T)
rm(dysplasia.df,coefs,plots,labels)

select.alpha = '0.9'
fitV = models[[select.alpha]]
lambda.opt = performance.at.1se[[select.alpha]][, 'lambda']

sampleNames = sub('_raw.txt', '', basename(files))
dspred = data.frame(matrix(ncol=2, nrow=length(samples), dimnames=list(samples,c('Prob','RR'))))

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
  
  dspred[n,] = c(preds[,1],RR[,1])
  print(dspred[1:n,])
}

save(dspred, file=paste(dir, 'predictions.Rdata', sep='/'))

range(dspred$RR)

#load(paste(dir, 'predictions.Rdata', sep='/'))

file = '~/Data/Ellie/Analysis/5e6_arms_disc/loo.Rdata'
if (!file.exists(file))
  stop(paste("Missing data file", file))
load(file, verbose=T)

pg.samp = do.call(rbind, pg.samp)

riskPal = rev(RColorBrewer::brewer.pal(11, 'RdYlBu'))
unadjRisk = ggplot(pg.samp, aes(OR)) + geom_histogram(aes(fill=..x..), bins=10, show.legend = T) +
  scale_fill_gradientn(colors = riskPal,  name='RR') + 
  labs(y='n Samples', x='Relative Risk', title='Unadjusted relative risk') + theme_light(base_size = 14)
unadjRisk

unadjRisk + geom_point(data=dspred, aes(x=RR, y=1)) + xlim(-15, 15)

ggplot(pg.samp, aes(Prediction)) + geom_histogram(aes(fill=..x..), bins=10, show.legend = T) +
  scale_fill_gradientn(colors = riskPal,  name='RR') + 
  labs(y='n Samples', x='Relative Risk', title=' relative risk') + theme_light(base_size = 14)
unadjRisk



preds = ggplot(pg.samp, aes(Prediction)) + geom_histogram(aes(fill=..x..), breaks=seq(0,1,0.1), show.legend = F) + 
  scale_fill_distiller(palette = 'RdYlBu', name='P(P)') +
  labs(title='Sample predictions', y='n Samples', x='Probability') + theme_light(base_size = 14)
preds

preds + geom_point(data=dspred, aes(x=Prob, y=(1:nrow(dspred))))
