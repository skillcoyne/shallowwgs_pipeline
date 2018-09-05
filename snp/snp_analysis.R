
source('lib/data_func.R')

ptdirs = list.dirs('/Volumes/fh/fast/reid_b/collab/Killcoyne/Data/PerPatient', full.names=T, recursive=F)


mergedSegs = NULL
mergedArms = NULL
for (pt in ptdirs) {
  print(pt)
  segvals = as.data.frame(data.table::fread(list.files(pt, '*wins_tiled.txt', full.names=T)))
  armvals = as.data.frame(data.table::fread(list.files(pt, '*arms_tiled.txt', full.names=T)))
  
  segvals = segment.matrix(segvals)
  armvals = segment.matrix(armvals)

  if (is.null(segvals) | is.null(armvals))
    stop(paste("Missing values in ", pt))
  
  if (is.null(mergedSegs)) {
    mergedSegs = segvals
    mergedArms = armvals
  } else {
    mergedSegs = rbind(mergedSegs, segvals)    
    mergedArms = rbind(mergedArms, armvals)    
  }
print(nrow(mergedSegs))
}

nrow(mergedSegs) == nrow(mergedArms)
rownames(mergedSegs) == rownames(mergedArms)

load('~/Data/Ellie/Analysis/5e6_arms_all_logR/model_data.Rdata', verbose = T)
load('~/Data/Ellie/Analysis/5e6_arms_all_logR/all.pt.alpha.Rdata', verbose = T)
fitV = models$`0.9`
l = performance.at.1se$`0.9`$lambda
rm(plots,performance.at.1se,dysplasia.df,models,cvs,labels)

#x = melt(mergedSegs)
#ggplot(x, aes(Var2, value, color=Var1)) + geom_point() + theme(legend.position = 'none') 

for (i in 1:ncol(mergedSegs))
  mergedSegs[,i] = unit.var(mergedSegs[,i], z.mean[i], z.sd[i])
for (i in 1:ncol(mergedArms))
  mergedArms[,i] = unit.var(mergedArms[,i], z.arms.mean[i], z.arms.sd[i])

cx = score.cx(mergedSegs, 1)
hist(cx)


arrayDf = subtract.arms(mergedSegs, mergedArms)
arrayDf = cbind(arrayDf, 'cx'=unit.var(cx, mn.cx, sd.cx))

prob = predict(fitV, newx=arrayDf, s=l, type='response')
rr = predict(fitV, newx=arrayDf, s=l, type='link')

offset = log(0.0225)
preds = cbind.data.frame(prob, rr)
colnames(preds) = c('Prob', 'RR')
preds = preds %>% dplyr::mutate( 
  'Adj. RR'=RR+offset,  
  'Adj. Prob'=1/(1+exp(-RR+abs(offset)) )
)
rownames(preds) = rownames(prob)

patients = basename(ptdirs)


normals = preds[grep('BLD|gastric', rownames(preds), ignore.case=T),]
hist(normals$Prob)
length(which(normals$Prob < 0.5))/nrow(normals)

bld = preds[grep('BLD', rownames(preds), ignore.case=T),]
hist(bld$Prob)
length(which(bld$Prob < 0.5))/nrow(bld)
length(which(bld$`Adj. Prob` < 0.5))/nrow(bld)

be = preds[grep('BLD|gastric', rownames(preds), ignore.case=T, invert = T),]

length(unique(sapply(rownames(be[which(be$Prob > 0.5),]), function(x) unlist(strsplit(x, '_'))[1])))
length(unique(sapply(rownames(be[which(be$`Adj. Prob` > 0.5),]), function(x) unlist(strsplit(x, '_'))[1])))

length(unique(sapply(rownames(be), function(x) unlist(strsplit(x, '_'))[1])))

length(unique(sapply(rownames(be[which(be$Prob < 0.5),]), function(x) unlist(strsplit(x, '_'))[1])))


myPal = rev(RColorBrewer::brewer.pal(11, 'RdYlBu'))
ggplot(be, aes(RR)) + geom_histogram(aes(fill=..x..), bins=15, show.legend = F) +
  scale_fill_gradientn(colors = myPal,  name='') + 
  labs(y='n Samples', x='Relative Risk', title='Unadjusted relative risk') 

ggplot(be, aes(`Adj. RR`)) + geom_histogram(aes(fill=..x..), bins=15, show.legend = F) +
  scale_fill_gradientn(colors = myPal,  name='') + 
  labs(y='n Samples', x='Relative Risk', title='Adjusted relative risk') 




be = be %>% dplyr::mutate( 'Adj. RR'=RR+offset )
be = be %>% dplyr::mutate( 'Adj. Prob'=1/(1+exp(-RR+abs(offset)) ))

288*2
length(which(be$`Adj. Prob` <= 0.3))
length(which(be$`Adj. Prob` >= 0.7))



ggplot(be, aes(Prob)) + geom_histogram(aes(fill=..x..), breaks=seq(0,1,0.1) , show.legend = F) +
  scale_fill_gradientn(colors = myPal,  name='') + 
  labs(title='Predictions, all samples model', y='n Samples', x='Unadjusted Probability') 

ggplot(be, aes(`Adj. Prob`)) + geom_histogram(aes(fill=..x..), breaks=seq(0,1,0.1) , show.legend = F) +
  scale_fill_gradientn(colors = myPal,  name='') + 
  labs(title='Predictions, all samples model', y='n Samples', x='Adjusted Probability') 


