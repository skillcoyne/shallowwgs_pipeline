
args = commandArgs(trailingOnly=TRUE)

if (length(args) <2)
  stop("Missing required arguments: <model dir> <split dir>")

library(pROC)
library(tidyverse)


cache.dir = args[1] # paste(data, 'Analysis/5e6_arms_all', sep='/')
split.dir = args[2] # paste(data, 'Analysis/5e6_arms_splits', sep='/')


  rocs = NULL
  performance = NULL
  select.alpha = '0.9' 
  split.dirs = list.files(split.dir, full.names = T)
  for (split in 1:length(split.dirs)) {
    print(split)  
    models = list.files(split.dirs[split],full.names = T, pattern='all.*Rdata')
    loo = list.files(split.dirs[split],full.names = T, pattern='loo.*Rdata')
    
    if (length(models) < 1 || length(loo) < 1) {
      warning(paste("Missing model or LOO Rdata file in split", split.dirs[split]))
      next
    }

    load(list.files(split.dirs[split],full.names = T, pattern='all.*Rdata'), verbose=T)
    ncoefs = nrow(coefs[[select.alpha]])
    performance = rbind(performance, cbind(performance.at.1se[[select.alpha]],ncoefs))
    

    rm(plots,coefs,performance.at.1se,models,cvs)
    
    load(list.files(split.dirs[split],full.names = T, pattern='loo.*Rdata'), verbose=T)
    rm(plots, performance.at.1se,coefs,nzcoefs,fits)
  
    predictions = do.call(rbind, pg.samp)
 head(pg.samp) 
  preds = pg.samp %>% dplyr::select(Status, Prediction)
  roc = pROC::roc(Status ~ Prediction, data=preds, auc=T, ci=T, of='thresholds', transpose = T)
    #rocs = rbind(rocs,c(ci.auc(roc), coords(roc, 'best') ))
    roc$model = 'in-silico'
    roc$is = split
    rocs[[as.character(basename(split.dirs[split]))]] = roc
  
    rm(dysplasia.df, labels)
  }
  #colnames(rocs) = c('auc.ci.min','auc','auc.ci.max','threshold','specificity','sensitivity')
  
  insilico.rocs=rocs
  insilico.perf=performance
  save(insilico.rocs, insilico.perf, file=paste(cache.dir, 'insilico_splits.Rdata', sep='/'))

print("Finished")
