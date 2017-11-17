#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(mice))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(GenomicRanges))

source('lib/load_patient_metadata.R')
source('lib/cv-pt-glm.R')

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3)
  stop("Missing required params: <data dir> <tile size/list> <demo:T/F>")

data = args[1]
tile.w = as.numeric(args[2])
if (is.na(tile.w)) tile.w = unlist(strsplit(args[2],','))

includeDemo = args[3]

print(args)
#data = '~/Data/Ellie'

#tile.w=c(5e+06, 'arms')
tileFile = paste('tile_patients_All_',as.character(tile.w),'.Rdata',sep='')
message( paste(tileFile, collapse=' '))
#includeDemo = F # demographics

data.files = list.files(paste(data, 'QDNAseq',sep='/'), full.names=T)
analysis.files = list.files(paste(data, 'Analysis', sep='/'), full.names=T)

tile.files = sapply(tileFile, function(f) grep(f, analysis.files, value=T,fixed=T))

if (length(tile.files) <= 0)
  stop(paste("No file ", tileFile, " found in ", data, '/Analysis', sep=''))

load(grep('All_patients.Rdata', analysis.files, value=T), verbose=T)

## Hospital.Research.ID info file
patient.file = grep('All_patient_info.xlsx', data.files, value=T)
if (length(patient.file) != 1)
  stop(paste("Missing/too many patient info file(s) in", data))
demo.file = grep('Demographics_full.xlsx', data.files, value=T)

all.patient.info = read.patient.info(patient.file, demo.file, set='all')$info

patient.info = subset(all.patient.info, Set == 'Training')

patient.info = plyr::arrange(patient.info, Status, Hospital.Research.ID, Endoscopy.Year, Pathology)
sum.patient.data = summarise.patient.info(patient.info)

patient.data = patient.data[sum.patient.data$Hospital.Research.ID]
patient.data = patient.data[!is.na(names(patient.data))]
length(patient.data)

sum.patient.data = as.data.frame(subset(sum.patient.data, Hospital.Research.ID %in% names(patient.data))) ## one patient was removed from the study
nrow(sum.patient.data)

message(paste("Loading df from", tile.files))

if (length(tile.files) > 1) {
  load(grep('arms',tile.files,value=T), verbose=T)
  armsDf = mergedDf
  armsDf = armsDf[grep('X|Y', rownames(armsDf), invert=T),]
  armsDf = armsDf[, intersect(colnames(armsDf), patient.info$Samplename)]
  dim(armsDf)
}

load(grep(tile.w[which(is.na(sapply(tile.w, function(x) as.numeric(x))))], tile.files, invert=T, value=T), verbose=T)
mergedDf = mergedDf[grep('X|Y', rownames(mergedDf), invert=T),]
mergedDf = mergedDf[, intersect(colnames(mergedDf), patient.info$Samplename)]
dim(mergedDf)


cache.dir = paste(data, 'Analysis',sub('\\+', '', paste(tile.w, collapse='_')), sep='/')
if (includeDemo)
  cache.dir = paste(cache.dir, '_with-demo', sep='')
message(paste("Writing to ", cache.dir))

if (!dir.exists(cache.dir)) dir.create(cache.dir)

## binomial: prog 1, np 0
labels = unlist(sapply(patient.data, function(df) {
  df$info = plyr::arrange(df$info, Endoscopy.Year, Pathology)
  label = as.integer(df$info$Status == 'P') #as.integer(df$info$Pathology %in% c('HGD', 'IMC'))
  names(label) = df$info$Samplename
  return(label)
}))
names(labels) = sub('.*\\.', '',  names(labels))

pts = do.call(rbind, lapply(patient.data, function(df) {
  df$info = arrange(df$info, Endoscopy.Year, Pathology)
  cbind(df$info[,c('Hospital.Research.ID','Samplename')])
}))
rownames(pts) = 1:nrow(pts)

# sort in label order
if (length(setdiff(colnames(mergedDf), names(labels))) > 3)
  warning("Labels vector is missing samples")

dysplasia.df = t(mergedDf[,intersect(names(labels), colnames(mergedDf))])
labels = labels[intersect(names(labels), colnames(mergedDf))]
dim(dysplasia.df)
dysplasia.df = apply(dysplasia.df, 2, unit.var)

get.loc<-function(df) {
  locs = do.call(rbind.data.frame, lapply(colnames(df), function(x) unlist(strsplit( x, ':|-'))))
  colnames(locs) = c('chr','start','end')
  locs[c('start','end')] = lapply(locs[c('start','end')], function(x) as.numeric(as.character(x)))
  locs$chr = factor(locs$chr, levels=c(1:22), ordered=T)
  locs
}


# Merge and subtract the average value of the arms from the segments
if (exists('armsDf')) {
  arms.df = t(armsDf[,intersect(names(labels), colnames(armsDf))])
  labels = labels[intersect(names(labels), colnames(armsDf))]
  dim(arms.df)
  arms.df = apply(arms.df, 2, unit.var)
  
  locs = get.loc(dysplasia.df)
  arm.locs = get.loc(arms.df)
  
  arms = makeGRangesFromDataFrame(arm.locs)
  locs = makeGRangesFromDataFrame(locs)
  
  tmp = dysplasia.df
  # subtract arms from 5e6 and merge both
  ov = findOverlaps(arms, locs)
  for (hit in unique(queryHits(ov))) {
    print(hit)
    cols = subjectHits(ov)[which(queryHits(ov) == hit)]
    for (i in 1:nrow(tmp)) {
      tmp[i,cols] = tmp[i,cols] - arms.df[i,hit]
    }
  }
  dysplasia.df = cbind(tmp, arms.df)
  
  
  loc = get.loc(dysplasia.df)
  dysplasia.df = dysplasia.df[,order(loc$chr, loc$start)]

  dim(dysplasia.df)
}

## TODO score complexity on the seg.vals directly rather than merged data frame?
dysplasia.df = score.cx(patient.data, dysplasia.df)

if (includeDemo) {
  ## Add in demographics
  dysplasia.df = add.demo.tocv(patient.info, dysplasia.df)
}


nl = 1000;folds = 10; splits = 5 

file = paste(data, 'Analysis/patient_folds.tsv', sep='/')
if (file.exists(file)) {
  message(paste("Reading folds file", file))
  sets = read.table(file, header=T, sep='\t')
} else {
  sets = create.patient.sets(pts, folds, splits, 0.2)  
  write.table(sets, quote=F, sep='\t', row.names=F, file=file)
  stop("Missing patient folds file")
}


alpha.values = c(0, 0.5,0.7,0.8,0.9,1)
## ----- All ----- ##
coefs = list(); plots = list(); performance.at.1se = list(); models = list(); cvs = list()
file = paste(cache.dir, 'all.pt.alpha.Rdata', sep='/')
if (file.exists(file)) {
  message(paste("loading", file))
  load(file, verbose=T)
} else {
  for (a in alpha.values) {
    fit0 <- glmnet(dysplasia.df, labels, alpha=a, nlambda=nl, family='binomial', standardize=F) # all patients
    
    l = fit0$lambda
    #l = more.l(fit0$lambda)
    
    cv.patient = crossvalidate.by.patient(x=dysplasia.df, y=labels, lambda=l, pts=sets, a=a, nfolds=folds, splits=splits, fit=fit0, select='deviance', opt=-1, standardize=F)
    
    lambda.opt = cv.patient$lambda.1se
    
    coef.opt = as.data.frame(non.zero.coef(fit0, lambda.opt))
    coefs[[as.character(a)]] = coef.stability(coef.opt, cv.patient$non.zero.cf)
    
    plots[[as.character(a)]] = arrangeGrob(cv.patient$plot+ggtitle('Classification'), cv.patient$deviance.plot+ggtitle('Binomial Deviance'), top=paste('alpha=',a,sep=''), ncol=2)
    
    performance.at.1se[[as.character(a)]] = subset(cv.patient$lambdas, `lambda-1se`)
    models[[as.character(a)]] = fit0
    cvs[[as.character(a)]] = cv.patient
  }
  save(plots, coefs, performance.at.1se, dysplasia.df, models, cv.patient, labels, file=file)
  p = do.call(grid.arrange, c(plots[ as.character(alpha.values) ], top='All samples, 10fold, 5 splits'))
  ggsave(paste(cache.dir, '/', 'all_samples_cv.png',sep=''), plot = p, scale = 1, width = 10, height = 10, units = "in", dpi = 300)
}
all.coefs = coefs
## --------------------- ##

## --------- No HGD --------- ##
`%nin%` <- Negate(`%in%`)

info = do.call(rbind, lapply(patient.data, function(df) df$info))

no.hgd.plots = list(); coefs = list(); performance.at.1se = list()
message("No HGD/IMC")
file = paste(cache.dir, 'nohgd.Rdata', sep='/')
if (file.exists(file)) {
  message(paste("loading file", file))
  load(file, verbose=T)
} else {
  # No HGD/IMC
  samples = intersect(rownames(dysplasia.df), subset(info, Pathology %nin% c('HGD', 'IMC'))$Samplename)
  for (a in alpha.values) {
    # all patients
    fitNoHGD <- glmnet(dysplasia.df[samples,], labels[samples], alpha=a, family='binomial', nlambda=nl, standardize = F) 
    cv.nohgd = crossvalidate.by.patient(x=dysplasia.df[samples,], y=labels[samples], lambda=fitNoHGD$lambda, pts=subset(pts, Samplename %in% samples), a=a, nfolds=folds, splits=splits, fit=fitNoHGD, standardize = F)
    
    no.hgd.plots[[as.character(a)]] = arrangeGrob(cv.nohgd$plot+ggtitle('Classification'), cv.nohgd$deviance.plot+ggtitle('Binomial Deviance'), top=paste('alpha=',a,sep=''), ncol=2)
    
    coef.1se = as.data.frame(non.zero.coef(fitNoHGD, cv.nohgd$lambda.1se))
    coefs[[as.character(a)]] = coef.stability(coef.1se, cv.nohgd$non.zero.cf)
    
    performance.at.1se[[as.character(a)]] = subset(cv.nohgd$lambdas, lambda == cv.nohgd$lambda.1se)
  }
  save(no.hgd.plots, coefs, performance.at.1se, file=file)
  p = do.call(grid.arrange, c(no.hgd.plots, top='No HGD/IMC samples, 10fold, 5 splits'))
  ggsave(paste(cache.dir, '/', 'nohgd_samples_cv.png',sep=''), plot = p, scale = 1, width = 10, height = 10, units = "in", dpi = 300)
}
# ----------------- #

## --------- No LGD ------------ ##
file = paste(cache.dir, 'nolgd.Rdata', sep='/')
message("No LGD")
nolgd.plots = list(); coefs = list(); performance.at.1se = list()
if (file.exists(file)) {
  load(file, verbose=T)
} else {
  # No LGD
  samples = intersect(rownames(dysplasia.df), c(subset(info, Status == 'NP')$Samplename, subset(info, Pathology %nin% c('HGD', 'IMC', 'LGD') & Status == 'P')$Samplename))
  
  # No HGD/IMC/LGD in all patients
  samples = intersect(rownames(dysplasia.df), subset(info, Pathology %nin% c('HGD', 'IMC', 'LGD'))$Samplename)
  
  for (a in alpha.values) {
    fitNoLGD <- glmnet(dysplasia.df[samples,], labels[samples], alpha=a, family='binomial', nlambda=nl, standardize = F) # all patients
    
    cv.nolgd = crossvalidate.by.patient(x=dysplasia.df[samples,], y=labels[samples], lambda=fitNoLGD$lambda, pts=subset(pts, Samplename %in% samples), a=a, nfolds=folds, splits=splits, fit=fitNoLGD, standardize = F)
    
    nolgd.plots[[as.character(a)]] = arrangeGrob(cv.nolgd$plot+ggtitle('Classification'), cv.nolgd$deviance.plot+ggtitle('Binomial Deviance'), top=paste('alpha=',a,sep=''), ncol=2)
    
    coef.1se = as.data.frame(non.zero.coef(fitNoLGD, cv.nolgd$lambda.1se))
    coefs[[as.character(a)]] = coef.stability(coef.1se, cv.nolgd$non.zero.cf)
    
    performance.at.1se[[as.character(a)]] = subset(cv.nolgd$lambdas, lambda == cv.nolgd$lambda.1se)
  }
  save(nolgd.plots, coefs, performance.at.1se, file=file)
  p = do.call(grid.arrange, c(nolgd.plots, top='No HGD/IMC/LGD samples, 10fold, 5 splits'))
  ggsave(paste(cache.dir, '/', 'nolgd_samples_cv.png',sep=''), plot = p, scale = 1, width = 10, height = 10, units = "in", dpi = 300)
}
## --------------------- ##

## --------- LOO --------- ##
message("LOO")
info = do.call(rbind, lapply(patient.data, function(df) df$info))
pg.samp = lapply(patient.data, function(pt) {
  info = pt$info
  info$SampleSD = NA
  info$SampleMEAN = NA
  if (length(info$Samplename) > 1) {
    info$SampleSD = apply(pt$seg.vals[,info$Samplename], 2, sd)
    info$SampleMEAN =  apply(pt$seg.vals[,info$Samplename], 2, mean)
  }
  info$Prediction = NA
  info$Prediction.Dev.Resid = NA
  info$PID = unlist(lapply(info$Path.ID, function(x) unlist(strsplit(x, 'B'))[1]))
  return(info)
})

select.alpha = 0.9


file = paste(cache.dir, 'loo.Rdata', sep='/')
if (file.exists(file)) {
  load(file, verbose=T)
  # Fix info
  pg.samp = lapply(pg.samp, function(df) {
    merge(df[,c('Samplename',setdiff(colnames(df), colnames(patient.info))) ], subset(patient.info, Patient == unique(df$Patient)), by='Samplename')
  })
} else {
  secf = all.coefs[[select.alpha]]
  a = select.alpha
  
  performance.at.1se = c(); coefs = list(); plots = list(); fits = list(); nzcoefs = list()
  # Remove each patient (LOO)
  for (pt in names(pg.samp)) {
    print(pt)
    tmp.patient.data = patient.data[subset(sum.patient.data, Hospital.Research.ID != pt)$Hospital.Research.ID]
    samples = as.vector(unlist(sapply(tmp.patient.data, function(df) df$info$Samplename )))
    
    train.rows = which(rownames(dysplasia.df) %in% samples)
    training = dysplasia.df[train.rows,]
    test = as.matrix(dysplasia.df[-train.rows,])
    #if (ncol(test) <= 1) next
    if ( nrow(test) == ncol(dysplasia.df) ) test = t(test)
    
    # Predict function giving me difficulty when I have only a single sample, this ensures the dimensions are the same
    sparsed_test_data <- Matrix(data=0, nrow=ifelse(length(pg.samp[[pt]]$Samplename) > 1, nrow(test), 1),  ncol=ncol(training),
                                dimnames=list(pg.samp[[pt]]$Samplename,colnames(training)), sparse=T)
    for(i in colnames(dysplasia.df)) sparsed_test_data[,i] = test[,i]
    
    # Fit generated on all samples, including HGD
    fitLOO <- glmnet(training, labels[train.rows], alpha=a, family='binomial', nlambda=nl, standardize=F) # all patients
    l = fitLOO$lambda

    #pf = rep(1, ncol(dysplasia.df))
    #pf[which(colnames(dysplasia.df) %in% rownames(secf))] = 0.01
    
    cv = crossvalidate.by.patient(x=training, y=labels[train.rows], lambda=l, a=a, nfolds=folds, splits=splits,
                                  pts=subset(pts, Samplename %in% samples), fit=fitLOO, standardize=F)


    plots[[pt]] = arrangeGrob(cv$plot+ggtitle('Classification'), cv$deviance.plot+ggtitle('Binomial Deviance'), top=pt, ncol=2)
    
    fits[[pt]] = cv  
    
    if ( length(cv$lambda.1se) > 0 ) {
      performance.at.1se = c(performance.at.1se, subset(cv$lambdas, lambda == cv$lambda.1se)$mean)
      
      #coef.1se = coef(fitLOO, cv$lambda.1se)[rownames(secf),]
        
      nzcoefs[[pt]] = as.data.frame(non.zero.coef(fitLOO, cv$lambda.1se))

      coefs[[pt]] = coef(fitLOO, cv$lambda.1se)[rownames(secf),]
        #coef.stability(coef.1se, cv$non.zero.cf)
      
      logit <- function(p){log(p/(1-p))}
      inverse.logit <- function(or){1/(1 + exp(-or))}
      
      pm = predict(fitLOO, newx=sparsed_test_data, s=cv$lambda.1se, type='response')
      sy = as.matrix(sqrt(binomial.deviance(pm, labels[pg.samp[[pt]]$Samplename])))
      
      pg.samp[[pt]]$Prediction = pm[,1]
      pg.samp[[pt]]$Prediction.Dev.Resid = sy[,1] 
      
    } else {
      warning(paste("Hospital.Research.ID", pt, "did not have a 1se"))
    }
  }
  save(plots, performance.at.1se, coefs, nzcoefs, fits, pg.samp, file=file)
}

message("Finished")