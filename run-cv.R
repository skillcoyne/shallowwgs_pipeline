library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)
library(gridExtra)
library(GenomicRanges)
library(glmnet)
library(CoxHD)
library(mice)

source('lib/load_patient_metadata.R')
source('lib/cv-pt-glm.R')


args = commandArgs(trailingOnly=TRUE)

data = args[1]
tile.w = as.character(args[2])
includeDemo = args[3]


data = '~/Data/Ellie'
tile.w=10e6
includeDemo = F # demographics

data.files = list.files(paste(data, 'QDNAseq',sep='/'), full.names=T)
analysis.files = list.files(paste(data, 'Analysis', sep='/'), full.names=T)

load(grep('All_patients.Rdata', analysis.files, value=T), verbose=T)

## Hospital.Research.ID info file
patient.file = grep('All_patient_info.xlsx', data.files, value=T)
if (length(patient.file) != 1)
  stop(paste("Missing/too many patient info file(s) in", data))
demo.file = grep('Demographics_full.xlsx', data.files, value=T)


all.patient.info = read.patient.info(patient.file, demo.file, set='all')$info

patient.info = subset(all.patient.info, Set == 'Training')

validation.patient.info =  subset(all.patient.info, Set == 'Test')
validation.patient.info = arrange(validation.patient.info, Status, Hospital.Research.ID, Endoscopy.Year, Pathology)
sum.validation.info = summarise.patient.info(validation.patient.info)

patient.info = arrange(patient.info, Status, Hospital.Research.ID, Endoscopy.Year, Pathology)
sum.patient.data = summarise.patient.info(patient.info)

validation.patient.data = patient.data[sum.validation.info$Hospital.Research.ID]
length(validation.patient.data)

patient.data = patient.data[sum.patient.data$Hospital.Research.ID]
patient.data = patient.data[!is.na(names(patient.data))]
length(patient.data)

sum.patient.data = as.data.frame(subset(sum.patient.data, Hospital.Research.ID %in% names(patient.data))) ## one patient was removed from the study
nrow(sum.patient.data)


tile.files = grep(as.character(tile.w), analysis.files, value=T, fixed=T)
load(tile.files[grep('tile', tile.files)], verbose=T)

mergedDf = mergedDf[grep('X|Y', rownames(mergedDf), invert=T),]


mergedDfValidation = mergedDf[,validation.patient.info$Samplename]
dim(mergedDfValidation)

mergedDf = mergedDf[, intersect(colnames(mergedDf), patient.info$Samplename)]
dim(mergedDf)

cache.dir = paste(data, 'Analysis',sub('\\+', '', tile.w),   sep='/')

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

unit.var <- function(x) {
  if (sd(x, na.rm=T) == 0) return(x)
  return((x-mean(x,na.rm=T))/sd(x,na.rm=T) )
}

dysplasia.df = apply(dysplasia.df, 2, unit.var)


score.cx <- function(pt.d, df) {
  get.length.variance <- function(pd) {
    samples = grep('^D', colnames(pd$seg.vals), value=T)
    len.var = var(pd$seg.vals$end.pos - pd$seg.vals$start.pos) # length won't vary between samples
    
    sample.var = apply(as.matrix(pd$seg.vals[, samples]), 2, var)
    t(cbind(len.var, sample.var))
  }
  
  complexity.measures = lapply(pt.d, get.length.variance)
  
  df = cbind(df, 'cx' = apply(df, 1, function(x) length(which(x >= sd(x)*2 | x <= -sd(x)*2))))
  df[,'cx'] = unit.var(df[,'cx'])
  
  df = cbind(df, 'segment.length'=0, 'sample.variance'=0)
  for (pt in names(complexity.measures)) {
    for (sample in colnames(complexity.measures[[pt]])) {
      df[sample,c('segment.length', 'sample.variance')] = as.numeric(complexity.measures[[pt]][,sample])
    }
  }
  
  df[, 'segment.length'] = unit.var(df[, 'segment.length'])
  df[, 'sample.variance'] = unit.var(df[, 'sample.variance'])
  
  return(df)
}

## TODO score complexity on the seg.vals directly rather than merged data frame?
dysplasia.df = score.cx(patient.data, dysplasia.df)

if (includeDemo) {
  
  ## Add in demographics
  
  # sort by samples in the order of the matrix
  patient.info = patient.info[match(rownames(dysplasia.df), patient.info$Samplename),]
  
  patient.info[which(patient.info$BMI > 100), 'BMI'] = NA # I think this one is wrong
  
  demo.data = patient.info[,c('Sex','Circumference','Maximal', 'Age.at.diagnosis', 'BMI', 'p53.Status', 'Smoking')]
  demo.data$Sex = as.integer(demo.data$Sex)-1
  demo.data$Smoking = as.integer(demo.data$Smoking)-1
  demo.data$p53.Status = as.integer(demo.data$p53.Status)-1
  demo.cols = c('Sex','C','M','Age', 'BMI', 'p53.Status', 'Smoking')
  
  encode.age <- function(x) {
    if (x > 70) { 
      (x - 70)/6 
    } else if (x < 50) {
      (x -50)/6
    } else {
      0
    }
  }
  #demo.data$Age.at.diagnosis = scale(demo.data$Age.at.diagnosis)
  
  # Encode based on presumed distributions or unit normalize
  demo.data = demo.data %>% rowwise() %>% mutate( 'age.encoded'= encode.age(Age.at.diagnosis))
  demo.data$Age.at.diagnosis = demo.data$age.encoded
  demo.data$age.encoded = NULL
  demo.data$BMI = unit.var(demo.data$BMI)
  
  demo.data[,c('Circumference','Maximal')] = apply(demo.data[,c('Circumference','Maximal')], 2, unit.var)
  
  ## Impute missing data
  imp = mice(as.matrix(demo.data), diagnostics = F)
  
  dysplasia.df = cbind(dysplasia.df, 'Sex'=NA,'C'=NA,'M'=NA,'Age'=NA, 'BMI'=NA, 'p53.Status'=NA, 'Smoking'=NA)
  dysplasia.df[,demo.cols] = as.matrix(complete(imp))
  
}


nl = 1000;folds = 10; splits = 5 

coefs = list(); plots = list(); performance.at.1se = list(); models = list(); cvs = list()
alpha.values = c(0, 0.5,0.7,0.8,0.9,1)
file = paste(data, 'Analysis/patient_folds.tsv', sep='/')
if (file.exists(file)) {
  message(paste("Reading folds file", file))
  sets = read.table(file, header=T, sep='\t')
} else {
  sets = create.patient.sets(pts, folds, splits, 0.2)  
  write.table(sets, quote=F, sep='\t', row.names=F, file=file)
  stop("Missing patient folds file")
}

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
  ggsave(paste(cache.dir, '/', as.character(tile.w),'_all_samples_cv.png',sep=''), plot = p, scale = 1, width = 10, height = 10, units = "in", dpi = 300)
}


## No HGD
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
  ggsave(paste(cache.dir, '/', as.character(tile.w),'_nohgd_samples_cv.png',sep=''), plot = p, scale = 1, width = 10, height = 10, units = "in", dpi = 300)
}

# ----------------- #

## No LGD
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
  ggsave(paste(cache.dir, '/', as.character(tile.w),'_nolgd_samples_cv.png',sep=''), plot = p, scale = 1, width = 10, height = 10, units = "in", dpi = 300)
}

# LOO
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

file = paste(cache.dir, 'loo.Rdata', sep='/')
select.alpha = 0.9
if (file.exists(file)) {
  load(file, verbose=T)
  # Fix info
  pg.samp = lapply(pg.samp, function(df) {
    merge(df[,c('Samplename',setdiff(colnames(df), colnames(patient.info))) ], subset(patient.info, Patient == unique(df$Patient)), by='Samplename')
  })
  
} else {
  performance.at.1se = c(); coefs = list(); plots = list(); fits = list()
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
    a = select.alpha
    fitLOO <- glmnet(training, labels[train.rows], alpha=a, family='binomial', nlambda=nl, standardize=F) # all patients
    l = fitLOO$lambda
    cv = crossvalidate.by.patient(x=training, y=labels[train.rows], lambda=l, 
                                  pts=subset(pts, Samplename %in% samples), a=a, nfolds=folds, splits=splits, fit=fitLOO, standardize=F)

    plots[[pt]] = arrangeGrob(cv$plot+ggtitle('Classification'), cv$deviance.plot+ggtitle('Binomial Deviance'), top=pt, ncol=2)
    
    fits[[pt]] = cv  
    
    if ( length(cv$lambda.1se) > 0 ) {
      performance.at.1se = c(performance.at.1se, subset(cv$lambdas, lambda == cv$lambda.1se)$mean)
      
      coef.1se = as.data.frame(non.zero.coef(fitLOO, cv$lambda.1se))
      coefs[[pt]] = coef.stability(coef.1se, cv$non.zero.cf)
      
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
  save(plots, performance.at.1se, coefs, fits, pg.samp, file=file)
}

message("Finished")