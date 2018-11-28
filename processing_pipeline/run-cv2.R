args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3)
  stop("Missing required params: <data dir> <outdir> <info file dir> <log transform DEF: FALSE>")


suppressPackageStartupMessages( library(BarrettsProgressionRisk) )
suppressPackageStartupMessages( library(tidyverse) )
suppressPackageStartupMessages(library(glmnet))

suppressPackageStartupMessages(source('~/workspace/shallowwgs_pipeline/lib/data_func.R'))
suppressPackageStartupMessages(source('~/workspace/shallowwgs_pipeline/lib/load_patient_metadata.R'))
suppressPackageStartupMessages(source('~/workspace/shallowwgs_pipeline/lib/cv-pt-glm.R'))

data = args[1]
# data = '~/Data/BarrettsProgressionRisk/Analysis/multipcf_perPatient/'
outdir = args[2]
# outdir = '~/Data/BarrettsProgressionRisk/Analysis'
infodir = args[3]
# infodir = '~/Data/BarrettsProgressionRisk/QDNAseq'
logT = F
if (length(args) == 4)
  logT = as.logical(args[4])

cache.dir = paste(outdir, '5e6_arms_all', sep='/')
if (logT) cache.dir = paste(cache.dir, '_logR', sep='')
dir.create(cache.dir, recursive=T, showWarnings=F)

## Hospital.Research.ID info file
patient.file = list.files(infodir, pattern='All_patient_info.xlsx', recursive=T, full.names=T)
demo.file = list.files(infodir, pattern='Demographics_full.xlsx', recursive=T, full.names=T)

if ( length(patient.file) != 1 | length(demo.file) != 1)
  stop("Missing files in info dir: All_patient_info.xlsx and Demographics_full.xlsx")

patient.info = read.patient.info(patient.file, demo.file, set='All')$info %>% dplyr::arrange(Status, Patient, Endoscopy.Year, Pathology)
sum.patient.data = summarise.patient.info(patient.info)
#nrow(sum.patient.data)


cleaned = list.files(path=data, pattern='tiled_segvals', full.names=T, recursive=T)
#raw = list.files(path=data, pattern='raw', full.names=T, recursive=T)

cleaned = grep(paste(sum.patient.data$Hospital.Research.ID, collapse = '|'), cleaned, value=T)
#raw = grep(paste(sum.patient.data$Hospital.Research.ID, collapse = '|'), raw, value=T)

arms = grep('arms', cleaned, value=T)
segs = grep('arms', cleaned, invert=T, value=T)

length(segs) == length(arms)

seg.tiles = do.call(bind_rows, lapply(segs, function(x) readr::read_tsv(x, col_types=cols(.default=col_double(), X1=col_character()))))
colnames(seg.tiles)[1] = 'sample'
samples = seg.tiles$sample

seg.tiles = as.matrix(seg.tiles[,-1])
rownames(seg.tiles) = samples

if (logT) seg.tiles = t(apply(seg.tiles, 1, logTransform))

segsList = prep.matrix(seg.tiles)

z.mean = segsList$z.mean
z.sd = segsList$z.sd
segs = segsList$matrix

# Complexity score
cx.score = BarrettsProgressionRisk::scoreCX(segs,1)
mn.cx = mean(cx.score)
sd.cx = sd(cx.score)

# After mean centering set all NA values to 0
segs = segsList$matrix
z.mean = segsList$z.mean
z.sd = segsList$z.sd

# Complexity score
cx.score = BarrettsProgressionRisk::scoreCX(segs,1)
mn.cx = mean(cx.score)
sd.cx = sd(cx.score)

# Load arm files  
arm.tiles = do.call(bind_rows, lapply(arms, function(x) readr::read_tsv(x, col_types=cols(.default=col_double(), X1=col_character()))))
colnames(arm.tiles)[1] = 'sample'

samples = arm.tiles$sample
arm.tiles = as.matrix(arm.tiles[,-1])
rownames(arm.tiles) = samples
if (logT) arm.tiles = t(apply(arm.tiles, 1, logTransform))

armsList = prep.matrix(arm.tiles)
arms = armsList$matrix

z.arms.mean = armsList$z.mean
z.arms.sd = armsList$z.sd

allDf = BarrettsProgressionRisk::subtractArms(segs, arms)
allDf = cbind(allDf, 'cx'=BarrettsProgressionRisk:::unit.var(cx.score, mn.cx, sd.cx))

## labels: binomial: prog 1, np 0
sampleStatus = subset(patient.info, Samplename %in% rownames(allDf), select=c('Samplename','Status'))
sampleStatus$Status = as.integer(sampleStatus$Status)-1
labels = sampleStatus$Status
names(labels) = sampleStatus$Samplename

dim(allDf)
dysplasia.df = as.matrix(allDf[sampleStatus$Samplename,])
dim(dysplasia.df)

#raw.segs = seg.tiles[sampleStatus$Samplename,]; raw.arms = arm.tiles[sampleStatus$Samplename,]

save(dysplasia.df, labels, mn.cx, sd.cx, z.mean, z.sd, z.arms.mean, z.arms.sd, file=paste(cache.dir, 'model_data.Rdata', sep='/'))
rm(raw.segs, raw.arms)


nl = 1000;folds = 10; splits = 5 

#file = list.files(data, pattern='patient_folds.tsv', recursive=T, full.names=T)
#if (length(file) == 1 && !allPts) {
#  message(paste("Reading folds file", file))
#  sets = read.table(file, header=T, sep='\t')
#} else {
sets = create.patient.sets(patient.info[c('Hospital.Research.ID','Samplename','Status')], folds, splits, 0.2)  
#}

alpha.values = c(0, 0.5,0.7,0.8,0.9,1)
## ----- All ----- ##
coefs = list(); plots = list(); performance.at.1se = list(); models = list(); cvs = list()
file = paste(cache.dir, 'all.pt.alpha.Rdata', sep='/')
if (file.exists(file)) {
  message(paste("loading", file))
  load(file, verbose=T)
} else {
  for (a in alpha.values) {
    fit0 <- glmnet(dysplasia.df, labels, alpha=a, nlambda=nl, family='binomial', standardize=F)    
    l = fit0$lambda

    cv.patient = crossvalidate.by.patient(x=dysplasia.df, y=labels, lambda=l, pts=sets, a=a, nfolds=folds, splits=splits, fit=fit0, select='deviance', opt=-1, standardize=F)
    
    lambda.opt = cv.patient$lambda.1se
    
    coef.opt = as.data.frame(non.zero.coef(fit0, lambda.opt))
    coefs[[as.character(a)]] = coef.stability(coef.opt, cv.patient$non.zero.cf)
    
    plots[[as.character(a)]] = arrangeGrob(cv.patient$plot+ggtitle('Classification'), cv.patient$deviance.plot+ggtitle('Binomial Deviance'), top=paste('alpha=',a,sep=''), ncol=2)
    
    performance.at.1se[[as.character(a)]] = subset(cv.patient$lambdas, `lambda-1se`)
    models[[as.character(a)]] = fit0
    cvs[[as.character(a)]] = cv.patient
  }
  save(plots, coefs, performance.at.1se, dysplasia.df, models, cvs, labels, file=file)
  p = do.call(grid.arrange, c(plots[ as.character(alpha.values) ], top='All samples, 10fold, 5 splits'))
  ggsave(paste(cache.dir, '/', 'all_samples_cv.png',sep=''), plot = p, scale = 2, width = 12, height = 10, units = "in", dpi = 300)
}
all.coefs = coefs
# ----------------- #


## --------- No HGD --------- ##
`%nin%` <- Negate(`%in%`)

no.hgd.plots = list(); coefs = list(); performance.at.1se = list(); models = list()
message("No HGD/IMC")
file = paste(cache.dir, 'nohgd.Rdata', sep='/')
if (file.exists(file)) {
  message(paste("loading file", file))
  load(file, verbose=T)
} else {
  # No HGD/IMC
  samples = intersect(rownames(dysplasia.df), subset(patient.info, Pathology %nin% c('HGD', 'IMC'))$Samplename)
  for (a in alpha.values) {
    # all patients
    fitNoHGD <- glmnet(dysplasia.df[samples,], labels[samples], alpha=a, family='binomial', nlambda=nl, standardize = F) 
    
    l = fitNoHGD$lambda
    #if (a == 0)
    #  l = more.l(fitNoHGD$lambda)
    
    cv.nohgd = crossvalidate.by.patient(x=dysplasia.df[samples,], y=labels[samples], lambda=l, pts=subset(sets, Samplename %in% samples), a=a, nfolds=folds, splits=splits, fit=fitNoHGD, standardize=F)
    
    no.hgd.plots[[as.character(a)]] = arrangeGrob(cv.nohgd$plot+ggtitle('Classification'), cv.nohgd$deviance.plot+ggtitle('Binomial Deviance'), top=paste('alpha=',a,sep=''), ncol=2)
    
    coef.1se = as.data.frame(non.zero.coef(fitNoHGD, cv.nohgd$lambda.1se))
    coefs[[as.character(a)]] = coef.stability(coef.1se, cv.nohgd$non.zero.cf)
    
    performance.at.1se[[as.character(a)]] = subset(cv.nohgd$lambdas, lambda == cv.nohgd$lambda.1se)
    
    models[[as.character(a)]] = fitNoHGD
  }
  save(no.hgd.plots, coefs, performance.at.1se, models, file=file)
  p = do.call(grid.arrange, c(no.hgd.plots, top='No HGD/IMC samples'))
  ggsave(paste(cache.dir, '/', 'nohgd_samples_cv.png',sep=''),  plot = p, scale = 2, width = 12, height = 10, units = "in", dpi = 300)
}
# ----------------- #

## --------- No LGD ------------ ##
file = paste(cache.dir, 'nolgd.Rdata', sep='/')
message("No LGD")
nolgd.plots = list(); coefs = list(); performance.at.1se = list()
if (file.exists(file)) {
  message(paste("loading file", file))
  load(file, verbose=T)
} else {
  # No HGD/IMC/LGD in all patients
  samples = intersect(rownames(dysplasia.df), subset(patient.info, Pathology %nin% c('HGD', 'IMC', 'LGD'))$Samplename)
  
  for (a in alpha.values) {
    fitNoLGD <- glmnet(dysplasia.df[samples,], labels[samples], alpha=a, family='binomial', nlambda=nl, standardize = F) # all patients
    
    l = fitNoLGD$lambda

    cv.nolgd = crossvalidate.by.patient(x=dysplasia.df[samples,], y=labels[samples], lambda=l, pts=subset(sets, Samplename %in% samples), a=a, nfolds=folds, splits=splits, fit=fitNoLGD, standardize = F)
    
    nolgd.plots[[as.character(a)]] = arrangeGrob(cv.nolgd$plot+ggtitle('Classification'), cv.nolgd$deviance.plot+ggtitle('Binomial Deviance'), top=paste('alpha=',a,sep=''), ncol=2)
    
    coef.1se = as.data.frame(non.zero.coef(fitNoLGD, cv.nolgd$lambda.1se))
    coefs[[as.character(a)]] = coef.stability(coef.1se, cv.nolgd$non.zero.cf)
    
    performance.at.1se[[as.character(a)]] = subset(cv.nolgd$lambdas, lambda == cv.nolgd$lambda.1se)
  }
  save(nolgd.plots, coefs, performance.at.1se, file=file)
  p = do.call(grid.arrange, c(nolgd.plots, top='No LGD/HGD/IMC samples, 10fold, 5 splits'))
  ggsave(paste(cache.dir, '/', 'nolgd_samples_cv.png',sep=''),   plot = p, scale = 2, width = 12, height = 10, units = "in", dpi = 300)
}
## --------------------- ##

## --------- LOO --------- ##
message("LOO")


pg.samp = patient.info %>% rowwise %>% dplyr::mutate(
  PID = sub('_$', '', unlist(strsplit(Path.ID, 'B'))[1])
) %>% filter(Samplename %in% rownames(dysplasia.df))

select.alpha = 0.9

file = paste(cache.dir, 'loo.Rdata', sep='/')
if (file.exists(file)) {
  message(paste("loading file", file))
  load(file, verbose=T)
} else {
  secf = all.coefs[[select.alpha]]
  a = select.alpha
  
  performance.at.1se = c(); coefs = list(); plots = list(); fits = list(); nzcoefs = list()
  # Remove each patient (LOO)
  for (pt in unique(pg.samp$Hospital.Research.ID)) {
  #for (pt in names(pg.samp)) {
    print(pt)
    #samples = subset(patient.info, Hospital.Research.ID != pt)$Samplename
    samples = subset(pg.samp, Hospital.Research.ID != pt)$Samplename
    
    train.rows = which(rownames(dysplasia.df) %in% samples)
    training = dysplasia.df[train.rows,]
    test = as.matrix(dysplasia.df[-train.rows,])
    #if (ncol(test) <= 1) next
    if ( nrow(test) == ncol(dysplasia.df) ) test = t(test)
    
    
    patient.samples = pg.samp %>% filter(Hospital.Research.ID == pt) %>% select(Samplename)
    # Predict function giving me difficulty when I have only a single sample, this ensures the dimensions are the same
    sparsed_test_data <- Matrix(data=0, nrow=ifelse(nrow(patient.samples) > 1, nrow(test), 1),  ncol=ncol(training),
                                dimnames=list(rownames(test),colnames(training)), sparse=T)
    
    for(i in colnames(dysplasia.df)) sparsed_test_data[,i] = test[,i]
    
    # Fit generated on all samples, including HGD
    fitLOO <- glmnet(training, labels[train.rows], alpha=a, family='binomial', nlambda=nl, standardize=F) # all patients
    l = fitLOO$lambda
    
    cv = crossvalidate.by.patient(x=training, y=labels[train.rows], lambda=l, a=a, nfolds=folds, splits=splits,
                                  pts=subset(sets, Samplename %in% samples), fit=fitLOO, standardize=F)
    
    plots[[pt]] = arrangeGrob(cv$plot+ggtitle('Classification'), cv$deviance.plot+ggtitle('Binomial Deviance'), top=pt, ncol=2)
    
    fits[[pt]] = cv  
    
    if ( length(cv$lambda.1se) > 0 ) {
      performance.at.1se = c(performance.at.1se, subset(cv$lambdas, lambda == cv$lambda.1se)$mean)
      
      nzcoefs[[pt]] = as.data.frame(non.zero.coef(fitLOO, cv$lambda.1se))
      
      coefs[[pt]] = coef(fitLOO, cv$lambda.1se)[rownames(secf),]

      logit <- function(p){log(p/(1-p))}
      inverse.logit <- function(or){1/(1 + exp(-or))}
      
      pm = predict(fitLOO, newx=sparsed_test_data, s=cv$lambda.1se, type='response')
      colnames(pm) = 'Prediction'
      
      rr = predict(fitLOO, newx=sparsed_test_data, s=cv$lambda.1se, type='link')
      colnames(rr) = 'RR'
      
      sy = as.matrix(sqrt(binomial.deviance(pm,labels[intersect(patient.samples$Samplename, names(labels))])))
      colnames(sy) = 'Prediction.Dev.Resid'
      
      pg.samp = pg.samp %>% 
        full_join(as_tibble(pm, rownames='Samplename'), by='Samplename') %>% 
        full_join(as_tibble(rr, rownames='Samplename'), by='Samplename') %>% 
        full_join(as_tibble(sy, rownames='Samplename'), by='Samplename')
      
    } else {
      warning(paste("Hospital.Research.ID", pt, "did not have a 1se"))
    }
  }
  save(plots, performance.at.1se, coefs, nzcoefs, fits, pg.samp, file=file)
}

## --------- LOO NO HGD--------- ##
message("LOO NO HGD")

pg.sampNOHGD = patient.info %>% rowwise %>% dplyr::mutate(
  PID = sub('_$', '', unlist(strsplit(Path.ID, 'B'))[1])
) %>% filter(Samplename %in% rownames(dysplasia.df) & Pathology %nin% c('HGD','IMC'))

select.alpha = 0.9

file = paste(cache.dir, 'loonohgd.Rdata', sep='/')
if (file.exists(file)) {
  load(file, verbose=T)
} else {
  dysplasia.dfNOHGD = dysplasia.df[intersect(rownames(dysplasia.df),subset(patient.info, Pathology %nin% c('HGD','IMC'))$Samplename),]
  
  secf = all.coefs[[select.alpha]]
  a = select.alpha
  
  performance.at.1se = c(); coefs = list(); plots = list(); fits = list(); nzcoefs = list()
  # Remove each patient (LOO)
  for (pt in unique(pg.sampNOHGD$Hospital.Research.ID)) {
  #for (pt in names(pg.sampNOHGD)) {
    print(pt)

    samples = subset(pg.sampNOHGD, Hospital.Research.ID != pt)$Samplename

    train.rows = which(rownames(dysplasia.dfNOHGD) %in% samples)
    training = dysplasia.df[train.rows,]
    test = as.matrix(dysplasia.dfNOHGD[-train.rows,])
    
    if (nrow(test) <= 0) next
    
    #if (ncol(test) <= 1) next
    if ( nrow(test) == ncol(dysplasia.df) ) test = t(test)
    
    patient.samples = pg.sampNOHGD %>% filter(Hospital.Research.ID == pt) %>% select(Samplename)
    # Predict function giving me difficulty when I have only a single sample, this ensures the dimensions are the same
    sparsed_test_data <- Matrix(data=0, nrow=ifelse(nrow(patient.samples) > 1, nrow(test), 1),  ncol=ncol(training),
                                dimnames=list(rownames(test),colnames(training)), sparse=T)

    for(i in colnames(dysplasia.dfNOHGD)) sparsed_test_data[,i] = test[,i]
    
    # Fit generated on all samples, including HGD
    fitLOO <- glmnet(training, labels[train.rows], alpha=a, family='binomial', nlambda=nl, standardize=F) # all patients
    l = fitLOO$lambda
    
    cv = crossvalidate.by.patient(x=training, y=labels[train.rows], lambda=l, a=a, nfolds=folds, splits=splits,
                                  pts=subset(sets, Samplename %in% samples), fit=fitLOO, standardize=F)
    
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
      or = predict(fitLOO, newx=sparsed_test_data, s=cv$lambda.1se, type='link')
      sy = as.matrix(sqrt(binomial.deviance(pm,labels[intersect(pg.sampNOHGD[[pt]]$Samplename, names(labels))])))

      pm = predict(fitLOO, newx=sparsed_test_data, s=cv$lambda.1se, type='response')
      colnames(pm) = 'Prediction'
      
      rr = predict(fitLOO, newx=sparsed_test_data, s=cv$lambda.1se, type='link')
      colnames(rr) = 'RR'
      
      sy = as.matrix(sqrt(binomial.deviance(pm,labels[intersect(patient.samples$Samplename, names(labels))])))
      colnames(sy) = 'Prediction.Dev.Resid'
      
      pg.sampNOHGD = pg.sampNOHGD %>% 
        full_join(as_tibble(pm, rownames='Samplename'), by='Samplename') %>% 
        full_join(as_tibble(rr, rownames='Samplename'), by='Samplename') %>% 
        full_join(as_tibble(sy, rownames='Samplename'), by='Samplename')
      
    } else {
      warning(paste("Hospital.Research.ID", pt, "did not have a 1se"))
    }
  }
  save(plots, performance.at.1se, coefs, nzcoefs, fits, pg.sampNOHGD, file=file)
}

message("Finished")

