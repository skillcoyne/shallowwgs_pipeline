args = commandArgs(trailingOnly=TRUE)

if (length(args) < 2)
  stop("Missing required params:<outdir> <info file dir>")

suppressPackageStartupMessages( library(BarrettsProgressionRisk) )
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(glmnet))


suppressPackageStartupMessages(source('~/workspace/shallowwgs_pipeline/lib/load_patient_metadata.R'))
suppressPackageStartupMessages(source('~/workspace/shallowwgs_pipeline/lib/cv-pt-glm.R'))
suppressPackageStartupMessages(source('~/workspace/shallowwgs_pipeline/lib/data_func.R'))


outdir = args[1]
# outdir = '~/Data/Ellie/Analysis'
infodir = args[2]
# infodir = '~/Data/Ellie/QDNAseq/training'
logT = F
if (length(args) == 4)
  logT = as.logical(args[4])

x = (unclass(Sys.time()) + sample(1:1000, 1))
cache.dir = paste(outdir, '5e6_arms_splits', sep='/')
if (logT) cache.dir = paste(cache.dir, '_logR', sep='')
cache.dir = paste(cache.dir, as.character(x) , sep='/')

dir.create(cache.dir, showWarnings = F, recursive = T)

## Hospital.Research.ID info file
patient.file = list.files(infodir, pattern='All_patient_info.xlsx', recursive=T, full.names=T)
demo.file = list.files(infodir, pattern='Demographics_full.xlsx', recursive=T, full.names=T)

if ( length(patient.file) < 1 | length(demo.file) < 1)
  stop("Missing files in info dir: All_patient_info.xlsx and Demographics_full.xlsx")

all.patient.info = read.patient.info(patient.file, demo.file, set='All')$info %>% arrange(Status, Patient, Endoscopy.Year, Pathology)
sum.patient.data = summarise.patient.info(all.patient.info)


modeldir = paste0(outdir, '/5e6_arms_all')
if (!dir.exists(modeldir))
  stop(paste0("Missing trained model in ", modeldir))
load(paste0(modeldir, '/all.pt.alpha.Rdata'),verbose=T)
rm(plots,coefs,performance.at.1se, cvs, models)

# By patient
leaveout = sample(sum.patient.data$Patient, round(nrow(sum.patient.data)*.2))
leaveoutSamples = patient.info %>% filter(Patient %in% leaveout & Samplename %in% rownames(dysplasia.df)) 

patient.info = all.patient.info %>% filter(!Patient %in% leaveout)
sum.patient.data = summarise.patient.info(patient.info)
nrow(sum.patient.data)

leaveoutDf = dysplasia.df[ intersect(leaveoutSamples$Samplename,rownames(dysplasia.df)), ]
dysplasia.df = dysplasia.df[intersect(patient.info$Samplename,rownames(dysplasia.df)), ]

## labels: binomial: prog 1, np 0
labels = labels[rownames(dysplasia.df)]
dim(dysplasia.df)

save(dysplasia.df, labels, file=paste(cache.dir, 'model_data.Rdata', sep='/'))


nl = 1000; folds = 10; splits = 5 
sets = create.patient.sets(patient.info[c('Hospital.Research.ID','Samplename','Status')], folds, splits, 0.2) 
#alpha.values = c(0, 0.5,0.7,0.8,0.9,1)
alpha.values = c(0.9)

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
    #l = more.l(l)
    
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
  ggsave(paste(cache.dir, '/', 'all_samples_cv.png',sep=''), plot = p, scale = 2, width = 15, height = 10, units = "in", dpi = 300)
}
all.coefs = coefs
# ----------------- #

# Predict leave out samples
select.alpha = 0.9

#file = paste(cache.dir, 'all.pt.alpha.Rdata', sep='/')
if (!file.exists(file))
  stop(paste("Missing data file", file))
load(file, verbose=T)

fitV = models[[as.character(select.alpha)]]
lambda.opt = performance.at.1se[[as.character(select.alpha)]][, 'lambda']

pred = predict(fitV, newx=leaveoutDf, s=lambda.opt, type='response')
colnames(pred) = 'Prediction'

rr = predict(fitV, newx=leaveoutDf, s=lambda.opt, type='link')
colnames(rr) = 'RR'

vpd = all.patient.info %>% dplyr::filter(Samplename %in% rownames(leaveoutDf)) %>%
  dplyr::mutate(  PID = sub('_$', '', unlist(strsplit(Path.ID, 'B'))[1])) %>% 
  full_join(as_tibble(pred, rownames='Samplename'), by='Samplename') %>% 
  full_join(as_tibble(rr, rownames='Samplename'), by='Samplename')

save(vpd, file=paste(cache.dir, 'leaveout_predictions.Rdata', sep='/'))

## --------- LOO --------- ##
message("LOO")

pg.samp = all.patient.info %>% rowwise %>% dplyr::mutate(
  PID = sub('_$', '', unlist(strsplit(Path.ID, 'B'))[1])
) %>% filter(Samplename %in% rownames(dysplasia.df))


file = paste(cache.dir, 'loo.Rdata', sep='/')
if (file.exists(file)) {
  load(file, verbose=T)
} else {
  secf = all.coefs[[select.alpha]]
  a = select.alpha
  
  performance.at.1se = c(); coefs = list(); plots = list(); fits = list(); nzcoefs = list()
  # Remove each patient (LOO)
  for (pt in unique(pg.samp$Hospital.Research.ID)) {
    print(pt)
    
    samples = all.patient.info %>% filter(Hospital.Research.ID != pt & Samplename %in% rownames(dysplasia.df)) %>% select(Samplename)
    
    train.rows = which(rownames(dysplasia.df) %in% samples$Samplename)
    training = dysplasia.df[train.rows,,drop=F]
    test = as.matrix(dysplasia.df[-train.rows,,drop=F])

    patient.samples = all.patient.info %>% filter(Hospital.Research.ID == pt & Samplename %in% rownames(dysplasia.df)) %>% select(Samplename)
    # Predict function giving me difficulty when I have only a single sample, this ensures the dimensions are the same
    sparsed_test_data <- Matrix(data=0, nrow=ifelse(nrow(patient.samples) > 1, nrow(test), 1),  ncol=ncol(training),
                                dimnames=list(rownames(test),colnames(training)), sparse=T)
    for(i in colnames(dysplasia.df)) sparsed_test_data[,i] = test[,i]
    
    # Fit generated on all samples, including HGD
    fitLOO <- glmnet(training, labels[train.rows], alpha=a, family='binomial', nlambda=nl, standardize=F) # all patients
    l = fitLOO$lambda
    
    cv = crossvalidate.by.patient(x=training, y=labels[train.rows], lambda=l, a=a, nfolds=folds, splits=splits,
                                  pts=subset(sets, Samplename %in% samples$Samplename), fit=fitLOO, standardize=F)
    
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




message("Finished")

