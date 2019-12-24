args = commandArgs(trailingOnly=TRUE)

set.seed(1234)

if (length(args) < 3)
  stop("Missing required params: <data dir> <outdir> <info file dir> <Set: All,Training>")


suppressPackageStartupMessages( library(BarrettsProgressionRisk) )
suppressPackageStartupMessages( library(tidyverse) )
suppressPackageStartupMessages( library(glmnet) )
suppressPackageStartupMessages( library(gridExtra) )

suppressPackageStartupMessages(source('~/workspace/shallowwgs_pipeline/lib/data_func.R'))
suppressPackageStartupMessages(source('~/workspace/shallowwgs_pipeline/lib/load_patient_metadata.R'))
suppressPackageStartupMessages(source('~/workspace/shallowwgs_pipeline/lib/cv-pt-glm.R'))


data = args[1]
# data = '~/Data/BarrettsProgressionRisk/Analysis/pcf_perPatient/50kb/'
outdir = args[2]
# outdir = '~/Data/BarrettsProgressionRisk/Analysis/cox_5e6_all/50kb/'
infodir = args[3]
# infodir = '~/Data/BarrettsProgressionRisk/QDNAseq'
set = 'All'
#set = 'Training'
#if (length(args) == 4) set = args[4]

select.alpha = 0.9
 if (length(args) == 4)
   select.alpha = as.numeric(args[4])

cache.dir = outdir

dir.create(cache.dir, recursive=T, showWarnings=F)

## Hospital.Research.ID info file
patient.file = list.files(infodir, pattern='All_patient_info.xlsx', recursive=T, full.names=T)
demo.file = list.files(infodir, pattern='Demographics_full.xlsx', recursive=T, full.names=T)

if ( length(patient.file) != 1 | length(demo.file) != 1)
  stop("Missing files in info dir: All_patient_info.xlsx and Demographics_full.xlsx")

all.patient.info = read.patient.info(patient.file, demo.file, set='All')$info %>% dplyr::arrange(Status, Patient, Endoscopy.Year, Pathology)

if (set != 'All') {
  patient.info = all.patient.info %>% filter(Set == set)
} else {
  patient.info = all.patient.info
}
length(unique(patient.info$Hospital.Research.ID))

sum.patient.data = as_tibble(summarise.patient.info(patient.info))

cleaned = list.files(path=data, pattern='tiled_segvals', full.names=T, recursive=T)
cleaned = grep(paste(sum.patient.data$Hospital.Research.ID, collapse = '|'), cleaned, value=T)

arm.files = c(grep('arm', cleaned, value=T))
seg.files = c(grep('arm', cleaned, value=T, invert=T))

if (length(seg.files) != length(arm.files)) stop("Tiled files do not match between short segments and arms.")

seg.tiles = do.call(bind_rows, lapply(seg.files, function(x) readr::read_tsv(x, col_types=cols(.default=col_double(), X1=col_character()))))
colnames(seg.tiles)[1] = 'sample'
samples = seg.tiles$sample

seg.tiles = as.matrix(seg.tiles[,-1])
rownames(seg.tiles) = samples
dim(seg.tiles)

per.samp.sd = apply(seg.tiles, 1, sd)
per.samp.mean = apply(seg.tiles,1,mean)

# After mean centering set all NA values to 0
segsList = prep.matrix(seg.tiles,scale=T, MARGIN=2)
#segsList = prep.matrix(seg.tiles,scale=F)
z.mean = segsList$z.mean
z.sd = segsList$z.sd
segs = segsList$matrix
dim(segs)

# Complexity score
cx.score = BarrettsProgressionRisk::scoreCX(segs,1)
mn.cx = mean(cx.score)
sd.cx = sd(cx.score)

# Load arm files  
arm.tiles = do.call(bind_rows, lapply(arm.files, function(x) readr::read_tsv(x, col_types=cols(.default=col_double(), X1=col_character())) %>% dplyr::filter(X1!='')))
colnames(arm.tiles)[1] = 'sample'
dim(arm.tiles)
samples = arm.tiles$sample
arm.tiles = as.matrix(arm.tiles[,-1])
rownames(arm.tiles) = samples

armsList = prep.matrix(arm.tiles,scale=T,MARGIN=2)
arms = armsList$matrix
z.arms.mean = armsList$z.mean
z.arms.sd = armsList$z.sd

# Merge segments and arms, subtract arms
allDf = BarrettsProgressionRisk::subtractArms(segs, arms)
mn.cx = sqrt(mean(cx.score^2))

# Add complexity score
allDf = cbind(allDf, 'cx'=BarrettsProgressionRisk:::unit.var(cx.score, mn.cx, sd.cx))

## labels: binomial: prog 1, np 0
sampleStatus = patient.info %>% filter(Samplename %in% rownames(allDf)) %>% dplyr::select(Samplename,Status) %>% 
  mutate(Status = ifelse(Status == 'OAC', 'P', Status)) %>% mutate(Status = factor(Status)) 

sampleStatus$Status = as.integer(sampleStatus$Status)-1
labels = sampleStatus$Status
names(labels) = sampleStatus$Samplename

dim(allDf)
dysplasia.df = as.matrix(allDf[sampleStatus$Samplename,])
dim(dysplasia.df)

save(dysplasia.df, labels, mn.cx, sd.cx, z.mean, z.sd, z.arms.mean, z.arms.sd, file=paste(cache.dir, 'model_data.Rdata', sep='/'))


nl = 1000;folds = 10; splits = 5 
sets = create.patient.sets(patient.info[c('Hospital.Research.ID','Samplename','Status')], folds, splits, 0.2)  

alpha.values = c(0, 0.5,0.7,0.8,0.9,1)
## ----- All ----- ##
coefs = list(); plots = list(); performance.at.1se = list(); models = list(); cvs = list()
file = paste(cache.dir, 'all.pt.alpha.Rdata', sep='/')
if (file.exists(file)) {
  message(paste("loading", file))
  load(file, verbose=T)
} else {
  for (a in alpha.values) {
    print(a)
    te = dplyr::select(patient.info, Samplename, Endoscopy.Year, Final.Endoscopy, Status) %>% group_by(Samplename, Status) %>% 
      dplyr::mutate(time = Final.Endoscopy+1 - Endoscopy.Year, status = as.integer(Status)-1) %>% filter( Samplename %in% rownames(dysplasia.df))
    
    dysplasia.df = dysplasia.df[te$Samplename,]
    
    fit0 <- glmnet(dysplasia.df, Surv(te$time, te$status), alpha = a, nlambda = nl, family = 'cox', standardize = F)
    #autoplot(fit0) + theme(legend.position = 'none')

    cv = cv.glmnet(dysplasia.df, Surv(te$time, te$status), type.measure = 'deviance', nfolds = folds, family = 'cox')
    #plot(cv)
    
    lambda.opt = cv$lambda.1se
    
    coef.opt = as.data.frame(non.zero.coef(fit0, lambda.opt))
    coefs[[as.character(a)]] = coef.opt
    
    plots[[as.character(a)]] = autoplot(cv) + labs(title='Cox', subtitle=paste0('alpha = ',a))
    #plots[[as.character(a)]] = arrangeGrob(cv.patient$plot+ggtitle('Classification'), cv.patient$deviance.plot+ggtitle('Binomial Deviance'), top=paste('alpha=',a,sep=''), ncol=2)
    
    #performance.at.1se[[as.character(a)]] = subset(cv.patient$lambdas, `lambda-1se`)
    models[[as.character(a)]] = fit0
    cvs[[as.character(a)]] = cv
  }
  save(plots, coefs, performance.at.1se, dysplasia.df, models, cvs, labels, file=file)
  p = do.call(grid.arrange, c(plots[ as.character(alpha.values) ], top='All samples, 10fold, 5 splits'))
  ggsave(paste(cache.dir, '/', 'all_samples_cv.png',sep=''), plot = p, scale = 2, width = 12, height = 10, units = "in", dpi = 300)
}
all.coefs = coefs
# ----------------- #


## --------- LOO --------- ##
message("LOO")

pg.samp = patient.info %>% rowwise %>% dplyr::mutate(
  PID = sub('_$', '', unlist(strsplit(Path.ID, 'B'))[1]),
  FittedRR = NA,
  LinearPred = NA,
) %>% filter(Samplename %in% rownames(dysplasia.df))

file = paste(cache.dir, paste0('loo_',select.alpha,'.Rdata'), sep='/')
if (file.exists(file)) {
  message(paste("loading file", file))
  load(file, verbose=T)
} else {
  secf = all.coefs[[select.alpha]]
  a = select.alpha
  
  te = dplyr::select(patient.info, Samplename, Endoscopy.Year, Final.Endoscopy, Status) %>% group_by(Samplename, Status) %>% 
    dplyr::mutate(time = Final.Endoscopy+1 - Endoscopy.Year, status = as.integer(Status)-1) %>% filter( Samplename %in% rownames(dysplasia.df))
  
  dysplasia.df = dysplasia.df[te$Samplename,]
  
  performance.at.1se = c(); coefs = list(); plots = list(); fits = list(); nzcoefs = list()
  # Remove each patient (LOO)
  for (pt in unique(pg.samp$Hospital.Research.ID)) {
    # pt.path = paste0(cache.dir, '/plots/', pt)
    # dir.create(pt.path, recursive = T, showWarnings = F)
    
  #for (pt in names(pg.samp)) {
    print(pt)
    samples = subset(pg.samp, Hospital.Research.ID != pt)$Samplename
    
    train.rows = which(rownames(dysplasia.df) %in% samples)
    training = dysplasia.df[train.rows,,drop=F]
    test = as.matrix(dysplasia.df[-train.rows,,drop=F])
    
    yt = filter(te, Samplename %in% rownames(training)) 
    
    if (ncol(test) <= 0) next # shouldn't be any but...

    patient.samples = pg.samp %>% filter(Hospital.Research.ID == pt) %>% dplyr::select(Samplename)
    # Predict function giving me difficulty when I have only a single sample, this ensures the dimensions are the same
    sparsed_test_data <- Matrix(data=0, nrow=ifelse(nrow(patient.samples) > 1, nrow(test), 1),  ncol=ncol(training),
                                dimnames=list(rownames(test),colnames(training)), sparse=T)
    
    for(i in colnames(dysplasia.df)) sparsed_test_data[,i] = test[,i]
    
    # Fit generated on all samples, including HGD
    fitLOO <- glmnet(training, Surv(yt$time, yt$status), alpha=a, family='cox', nlambda=nl, standardize=F) # all patients
    l = fitLOO$lambda
    
    cv = cv.glmnet(training, Surv(yt$time, yt$status), type.measure = 'deviance', nfolds = folds, family = 'cox')

    lambda.opt = cv$lambda.1se
    
    plots[[pt]] = autoplot(cv) + labs(title='Cox', subtitle=paste0('alpha = ',a))

    fits[[pt]] = cv

    if ( length(cv$lambda.1se) > 0 ) {
      #performance.at.1se = c(performance.at.1se, subset(cv$lambdas, lambda == cv$lambda.1se)$mean)
      
      nzcoefs[[pt]] = as.data.frame(non.zero.coef(fitLOO, cv$lambda.1se))
      coef.opt = as.data.frame(non.zero.coef(fitLOO, lambda.opt))
      nzcoefs[[pt]] = coef.opt

      coefs[[pt]] = coef(fitLOO, cv$lambda.1se)[rownames(secf),]

      pm = predict(fitLOO, newx=sparsed_test_data, s=cv$lambda.1se, type='response')
      colnames(pm) = 'FittedRR'
      
      rr = predict(fitLOO, newx=sparsed_test_data, s=cv$lambda.1se, type='link')
      colnames(rr) = 'LinearPred'

      patient = pg.samp %>% filter(Hospital.Research.ID == pt) %>% arrange(Samplename) 
      
      patient$FittedRR = pm[patient$Samplename,]
      patient$LinearPred = rr[patient$Samplename,]

      pg.samp[which(pg.samp$Hospital.Research.ID == pt),] = patient

    } else {
      warning(paste("Hospital.Research.ID", pt, "did not have a 1se"))
    }
  }
  save(plots, performance.at.1se, coefs, nzcoefs, fits, pg.samp, file=file)
}


message("Finished")

