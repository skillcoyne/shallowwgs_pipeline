library(ggplot2)
library(gridExtra)
library(glmnet)
library(readxl)
library(pROC)
library(ggfortify)
library(preprocessCore)

suppressPackageStartupMessages(source('~/workspace/shallowwgs_pipeline/lib/data_func.R'))
suppressPackageStartupMessages(source('~/workspace/shallowwgs_pipeline/lib/common_plots.R'))
suppressPackageStartupMessages(source('~/workspace/shallowwgs_pipeline/lib/cv-pt-glm.R'))

adjust.cols<-function(mtx, means=NULL, sds=NULL, na.replace=0) {
  if (!is.null(means)) {
    if (ncol(mtx) != length(means))
      stop("Vector of means needs to match the columns in the matrix")
    for (i in 1:ncol(mtx)) 
      mtx[,i] = unit.var(mtx[,i], means[i], sds[i])
  } else {
    for (i in 1:ncol(mtx)) 
      mtx[,i] = unit.var(mtx[,i])
  }
  mtx[is.na(mtx)] = na.replace
  
  return(mtx)
}

chr.info = get.chr.lengths(file='/tmp/hg19_info.txt')

load('~/Data/Ellie/Analysis/5e6_arms_all/model_data.Rdata', verbose = T)
rm(labels, dysplasia.df,mn.cx,sd.cx,z.mean,z.sd,z.arms.mean,z.arms.sd)

load('~/Data/Reid_SNP/PerPatient/allpts_ascat.Rdata', verbose=T)
qcdata = qcdata %>% rowwise %>% dplyr::mutate(
  PatientID = unlist(strsplit(Samplename, '_'))[1],
  SampleID = unlist(strsplit(Samplename, '_'))[2],
  EndoID = unlist(strsplit(Samplename, '_'))[3],
  Level = unlist(strsplit(Samplename, '_'))[4]
)
segments.list = lapply(segments.list, function(pt) {
  pt$sample = sub('\\.LogR','',pt$sample)
  pt
})

qcdata$ASCAT.SCA.ratio = apply(qcdata,1,function(s) {
  smp = subset(segments.list[[ s[['PatientID']] ]], sample == s[['Samplename']] & chr %in% c(1:22))
  smp = smp %>% rowwise %>% dplyr::mutate(
    #'Total' = nAraw + nBraw,
    'Total' = nMajor + nMinor,
    'CNV' = round(Total) - round(as.numeric(s[['Ploidy']])) )
  
  x = subset(smp, CNV != 0 & chr %in% c(1:22))
  sum(x$endpos - x$startpos) / chr.info[22,'genome.length']
})

qcdata$SampleType = 'BE'
qcdata$SampleType[grep('BLD',qcdata$Level, ignore.case=T)] = 'Blood Normal'
qcdata$SampleType[grep('gastric',qcdata$Level, ignore.case=T)] = 'Gastric Normal'

patient.info = as.data.frame(read_xlsx('~/Data/Ellie/Analysis/SNP/metadata_T1T2.xlsx'))
patient.info$UniqueSampleID = paste(patient.info$PatientID, patient.info$`Timepoint Code`, sep='_')
patient.info$Path.Status = patient.info$Status
patient.info[patient.info$PatientID %in% subset(patient.info, Pathology %in% c('IMC','HGD'), select='PatientID')[,1], 'Path.Status'] = 'P'

sample.list = read_xlsx('~/Data/Ellie/Analysis/SNP/20180604_Reid_1M_SampleList.xlsx', sheet = 2)
nrow(sample.list)
colnames(sample.list)[3] = 'Total.SCA'

sample.list$SCA.Ratio = sample.list$Total.SCA/max(chr.info$genome.length)

sample.info = qcdata[,c('PatientID','SampleID','EndoID','Level','Samplename','ASCAT.SCA.ratio')]
sample.info = base::merge(sample.info, patient.info[,c('PatientID','UniqueSampleID','Timepoint','Timepoint Code','Status','Pathology')], by.x=c('PatientID','EndoID'), by.y=c('PatientID','Timepoint Code'))
head(sample.info)

message(paste(length(unique(sample.info$PatientID)), 'patients listed in metadata file'))
message(paste(nrow(qcdata), 'samples available'))

get.ratio<-function(PID, SID) { subset(sample.list, PatientID == PID & SampleNum == SID)$SCA.Ratio }
sample.info = sample.info %>% rowwise() %>% dplyr::mutate(
  SCA.Ratio = ifelse(Level %in% c('BLD','GASTRIC'), 0, get.ratio(PatientID, SampleID))
)

patient.info = patient.info %>% rowwise() %>% dplyr::mutate( Endoscopy = ifelse(Timepoint == 'T1', 1, 2) )

qcdata = base::merge(qcdata, sample.info[,c('Samplename','Status','UniqueSampleID','Timepoint','SCA.Ratio','Pathology')], by='Samplename')

## PER ENDOSCOPY
qcdata.samples = qcdata %>% group_by(PatientID, EndoID, UniqueSampleID, Status, Timepoint, SampleType, Pathology) %>% dplyr::summarise_if(is.numeric, c('mean','max','min','sd')) %>% mutate_if(is.numeric, round, digits=3)

info = qcdata.samples %>% group_by(Status, PatientID, UniqueSampleID) %>% dplyr::mutate(
  exclude=(Status == 'NP' & (length(which(Pathology == 'HGD')) > 0)),
  wgd = (Status == 'NP' & Ploidy_max > 2.7 ),
  lowsca = (ASCAT.SCA.ratio_mean < 0.02 & Purity_mean > 0.95)
)

table(unique(info[,c('PatientID','Status')])$Status)
table(unique(subset(info, !exclude & !wgd & !lowsca, select=c('PatientID','Status')))$Status)

#info = subset(info, !exclude & !wgd & !lowsca & SampleType == 'BE')
#info = subset(info, !exclude & !wgd & SampleType == 'BE')
info = subset(info, SampleType == 'BE')
nrow(info)

if (file.exists('~/Data/Reid_SNP/PerPatient/tmp_seg_pt.Rdata')) {
  load('~/Data/Reid_SNP/PerPatient/tmp_seg_pt.Rdata', verbose=T) 
} else {

  mergedSegs = NULL
  mergedArms = NULL
  
  length(ptdirs)
  
  for (pt in ptdirs) {
    print(pt)
    if (length(list.files(pt, '*wins_tiled.txt', full.names=T)) <= 0) {
      message(paste("No tiled files for",pt))
      next
    }
    
    segvals = as.data.frame(data.table::fread(list.files(pt, '*wins_tiled.txt', full.names=T)))
    armvals = as.data.frame(data.table::fread(list.files(pt, '*arms_tiled.txt', full.names=T)))
    
    segvals = segment.matrix(segvals)
    segvals[is.na(segvals)] = mean(segvals,na.rm=T)
    
    armvals = segment.matrix(armvals)
    armvals[is.na(armvals)] = mean(armvals,na.rm=T)
    
    if (is.null(segvals) | is.null(armvals))
      stop(paste("Missing values in ", pt))
    
    if (is.null(mergedSegs)) {
      mergedSegs = segvals
      mergedArms = armvals
    } else {
      mergedSegs = rbind(mergedSegs, segvals)    
      mergedArms = rbind(mergedArms, armvals)    
    }
  }
  nrow(mergedSegs) == nrow(mergedArms)
  #setdiff(rownames(mergedSegs), rownames(mergedArms))
  dim(mergedSegs)
  
  save(mergedSegs, mergedArms, file='tmp_seg_pt.Rdata')
}

dim(mergedSegs)
dim(mergedArms)

## Select samples
mergedSegs = mergedSegs[unique(info$UniqueSampleID),]
mergedArms = mergedArms[unique(info$UniqueSampleID),]

dim(mergedSegs)
dim(mergedArms)


# Adjust for ploidy
for (i in 1:nrow(mergedSegs)) {
  sample = rownames(mergedSegs)[i]
  x = unlist(strsplit(sample, '_'))
  ploidy = round(subset(info, UniqueSampleID == sample)$Ploidy_mean)
  if (grepl('BLD|gastric', sample, ignore.case = T)) ploidy = 2
  mergedSegs[i,] = mergedSegs[i,]-(ploidy-1)
  mergedArms[i,] = mergedArms[i,]-(ploidy-1)
}
segmentVariance = apply(mergedSegs, 1, var)

# Q-norm
q.norm.by<-function(df, target) {
  dfq = do.call(cbind.data.frame,lapply(1:ncol(df), function(i) {
    as.data.frame(normalize.quantiles.use.target(df[,i,drop=F], as.vector(target[,i])))
  }))
  
  colnames(dfq) = colnames(df)
  rownames(dfq) = rownames(df)
  
  return(as.matrix(dfq))
}

copySegs = q.norm.by(mergedSegs, raw.segs)
copyArms = q.norm.by(mergedArms, raw.arms)

#copySegs = mergedSegs
#copyArms = mergedArms

# grid.arrange(
#   ggplot( melt(raw.segs), aes(sample=value)) + stat_qq() + labs(title='Normal Q-Q plot, sWGS bins (no arms)'),
#   ggplot( melt(mergedSegs), aes(sample=value)) + stat_qq() + labs(title='Normal Q-Q plot, SNP bins (no arms)'),
#   ggplot( melt(copySegs), aes(sample=value)) + stat_qq() + labs(title='Normal Q-Q plot, SNP bins q-normed')
# )

ms = adjust.cols(copySegs)
range(ms)
dim(ms)

ma = adjust.cols(mergedArms)
range(ma)
dim(ma)

cx = score.cx(ms, 1)
arrayDf = subtract.arms(ms, ma)
arrayDf = cbind(arrayDf, 'cx'=unit.var(cx))


slabels = as.data.frame(subset(info, !exclude & !wgd & !lowsca & SampleType == 'BE', select=c('PatientID','UniqueSampleID','Status')))
slabels$Status = factor(slabels$Status)
#slabels = as.data.frame(pred.qc, select=c('PatientID','Samplename','Status'))
rownames(slabels) = slabels$UniqueSampleID

df = arrayDf[which(rownames(arrayDf) %in% slabels$UniqueSampleID),]

status = as.integer(slabels$Status)-1
names(status) = slabels$UniqueSampleID

nl = 1000;folds = 10; splits = 5 
sets = create.patient.sets(slabels, folds, splits, 0.15)  
## ----- All ----- ##
coefs = list(); plots = list(); performance.at.1se = list(); models = list(); cvs = list()
cache.dir = '~/Data/Ellie/Analysis/SNP/no-low'
dir.create(cache.dir, recursive = T, showWarnings = F)
file = paste(cache.dir, 'all.pt.alpha.Rdata', sep='/')
if (file.exists(file)) {
 message(paste("loading", file))
 load(file, verbose=T)
} else {
  alpha.values = c(0,0.5,0.7,0.9,1)

  for (a in alpha.values) {
    message(paste('alpha=',a))
    fit0 <- glmnet(df, status, alpha=a, nlambda=nl, family='binomial', standardize=F)    
    autoplot(fit0) + theme(legend.position="none")
    l = fit0$lambda
    if (a > 0) l = more.l(l)
    
    cv.patient = crossvalidate.by.patient(x=df, y=status, lambda=l, pts=sets, a=a, nfolds=folds, splits=splits, fit=fit0, select='deviance', opt=-1, standardize=F)
    
    lambda.opt = cv.patient$lambda.1se
    
    coef.opt = as.data.frame(non.zero.coef(fit0, lambda.opt))
    coefs[[as.character(a)]] = coef.stability(coef.opt, cv.patient$non.zero.cf)
    
    plots[[as.character(a)]] = arrangeGrob(cv.patient$plot+ggtitle('Classification'), cv.patient$deviance.plot+ggtitle('Binomial Deviance'), top=paste('alpha=',a,sep=''), ncol=2)
    
    performance.at.1se[[as.character(a)]] = subset(cv.patient$lambdas, `lambda-1se`)
    models[[as.character(a)]] = fit0
    cvs[[as.character(a)]] = cv.patient
  }
  save(plots, coefs, performance.at.1se, df, models, cvs, slabels, status, info, file=file)
  p = do.call(grid.arrange, c(plots[ as.character(alpha.values) ], top='All samples, 10fold, 5 splits'))
  ggsave(paste(cache.dir, '/', 'all_samples_cv.png',sep=''), plot = p, scale = 2, width = 12, height = 10, units = "in", dpi = 300)
}

#do.call(grid.arrange, c(plots))

#model.performance(performance.at.1se, coefs)$plot + labs(title='SNP CV models')


## Leave one out

slabels$Prediction = NA
slabels$RR = NA

colnames(sets)[2] = 'Samplename'

file = paste(cache.dir, 'loo.Rdata', sep='/')
# if (file.exists(file)) {
#   load(file, verbose=T)
# } else {
  select.alpha = 1
  secf = coefs[[select.alpha]]
  a = select.alpha
  
  performance.at.1se = c(); plots = list(); fits = list(); nzcoefs = list()
  # Remove each patient (LOO)
  for (pt in unique(slabels$PatientID)) {
    print(pt)
    samples = subset(slabels, PatientID != pt)[,2]
    
    train.rows = which(rownames(df) %in% samples)
    training = df[train.rows,]
    test = as.matrix(df[-train.rows,])
    if ( nrow(test) == ncol(df) ) test = t(test)
    
    # Predict function giving me difficulty when I have only a single sample, this ensures the dimensions are the same
    sparsed_test_data <- Matrix(data=0, nrow=nrow(test),  ncol=ncol(training),
                                dimnames=list(rownames(test),colnames(training)), sparse=T)
    for(i in colnames(arrayDf)) sparsed_test_data[,i] = test[,i]
    
    # Fit generated on all samples, including HGD
    fitLOO <- glmnet(training, status[train.rows], alpha=a, family='binomial', nlambda=nl, standardize=F) # all patients
    l = fitLOO$lambda
    
    cv = crossvalidate.by.patient(x=training, y=status[train.rows], lambda=l, a=a, nfolds=folds, splits=splits,
                                  pts=subset(sets, Samplename %in% samples), fit=fitLOO, standardize=F)
    
    plots[[pt]] = arrangeGrob(cv$plot+ggtitle('Classification'), cv$deviance.plot+ggtitle('Binomial Deviance'), top=pt, ncol=2)
    fits[[pt]] = cv  
    
    if ( length(cv$lambda.1se) > 0 ) {
      performance.at.1se = c(performance.at.1se, subset(cv$lambdas, lambda == cv$lambda.1se)$mean)
      
      coef.1se = coef(fitLOO, cv$lambda.1se)[rownames(secf),]
      
      nzcoefs[[pt]] = as.data.frame(non.zero.coef(fitLOO, cv$lambda.1se))
      
      #coefs[[pt]] = coef(fitLOO, cv$lambda.1se)[rownames(secf),]
      #coef.stability(coef.1se, cv$non.zero.cf)
      
      logit <- function(p){log(p/(1-p))}
      inverse.logit <- function(or){1/(1 + exp(-or))}
      
      pm = predict(fitLOO, newx=sparsed_test_data, s=cv$lambda.1se, type='response')
      or = predict(fitLOO, newx=sparsed_test_data, s=cv$lambda.1se, type='link')
      
      slabels[which(slabels$UniqueSampleID %in% rownames(pm)),'Prediction'] = pm[,1]
      slabels[which(slabels$UniqueSampleID %in% rownames(pm)),'RR'] = or[,1]
      
    } else {
      warning(paste(pt, "did not have a 1se"))
    }
  }
  save(plots, performance.at.1se, nzcoefs, fits, slabels, file=file)
#}
message('Finished')

