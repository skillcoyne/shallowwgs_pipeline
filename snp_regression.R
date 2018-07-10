library(ggplot2)
library(gridExtra)
library(glmnet)
library(pander)
library(readxl)
library(pROC)
library(ggfortify)

source('lib/cv-pt-glm.R')
source('lib/load_patient_metadata.R')
source('lib/cv-pt-glm.R')
source('lib/data_func.R')
source('lib/common_plots.R')

chr.info = get.chr.lengths(file='hg19_info.txt')

ptdirs = list.dirs('/Volumes/fh/fast/reid_b/collab/Killcoyne/Data/PerPatient', full.names=T, recursive=F)
patients = basename(ptdirs)

load('~/Data/Ellie/Analysis/5e6_arms_all_logR/model_data.Rdata', verbose = T)
load('~/Data/Ellie/Analysis/5e6_arms_all_logR/all.pt.alpha.Rdata', verbose = T)
fitV = models$`0.9`
l = performance.at.1se$`0.9`$lambda

swgs_labels = labels

rm(plots,performance.at.1se,models,cvs,labels)

pv = var(dysplasia.df[swgs_labels == 1])  # P
npv = var(dysplasia.df[swgs_labels == 0]) # NP

range(dysplasia.df[swgs_labels == 1])
range(dysplasia.df[swgs_labels == 0])


#load('/Volumes/fh/fast/reid_b/collab/Killcoyne/SNP_R/allpts_ascat.Rdata', verbose=T)
load('~/Data/Ellie/Analysis/SNP/allpts_ascat.Rdata', verbose=T)
#patient.info = as.data.frame(read_xlsx('/Volumes/fh/fast/reid_b/collab/Killcoyne/SNP_Project/metadata_T1T2.xlsx'))
patient.info = as.data.frame(read_xlsx('~/Data/Ellie/Analysis/SNP/metadata_T1T2.xlsx'))
patient.info$UniqueSampleID = paste(patient.info$PatientID, patient.info$`Timepoint Code`, sep='_')
patient.info$Path.Status = patient.info$Status
patient.info[patient.info$PatientID %in% subset(patient.info, Pathology %in% c('IMC','HGD'), select='PatientID')[,1], 'Path.Status'] = 'P'



#sample.list = read_xlsx('/Volumes/fh/fast/reid_b/collab/Killcoyne/Data/20180604_Reid_1M_SampleList.xlsx', sheet = 2)
#nrow(sample.list)
#colnames(sample.list)[3] = 'Total.SCA'

#sample.list$SCA.Ratio = sample.list$Total.SCA/max(chr.info$genome.length)

sample.info  = do.call(rbind.data.frame, sapply(qcdata$Samplename, strsplit, '_'))
colnames(sample.info) = c('PatientID','SampleID','EndoID','Level')
sample.info[] = lapply(sample.info[], as.character)
sample.info[which(sample.info$PatientID == 524 & sample.info$Level == 524), 'Level'] = ''
sample.info$Samplename = qcdata$Samplename

message(paste(length(unique(sample.info$PatientID)), 'patients listed in metadata file'))
message(paste(nrow(qcdata), 'samples available'))

sample.info = base::merge(sample.info, patient.info[,c('PatientID','Timepoint','Timepoint Code','Status','Pathology')], by.x=c('PatientID','EndoID'), by.y=c('PatientID','Timepoint Code'))


if (file.exists('tmp_seg_pt.Rdata')) {
  load('tmp_seg_pt.Rdata', verbose=T) 
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

rownames(mergedSegs) = sub('\\.LogR','', rownames(mergedSegs))
rownames(mergedArms) = sub('\\.LogR','', rownames(mergedArms))


rows = grep('BLD|gastric', rownames(mergedSegs), ignore.case=T)

mergedSegs = mergedSegs[-rows,]
mergedArms = mergedArms[-rows,]

#ggplot(melt(mergedSegs), aes(Var2, value)) + geom_point() 
#ggplot(melt(mergedArms), aes(Var2, value)) + geom_point()


ms = mergedSegs
for (i in 1:ncol(ms)) {
  #ms[,i] = unit.var(ms[,i])
  ms[,i] = unit.var(ms[,i], z.mean[i], z.sd[i])
}

ggplot(melt(ms), aes(Var2, value)) + geom_point()
range(ms)
dim(ms)
# ma  = t(apply(mergedArms, 1, function(sample) {
#   sapply(1:ncol(mergedArms), function(i){
#     unit.var(sample[i], z.arms.mean[i], z.arms.sd[i])
#   })
# }))

ma = mergedArms
for (i in 1:ncol(ma)) {
  #ma[,i] = unit.var(ma[,i])
  ma[,i] = unit.var(ma[,i], z.arms.mean[i], z.arms.sd[i])
}
range(ma)
dim(ma)

nm = melt(ms[grep('BLD|gastric', rownames(ms), ignore.case=T),])
be = melt(ms[grep('BLD|gastric', rownames(ms), ignore.case=T, invert=T),])

grid.arrange(
  ggplot(be, aes(Var2, value, group=Var2)) + geom_point() + labs(title='BE samples'),
  ggplot(nm, aes(Var2, value, group=Var2)) + geom_point() + labs(titls='Blood/gastric')
)

cx = score.cx(ms, 1)
arrayDf = subtract.arms(ms, ma)
arrayDf = cbind(arrayDf, 'cx'=unit.var(cx, mn.cx, sd.cx))

ggplot(melt(arrayDf), aes(Var2, value, group=Var2)) + geom_point() + labs(title='All, scaled')


slabels = sample.info[,c('PatientID','Status', 'Samplename')]

missing = setdiff(sample.info$Samplename, rownames(arrayDf))
rows = grep('BLD|gastric', rownames(arrayDf), ignore.case=T)
bld_gastric = arrayDf[rows,]

df = arrayDf[-rows,]

slabels = subset(slabels, Samplename %in% rownames(df))

status = as.integer(as.factor(slabels$Status))-1
names(status) = slabels$Samplename

df = df[-grep(setdiff(rownames(df),names(status)), rownames(df)),]

#df = arrayDf[which(!rownames(arrayDf) %in% c(lowsca$Samplename,highPloidyCN$Samplename)),]



nl = 1000;folds = 10; splits = 5 
sets = create.patient.sets(slabels, folds, splits, 0.2)  

#sets = create.patient.sets(subset(sample.info, !Samplename %in% c(lowsca$Samplename,highPloidyCN$Samplename)) [c('PatientID','Samplename','Status')], folds, splits, 0.2)  

alpha.values = c(0,0.5,0.9,1)
## ----- All ----- ##
coefs = list(); plots = list(); performance.at.1se = list(); models = list(); cvs = list()
cache.dir = '~/Data/Ellie/Analysis/SNP/'
dir.create(cache.dir, recursive = T, showWarnings = F)
file = paste(cache.dir, 'all.pt.alpha.Rdata', sep='/')
#if (file.exists(file)) {
#  message(paste("loading", file))
#  load(file, verbose=T)
#} else 
{
  for (a in alpha.values) {
    fit0 <- glmnet(df, status, alpha=a, nlambda=nl, family='binomial', standardize=F)    
    autoplot(fit0) + theme(legend.position="none")
    l = fit0$lambda
    
    cv.patient = crossvalidate.by.patient(x=df, y=status, lambda=l, pts=sets, a=a, nfolds=folds, splits=splits, fit=fit0, select='deviance', opt=-1, standardize=F)
    
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
save(df,status,patient.info, file=paste(cache.dir,'snp.Rdata',sep='/'))
