

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 2)
  stop("Missing parameters: <patient ID> <model dir>")

library(ggplot2)
library(gridExtra)
library(glmnet)
library(readxl)
library(pROC)
library(ggfortify)
library(preprocessCore)

suppressPackageStartupMessages(source('lib/data_func.R'))
suppressPackageStartupMessages(source('lib/common_plots.R'))
suppressPackageStartupMessages(source('lib/cv-pt-glm.R'))



pt = args[1]
#pt = '512'

cache.dir = args[2]
#cache.dir = '~/Data/Ellie/Analysis/SNP'

dir.create(cache.dir, recursive = T, showWarnings = F)
file = paste(cache.dir, 'all.pt.alpha.Rdata', sep='/')
if (!file.exists(file))
  stop("Missing full model file, run snp-regression.R first")
load(file, verbose = T)

nl = 1000;folds = 10; splits = 5 
sets = create.patient.sets(slabels, folds, splits, 0.15)  

slabels$Prediction = NA
slabels$RR = NA

colnames(sets)[2] = 'Samplename'

file = paste(cache.dir, 'loo.Rdata', sep='/')

select.alpha = 0.9
secf = coefs[[select.alpha]]
a = select.alpha

performance.at.1se = c(); plots = list(); fit = c(); nzcoefs = c() 
preds = matrix(ncol=2, nrow=length(subset(slabels, PatientID == pt)[,2]), dimnames=list(subset(slabels, PatientID == pt)[,2],c('Prob','RR')))
# Remove each patient (LOO)
#for (pt in unique(slabels$PatientID)) {
print(pt)
samples = subset(slabels, PatientID != pt)[,2]

train.rows = which(rownames(df) %in% samples)
training = df[train.rows,]
test = as.matrix(df[-train.rows,])
if ( nrow(test) == ncol(df) ) test = t(test)

# Predict function giving me difficulty when I have only a single sample, this ensures the dimensions are the same
sparsed_test_data <- Matrix(data=0, nrow=nrow(test),  ncol=ncol(training),
                            dimnames=list(rownames(test),colnames(training)), sparse=T)
for(i in colnames(df)) sparsed_test_data[,i] = test[,i]

# Fit generated on all samples, including HGD
fitLOO <- glmnet(training, status[train.rows], alpha=a, family='binomial', nlambda=nl, standardize=F) # all patients
l = fitLOO$lambda
l = more.l(l)

cv = crossvalidate.by.patient(x=training, y=status[train.rows], lambda=l, a=a, nfolds=folds, splits=splits,
                              pts=subset(sets, Samplename %in% samples), fit=fitLOO, standardize=F)

plots = arrangeGrob(cv$plot+ggtitle('Classification'), cv$deviance.plot+ggtitle('Binomial Deviance'), top=pt, ncol=2)
fit = cv  

if ( length(cv$lambda.1se) > 0 ) {
  performance.at.1se = subset(cv$lambdas, lambda == cv$lambda.1se)$mean
  
  coef.1se = coef(fitLOO, cv$lambda.1se)[rownames(secf),]
  
  nzcoefs = as.data.frame(non.zero.coef(fitLOO, cv$lambda.1se))
  
  pm = predict(fitLOO, newx=sparsed_test_data, s=cv$lambda.1se, type='response')
  or = predict(fitLOO, newx=sparsed_test_data, s=cv$lambda.1se, type='link')
  
  preds[rownames(pm),'Prob'] = pm[,1]
  preds[rownames(pm),'RR'] = or[,1]
  
} else {
  warning(paste(pt, "did not have a 1se"))
}

save(plots, performance.at.1se, nzcoefs, fit, preds, file=file)
message("Finished")

