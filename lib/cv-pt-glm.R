require(mice)
require(glmnet)
require(ggplot2)
require(plyr)
require(GenomicRanges)


set.up.dataX<-function(tile.files, samplenames, labels) {
  get.loc<-function(df) {
    locs = do.call(rbind.data.frame, lapply(colnames(df), function(x) unlist(strsplit( x, ':|-'))))
    colnames(locs) = c('chr','start','end')
    locs[c('start','end')] = lapply(locs[c('start','end')], function(x) as.numeric(as.character(x)))
    locs$chr = factor(locs$chr, levels=c(1:22), ordered=T)
    locs
  }

  if (length(tile.files) > 1) {
    message('Loading arms')
    load(grep('arms',tile.files,value=T), verbose=F)
    armsDf = mergedDf
    armsDf = armsDf[grep('X|Y', rownames(armsDf), invert=T),]
    armsDf = armsDf[, intersect(colnames(armsDf), samplenames)]
    print(dim(armsDf))
  }

  # tile
  message("Loading tile")
  load(grep('arms', tile.files, value=T, invert=T), verbose=F)
  mergedDf = mergedDf[grep('X|Y', rownames(mergedDf), invert=T),]
  mergedDf = mergedDf[, intersect(colnames(mergedDf), samplenames)]
  print(dim(mergedDf))
  
  dysplasia.df = t(mergedDf[,intersect(names(labels), colnames(mergedDf))])
  labels = labels[intersect(names(labels), colnames(mergedDf))]
  
  dysplasia.df = apply(dysplasia.df, 2, unit.var)
  
  if (exists('armsDf')) {
    arms.df = t(armsDf[,intersect(names(labels), colnames(armsDf))])
    labels = labels[intersect(names(labels), colnames(armsDf))]
    arms.df = apply(arms.df, 2, unit.var)
    
    locs = get.loc(dysplasia.df)
    arm.locs = get.loc(arms.df)
    
    arms = makeGRangesFromDataFrame(arm.locs)
    locs = makeGRangesFromDataFrame(locs)
    
    tmp = dysplasia.df
    # subtract arms from 5e6 and merge both
    ov = findOverlaps(arms, locs)
    for (hit in unique(queryHits(ov))) {
      cols = subjectHits(ov)[which(queryHits(ov) == hit)]
      for (i in 1:nrow(tmp)) {
        tmp[i,cols] = tmp[i,cols] - arms.df[i,hit]
      }
    }
    dysplasia.df = cbind(tmp, arms.df)

    loc = get.loc(dysplasia.df)
    dysplasia.df = dysplasia.df[,order(loc$chr, loc$start)]
  }

  dysplasia.df = score.cx(patient.data, dysplasia.df)
  
  return(dysplasia.df)
}

add.demo.tocv <- function(patient.info, df) {

  patient.info = patient.info[match(rownames(df), patient.info$Samplename),]
  
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
  
  df = cbind(df, 'Sex'=NA,'C'=NA,'M'=NA,'Age'=NA, 'BMI'=NA, 'p53.Status'=NA, 'Smoking'=NA)
  df[,demo.cols] = as.matrix(complete(imp))
  
  return(df)
}

unit.varX <- function(x) {
  if (sd(x, na.rm=T) == 0) return(x)
  return((x-mean(x,na.rm=T))/sd(x,na.rm=T) )
}

score.cxX <- function(pt.d, df) {
  get.length.variance <- function(pd) { # per patient
    samples = grep('^D', colnames(pd$seg.vals), value=T)
    len.var = var(pd$seg.vals$end.pos - pd$seg.vals$start.pos) # length won't vary between samples
    
    sample.var = apply(as.matrix(pd$seg.vals[, samples]), 2, var)
    t(cbind(len.var, sample.var))
  }
  
  complexity.measures = lapply(pt.d, get.length.variance)
  
  df = cbind(df, 'cx' = apply(df, 1, function(x) length(which(x >= sd(x)*2 | x <= -sd(x)*2))))
  df[,'cx'] = unit.var(df[,'cx'])
  
  #df = cbind(df, 'segment.length'=0, 'sample.variance'=0)
  # for (pt in names(complexity.measures)) {
  #   for (sample in colnames(complexity.measures[[pt]])) {
  #     df[sample,c('segment.length', 'sample.variance')] = as.numeric(complexity.measures[[pt]][,sample])
  #   }
  # }
  
  #df[, 'segment.length'] = unit.var(df[, 'segment.length'])
  #df[, 'sample.variance'] = unit.var(df[, 'sample.variance'])
  
  return(df)
}


pi.hat<-function(x) exp(x)/(1+exp(x))

non.zero.coef<-function(fit, s) {
  cf =  as.matrix(coef(fit, s))
  cf[which(cf != 0),][-1]
}

coef.stability<-function(opt, nz.list) {
  for (i in 1:length(nz.list)) {
    df = as.data.frame(nz.list[[i]])
    colnames(df) = i
    opt = merge(opt, df, by='row.names', all.x=T)
    rownames(opt) = opt$Row.names
    opt$Row.names = NULL
    opt[,(i+1)] = as.integer(!is.na(opt[,(i+1)]))
  }
  opt = opt[order(sapply(rownames(opt), function(x) as.numeric(unlist(strsplit(x, ':'))[1])  )),]
  return(opt)
}

create.sample.sets<-function(data, uniqueID, n, splits, minR=0.2) {
  data = data %>% mutate(label = ifelse(Status == 'P', 1, 0))
  
  check.sets<-function(df, grpCol, min) {
    sets = table(cbind.data.frame('set'=df[[grpCol]], 'labels'=df$label))
    while ( (length(which(sets/rowSums(sets) < minR) ) >= 2 | length(which(sets/rowSums(sets) == 0)) > 0) ) {
      print(sets/rowSums(sets))
      s = sample(rep(seq(5), length = length(unique(df[[uniqueID]]))))
      df2 = merge(df, cbind(unique(df[[uniqueID]]), 'tmpgrp'=s), by=1)
      df[[grpCol]] = df2$tmpgrp
      sets = table(cbind.data.frame('set'=df[[grpCol]], 'labels'=df$label))
    }
    return(df[[grpCol]])
  }
  
  s = sample(rep(seq(splits), length = length(unique(data[[uniqueID]]))))
  data = merge(data %>% dplyr::select(matches(uniqueID), everything()), cbind(unique(data[[uniqueID]]), 'group'=s), by=1)  
  colnames(data)[ncol(data)] = c('fold.1')
  data$fold.1 = check.sets(data, 'fold.1', minR)
  
  for (i in 2:n) {
    s = sample(rep(seq(splits), length = length(unique(data[[uniqueID]]))))
    data = merge(data %>% dplyr::select(matches(uniqueID), everything()), cbind(unique(data[[uniqueID]]), 'group'=s), by=1)  
    foldcol = grep('group',colnames(data))
    colnames(data)[foldcol] = paste('fold',i,sep='.')
    data[[paste('fold',i,sep='.')]] = check.sets(data, paste('fold',i,sep='.'))
  }
  return(data)
  
  
}

create.patient.sets<-function(pts, n, splits, minR=0.2) {
  # This function just makes sure the sets don't become too unbalanced with regards to the labels.
  pts$label = 0
  pts[which(pts$Status == 'P'), 'label'] = 1

  uniquePtID = colnames(pts)[1]

  check.sets<-function(df, grpCol, min) {
    sets = table(cbind.data.frame('set'=df[[grpCol]], 'labels'=df$label))
    while ( (length(which(sets/rowSums(sets) < minR) ) >= 2 | length(which(sets/rowSums(sets) == 0)) > 0) ) {
      print(sets/rowSums(sets))
      s = sample(rep(seq(5), length = length(unique(df[[uniquePtID]]))))
      df2 = merge(df, cbind(unique(df[[uniquePtID]]), 'tmpgrp'=s), by=1)
      df[[grpCol]] = df2$tmpgrp
      sets = table(cbind.data.frame('set'=df[[grpCol]], 'labels'=df$label))
    }
    return(df[[grpCol]])
  }
  
  s = sample(rep(seq(splits), length = length(unique(pts[[uniquePtID]]))))
  patients = merge(pts, cbind(unique(pts[[uniquePtID]]), 'group'=s), by=1)  
  colnames(patients)[5] = c('fold.1')
  patients$fold.1 = check.sets(patients, 'fold.1', minR)
  
  for (i in 2:n) {
    s = sample(rep(seq(5), length = length(unique(pts[[uniquePtID]]))))
    patients = merge(patients, cbind(unique(patients[[uniquePtID]]), 'group'=s), by=1)  
    foldcol = grep('group',colnames(patients))
    colnames(patients)[foldcol] = paste('fold',i,sep='.')
    patients[[paste('fold',i,sep='.')]] = check.sets(patients, paste('fold',i,sep='.'))
  }
  return(patients)
}

binomial.deviance<-function(pmat, y) {
  # Binomial deviance, lifted from cv.lognet
  prob_min = 1e-05; prob_max = 1 - prob_min
  pmat = pmin(pmax(pmat, prob_min), prob_max)
  dev = apply(pmat, 2, function(p)  -2*((y==1)*log(p)+(y==0)*log(1-p)) )
  return(dev)
}

crossvalidate.by.patient<-function(x,y,lambda,pts,a=1, sampleID = 1, nfolds=10, splits=5, fit=NULL, minR=0.2, select='deviance', opt=-1, ...) {
  if (nfolds > 5) minR = 0.1
  x = as.matrix(x)
  message(paste("Running", splits, "splits",nfolds,"times on", paste(dim(x), collapse=':'), 'alpha=',a ))
  fit.e = list()
  
  if (length(grep('fold',colnames(pts))) == nfolds) {
    tpts = pts %>% data.frame
  } else {
    tpts = create.sample.sets(pts, sampleID, nfolds, splits, minR) %>% data.frame
  }
  
  cv.pred = (matrix(nrow=0, ncol=length(lambda)))
  cv.binomial.deviance = (matrix(nrow=0, ncol=length(lambda)))
  for (n in 1:nfolds) {  
    message(paste("Fold ", n))
    setCol = grep(paste('^fold.',n,'$',sep=''), colnames(tpts))
    
    cv.class = matrix(nrow=splits, ncol=length(lambda))
    deviance = matrix(nrow=splits, ncol=length(lambda))
    for (i in 1:splits) { 
      message(paste0('  ', i, " split"))
      test.rows = which(rownames(x) %in% tpts[which(tpts[,setCol] == i), sampleID])
      test = x[test.rows,]
      training = x[-test.rows,]
      # pre-spec lambda seq
      fitCV = glmnet(training, y[-test.rows], lambda=lambda, alpha=a, family='binomial',  standardize=F) 
      # autoplot(fit) + theme(legend.position="none")
      
      # Confusion matrix: quantitiative, pos results + neg results / number of test rows
      pred <- pi.hat(predict(fitCV, test, type='link'))  # pi.hat turns these into probabilities from logit
      cv.class[i,] = apply(pred, 2, function(p) {
        (p%*%y[test.rows] + (1-p) %*% (1-y[test.rows]))/length(test.rows)
      })

      # Binomial deviance, lifted from cv.lognet
      dev = binomial.deviance(predict(fitCV, test, type='response'), as.factor(y[test.rows]))
      deviance[i,] = apply(dev, 2, weighted.mean, w=rep(1, nrow(dev)), na.rm=T)
      
      fit.e[[length(fit.e)+1]] = fitCV
    }
    cv.pred = rbind(cv.pred, cv.class)    
    cv.binomial.deviance = rbind(cv.binomial.deviance, deviance)
  }
  
  # Not really sure what to do with the binomial deviance now...
  lambdas = cbind.data.frame('lambda.at'=1:ncol(cv.pred),
                              'mean'= colMeans(cv.pred), 
                              'sme'= apply(cv.pred, 2, sd)/sqrt(nrow(cv.pred)), 
                              'sd'= apply(cv.pred, 2, sd), 
                              'lambda'= lambda,
                              'log.lambda' = log(lambda),
                              'mean.b.dev' = colMeans(cv.binomial.deviance),
                              'sd.b.dev' = apply(cv.binomial.deviance, 2,sd),
                              'sme.b.dev' = apply(cv.binomial.deviance, 2, sd)/sqrt(nrow(cv.binomial.deviance)),
                              'lambda.min' = F, 'lambda+1se' = F, 'lambda-1se' = F) 
  
  if (select == 'classification') {
    # Minimize the classification mean error & select min lambda
    lambdas[which.max(subset(lambdas, sme < median(sme) & sme > min(sme))$mean), 'lambda.min'] = T
    
    lmin = subset(lambdas, lambda.min)
    se1 = lmin$mean + (sd(lambdas$mean)/sqrt(nrow(lambdas)))
    
    
    lambda.1se.search = subset(lambdas, mean >= se1 & !lambda.min & log.lambda > lmin$log.lambda) 
    if (nrow(lambda.1se.search) > 0)
      lambdas[arrange(lambda.1se.search, mean)[1,'lambda.at'], 'lambda+1se'] = T
    
    lambda.1se.search = subset(lambdas, mean <= se1 & !lambda.min & log.lambda < lmin$log.lambda) 
    if (nrow(lambda.1se.search) > 0)
      lambdas[arrange(lambda.1se.search, -mean)[1,'lambda.at'], 'lambda-1se'] = T
    
  } else if (select == 'deviance') {
    lambdas[which.min(lambdas$mean.b.dev), 'lambda.min'] = T
    lmin = subset(lambdas, lambda.min)
    
    se1 = lmin$mean.b.dev + (sd(lambdas$mean.b.dev)/sqrt(nrow(lambdas)))
    
    lambda.1se.search = subset(lambdas, mean.b.dev >= se1 & !lambda.min & log.lambda > lmin$log.lambda) 
    if (nrow(lambda.1se.search) > 0)
      lambdas[arrange(lambda.1se.search, mean.b.dev)[1,'lambda.at'], 'lambda+1se'] = T
    
    lambda.1se.search = subset(lambdas, mean.b.dev <= se1 & !lambda.min & log.lambda < lmin$log.lambda) 
    if (nrow(lambda.1se.search) > 0)
      lambdas[arrange(lambda.1se.search, -mean.b.dev)[1,'lambda.at'], 'lambda-1se'] = T
  }
  
  opt = ifelse (opt > 0, 'lambda+1se', 'lambda-1se')
  
  
  lambda.min = subset(lambdas, lambda.min)$lambda
  lambda.1se = lambdas[which(lambdas[[opt]]), 'lambda']
    
  nzcf = lapply(fit.e, non.zero.coef, s=lambda.1se)
  
  plots = plot.patient.cv(lambdas, fit)
  plots$performance = plots$performance + theme(legend.position='bottom') + scale_colour_discrete(name = "") +
    labs(title=paste(splits,' splits, ',nfolds,' folds, alpha=',a, sep=''), y='mean pred.', x='log(lambda)') 
  plots$deviance = plots$deviance + theme(legend.position='bottom') + scale_colour_discrete(name="") + 
    labs(title=paste(splits,' splits, ',nfolds,' folds, alpha=',a, sep=''), y='mean Binomial Deviance', x='log(lambda)') 
  
  return(list('max.cm'=lambdas[lambdas$`lambda+1se` ==T,'mean'], 
              'lambda.min'=lambda.min, 'lambda.1se'=lambda.1se, 'lambdas'=lambdas, 'plot'=plots$performance, 'deviance.plot'=plots$deviance, 'non.zero.cf'=nzcf))
}


plot.patient.cv<-function(df, fit=NULL) {
  gp = ggplot(df, aes(y=mean,x=log.lambda)) + geom_point() + geom_errorbar(aes(ymin=mean-sme, ymax=mean+sme)) + 
    geom_point(data=subset(df, lambda.min == T), aes(y=mean, x=log.lambda, colour="lambda.min"), size=2 ) +  
    geom_vline(xintercept = subset(df, lambda.min == T)$log.lambda, colour='grey') +
    annotate("text", x=subset(df, lambda.min == T)$log.lambda, y=min(df$mean)+sd(df$mean), 
             label=paste(round(df[subset(df, lambda.min == T)$lambda.at,c('mean', 'sme')],3), collapse='\n+/-')) 
  
  gpD = ggplot(df, aes(y=mean.b.dev,x=log.lambda)) + geom_point() +
    geom_errorbar(aes(ymin=mean.b.dev-sme.b.dev, ymax=mean.b.dev+sme.b.dev)) + 
    geom_point(data=subset(df, lambda.min == T), aes(y=mean.b.dev, x=log.lambda, colour="lambda.min"), size=2 ) +          
    geom_vline(xintercept = subset(df, lambda.min == T)$log.lambda, colour='grey') +
    annotate("text", x=subset(df, lambda.min == T)$log.lambda, y=min(df$mean)+sd(df$mean), 
             label=paste(round(df[subset(df, lambda.min == T)$lambda.at,c('mean', 'sme')],3), collapse='\n+/-')) 
  
  
  if (nrow(subset(df, `lambda+1se` == T)) > 0) {
    gp = gp + geom_point(data=subset(df, `lambda+1se` == T), aes(y=mean, x=log.lambda, colour="lambda+1se"), size=2 ) +
      geom_vline(xintercept = subset(df, `lambda+1se` == T)$log.lambda, colour='grey') +
      annotate("text", x=subset(df, `lambda+1se` == T)$log.lambda, y=min(df$mean)+sd(df$mean), 
               label=paste(round(df[subset(df, `lambda+1se` == T)$lambda.at,c('mean', 'sme')],3), collapse='\n+/-')) 
    
    gpD = gpD + geom_point(data=subset(df, `lambda+1se` == T), aes(y=mean.b.dev, x=log.lambda, colour="lambda+1se"), size=2 ) +
      geom_vline(xintercept = subset(df, `lambda+1se` == T)$log.lambda, colour='grey') +
      annotate("text", x=subset(df, `lambda+1se` == T)$log.lambda, y=min(df$mean)+sd(df$mean), 
               label=paste(round(df[subset(df, `lambda+1se` == T)$lambda.at,c('mean', 'sme')],3), collapse='\n+/-')) 
  }

  if (nrow(subset(df, `lambda-1se` == T)) > 0) {
    gp = gp + geom_point(data=subset(df,  `lambda-1se` == T), aes(y=mean, x=log.lambda, colour="lambda-1se"), size=2 ) +
      geom_vline(xintercept = subset(df,  `lambda-1se` == T)$log.lambda, colour='grey') +
      annotate("text", x=subset(df,  `lambda-1se` == T)$log.lambda, y=min(df$mean)+sd(df$mean), 
               label=paste(round(df[subset(df,  `lambda-1se` == T)$lambda.at,c('mean', 'sme')],3), collapse='\n+/-')) 
    
    gpD = gpD + geom_point(data=subset(df,  `lambda-1se` == T), aes(y=mean.b.dev, x=log.lambda, colour="lambda-1se"), size=2 ) +
      geom_vline(xintercept = subset(df,  `lambda-1se` == T)$log.lambda, colour='grey') +
      annotate("text", x=subset(df,  `lambda-1se` == T)$log.lambda, y=min(df$mean)+sd(df$mean), 
               label=paste(round(df[subset(df,  `lambda-1se` == T)$lambda.at,c('mean', 'sme')],3), collapse='\n+/-')) 
  }

  
  gp = gp + theme(legend.position='bottom') + scale_colour_discrete(name = "") 
  gpD = gpD + theme(legend.position='bottom') + scale_colour_discrete(name="") 
  
  if (!is.null(fit)) {
    df$nzcoef = sapply(df$lambda, function(l) length(non.zero.coef(fit, l)))
    
    coef.min = subset(df, lambda.min == T)$nzcoef
    coefP1se = subset(df, `lambda+1se` == T)$nzcoef
    coefM1se = subset(df, `lambda-1se` == T)$nzcoef
    
    coef.text = arrange(subset(df,nzcoef %in% c(0,  coef.min, coefP1se, coefM1se, max(nzcoef) )), -lambda.min, -`lambda+1se`, -`lambda-1se`)
    coef.text = arrange(coef.text[!duplicated(coef.text$nzcoef),], -nzcoef)
    
    gp = gp + annotate("text", y=min(df$mean), x=coef.text$log.lambda, label=coef.text$nzcoef, color='darkblue')
    gpD = gpD + annotate("text", y=max(df$mean.b.dev), x=coef.text$log.lambda, label=coef.text$nzcoef, color='darkblue')
  }
  return(list('performance'=gp, 'deviance'=gpD))
}

plot.cv<-function(cv, title='') {
  grid.arrange( arrangeGrob(cv$plot+ggtitle('Classification'), cv$deviance.plot+ggtitle('Binomial Deviance'), top=title, ncol=2) )
}


more.l<-function(lambda) {
  if (min(lambda) > exp(-5))
    lambda = c(lambda, seq( min(lambda), exp(-5), -1e-3))
  l = unique(sort(c(lambda, seq(exp(-5), exp(-10), -1e-6),
             seq(exp(-10), exp(-15), -1e-8),  seq(exp(-15), exp(-20), -1e-9)
             ), decreasing=T))
  return(l)
}


run.fit<-function(df, labels, a=1, more.l=F, opt.lambda=-1) {
  print(a)
  fit0 <- glmnet(df, labels, alpha=a, nlambda=nl, family='binomial', standardize=F) # all patients
  
  l = fit0$lambda
  if (more.l) l = more.l(fit0$lambda)
  
  cv = crossvalidate.by.patient(x=df, y=labels, lambda=l, pts=sets, a=a, nfolds=folds, splits=splits, fit=fit0, select='deviance', opt=opt.lambda, standardize=F)
  
  opt.lambda = ifelse(opt.lambda > 0, 'lambda+1se', 'lambda-1se')
  
  return( list('cv'=cv, 'perf'=cv$lambdas[which(cv$lambdas[[opt.lambda]]),], 'coefs'=coef.stability(as.data.frame(non.zero.coef(fit0, cv$lambda.1se)), cv$non.zero.cf), 'fit0'=fit0 ) )
}


