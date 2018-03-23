predict.progression <- function(rs,beObj, fitV, lambda, hdr=T, model.norm=T) {
  if (is.null(beObj) | is.null(fitV) | is.null(lambda))
    stop("No computation possible without model object or glmnet fit")

  if (ncol(rs) < 4)
    stop("File needs to contain 4 columns: chr, start, end, sample_value")
  #head(rs)
  
  cnCol = grep('Total_CN', colnames(rs))
  if (length(cnCol) >= 1)
    rs = rs[,c(1:cnCol)]
  
  segs <- tile.segmented.data(rs, size=5e6 )
  segM = as.matrix(segs[,-c(1:3)])
  rownames(segM) = paste(segs$chr, ':', segs$start, '-', segs$end, sep='')
  segM = t(segM)
  
  if (model.norm) {
    for (i in 1:ncol(segM)) 
      segM[,i] = unit.var(segM[,i], beObj@z.mean[i], beObj@z.sd[i])
  } else {
    segM = unit.var(segM)
  }
  
  arms <- tile.segmented.data(rs, size='arms')
  armsM = as.matrix(arms[,-c(1:3)])
  rownames(armsM) = paste(arms$chr, ':', arms$start, '-', arms$end, sep='')
  armsM = t(armsM)
  
  if (model.norm) {
    for (i in 1:ncol(armsM)) 
      armsM[,i] = unit.var(armsM[,i], beObj@z.arms.mean[i], beObj@z.arms.sd[i])
  } else {
    armsM = unit.var(armsM)
  }

  cx.score = score.cx(segM,1)
  
  mergedDf = subtract.arms(segM, armsM)
  mergedDf = cbind(mergedDf, 'cx' = unit.var(cx.score, beObj@mn.cx, beObj@sd.cx))
  
  # Predict function giving me difficulty when I have only a single sample, this ensures the dimensions are the same
  sparsed_test_data <- Matrix(data=0, nrow=nrow(mergedDf),  ncol=ncol(mergedDf),
                              dimnames=list(rownames(mergedDf),colnames(mergedDf)), sparse=T)
  for(i in colnames(mergedDf)) sparsed_test_data[,i] = mergedDf[,i]
  
  preds = predict(fitV, newx=sparsed_test_data, s=lambda, type='response')
  RR = predict(fitV, newx=sparsed_test_data, s=lambda, type='link')
  
  return(list('prob'=preds, 'rel.risk'=RR))
}

risk<-function(p) {
  if (!is.numeric(p) | (p > 1 | p < 0) ) stop("Numeric probability between 0-1 required")
  risk = ''
  if (p < 0.3 ) {
    risk = 'Low'
  } else if (p >= 0.3 & p < 0.7) {
    risk = 'Moderate'
  } else if ( p >= 0.7 ) {
    risk = 'High'
  } 
  return(risk)
}

rx<-function(pR) {
  pR$Risk = factor(pR$Risk, levels=c('Low','Moderate','High'), ordered=T)
  pR$Sample = as.character(pR$Sample)
  p53Col = grep('p53', colnames(pR), value=T, ignore.case=T)
  pathCol = grep('path', colnames(pR), value=T, ignore.case=T)
  
  rules = as.data.frame(matrix(ncol=3, nrow=nrow(pR), dimnames=list(c(), c('Time 1','Time 2','rule'))))
  # Consecutive
  for (i in 1:nrow(pR)) {
    risks = table(pR$Risk[i:(i+1)])
    p53 = NULL
    if (length(p53Col) > 0)
      p53 = table(pR[[p53Col]][i:(i+1)])
    
    rule = 'None'
    if ( risks['High'] == 2 || (length(pathCol) > 0 && length(which(grepl('HGD|IMC', pR[[pathCol]][i:(i+1)]))) > 0) ) {
      rule = 1
    } else if ( risks['High'] == 1 || (!is.null(p53) && p53['1'] > 0) ) {
      rule = 2
    } else if ( risks['Moderate'] > 0 ) {
      rule = 3
    } else if ( risks['Low'] == 2) {
      rule = 4
    }
    rules[i,] = cbind( pR$Sample[i], pR$Sample[(i+1)], as.integer(rule) )
  }  
  rules$rule = as.integer(rules$rule)
  
  rules$Rx = sapply(rules$rule, rule.rx)
  return(rules)
}

rule.rx<-function(n) {
  rr = c('Immediate RFA', 'Recheck 6-12 months',
         'Recheck 12-24 months','Regular surveillance 3-5 years')
  return(rr[n])
}