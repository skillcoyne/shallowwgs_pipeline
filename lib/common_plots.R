
myPal = rev(RColorBrewer::brewer.pal(11, 'RdYlBu'))

plot.theme = theme(text=element_text(size=12), panel.background=element_blank(), strip.background =element_rect(fill="white"),  
                   strip.text = element_text(size=12), 
                   axis.line=element_line(color='black'), panel.grid.major=element_line(color='grey90'),
                   panel.border = element_rect(color="grey", fill=NA, size=0.5), panel.spacing = unit(0.1, 'lines')  ) 

get.roc.stat<-function(roc, x='best') {
  auc = c( pROC::auc(roc), range(pROC::ci.auc(roc))  )
  sens = c(pROC::coords(roc, x, transpose=T)[['sensitivity']], pROC::ci.se(roc,pROC::coords(roc, x, transpose=T)[['specificity']], conf.level=0.95))
  spec = c(pROC::coords(roc, x, transpose=T)[['specificity']], pROC::ci.sp(roc,pROC::coords(roc, x, transpose=T)[['sensitivity']], conf.level=0.95))
  
  stat = rbind('AUC'=auc, 'Sensitivity'=sens[c(1,2,4)], 'Specificity'=spec[c(1,2,4)])
  colnames(stat) = c ('value','CI.min','CI.max')
  stat = as.data.frame(stat)
  stat$model = roc$model
  
  if (!is.null(roc$n) & is.numeric(roc$n)) stat$n = roc$n
  
  stat$Measure = rownames(stat)
  stat
}

get.auc.at<-function(roc, foc='spec', p=c(1,0.99)) {
  # To find the auc at a given spec or sens
  pROC::auc(roc, partial.auc=p, partial.auc.focus=foc, partial.auc.correct=T)
}


get.legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

cn.mtn.plot<-function(df, label, type='bar') {
  df = arrange(df, chr, start)
  df$mean.value = df$value
  rows = which(df$mean.value > 1 | df$mean.value < -1)
  
  df$cuts = cut(df$mean.value, breaks=c(min(df$mean.value), -1,0,1, max(df$mean.value)), include.lowest = T, ordered_result = T)
  medDf = df[-rows,]
  rangeDf = df[rows,]
  
  medDf$cuts = factor('-1-1')
  rangeDf$cuts = factor(ifelse(rangeDf$mean.value > 1, '> 1', '< -1'))

  p = ggplot(df, aes(x=chr.length)) + facet_grid(Status~chr, scales='free_x', space='free_x', labeller=label) 
  
  if (type == 'seg') {
    riskCols = brewer.pal(10, "PRGn")[c(1, 5, 10)]
    p = p + geom_segment(aes(x=start, xend=end, y=mean.value, yend=mean.value, color=cn), show.legend=T, size=1) + 
      scale_color_manual(values=riskCols) 
  }
  if (type == 'bar') {
    highCol = rev(brewer.pal(6, "Purples"))[c(1,3)]
    lowCol = rev(brewer.pal(6, "Greens"))[c(3,1)]
    # grey = brewer.pal(3,'Greys')[1]
    
    colors = c('#74C476', 'grey44','#9E9AC8')
    
    p = p + 
      geom_rect(data=rangeDf, aes(xmin=start, xmax=end, ymin=0, ymax=mean.value, fill=cuts), show.legend=T ) + 
      geom_rect(data=medDf, aes(xmin=start, xmax=end, ymin=0, ymax=mean.value, fill=cuts), show.legend=F ) + 
      scale_fill_manual(values=colors, limits=c( levels(rangeDf$cuts)[1], levels(medDf$cuts), levels(rangeDf$cuts)[2]), name='') 
  }
  p + labs(x='', y='Mean adjusted segmentation value', title='All samples') + theme_minimal() +
    theme(text=element_text(size=12), axis.text.x=element_blank(), legend.position='right', 
          panel.grid = element_blank(), panel.background = element_blank(),   
          panel.border = element_rect(color="grey88", fill=NA, size=0.5), panel.spacing.x = unit(0, 'lines') ) 
}

get.cn.state = function(x) {
  str = 'neutral'
  if(x <= -1) {
    str = 'loss'
  } else if (x >= 1) {
    str = 'gain'
  }
  return(str)
}

model.performance<-function(p1se, coeffs, folds, splits, lims=NULL) {
  performance = do.call(rbind.data.frame, p1se)
  performance$alpha = rownames(performance)
  
  if (grepl('chr', performance$alpha[1])) {
    performance$alpha = factor(performance$alpha, levels=c(paste('chr',c(1:22), sep=''), setdiff(rownames(performance), paste('chr',c(1:22),sep=''))) )
  }
  
  ## check how often the feature(s) is selected at that lamda in each split. "stability selection"
  coef.stable = lapply( coeffs, function(cf) {
    sort(rowSums(cf[,-1]), decreasing=T)
  })
  coef.stable = coef.stable[names(which(sapply(coef.stable, length) > 0))]
  coef.stable = lapply(coef.stable, function(cf) cf/(folds*splits))
  
  performance = cbind(performance, do.call(rbind,  lapply(coef.stable, function(x) {
    cbind('n.Coef'=length(x), '25%'=length(which(x>=.25)), 
          '50%'=length(which(x>=.5)), '75%'=length(which(x>=.75 )), 
          '100%'=length(which(x==1)) )
  })))
  
  p = ggplot(performance, aes(alpha, mean, fill=(alpha != 0.9))) + ylim(0,1) + geom_col() +
    geom_text( aes(label=round(mean, 3)), vjust=3) +
    geom_errorbar(aes(ymin=mean-sme, ymax=mean+sme), width=0.3, size=1, color='grey39') + 
    geom_text(aes(y=0.1, label=n.Coef), ) +
    geom_text(aes(label=paste('(',`75%`,')',sep=''), y=0.05)) + 
    scale_fill_brewer(palette='Paired') +
    #scale_fill_manual(values=c('red', 'grey39')) + 
    plot.theme + theme(legend.position = 'none') +
    labs(x='Elasticnet penalty value: ridge <-> lasso', y='Model classification at lambda-1se')
  
  
  if (!is.null(lims))
    p = p+ylim(lims)
  
  cfs = as.data.frame(matrix(data=0,nrow=length(unique(names(table(unlist(lapply(coef.stable, names)))))), ncol=length(names(coef.stable)), 
                             dimnames=list( unique(names(table(unlist(lapply(coef.stable, names))))), names(coef.stable) )))
  for (i in names(coef.stable)) 
    cfs[intersect(rownames(cfs), names(coef.stable[[i]])), i] = 1
  return(list('performance'=performance, 'stable.coeffients'=coef.stable, 'cf.per.feature'=cfs,'plot'=p))
}

roc.plot <- function(roc,title="", incS = T, threshold = 'best') {
  df = cbind.data.frame('specificity'=rev(roc$specificities), 'sensitivity'=rev(roc$sensitivities))
  
  best = as_tibble(round(pROC::coords(roc, threshold, transpose=F),2)) 
  label = paste0('AUC:',round(roc$auc, 2))
  if (incS) {
    label = paste0('(FPR ', round((1-best$specificity),2),' TPR ', round(best$sensitivity,2),')\n', label )
  }
  
  st = paste0('95% CI: ', paste(round(range(pROC::ci.auc(roc$auc)),2), collapse='-'))
  
  ggplot(df, aes(specificity, sensitivity)) + geom_line() + scale_x_reverse() + 
    geom_label(data=as.data.frame(t(round(pROC::coords(roc, threshold, transpose=T),2))), 
               aes(x=specificity, y=sensitivity,label=label), nudge_y=ifelse(incS,0.1,0), nudge_x = ifelse(incS, 0.15, 0)) +
    labs(title=title, subtitle=st, x='Specificity (1-FPR)', y='Sensitivity (TPR)') + plot.theme
}

multi.roc.plot <- function(rocList, threshold='best', title="ROC", palette=NULL, colors=NULL) {
  aucs = do.call(rbind, lapply(rocList,function(r) 
    cbind.data.frame('AUC'=r$auc,'model'=r$model,'Specificity'=pROC::coords(r, threshold,transpose=T)[['specificity']], 'Sensitivity'=pROC::coords(r,threshold,transpose=T)[['sensitivity']])))
  aucs$best.y = aucs$Sensitivity
  aucs$best.x = aucs$Specificity

  if (nrow(aucs) > 2) {
    aucs$y = rev(seq(0.1,1, 0.1))[1:nrow(aucs)]
    aucs$x = 0.1
  }
  
  aucs$label=with(aucs, paste('(TPR=', round(Sensitivity,2), ' FPR=', round((1-Specificity),2), ')\nAUC ', round(AUC,2), sep=''))

  #aucs$label=with(aucs, paste(model, '\nAUC ', round(AUC,2), sep=''))
  
    
  df = do.call(rbind, lapply(rocList, function(r) 
    cbind.data.frame('Specificity'=rev(r$specificities), 'Sensitivity'=rev(r$sensitivities),'model'=r$model) 
  ))

  p = ggplot(df, aes(Specificity, Sensitivity,color=model,fill=model)) +
        geom_line() + scale_x_reverse() +
        geom_point(data=aucs,aes(best.x,best.y), color='grey39', shape='*', size=8, show.legend = F) +
        labs(title=title, x='Specificity (1-FPR)', y='Sensitivity (TPR)')  + plot.theme
  
  if (nrow(aucs) <= 2) {
    p = p + geom_label_repel(data=aucs,aes(best.x,best.y,label=label), color='white', show.legend = F, segment.size = 0.2, nudge_x=0.5, segment.colour = 'grey39') 
  } else {
    p = p + geom_label_repel(data=aucs,aes(x,y,label=label), color='white', show.legend = F) 
  }
  
  if (is.null(colors) & !is.null(palette)) {
    colors = RColorBrewer::brewer.pal(length(rocList)+1, palette)
  } else if (is.null(colors) & is.null(palette)) {
    colors = RColorBrewer::brewer.pal(length(rocList)+1,'Set1')
  }
  p = p + scale_color_manual(name='', values=colors) + scale_fill_manual(name = '', values=colors) + theme(legend.position = 'bottom')
 
  return(p)
}

pred.tiles <- function(df, probs=F, OR=F, path=5, colors=c('green3','orange','red3'), ylim=NULL, ...) {
  #if (is.null(ylim)) ylim = 1:max(as.integer(levels(droplevels(df$Block))))
  if (is.null(ylim)) ylim = 1:max(df$Block)

  df$months.before.final = factor(df$months.before.final, ordered = T)
  df$months.before.final = factor(df$months.before.final, levels=rev(levels(df$months.before.final)))
  
  p = ggplot(df, aes(months.before.final, Block)) # + scale_y_discrete(limits=ylim)

  if (OR) {
    p = p + geom_tile(aes(fill=OR), color='white') + scale_fill_gradientn(colors = myPal, limits=c(-9.3,15), name='RR') 
  } else if (probs) {
    p = p + geom_tile(aes(fill=Prediction), color='white') + scale_fill_gradientn(colors = myPal, limits=c(0,1), name='P(P)') 
  } else {
    p = p + geom_tile(aes(fill=Risk), color='white', size=2) + scale_fill_manual(values=colors, limits=levels(df$Risk), name='Progression')
  }
  
  if (!is.null(path)) {
    p = p + geom_point(aes(shape=Pathology), fill='white', color='white', size=path) + 
        scale_shape_manual(values=c(1,0,15,24,25), limits=levels(df$Pathology), labels=c('NDBE','ID','LGD','HGD','IMC'), guide=guide_legend(override.aes=list(fill='white', color='white')))
  }
  p + facet_wrap(~Patient, ncol=1) +
    labs(y='Oesophageal Location (OGJ..)',x='Endoscopy (Base...End)', title='') + scale_y_continuous(expand=c(0,0), breaks = ylim) + 
    plot.theme + theme(axis.text = element_blank(), legend.key=element_rect(fill='grey39'), panel.background=element_rect(colour = 'black'), panel.grid.major=element_blank(), panel.spacing = unit(0.2, 'lines'), panel.border = element_rect(color="black", fill=NA, size=0.5)  ) 
}

plot.cn.chr<-function(df, chrom='17', info=NULL, haz=NULL, genes=NULL, label=NULL) {
  #df = subset(df, chr == chrom)
  df$mean.value = df$value  
  df = subset(df, chr == chrom)# & seg.type != 'arms')  
  
  df$cuts = cut(df$mean.value, breaks=c(min(df$mean.value), -1,0,1, max(df$mean.value)), include.lowest = T, ordered_result = T)
  rows = which(df$mean.value > 1 | df$mean.value < -1)
  medDf = df[-rows,]
  rangeDf = df[rows,]
  arms = which(rangeDf$seg.type == 'arms')
  
  medDf$cuts = factor('-1-1')
  rangeDf$cuts = factor(ifelse(rangeDf$mean.value > 1, '> 1', '< -1'))

  colors = c('#74C476', 'darkblue','#9E9AC8')
  
  if (is.null(label)) label = labeller(chr=label_both)
  
  p = ggplot(df, aes(x=chr.length)) + facet_grid(Status~chr, scales='free_x', labeller=label)
  
  if (!is.null(info)) {
    info = info[which(info$chrom == chrom),]
    p = p + geom_vline(data=info, aes(xintercept=chr.cent-cent.gap), color='grey78', size=2 ) 
  }
  
  p = p +
    geom_rect(data=rangeDf[arms,], aes(xmin=start, xmax=end, ymin=0, ymax=mean.value, fill=cuts), show.legend=F, alpha=0.3 ) + 
    geom_rect(data=rangeDf[-arms,], aes(xmin=start, xmax=end, ymin=0, ymax=mean.value, fill=cuts), show.legend=T, alpha=0.8 ) + 
    geom_rect(data=medDf, aes(xmin=start, xmax=end, ymin=0, ymax=mean.value, fill=cuts), show.legend=F ) + 
    scale_fill_manual(values=colors, limits=c( levels(rangeDf$cuts)[1], levels(medDf$cuts), levels(rangeDf$cuts)[2]), name='') 
  
  if (!is.null(haz)) {
    haz = subset(haz, chr == chrom)
    if (nrow(haz) > 0 & length(grep('loss',unique(haz$gl))) > 0 ) {
      p = p + geom_rect(data=subset(haz, gl == 'loss'), aes(xmin=start, xmax=end, ymin=0, ymax=max), color='blue3', alpha=0, show.legend=F, linetype='dotted',size=0.5) 
    } 
    
    if (nrow(haz) > 0 & length(grep('gain',unique(haz$gl))) > 0) {
      p = p + geom_rect(data=subset(haz, gl == 'gain'), aes(xmin=start, xmax=end, ymin=0, ymax=max), color='goldenrod2', alpha=0, show.legend=F, linetype='dotted', size=0.5) 
    }
  }
  
  if (!is.null(genes)) {
    genes = subset(genes, chr == chrom)
    if (nrow(genes) > 0) p = p + geom_point(data=genes, aes(x=start,  y=value), color='darkblue', size=3) +
        geom_text_repel(data=genes, aes(x=start, y=value, label=`Gene Symbol`), color='darkblue', fontface='italic') 
  }
  
  if (!is.null(info)) {
    info = info[which(info$chrom == chrom),]
    p = p + #geom_vline(data=info, aes(xintercept=chr.cent-cent.gap), color='grey45', size=2, alpha=0.5 ) +
      geom_label(data=info, aes(x=chr.cent-cent.gap*3, y=max(df$mean.value)), label='p-arm') +
      geom_label(data=info, aes(x=chr.cent+cent.gap*1.5, y=max(df$mean.value)), label='q-arm')
  }
  
  p + plot.theme + 
    labs(x='', y='Mean adj. segmentation value', title='') +
    theme(axis.text.x=element_blank(), legend.position='none', panel.grid=element_blank(), panel.background=element_blank())  
}


plot.endo<-function(df, minM=NULL, maxM=NULL, bin=F) {
  if (is.null(minM)) minM = min(df$months.before.final)
  if (is.null(maxM)) maxM = max(df$months.before.final)
  
  colors = c('green4','red4')
  
  df = arrange(df %>% group_by(Patient) %>% mutate('overall.pred'=sum(as.integer(pred > 0))/length(Patient)), Patient)
  
  df = transform(df, Patient=reorder(Patient, -overall.pred) ) 
  
  if (bin) {
    df$endo.pred = ifelse(df$pred > 0, 'Pos','Neg')
    ggplot(df, aes(months.before.final, endo.pred, group=1, color=endo.pred)) + facet_grid(Patient ~.) + 
      scale_x_continuous(breaks=seq(minM, maxM, by = 12)) +
      geom_line( color='grey39') + geom_point(aes(color=endo.pred)) +
      geom_label_repel(aes(label=max.path, fill=endo.pred), color='white', fontface='bold', size=4, label.size=0, label.r=unit(0.25, 'lines'), label.padding=unit(0.15, 'lines'), show.legend=F) + 
      scale_fill_manual(values=colors) + 
      scale_color_manual(values=colors, name='Endoscopy Prediction') +
      labs(x='Months before diagnosis', y='Predictive Endoscopy') + plot.theme + 
      theme(panel.border = element_rect(color = "grey50", fill = NA, size = 1), legend.position='bottom', strip.background =element_rect(fill="lightcyan2"), panel.grid.major.x=element_blank(), axis.text.x=element_text(angle=45, hjust=1))
    
  } else {
    ggplot(df, aes(months.before.final, pred.ratio, group=1, color=pred.ratio, fill=pred.ratio)) + facet_grid(Patient ~.) +
      ylim(0,1) + geom_abline(slope=0, intercept=seq(.25,.75,.25), color='grey', linetype='dashed') +
      geom_col(position='dodge', aes(y=pred.ratio), color='white') +# geom_line() +
      scale_x_continuous(breaks=seq(minM, maxM, by = 12)) + 
      geom_label(aes(y=ifelse(pred.ratio==1, 0.75, pred.ratio), label=max.path, fill=pred.ratio), color='white', fontface='bold', size=4, label.size=0, label.r=unit(0.25, 'lines'), label.padding=unit(0.15, 'lines'), nudge_y=0.1) +
      labs(x='Months before diagnosis', y='Sample prediction ratio') + 
      scale_fill_continuous(high='orangered4', low='green4', name='Sample Prediction Ratio') + plot.theme + 
      theme(panel.border = element_rect(color = "grey50", fill = NA, size = 1), legend.position='bottom',strip.background =element_rect(fill="lightcyan2"), panel.grid.major.x=element_blank(), axis.text.x=element_text(angle=45, hjust=1))
  }
}

vig.plot<-function(df, cols=BarrettsProgressionRisk::riskColors(), pathology=F, add.title=F) {
  scale = seq(min(df$Endoscopy.Year), max(df$Endoscopy.Year), 1)
  subtitle = with(df, paste(unique(Age.at.diagnosis),unique(Sex),': ', sep=''))
  if (!is.na(unique(df$Smoking)))                             
    subtitle = paste(subtitle, ifelse(unique(df$Smoking) == 'Y', 'smoker, ','non-smoker, '),sep='')
  
  subtitle = paste(subtitle, 'M', unique(df$Maximal), sep='')
  
  if (!is.na(unique(df$Circumference)))
    subtitle = paste(subtitle, 'C', df$Circumference, sep='')
  
  subtitle = paste(subtitle, " BE diagnosed ", unique(df$Initial.Endoscopy), sep='')
  
  df = df %>% mutate(Risk = factor(df$Risk, levels=c('Low','Moderate','High'), ordered = T), p53.Status = factor(p53.Status, levels=c(0,1), ordered = T))

  p = ggplot(df, aes(Endoscopy.Year, nl)) + geom_tile(aes(fill=Risk), color='white', size=1.5) + 
    scale_fill_manual(values=cols, limits=levels(df$Risk),name='Risk') 
  
  if (!pathology) {
    p = p + geom_point(aes(shape=Pathology), color='white', fill='white', size=5) + 
    scale_shape_manual(values=c(1,0,15,17,25), limits=levels(df$Pathology), guide=guide_legend(override.aes=list(fill='white', color='white'))) 
  } else {
    p = p + geom_text(aes(label=Pathology, color=Risk), show.legend=F) + scale_color_manual(values=c('white','black','white')) 
  }
    p = p + labs(y='Oesophageal Location',x='Endoscopies') + scale_x_continuous(breaks=scale) + 
        plot.theme + theme(axis.text.x=element_text(angle=45, hjust=1),legend.position='right',  legend.key = element_rect(fill='grey39', color=NA)) 
    if (add.title)
      p = p + labs(title=paste('Patient',unique(info$Patient)), subtitle=subtitle)
    
  return(p)   
}



min.theme = theme(text=element_text(size=12), panel.background=element_blank(), strip.background=element_rect(fill="grey97"), 
                  strip.text.y = element_text(size=6),
                  strip.text.x = element_text(size=9), 
                  axis.line=element_line(color='grey', size=0.1), panel.grid.major=element_blank(), 
                  panel.border = element_rect(color="grey", fill=NA, size=0.2), 
                  panel.spacing.x = unit(0, 'lines'), 
                  panel.spacing.y = unit(0.05, 'lines')   ) 


mtn.spatial<-function(df, b, limits, title='', pal=NULL, ylim=c(-5,5), axis.text.y=F) {
  df$mean.value = df$value
  
  if (is.null(ylim)) ylim = range(df$value)
  #print(ylim)
  
  if (is.null(pal))
    pal = RColorBrewer::brewer.pal(length(limits), 'Set2')
  
  df = left_join(df,  dplyr::select(b, 'Samplename','Endoscopy.Year','months.before.final','Block','ogj'), by=c('variable' = 'Samplename'))
  
  df$ogj = factor(df$ogj, ordered = T)
  df$ogj = factor(df$ogj, levels=rev(levels(df$ogj)), ordered = T)
  #pal = pal[which(limits %in% levels(df$nl))]
  
  altcolor = c(RColorBrewer::brewer.pal(9, 'Greys')[2],RColorBrewer::brewer.pal(9, 'Greys')[4])
  fills = c()
  for (yr in sort(unique(df$Endoscopy.Year))) {
    ac = (which(sort(unique(df$Endoscopy.Year)) == yr) %% 2) + 1
    fills = c(fills, rep(altcolor[ac], length(which(table(subset(df, Endoscopy.Year == yr)$ogj) > 0))))
  }
  
  blanklabel <- function(s) {
    return('')
  }
  
  p = ggplot(df, aes(x=chr.length)) + 
    facet_grid(Endoscopy.Year+ogj~chr, scales='free', space='free_x', labeller = labeller(.multi_line = F)) +
    geom_rect(aes(xmin=start, xmax=end, ymin=0, ymax=mean.value, fill=ogj), show.legend=T) + 
    scale_y_continuous(limits=ylim, breaks = c(0,ylim)) + 
    scale_fill_manual(values=pal, limits=limits,  name='OGJ Dist.') +
    labs(title=title, x='') + min.theme +
    theme(axis.line = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), strip.background=element_rect(fill="white", color=NA), 
          strip.text.x = element_text(size = 11),panel.spacing.y = unit(0.1, 'lines'), axis.text.y = element_text(size=6))    
  
  if (!axis.text.y)
    p = p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())    

  legend = get.legend(p)
  p = p + theme(legend.position = 'none')
  
  g = ggplot_gtable(ggplot_build(p))
  strips = grep('strip-r', g$layout$name)
  
  for (n in 1:length(strips)) {
    i = strips[n]
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[n]
  }
  
  return(list('grob'=g, 'legend'=legend))
}


plot.risk.by.path<-function(preds) {
  riskCols = BarrettsProgressionRisk:::riskColors()
  riskPath = preds %>% dplyr::group_by(Status, Pathology, Risk) %>% tally(name='Freq') %>% ungroup %>%
    mutate(Status = recode_factor(Status, NP = 'Non-Progressor', P = 'Progressor', .ordered = T))
  
  pathTotal = preds %>% dplyr::group_by(Status, Pathology) %>% tally(name='nPath') %>% ungroup %>%
    mutate(Status = recode_factor(Status, NP = 'Non-Progressor', P = 'Progressor', .ordered = T))
  
  riskPath = full_join(riskPath, pathTotal, by=c('Status','Pathology')) %>% 
    mutate(ratio = round(Freq/nPath, 3), Risk = factor(Risk, levels=c('Low','Moderate','High'), ordered = T))
  
  pathTotal = riskPath %>% dplyr::group_by(Status, Pathology) %>% dplyr::summarise('nPath'=unique(nPath))
  
  rp = ggplot(riskPath, aes(Pathology, ratio*100)) + geom_bar(aes(group=Risk, fill=Risk),stat='identity', position=position_stack()) + 
    scale_fill_manual(values=riskCols) + scale_color_manual(values=c('white','black','white')) +
    geom_text(data=pathTotal, aes(label=paste('n=',nPath,sep=''), x=Pathology, y=102), position=position_stack(), size=4) +
    geom_text(data=subset(riskPath, ratio>0), aes(group=Risk,label=paste(ratio*100,'%',sep=''), ), position=position_stack(vjust = 0.5), size=4, show.legend = F) +
    #geom_text(data=subset(riskPath, ratio>0), aes(group=Risk,label=paste(ratio*100,'%',sep=''), color=Risk), position=position_stack(vjust = 0.5), size=5, show.legend = F) +
    facet_grid(~Status, scales = 'free_x', space = 'free_x') + plot.theme + labs(title='Samples predicted by pathology', y='% Predicted', x='Pathology') + 
    theme(legend.position = 'bottom', text=element_text(size=14))
  
  return(rp)
}


transform.by.time<-function(pt) {
  pt = pt %>% dplyr::group_by(PID, Samplename) %>%  
    dplyr::mutate( 'months.before.final' = (Final.Endoscopy - Endoscopy.Year)*12) %>% 
    arrange(Pathology, PID, Endoscopy.Year) %>% ungroup 
  
  tb = with(pt, table(Endoscopy.Year, PID))
  tb[tb>0] = 1
  for (yr in names(which(rowSums(tb,na.rm=T) > 1))) {
    tmp = filter(pt, Endoscopy.Year == yr)
    tmp = tmp %>% dplyr::group_by(PID) %>% dplyr::summarise(max(Pathology))

    match = which(pt$Endoscopy.Year == yr & pt$PID == as.character(tmp[which.min(tmp$`max(Pathology)`), 'PID']))
    pt[match,'months.before.final'] = pt[match,'months.before.final'] + 6
  }
  pt$nl = (pt$Block-1)*2
  
  tb = table(pt[,c('months.before.final', 'PID')])
  if (nrow(tb) > 0) {
    for (i in 1:nrow(tb)) {
      pids = tb[i,]
      yrs = which(pids > 0)
      if (length(yrs) > 1) {
        for (j in 2:length(yrs)) 
          pt[which(pt$PID == names(yrs)[j]), 'nl'] = pt[which(pt$PID == names(yrs)[j]), 'nl'] + 1
      }
    }
  } else {
    pt$nl = 1
  }
  arrange(pt, -Endoscopy.Year, PID)
}
