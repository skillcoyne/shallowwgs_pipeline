

get.legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

plot.theme = theme(text=element_text(size=14, face='bold'), panel.background=element_blank(), strip.background =element_rect(fill="lightcyan2"),  
                   strip.text = element_text(size=14, face='bold'), 
                   axis.line=element_line(color='black'), panel.grid.major=element_line(color='grey90'),
                   panel.border = element_rect(color="grey", fill=NA, size=0.5), panel.spacing = unit(0.1, 'lines')  ) 

cn.mtn.plot<-function(df, label, type='bar') {
  low = c( min(df$mean.value), quantile(subset(df, mean.value < -1)$mean.value, probs=c(0.25, 0.75)), -1)
  high = c(1, quantile(subset(df, mean.value > 1)$mean.value, probs=c(0.25, 0.75)), max(df$mean.value))
  df$cuts = cut(df$mean.value, breaks=c(low,high), include.lowest = T, ordered_result = T)
  
  p = ggplot(df, aes(x=chr.length)) + facet_grid(Status~chr, scales='free_x', space='free_x', labeller=label) 
  
  if (type == 'seg') {
    riskCols = brewer.pal(10, "PRGn")[c(1, 5, 10)]
    p = p + geom_segment(aes(x=start, xend=end, y=mean.value, yend=mean.value, color=cn), show.legend=T, size=1) + scale_color_manual(values=riskCols) 
  }
  if (type == 'bar') {
    highCol = rev(brewer.pal(6, "Purples"))[c(3,1,3)]
    lowCol = rev(rev(brewer.pal(6, "Greens"))[c(3,1,3)])
    grey = brewer.pal(3,'Greys')[1]
    
    p = p + 
      geom_rect(aes(xmin=start, xmax=end, ymin=0, ymax=mean.value, fill=cuts), show.legend=T ) + 
      scale_fill_manual(values=c(lowCol,grey,highCol),labels=c('<25%','25-75%','75%','0','<25%','25-75%','>75%'), limits=levels(df$cuts), name='') 
  }
  p + labs(x='', y='Mean adjusted segmentation value', title='All samples') +
    theme(text=element_text(face='bold', size=14), axis.text.x=element_blank(), legend.position='right', 
          panel.grid = element_blank(), panel.background = element_blank(),strip.background =element_rect(fill="whitesmoke"),  
          panel.border = element_rect(color="grey", fill=NA, size=0.5), panel.spacing.x = unit(0, 'lines') )
}

cn = function(x) {
  str = 'neutral'
  if(x <= -1) {
    str = 'loss'
  } else if (x >= 1) {
    str = 'gain'
  }
  return(str)
}

model.performance<-function(p1se, coeffs, lims=NULL) {
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
  
  p = ggplot(performance, aes(alpha, mean, color=(alpha != 0.9))) + 
    geom_point(size=2) +
    geom_text( aes(label=round(mean, 3)), hjust=-0.5) +
    geom_errorbar(aes(ymin=mean-sme, ymax=mean+sme), width=0.3, size=1) + 
    geom_text(aes(label=n.Coef), nudge_y=0.03, fontface='bold') +
    geom_text(aes(label=paste('(',`75%`,')',sep='')), nudge_y=0.02) + 
    scale_colour_manual(values=c('red', 'grey39')) + plot.theme + theme(legend.position = 'none') +
    labs(x='Elasticnet penalty value: ridge <-> lasso', y='Model classification at lambda-1se')
  
  
  if (!is.null(lims))
    p = p+ylim(lims)
  
  cfs = as.data.frame(matrix(data=0,nrow=length(unique(names(table(unlist(lapply(coef.stable, names)))))), ncol=length(names(coef.stable)), 
                             dimnames=list( unique(names(table(unlist(lapply(coef.stable, names))))), names(coef.stable) )))
  for (i in names(coef.stable)) 
    cfs[intersect(rownames(cfs), names(coef.stable[[i]])), i] = 1
  return(list('performance'=performance, 'stable.coeffients'=coef.stable, 'cf.per.feature'=cfs,'plot'=p))
}

roc.plot <- function(roc,title="") {
  df = cbind.data.frame('specificity'=rev(roc$specificities), 'sensitivity'=rev(roc$sensitivities))

  ggplot(df, aes(specificity, sensitivity)) + geom_line() + scale_x_reverse() + 
    geom_label(data=as.data.frame(t(round(pROC::coords(roc, "best"),2))), 
               aes(x=specificity, y=sensitivity, 
                   label=paste(' (FPR ', round((1-specificity),2)*100, '%, TPR ', round(sensitivity,2)*100,'%)\nAUC:',round(roc$auc, 2)*100, '%', sep='')), nudge_x=.25) +
    labs(title=title, x='Specificity (FPR)', y='Sensitivity (TPR)')  + plot.theme
}

multi.roc.plot <- function(rocList, title="ROC", palette=NULL, colors=NULL) {
  aucs = do.call(rbind, lapply(rocList,function(r) 
    cbind.data.frame('AUC'=r$auc*100,'model'=r$model,'Specificity'=pROC::coords(r,'best')[['specificity']], 'Sensitivity'=pROC::coords(r,'best')[['sensitivity']])))
  aucs$y = aucs$Sensitivity
  aucs$x = aucs$Specificity

  if (nrow(aucs) > 3) {
    aucs$y = seq(0.1,1, 1/nrow(aucs))
    aucs$x = 0.5
  }
  
  aucs$label=with(aucs, paste(model, '\nAUC ', round(AUC), '%\nTPR=', round(Sensitivity*100), '% FPR=', round((1-Specificity)*100), '%', sep=''))
  
  df = do.call(rbind, lapply(rocList, function(r) 
    cbind.data.frame('Specificity'=rev(r$specificities), 'Sensitivity'=rev(r$sensitivities),'model'=r$model) 
  ))

  p = ggplot(df, aes(Specificity, Sensitivity,color=model,fill=model)) +
        geom_line(show.legend = F, size=1.5) + scale_x_reverse() +
        geom_label(data=aucs, aes(x=x, y=y, label=label), color='white', show.legend = F) +
        labs(title=title, x='Specificity (FPR)', y='Sensitivity (TPR)')  + plot.theme
  if (is.null(colors) & !is.null(palette)) {
    colors = RColorBrewer::brewer.pal(length(rocList)+1, palette)
  } else if (is.null(colors) & is.null(palette)) {
    colors = RColorBrewer::brewer.pal(length(rocList)+1,'Set1')
  }
  p = p + scale_color_manual(values=colors) + scale_fill_manual(values=colors)
 
  return(p)
}

pred.tiles <- function(df, probs=F, OR=F, cols=c('green3','orange','red3'), ...) {
  p = ggplot(df, aes(brks, nl))

  if (OR) {
    p = p + geom_tile(aes(fill=OR), color='white') + scale_fill_gradientn(colors = myPal, limits=c(-9.3,15), name='RR') 
  } else if (probs) {
    p = p + geom_tile(aes(fill=Prediction), color='white') + scale_fill_gradientn(colors = myPal, limits=c(0,1), name='P(P)') 
  } else {
    p = p + geom_tile(aes(fill=Risk)) + 
      scale_fill_manual(values=cols, limits=levels(df$Risk), name='Progression')
  }
  p + 
    geom_point(aes(shape=Pathology, color=Risk), fill='white', size=3) + scale_color_manual(values=c('white','black','white'), guide=F) +
    scale_shape_manual(values=c(1,0,15,24,25), limits=levels(df$Pathology), guide=guide_legend(override.aes=list(fill='black', color='black'))) +
    #scale_fill_manual(values=cols, limits=levels(df$Risk),name='Progression') +  #scale_color_manual(values=cols, limits=levels(df$Risk),name='Progression') +
    #geom_label(aes(label=Pathology), fill='white', show.legend=F) + 
    #scale_color_manual(values=c('white','black','white'), limits=levels(df$Risk)) +
    facet_wrap(~Patient, scales='free', labeller=label_both, ...) +
    labs(y='Oesophageal Location',x='Months Prior to Endpoint', title='') +
    plot.theme + theme(axis.text.x=element_text(angle=45, hjust=1)) 
}

plot.cn.chr<-function(df, chrom='17', info=NULL, haz=NULL, genes=NULL) {
  low = c( min(df$mean.value), quantile(subset(df, mean.value < -1)$mean.value, probs=c(0.25, 0.75)), -1)
  high = c(1, quantile(subset(df, mean.value > 1)$mean.value, probs=c(0.25, 0.75)), max(df$mean.value))
  
  df$cuts = cut(df$mean.value, breaks=c(low,high), include.lowest = T, ordered_result = T)
  
  highCol = rev(brewer.pal(6, "Purples"))[c(3,1,3)]
  lowCol = rev(rev(brewer.pal(6, "Greens"))[c(3,1,3)])
  cols = c(lowCol,'#F7F7F7', highCol)
  
  df = subset(df, chr == chrom) # & seg.type != 'arms')  
  p = ggplot(df, aes(x=chr.length)) + facet_grid(Status~chr, scales='free_x', labeller=labeller(chr=label_both)) +
    geom_rect(aes(xmin=start, xmax=end, ymin=0, ymax=mean.value, fill=cuts),show.legend = T) + 
    #geom_segment(aes(x=start, xend=end, y=mean.value, yend=mean.value, color=cuts), show.legend=T, size=1) + 
    #geom_segment(data=subset(df, seg.type == 'arms'), aes(x=start, xend=end, y=mean.value, yend=mean.value, color=cn), alpha=0.8, show.legend = F) + 
    scale_color_manual(values=cols) + scale_fill_manual(values=cols,labels=c('<25%','25-75%','75%','0','<25%','25-75%','>75%'), name='') 
  
  if (!is.null(info)) {
    info = info[which(info$chrom == chrom),]
    p = p + geom_vline(data=info, aes(xintercept=chr.cent-cent.gap), color='grey39', linetype='longdash' ) +
      geom_label(data=info, aes(x=chr.cent-cent.gap*3, y=max(df$mean.value)), label='p-arm')
  }
  
  if (!is.null(genes)) {
    genes = subset(genes, chr == chrom)
    if (nrow(genes) > 0) p = p + geom_point(data=genes, aes(x=start,  y=value), color='darkblue', size=3) +
        geom_text_repel(data=genes, aes(x=start, y=value, label=Gene.Symbol), color='darkblue') 
  }
  if (!is.null(haz)) {
    haz = subset(haz, chr == chrom)
    if (nrow(haz) > 0 & length(grep('loss',unique(haz$gl))) > 0 ) {
      p = p + geom_rect(data=subset(haz, gl == 'loss'), aes(xmin=start, xmax=end, ymin=0, ymax=max), color='blue3', alpha=0, show.legend=F, size=1) 
    } 
    if (nrow(haz) > 0 & length(grep('gain',unique(haz$gl))) > 0) {
      p = p + geom_rect(data=subset(haz, gl == 'gain'), aes(xmin=start, xmax=end, ymin=0, ymax=max), color='goldenrod2', alpha=0, show.legend=F, size=1) 
    }
  }
  
  p + theme(text=element_text(face='bold', size=14), axis.text.x=element_blank(), legend.position='bottom', 
            panel.grid = element_blank(), panel.background = element_blank(),strip.background =element_rect(fill="lightcyan2"),  
            panel.border = element_rect(color="grey", fill=NA, size=0.5), panel.spacing.x = unit(0.1, 'lines') ) + labs(x='', y='Mean adjusted segmentation value')
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

vig.plot<-function(df, info, cols=c('green3','orange','red3'), pathology=F) {
  scale = seq(min(info$Endoscopy.Year), max(info$Endoscopy.Year), 1)
  subtitle = with(info, paste(unique(Age.at.diagnosis),unique(Sex),': ', sep=''))
  if (!is.na(unique(info$Smoking)))                             
    subtitle = paste(subtitle, ifelse(unique(info$Smoking) == 'Y', 'smoker, ','non-smoker, '),sep='')
  
  subtitle = paste(subtitle, 'M', unique(info$Maximal), sep='')
  
  if (!is.na(unique(info$Circumference)))
    subtitle = paste(subtitle, 'C', info$Circumference, sep='')
  
  subtitle = paste(subtitle, " BE diagnosed ", unique(info$Initial.Endoscopy), sep='')
  
  df$Risk = factor(df$Risk, levels=c('Low','Moderate','High'), ordered = T)
  df$p53.Status = factor(df$p53.Status, levels=c(0,1), ordered = T)

  p = ggplot(df, aes(Endoscopy.Year, nl)) + 
    geom_tile(aes(fill=Risk), color='white', size=1.5) + scale_fill_manual(values=cols, limits=levels(df$Risk),name='Progression') 
  
  if (!pathology) {
    p = p + geom_point(aes(shape=Pathology, color=Risk), fill='white', size=5) + scale_color_manual(values=c('white','black','white'), guide=F) + 
    scale_shape_manual(values=c(1,0,15,24,25), limits=levels(df$Pathology), guide=guide_legend(override.aes=list(fill='black', color='black'))) 
    #scale_shape_manual(values=c('B','I','L','H','C'), limits=levels(df$Pathology), guide=guide_legend(override.aes=list(fill='black', color='black'))) +
  } else {
    p = p + geom_text(aes(label=Pathology, color=Risk), show.legend=F) + scale_color_manual(values=c('white','black','white')) 
  }
    p = p + labs(y='Oesophageal Location',x='Months Prior to Endpoint', title=paste('Patient',unique(info$Patient)), subtitle=subtitle) + scale_x_continuous(breaks=scale) + 
        plot.theme + theme(axis.text.x=element_text(angle=45, hjust=1),legend.position='right') 
    
}



