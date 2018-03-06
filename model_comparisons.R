library(ggplot2)

hazards = apply( exp(t(dysplasia.df[, rownames(coefs[[select.alpha]])]) *  coefs[[select.alpha]][,1]), 1, function(x) {
  sd(x)/mean(x)
})

calc.haz <- function(df, coefs) {
  apply( exp(t(df[, rownames(coefs)]) *  coefs[,1]), 1, function(x)  sd(x)/mean(x) ) 
}

load.runs<-function(dir) {
  file = paste(dir, 'all.pt.alpha.Rdata', sep='/')
  message(paste("loading", file))
  load(file, verbose=T)
  
  cfs = cbind(
    'coef'=coefs[['0.9']][,1], 
    'stability'=rowSums(coefs[['0.9']][,-1])/50, 
    'hazard'=calc.haz(dysplasia.df, coefs[['0.9']])
    )
  return(list('coefs'=cfs, 'df'=dysplasia.df, 'performance'=performance.at.1se[['0.9']], 'model'=models[['0.9']], 'plots'=plots[['0.9']]))
}


data = '~/Data/Ellie/Analysis'

pcols = c('mean','sme')
prows = c('genomic','genomic.cx', 'genomic.demo', 'genomic.cx.demo','cx','demo')
perf.df = data.frame(matrix(ncol=length(pcols), nrow=length(prows), dimnames=list(prows, pcols)))

genomic = load.runs(paste(data, '5e06-genomic', sep='/')) # genomic only
perf.df['genomic',] = genomic$performance[pcols]

cx.genomic = load.runs(paste(data, '5e06-with_cx', sep='/'))
perf.df['genomic.cx',] = cx.genomic$performance[pcols]

demo.genomic = load.runs(paste(data, '5e06-with_demo', sep='/'))
perf.df['genomic.demo',] = demo.genomic$performance[pcols]

cx.demo.gen = load.runs(paste(data, '5e06-with_cx_demo', sep='/'))
perf.df['genomic.cx.demo',] = cx.demo.gen$performance[pcols]

demo.only = load.runs(paste(data, '5e06-demo', sep='/'))
perf.df['demo',] = demo.only$performance[pcols]

cx.only = load.runs(paste(data, '5e06-cx', sep='/'))
perf.df['cx',] = cx.only$performance[pcols]

perf.df$data = prows

ggplot(perf.df, aes(data, mean, group=data, fill=data)) + geom_col() + geom_errorbar(aes(ymin=mean-sme, ymax=mean+sme), color='darkgrey', width=0.5) +
  geom_text(aes(label=round(mean, 3)), nudge_y = -0.03) + labs(x='', y='Classification', title='GLM dataset comparison')





cx.demo.gen$performance[pcols]
cx.genomic$performance[pcols]
demo.genomic$performance[pcols]

