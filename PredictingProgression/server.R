#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/

library(shiny)
library(ggplot2)
library(ggrepel)
library(glmnet)


source('lib/data_func.R')

chr.info = read.table('data/chr_hg19.txt', sep='\t', header=T)


load('data/model_data.Rdata', verbose=T)
rm(dysplasia.df, labels)

load('data/all.pt.alpha.Rdata', verbose=T)
select.alpha = '0.9'
fitV = models[[select.alpha]]
lambda.opt = performance.at.1se[[select.alpha]][, 'lambda']
rm(plots, coefs, dysplasia.df, labels, models, performance.at.1se)


load('data/loo.Rdata', verbose=T)
rm(plots, coefs, nzcoefs, fits)

pg.samp = do.call(rbind, pg.samp)
pg.samp$Hospital.Research.ID = NULL

#hist(pg.samp$OR)

riskPal = rev(RColorBrewer::brewer.pal(11, 'RdYlBu'))
unadjRR = ggplot(pg.samp, aes(OR)) + geom_histogram(aes(fill=..x..), bins=10, show.legend = T) +
  scale_fill_gradientn(colors = riskPal,  name='RR') + 
  labs(y='n Samples', x='Relative Risk', title='Unadjusted relative risk') + theme_light(base_size = 14)


probs = ggplot(pg.samp, aes(Prediction)) + geom_histogram(aes(fill=..x..), bins=10, show.legend = T) +
  scale_fill_gradientn(colors = riskPal,  name='RR') + 
  labs(y='n Samples', x='P(Progression)', title='Predictions Probabilities') + theme_light(base_size = 14)


# p = ggplot(df, aes(x=chr.length)) + facet_grid(Endoscopy.Year+ogj~chr, scales='free', space='free_x', labeller = labeller(.multi_line = F)) +
#   geom_rect(aes(xmin=start, xmax=end, ymin=0, ymax=mean.value, fill=ogj), show.legend=T) + 
#   scale_fill_manual(values=pal, limits=limits,  name='OGJ Dist.') +
#   labs(title=title, x='') + min.theme  


# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
  ## how to correctly deal with actual CN calls?  e.g. Total 2, Minor 1 etc
  
  predict.progression <- function(file, model.norm=T) {
    #rs <- read.table(file, header=T, sep="\t")
    rs <- as.data.frame(fread(file, header=T))
    if (ncol(rs) < 4)
      stop("File needs to contain 4 columns: chr, start, end, sample_value")
    #head(rs)
    
    cnCol = grep('Total_CN', colnames(rs))
    if (length(col) >= 1)
      rs = rs[,c(1:cnCol)]
    
    segs <- tile.segmented.data(rs, size=5e6 )
    segM = as.matrix(segs[,-c(1:3)])
    rownames(segM) = paste(segs$chr, ':', segs$start, '-', segs$end, sep='')
    segM = t(segM)
    
    if (model.norm) {
      for (i in 1:ncol(segM)) 
        segM[,i] = unit.var(segM[,i], z.mean[i], z.sd[i])
    } else {
      segM = unit.var(segM)
    }
    
    arms <- tile.segmented.data(rs, size='arms')
    armsM = as.matrix(arms[,-c(1:3)])
    rownames(armsM) = paste(arms$chr, ':', arms$start, '-', arms$end, sep='')
    armsM = t(armsM)
    
    if (model.norm) {
      for (i in 1:ncol(armsM)) 
        armsM[,i] = unit.var(armsM[,i], z.arms.mean[i], z.arms.sd[i])
    } else {
      armsM = unit.var(armsM)
    }
    
    
    #nrow(armsM) == nrow(segM)
    
    cx.score = score.cx(segM,1)
    
    mergedDf = subtract.arms(segM, armsM)
    mergedDf = cbind(mergedDf, 'cx' = unit.var(cx.score, mn.cx, sd.cx))
    
    # Predict function giving me difficulty when I have only a single sample, this ensures the dimensions are the same
    sparsed_test_data <- Matrix(data=0, nrow=nrow(mergedDf),  ncol=ncol(mergedDf),
                                dimnames=list(rownames(mergedDf),colnames(mergedDf)), sparse=T)
    for(i in colnames(mergedDf)) sparsed_test_data[,i] = mergedDf[,i]
    
    preds = predict(fitV, newx=sparsed_test_data, s=lambda.opt, type='response')
    RR = predict(fitV, newx=sparsed_test_data, s=lambda.opt, type='link')
    
    return(list('prob'=preds, 'rel.risk'=RR))
  }
  
  output$contents <- renderTable({
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, it will be a data frame with 'name',
    # 'size', 'type', and 'datapath' columns. The 'datapath'
    # column will contain the local filenames where the data can
    # be found.
    req(input$file)
    
    inFile <- input$file
    
    df <- data.table::fread(input$file$datapath,
                            header = input$header)
                   #sep = input$sep,
                   #quote = input$quote)
    
    
    return(head(df))
    
    #pp = predict.progression(inFile$datapath)
    
    #pp$prob
    
  })     
   output$plot <- renderPlot({
     inFile <- input$file
     
     norm = input$norm
     
    if (is.null(inFile))  
      return(unadjRR)
    
    pp = predict.progression(inFile$datapath, norm)
      
    RR = as.data.frame(pp$rel.risk)
    colnames(RR) = 'rel.risk'

    pR = as.data.frame(round(pp$prob, 3))
    colnames(pR) = 'Probability'
    
    
    plot = unadjRR + geom_point(data=RR, aes(x=rel.risk, y=3), size=2) +
      geom_text_repel(data=RR, aes(x=rel.risk, y=3, label=round(rel.risk, 2)))
    
    gridExtra::grid.arrange(plot, gridExtra::tableGrob(pR, theme=ttheme_minimal()), nrow=2, heights=c(4,1))

   })
  
   
   
})
