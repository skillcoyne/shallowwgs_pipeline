#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/

library(shiny)
library(ggplot2)
library(glmnet)

source('lib/data_func.R')

chr.info = read.table('data/chr_hg19.txt', sep='\t', header=T)


load('data/model_data.Rdata', verbose=T)
rm(dysplasia.df, labels)

load('data/all.pt.alpha.Rdata', verbose=F)
select.alpha = '0.9'
fitV = models[[select.alpha]]
lambda.opt = performance.at.1se[[select.alpha]][, 'lambda']
rm(plots, coefs, dysplasia.df, cv.patient, labels, models, performance.at.1se)


load('data/loo.Rdata', verbose=T)
rm(plots, coefs, nzcoefs, fits)

pg.samp = do.call(rbind, pg.samp)
pg.samp$Hospital.Research.ID = NULL

#hist(pg.samp$OR)

riskPal = rev(RColorBrewer::brewer.pal(11, 'RdYlBu'))
unadjRR = ggplot(pg.samp, aes(OR)) + geom_histogram(aes(fill=..x..), bins=10, show.legend = T) +
  scale_fill_gradientn(colors = riskPal,  name='RR') + 
  labs(y='n Samples', x='Relative Risk', title='Unadjusted relative risk') + theme_light(base_size = 14)



# p = ggplot(df, aes(x=chr.length)) + facet_grid(Endoscopy.Year+ogj~chr, scales='free', space='free_x', labeller = labeller(.multi_line = F)) +
#   geom_rect(aes(xmin=start, xmax=end, ymin=0, ymax=mean.value, fill=ogj), show.legend=T) + 
#   scale_fill_manual(values=pal, limits=limits,  name='OGJ Dist.') +
#   labs(title=title, x='') + min.theme  


# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
  output$contents <- renderTable({
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, it will be a data frame with 'name',
    # 'size', 'type', and 'datapath' columns. The 'datapath'
    # column will contain the local filenames where the data can
    # be found.
    inFile <- input$file
    
    if (is.null(inFile))
      return(NULL)
    
    rs <- read.table(inFile$datapath, header=T, sep="\t")
    if (ncol(rs) < 4)
      stop("File needs to contain 4 columns: chr, start, end, sample_value")
    
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    
    progress$set(message = "Setting up data...", value = 0)
    
    n <- 20
    updateProgress <- function(detail = NULL) {
      progress$inc(amount = 1/n, detail = detail)
    }
    
    segs <- tile.segmented.data(rs, size=5e6 )
    segM = as.matrix(segs[,-c(1:3)])
    rownames(segM) = paste(segs$chr, ':', segs$start, '-', segs$end, sep='')
    segM = t(segM)

    for (i in 1:ncol(segM)) 
      segM[,i] = unit.var(segM[,i], z.mean[i], z.sd[i])
    

    arms <- tile.segmented.data(rs, size='arms')
    armsM = as.matrix(arms[,-c(1:3)])
    rownames(armsM) = paste(arms$chr, ':', arms$start, '-', arms$end, sep='')
    armsM = t(armsM)

    for (i in 1:ncol(armsM)) 
      armsM[,i] = unit.var(armsM[,i], z.arms.mean[i], z.arms.sd[i])
    
    nrow(armsM) == nrow(segM)
        
    cx.score = score.cx(segM,1)

    mergedDf = subtract.arms(segM, armsM)
    mergedDf = cbind(mergedDf, 'cx' = unit.var(cx.score, mn.cx, sd.cx))
    
    # Predict function giving me difficulty when I have only a single sample, this ensures the dimensions are the same
    sparsed_test_data <- Matrix(data=0, nrow=nrow(mergedDf),  ncol=ncol(mergedDf),
                                dimnames=list(rownames(mergedDf),colnames(mergedDf)), sparse=T)
    for(i in colnames(mergedDf)) sparsed_test_data[,i] = mergedDf[,i]

    preds = predict(fitV, newx=sparsed_test_data, s=lambda.opt, type='response')
    RR = predict(fitV, newx=sparsed_test_data, s=lambda.opt, type='link')
    
    preds
    
  })     
   output$plot <- renderPlot({
  
    #if (is.null(inFile))  
    #  return(unadjRR)
     
    #unadjRR + geom_point(data=as.data.frame(RR), aes(x=`1`, y=3), size=2)
  
  #   
  #   # generate bins based on input$bins from ui.R
  #   x    <- faithful[, 2] 
  #   bins <- seq(min(x), max(x), length.out = input$bins + 1)
  #   
  #   # draw the histogram with the specified number of bins
  #   #hist(x, breaks = bins, col = 'darkgray', border = 'white')
  #   
  #   ggplot(faithful, aes(x=waiting)) + geom_histogram(breaks=bins, color='lightgrey') + theme_minimal()
  #   
   })
  
})
