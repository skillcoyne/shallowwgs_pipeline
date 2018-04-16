
library(shiny)
library(ggplot2)
library(ggrepel)
library(glmnet)


source('lib/data_func.R')
source('lib/risk_helper.R')

adjustRisk <- function(RR, offset, type='risk') {
  if (type == 'prob') {
    x = 1/(1+exp(-RR+abs(offset)))
  } else {
    x = RR+offset
  }
  return(x)
}

loadModelData <- function() {
  message("Loading model data.")
  chr.info = read.table('data/chr_hg19.txt', sep='\t', header=T)
  
  load('data/model_data.Rdata', verbose=T)
  rm(dysplasia.df, allDf, labels)
  
  load('data/all.pt.alpha.Rdata', verbose=T)
  select.alpha = '0.9'
  fitV = models[[select.alpha]]
  lambda.opt = performance.at.1se[[select.alpha]][, 'lambda']
  rm(plots, coefs, dysplasia.df, labels, models, performance.at.1se)
  
  load('data/loo.Rdata', verbose=T)
  rm(plots, coefs, nzcoefs, fits)
  
  pg.samp = do.call(rbind, pg.samp)
  pg.samp$Hospital.Research.ID = NULL
  
  status = unique(pg.samp[,c('Patient','Status')])$Status
  
  cases = table(status)['P']
  
  mn = round((table(status)['P']/(0.01*100))*100)
  m = round((table(status)['P']/(0.0225*100))*100)
  mx = round((table(status)['P']/(0.035*100))*100)
  
  offsetMin = log(cases/mn)
  offsetMean = log(cases/m)
  offsetMax = log(cases/mx)
  
  pg.samp = pg.samp %>% mutate(  
    adjRRmin =  adjustRisk(OR, offsetMin),
    adjRRmax = adjustRisk(OR, offsetMax),
    adjRR =  adjustRisk(OR, offsetMean),
    
    adjProb = adjustRisk(OR, offsetMean, 'prob'),
    adjProbMin = adjustRisk(OR, offsetMin, 'prob'),
    adjProbMax = adjustRisk(OR, offsetMax, 'prob'),
  )
  
  riskPal = rev(RColorBrewer::brewer.pal(11, 'RdYlBu'))
  adjRRplot = ggplot(pg.samp, aes(adjRR)) + geom_histogram(aes(fill=..x..), bins=10, show.legend = T) +
    scale_fill_gradientn(colors = riskPal,  name='RR') + 
    labs(y='n Samples', x='Relative Risk', title='Adjusted relative risk') + theme_light(base_size = 14)
  
  probs = ggplot(pg.samp, aes(adjProb)) + geom_histogram(aes(fill=..x..), bins=10, show.legend = T) +
    scale_fill_gradientn(colors = riskPal,  name='RR') + 
    labs(y='n Samples', x='P(Progression)', title='Adjusted Predictions Probabilities') + theme_light(base_size = 14)
  
  setOldClass("glmnet"); setOldClass("lognet")
  beModel = setClass('BEModel', representation('mn.cx'="numeric", 'sd.cx'="numeric", 'z.mean'="numeric",
                                               'z.sd'="numeric",'z.arms.mean'="numeric",'z.arms.sd'='numeric',
                                               'minOffset'="numeric", 'maxOffset'="numeric", 'meanOffset'="numeric"))
  beModel = new("BEModel", 'mn.cx'=mn.cx, 'sd.cx'=sd.cx, 'z.mean'=z.mean,'z.sd'=z.sd,'z.arms.mean'=z.arms.mean, 'z.arms.sd'=z.arms.sd, 
                'minOffset'=offsetMin, 'maxOffset'=offsetMax, 'meanOffset'=offsetMean)
  
  message("Model data loading completed.")
  return(list('be'=beModel, 'fit'=fitV, 'lambda.opt'=lambda.opt, 'plot'=adjRRplot))
}

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
  modelData = loadModelData()
  
  readDataFile <- reactive({
    req(input$file)
    dataFile <- input$file$datapath
    message(paste("Reading data file", dataFile))
    print(dataFile)
    df <- read.table(dataFile, header=T, sep='\t')
    samplenames = grep('chr|arm|start|end|probes', colnames(df), ignore.case=T, invert=T, value=T)
    return(list('data'=df, 'samples'=samplenames))
  })
  
  perSampleMeta <- reactive({
    if (!is.null(input$pathfile)) {
      fi = readDataFile()
      samplenames = fi$samples
      
      if (!is.null(input$pathfile)) {
        pf = read.table(input$pathfile$datapath, sep='\t', header=T, stringsAsFactors=F)
        if (length(samplenames) < nrow(pf)) {
          warning("Samples do not match in data and path files, ignoring path.")
          pf = NULL
        } else {
          print(input$sampleMeta)
        }
      }
    } else {
      pf = NULL
    }
    return(pf)
  })
  
  generatePred <- reactive({
    message("Generating predictions")
    fi = readDataFile()
    df = fi$data

    pf = perSampleMeta()
    print(paste("per sample meta:", is.null(pf)))
    pp = predict.progression(df, modelData$be, modelData$fit, modelData$lambda.opt, T, T)

    samplenames = rownames(pp$rel.risk)

    pR = cbind.data.frame(samplenames, sapply(pp$rel.risk, adjustRisk, offset=modelData$be@meanOffset, type='prob'), 
                            sapply(pp$rel.risk, adjustRisk, offset=modelData$be@meanOffset))
    colnames(pR) = c('Sample', 'Probability', 'Relative Risk')
    rownames(pR) = 1:nrow(pR)
      
    pR$Risk = sapply(pR$Probability, risk)
    
    recommendations = rx(pR)

    if (!is.null(pf)) {
      pR = merge(pf, pR, by='Sample')
      recommendations = rx(pR)
    }

    #recommendations = recommendations[complete.cases(recommendations),]
    return(list('pR'=pR, 'rx'=recommendations))
  })

  output$samples <- renderUI({
    fi = readDataFile()
    
    if (is.null(input$pathfile))  {
      samplenames = fi$samples

      myTabs = lapply(samplenames, function(x) {
        tabPanel(paste('Sample',x,sep=":"), 
                 radioButtons('path','Sample pathology', 
                              choices=list('Non-dysplastic BE'='ndbe','Indeterminate'='id','Low grade dysplasia'='lgd','High grade dysplasia'='hgd','Intermucosal carcinoma'='imc'), selected='ndbe'),
                 radioButtons('p53','p53 IHC',
                              choices=list('Not performed'=NA,'Normal'=0,'Aberrant'=1)))
      })
      do.call(tabsetPanel, c('sampleMeta',myTabs, type='pills'))
    }
    
  })

  ## how to correctly deal with actual CN calls?  e.g. Total 2, Minor 1 etc
  output$example <- renderTable({
    df <- read.table('data/example_pt_segvals.txt', header=T, sep='\t')
    return(head(df))
  })

  output$examplep53 <- renderTable({
    df <- read.table('data/path_p53_example.txt', header=T, sep='\t')
    return(head(df))
  })
  
  output$contents <- renderTable({
    req(input$file)
    inFile <- input$file
      
    df <- read.table(input$file$datapath, header=input$header, sep='\t')

    return(head(df))
  })     
  
  output$plot <- renderPlot({
    di = generatePred()
    pR = di$pR
    
    plot = modelData$plot + geom_point(data=pR, aes(x=`Relative Risk`, y=3), size=2) +
      geom_text_repel(data=pR, aes(x=`Relative Risk`, y=3, label=round(`Relative Risk`, 2)))

    riskCols = RColorBrewer::brewer.pal(11, "RdYlBu")[seq(1,11, 3)]
    #RColorBrewer::display.brewer.pal(11,'RdYlBu')

    #gridExtra::grid.arrange(table)
    
    gridExtra::grid.arrange(plot)

   })
  
  output$riskTable <- renderPlot({
    di = generatePred()
    rx = di$rx
    
    riskCols = RColorBrewer::brewer.pal(11, "RdYlBu")[seq(1,11, 3)]
    fonts = c('bold.italic','italic','plain','plain' )
    
    tt3 <- gridExtra::ttheme_minimal(
      core=list(bg_params = list(fill=riskCols[rx$rule], col=NA),
                fg_params=list(fontface=fonts[rx$rule]), col='white'),
      colhead=list(fg_params=list(col="black", fontface=4L)))
    
    table = gridExtra::tableGrob(format(rx[,c('Time 1', 'Time 2', 'Rx')], core.just='left'),theme=tt3  )
    gridExtra::grid.arrange(table, top='Recommendations per pair of samples')
  }) 
  
   
})
