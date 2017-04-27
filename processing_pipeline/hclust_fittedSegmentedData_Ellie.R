# Basically just plots per patient
library(ggplot2)
library(GGally)
library(plyr)
library(xlsx)
library(dplyr)

source('lib/load_patient_metadata.R')

data = '~/Data/Ellie'

data.files = list.files(paste(data, 'QDNAseq',sep='/'), full.names=T)
plot.dir = paste(data, 'Analysis/multipcf_plots_fitted_perPatient', sep='/')

if (length(list.files(plot.dir)) <= 0)
  stop(paste("No analysis files found in", plot.dir ))

## Patient info file
patient.file = grep('patient_info.xls', data.files, value=T)
if (length(patient.file) <= 0)
  stop(paste("Missing patient info file in", data))

patient.info = read.patient.info(patient.file)


patient.info$Patient = gsub("/", "_", patient.info$Patient)
head(patient.info)


#67 probes is equivalent to 1MB
min.probes=67

my_custom_cor <- function(data, mapping, color = I("grey10"), sizeRange = c(1, 7), ...) {
  # get the x and y data to use the other code
  x <- eval(mapping$x, data)
  y <- eval(mapping$y, data)
  
  ct <- cor.test(x,y)
  sig <- symnum(
    ct$p.value, corr = FALSE, na = FALSE,
    cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
    symbols = c("***", "**", "*", ".", " ")
  )
  
  r <- unname(abs(ct$estimate))
  rt <- format(r, digits=2)[1]
  
  # since we can't print it to get the strsize, just use the max size range
  cex <- max(sizeRange)
  
  # helper function to calculate a useable size
  percent_of_range <- function(percent, range) {
    percent * diff(range) + min(range, na.rm = TRUE)
  }
  
  sig_star_color = color
  if (ct$p.value <= 0.05)
    sig_star_color = I('firebrick')
  
  # plot the cor value
  p = ggally_text(
    label = as.character(rt), 
    mapping = aes(),
    xP = 0.5, yP = 0.5, 
    size = I(percent_of_range(cex * abs(r), sizeRange)),
    color = color,
    ... ) +
    
    # add the sig stars
    geom_text(
      aes_string( x = 0.8, y = 0.8 ),
      label = sig, 
      size = I( (percent_of_range(cex * abs(r), sizeRange)/2) ),
      #size = I(cex),
      color = sig_star_color,
      ... ) + 
    # remove all the background stuff and wrap it with a dashed line
    theme_classic() + 
    theme(
      panel.background = element_rect(
        color = color, 
        size = 2,
        linetype = "longdash"
      ), 
      axis.line = element_blank(), 
      axis.ticks = element_blank(), 
      axis.text.y = element_blank(), 
      axis.text.x = element_blank()
    )
}

panel.cor <- function(x, y, digits=2, prefix="", cex.cor) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, use= "pairwise.complete.obs"))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
  
  test <- cor.test(x,y, use= "pairwise.complete.obs")
  # borrowed from printCoefmat
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE,
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " "))
  
  text(0.5, 0.5, txt, cex = cex * r)
  text(.8, .8, Signif, cex=cex, col=2)
}

my_bars<-function(data, mapping, sizeRange=c(1,7),...) {
  # get the x and y data to use the other code
  x <- eval(mapping$x, data)
  y <- eval(mapping$y, data)
  
  name = as.character(mapping$x)
  ggplot(data=data,mapping=mapping) + geom_histogram(bins=10,col='green4',fill='palegreen',alpha=0.7) + 
    annotate("text",label=name, x=Inf, y=Inf, vjust=1,hjust=1,col=I('grey14'),fontface=2,size=max(sizeRange) )
}

panel.hist <- function(x, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}


for (patient.name in unique(patient.info$Patient)  )  {
  patient.plot.dir = paste(plot.dir, patient.name, sep='/')
  
  existing.plots = list.files(patient.plot.dir, '*.png', recursive=T )
  
  if (length(grep('_NormalizedLargeFittedSegments_correlationPlots_gamma', existing.plots, value=T)) == 7)
    next
  print(patient.name)
  
  for(gamma2 in c(250,5,10,25,50,100,500,1000)) { 
    print(gamma2)
    pt.gamma.plot.dir = paste(patient.plot.dir, paste('gamma2',gamma2,sep='_'),sep='/')
    print(pt.gamma.plot.dir)
    
    segvals = read.table(paste(patient.plot.dir, "/",patient.name,"_segmentedCoverage_fitted_gamma",gamma2,".txt",sep=""),sep="\t",stringsAsFactors=F,header=T)
    segvals = segvals[segvals$n.probes>=min.probes,]
    
    if(nrow(segvals)>0) {
      HC = hclust(dist(t(segvals[,-(1:5)])))
      
      png(paste(pt.gamma.plot.dir,"/hierarchicalClustering_fromFittedSegments_gamma",gamma2,".png",sep=""),width=1000)
      plot(HC, main=paste(patient.name,"clustering by fitted segment coverage"),xlab="",sub="")
      dev.off()

      # Normalize  (value-mean(value))/sd(value)
      normalised.segvals = segvals[,-(1:5)]
      for(c in 1:nrow(normalised.segvals)) {
        normalised.segvals[c,] = (normalised.segvals[c,]-mean(unlist(normalised.segvals[c,])))/sd(unlist(normalised.segvals[c,]))
      }
      
      HC = hclust(dist(t(normalised.segvals)))
      png(paste(pt.gamma.plot.dir,"/", patient.name,"_hierarchicalClustering_fromNormalizedLargeFittedSegments_gamma",gamma2,".png",sep=""),width=1000)
      plot(HC, main=paste(patient.name, "clustering by fitted segment coverage"),xlab="",sub="")
      dev.off()
      
      no.samples = ncol(normalised.segvals)
      no.rows = ceiling(no.samples/3)
      
      png(paste(pt.gamma.plot.dir,"/", patient.name,"_NormalizedLargeFittedSegments_histogram_gamma",gamma2,".png",sep=""),width=1500,height=600*no.rows)
      par(mfrow=c(no.rows,3))
      for(s in 1:no.samples) 
        hist(normalised.segvals[,s],col="blue",main=names(normalised.segvals)[s],xlab="normalised coverage")
      dev.off()
      
      if(nrow(segvals)>=3) {
        # order by date of endoscopy
        patient = subset(patient.info, Patient == patient.name)
        patient = arrange(patient, Endoscopy.Year, Pathology)
        
        ## TODO Mistake in earlier version of the patient file caused this will be fixed for the next run
        snames = grep('_1072(5|9)$', patient$Samplename)
        if ( length(snames > 0) ) {
          patient$Samplename[snames] =  paste( sub('-','_', patient$Plate.Index[snames]), '10725_10729', sep='_' )
        }
        
        sample.normalised.segvals = normalised.segvals[, intersect(colnames(normalised.segvals), patient$Samplename)]

        png(paste(pt.gamma.plot.dir,"/", patient.name,"_NormalizedLargeFittedSegments_correlationPlots_gamma",gamma2,".png",sep=""),width=600*no.rows,height=600*no.rows)			
        pairs = ggpairs( sample.normalised.segvals, 
                         lower=list(continuous=wrap("points", colour='darkblue', shape=20, alpha=0.3, size=3)),
                         diag=list(continuous=wrap(my_bars, sizeRange=c(1,5))), 
                         upper=list(continuous=my_custom_cor), 
                         axisLabels='show', title=patient.name, columnLabels=paste(paste(patient$Endoscopy.Year, ' (',patient$Pathology, ')', sep=''), patient$Plate.Index, sep='\n') )
        
        print(pairs  + theme( # this doesn't work...
          text = element_text(size=20, face='bold')
        )) 
        dev.off()
      }
    }
  }
}
