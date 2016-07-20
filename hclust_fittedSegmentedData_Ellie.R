# Basically just plots per patient
library(ggplot2)
library(GGally)

#setwd("/lustre/scratch112/sanger/cgppipe/PanCancerDownloads/workspace/dw9/Ellie")

#just use samples from one patient
patient.info = read.table("data/All_patient_info.txt",sep="\t",header=T,stringsAsFactors=F)
patient.info$SLX_ID[patient.info$SLX_ID=="SLX-10729"] = "SLX-10725_10729"
patient.info$SLX_ID[patient.info$SLX_ID=="SLX-10725"] = "SLX-10725_10729"
patient.info$samplename = paste(patient.info$Index,gsub("SLX-","",patient.info$SLX_ID),sep="_")

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

patient.info$Patient = gsub("/", "_", patient.info$Patient)
for (patient.name in c('AD0531','AD0591','AHM0363','PR1_HIN_044') ) { # unique(patient.info$Patient)  ) {
#for(patient.index in 1:10) {
#	patient.name = unique(patient.info$Patient)[patient.index]
	print(patient.name)
	#patient.name = gsub("/","_",patient.name)
	
	for(gamma2 in c(5,10,25,50,100,250,500,1000)){
		print(gamma2)
	  plotdir = paste("hierarchicalClustering",patient.name, paste("gamma2", gamma2, sep="_"), sep="/")
	  if (!dir.exists(plotdir)) 
	    dir.create(plotdir, recursive=T)
	  
		segvals = read.table(paste("multipcf_plots_fitted_perPatient/",patient.name,"_segmentedCoverage_fitted_gamma",gamma2,".txt",sep=""),sep="\t",stringsAsFactors=F,header=T)
		segvals = segvals[segvals$n.probes>=min.probes,]
		if(nrow(segvals)>0) {
			HC = hclust(dist(t(segvals[,-(1:5)])))
			
			png(paste(plotdir,"/",patient.name,"_hierarchicalClustering_fromFittedSegments_gamma",gamma2,".png",sep=""),width=1000)
			plot(HC, main="clustering by fitted segment coverage",xlab="",sub="")
			dev.off()

			normalised.segvals = segvals[,-(1:5)]
			for(c in 1:nrow(normalised.segvals)) {
				normalised.segvals[c,] = (normalised.segvals[c,]-mean(unlist(normalised.segvals[c,])))/sd(unlist(normalised.segvals[c,]))
			}
			
			HC = hclust(dist(t(normalised.segvals)))
			png(paste(plotdir,"/",patient.name,"_hierarchicalClustering_fromNormalizedLargeFittedSegments_gamma",gamma2,".png",sep=""),width=1000)
			plot(HC, main="clustering by fitted segment coverage",xlab="",sub="")
			dev.off()
			
			no.samples = ncol(normalised.segvals)
			no.rows = ceiling(no.samples/3)
			png(paste(plotdir,"/",patient.name,"_NormalizedLargeFittedSegments_histogram_gamma",gamma2,".png",sep=""),width=1500,height=600*no.rows)
			par(mfrow=c(no.rows,3))
			for(s in 1:no.samples) {
				hist(normalised.segvals[,s],col="blue",main=names(normalised.segvals)[s],xlab="normalised coverage")
			}
			dev.off()
			if(nrow(segvals)>=3) {
			  patient = subset(patient.info, Patient == patient.name)
			  patient = patient[order(patient$Endoscopy_Year),]
			  
			  # order by date of endoscopy
			  normalised.segvals = normalised.segvals[, patient$samplename]
        
			  png(paste(plotdir,"/",patient.name,"_NormalizedLargeFittedSegments_correlationPlots_gamma",gamma2,".png",sep=""),width=600*no.rows,height=600*no.rows)			
			  pairs = ggpairs( normalised.segvals, 
			           lower=list(continuous=wrap("points", colour='darkblue', shape=20, alpha=0.3, size=3)),
			           diag=list(continuous=wrap(my_bars, sizeRange=c(1,5))), 
			           upper=list(continuous=my_custom_cor), 
			           axisLabels='show', title=patient.name, verbose=T, columnLabels = patient$Endoscopy_Year )
			  
        print(pairs  + theme( # this doesn't work...
          text = element_text(size=20, face='bold')
        )) 
				dev.off()
				
 				png(paste(plotdir,"/",patient.name,"_NormalizedLargeFittedSegments_correlationPlots_gamma",gamma2,"_old.png",sep=""),width=600*no.rows,height=600*no.rows)			
 				pairs( normalised.segvals,upper.panel=panel.cor, diag.panel=panel.hist, labels=paste(patient$samplename, patient$Endoscopy_Year, sep="\n"))
 				dev.off()
			}
			
			#get variable regions
			#segvals$variable = sapply(1:nrow(segvals),function(i){max(normalised.segvals[i,])-median(normalised.segvals[i,])>=1 | median(normalised.segvals[i,])-min(normalised.segvals[i,])>=1})
			#write.table(segvals,paste("multipcf_plots_fitted_perPatient/",patient.name,"_segmentedCoverage_fitted_gamma",gamma2,".txt",sep=""),sep="\t",quote=F,row.names=F)
		}
	}
}

q(save="no")
