library(xlsx)
library(plyr)

strip.whitespace <- function (x) gsub("\\s+|\\s+", "", x)

read.patient.info<-function(file) {
  if (!file.exists(file))
    stop(paste(file, "doesn't exist or is not readable."))

  patient.info = NULL
  for (i in 1:length(getSheets(loadWorkbook(file)))) {
    print(i)
    if (names(getSheets(loadWorkbook(file))[i]) == 'Technical Repeats') next

    ws = read.xlsx2(file, sheetIndex=i, stringsAsFactors=F, header=T)
    ws = ws[which( ws$Patient.ID != ""), grep('Patient|Path\\.ID|Endoscopy.Year|Pathology$|Progressor|Plate.Index|SLX|cellularity|p53', colnames(ws), value=T, ignore.case=T)]
    head(ws)
    #head(patient.info)
    
    slx.cols = grep('SLX', colnames(ws), value=T)
    if (length(slx.cols) > 1) {
      slx.rows = ws[,slx.cols]
      ws[slx.cols] = lapply(ws[slx.cols], as.null)
      ws$SLX.ID = strip.whitespace(apply(slx.rows, 1, function(x) paste(x, collapse='_')))
    }
    
    ## Order columns
    ws = ws[,c(grep('Patient',colnames(ws), ignore.case=T),
               grep('Path\\.ID',colnames(ws), ignore.case=T),
               grep('Progressor',colnames(ws), ignore.case=T),
               grep('Endoscopy.year',colnames(ws), ignore.case=T),
               grep('Pathology',colnames(ws),ignore.case=T),
               grep('Plate.Index',colnames(ws),ignore.case=T),
               grep('SLX',colnames(ws),ignore.case=T),
               grep('cellularity',colnames(ws),ignore.case=T),
               grep('p53', colnames(ws),ignore.case=F)) ]
    
    
    colnames(ws) = c('Patient','Path.ID','Status','Endoscopy.Year','Pathology','Plate.Index','SLX.ID','Barretts.Cellularity', 'p53.Status')
    
    if (is.null(patient.info)) {
      patient.info = ws
    } else {
      patient.info = rbind(patient.info, ws)
    }
  }
  patient.info$SLX.ID = gsub('SLX-', '', strip.whitespace( patient.info$SLX.ID ) )

  patient.info$Status = factor(strip.whitespace(patient.info$Status))
  patient.info$Pathology = factor(strip.whitespace(patient.info$Pathology))
  patient.info$Endoscopy.Year = as.numeric(patient.info$Endoscopy.Year)
  patient.info$Plate.Index = strip.whitespace(patient.info$Plate.Index)
  patient.info$Barretts.Cellularity = as.integer(patient.info$Barretts.Cellularity)
  patient.info$p53.Status = as.integer(patient.info$p53.Status)
  
  patient.info$Samplename = gsub('-', '_', paste(patient.info$Plate.Index,strip.whitespace(patient.info$SLX.ID),sep="_"))
  ## TODO Mistake in earlier version of the patient file caused this will be fixed for the next run
  snames = grep('_1072(5|9)$', patient.info$Samplename)
  if ( length(snames > 0) ) {
    patient.info$Samplename[snames] =  paste( sub('-','_', patient.info$Plate.Index[snames]), '10725_10729', sep='_' )
  }

  head(patient.info)
  patient.info$Pathology
  
  # Remove 'normal' samples
  patient.info = subset(patient.info, Pathology != 'D2')
  patient.info = subset(patient.info, Pathology != 'Gastriccardia')
  patient.info$Pathology = droplevels(patient.info$Pathology)
  patient.info$Pathology = ordered( patient.info$Pathology, 
                                    levels=c("normal","?","BE","ID","LGD", "IMC","HGD" ))
  
  patient.info = arrange(patient.info, Status, Patient, Endoscopy.Year, Pathology)

  return(patient.info)
}
