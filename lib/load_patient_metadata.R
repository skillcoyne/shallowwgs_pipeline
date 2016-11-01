library(xlsx)
library(plyr)

strip.whitespace <- function (x) gsub("\\s+|\\s+", "", x)

read.patient.info<-function(file) {
  if (!file.exists(file))
    stop(paste(file, "doesn't exist or is not readable."))

  patient.info = NULL
  for (i in 1:length(getSheets(loadWorkbook(file)))) {
    print(i)
    ws = read.xlsx2(file, sheetIndex=i, stringsAsFactors=F, header=T)
    ws = ws[which( ws$Patient.ID != ""), grep('Patient|Endoscopy.Year|Pathology$|Progressor|Plate.Index|SLX', colnames(ws), value=T, ignore.case=T)]
    head(ws)
    head(patient.info)
    
    slx.cols = grep('SLX', colnames(ws), value=T)
    if (length(slx.cols) > 1) {
      slx.rows = ws[,slx.cols]
      ws[slx.cols] = lapply(ws[slx.cols], as.null)
      ws$SLX.ID = strip.whitespace(apply(slx.rows, 1, function(x) paste(x, collapse='_')))
    }
    
    ## Order columns
    ws = ws[,c(grep('Patient',colnames(ws), ignore.case=T),
               grep('Progressor',colnames(ws), ignore.case=T),
               grep('Endoscopy.year',colnames(ws), ignore.case=T),
               grep('Pathology',colnames(ws),ignore.case=T),
               grep('Plate.Index',colnames(ws),ignore.case=T),
               grep('SLX',colnames(ws),ignore.case=T))]
    colnames(ws) = c('Patient','Status','Endoscopy.Year','Pathology','Plate.Index','SLX.ID')
    
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
  patient.info = arrange(patient.info, Status, Patient, Endoscopy.Year)
  
  return(patient.info)
}
