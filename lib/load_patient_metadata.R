
suppressPackageStartupMessages( library(plyr) )

strip.whitespace <- function (x) gsub("\\s+|\\s+", "", x)

summarise.patient.info<-function(df) {
  sum.pt = (df %>% 
    group_by(Patient, Status) %>%
    summarise(
      start.year=range(Endoscopy.Year)[1],
      end.year=range(Endoscopy.Year)[2],
      total.samples=length(Samplename), 
      highest.path=sort(Pathology, decreasing=T)[1],
      first.batch=sort(Batch.Name)[1],
      med.cellularity=median(Barretts.Cellularity,na.rm=T),
      Initial.Analysis=(sort(Batch.Name)[1] %in% levels(patient.info$Batch.Name)[1:5])
    ))
  return(arrange(sum.pt, Status, total.samples, highest.path))
}

read.patient.info<-function(file) {
  if (!file.exists(file))
    stop(paste(file, "doesn't exist or is not readable."))

  if (grepl('\\.xlsx$', basename(file))) {
    library(xlsx)

    patient.info = NULL
    for (i in 1:length(getSheets(loadWorkbook(file)))) {
    print(i)
    if (names(getSheets(loadWorkbook(file))[i]) == 'Technical Repeats') next
    batch_name = gsub(' ', '_', (names(getSheets(loadWorkbook(file))[i])))
      
    ws = read.xlsx2(file, sheetIndex=i, stringsAsFactors=F, header=T)
    ws = ws[which( ws$Patient.ID != ""), grep('Patient|Path\\.ID|Endoscopy.Year|Pathology$|Progressor|Plate.Index|SLX|cellularity|p53|Number', colnames(ws), value=T, ignore.case=T)]
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
               grep('p53', colnames(ws),ignore.case=F),
               grep('Number.of', colnames(ws), ignore.case=F)) ]
    
    
    colnames(ws) = c('Patient','Path.ID','Status','Endoscopy.Year','Pathology','Plate.Index','SLX.ID','Barretts.Cellularity', 'p53.Status', 'Total.Reads')
    
    ws[ws$p53.Status == "", 'p53.Status'] = NA
    ws[ws$Barretts.Cellularity == "", 'Barretts.Cellularity'] = NA
    
    ws$Batch.Name = batch_name
    
    if (is.null(patient.info)) {
      patient.info = ws
    } else {
      patient.info = rbind(patient.info, ws)
    }
  }
  } else {
    patient.info = read.table(file, header=T, sep='\t', stringsAsFactors=F)
    
    slx.cols = grep('SLX', colnames(patient.info), value=T)
    if (length(slx.cols) > 1) {
      slx.rows = patient.info[,slx.cols]
      patient.info[slx.cols] = lapply(patient.info[slx.cols], as.null)
      patient.info$SLX.ID = strip.whitespace(apply(slx.rows, 1, function(x) paste(x, collapse='_')))
      patient.info$SLX.ID = sub('_$', '', patient.info$SLX.ID)
    }
    patient.info = patient.info[, c(grep('Patient',colnames(patient.info), ignore.case=T),
                                    grep('Path\\.ID',colnames(patient.info), ignore.case=T),
                                    grep('Progressor',colnames(patient.info), ignore.case=T),
                                    grep('Endoscopy.year',colnames(patient.info), ignore.case=T),
                                    grep('Pathology',colnames(patient.info),ignore.case=T),
                                    grep('Plate.Index',colnames(patient.info),ignore.case=T),
                                    grep('SLX',colnames(patient.info),ignore.case=T),
                                    grep('cellularity',colnames(patient.info),ignore.case=T),
                                    grep('p53', colnames(patient.info),ignore.case=F),
                                    grep('Number.of', colnames(patient.info), ignore.case=F)) ]
    colnames(patient.info) = c('Patient','Path.ID','Status','Endoscopy.Year','Pathology','Plate.Index','SLX.ID','Barretts.Cellularity', 'p53.Status', 'Total.Reads')
    
    
  }
  
  patient.info$SLX.ID = gsub('SLX-', '', strip.whitespace( patient.info$SLX.ID ) )
  patient.info$Plate.Index = strip.whitespace(patient.info$Plate.Index)
  
  patient.info[c('Status','Pathology','p53.Status')] = 
    lapply(patient.info[c('Status','Pathology','p53.Status')], function(x) factor(strip.whitespace(x)))
  
  patient.info[c('Endoscopy.Year','Barretts.Cellularity','Total.Reads')] = 
    lapply(patient.info[c('Endoscopy.Year','Barretts.Cellularity','Total.Reads')], as.numeric)

  
  patient.info$Samplename = gsub('-', '_', paste(patient.info$Plate.Index,strip.whitespace(patient.info$SLX.ID),sep="_"))
  ## TODO Mistake in earlier version of the patient file caused this will be fixed for the next run
  snames = grep('_1072(5|9)$', patient.info$Samplename)
  if ( length(snames > 0) ) {
    patient.info$Samplename[snames] =  paste( sub('-','_', patient.info$Plate.Index[snames]), '10725_10729', sep='_' )
  }

  # Remove 'normal' samples
  patient.info = subset(patient.info, Pathology != 'D2')
  patient.info = subset(patient.info, Pathology != 'Gastriccardia')
  patient.info$Pathology = droplevels(patient.info$Pathology)
  patient.info$Pathology = ordered( patient.info$Pathology, 
                                    levels=c("normal","?","BE","ID","LGD", "IMC","HGD" ))
  
  patient.info$Batch.Name = ordered(patient.info$Batch.Name, 
                                    levels=c("NP_Pilot_Study","Progressor_Pilot_Study", "Exome_subcohort", paste("Batch", c(1:10, '16S'), sep='_')))
  
  patient.info = arrange(patient.info, Status, Patient, Endoscopy.Year, Pathology)

  return(patient.info)
}
