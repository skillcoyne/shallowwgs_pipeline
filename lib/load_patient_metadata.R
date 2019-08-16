
strip.whitespace <- function (x) gsub("\\s+|\\s+", "", x)

summarise.patient.info<-function(df) {
  require(plyr)
  
  sum.pt = df %>% dplyr::group_by(Patient, Hospital.Research.ID, Status, Set) %>% 
    dplyr::summarise( 
        age.diagnosis = unique(Age.at.diagnosis),
        wt = unique(`Weight.(kg)`),
        ht = unique(`Height.(cm)`),
        BMI = unique(BMI),
        gender = unique(Sex),
        smoking = unique(Smoking),
        C = unique(Circumference),
        M = unique(Maximal),
        start.year=range(Endoscopy.Year)[1], 
        end.year=range(Endoscopy.Year)[2],
        total.samples=length(Samplename), 
        total.endos=length(unique(PID)),
        initial.endo = unique(Initial.Endoscopy),
        final.endo = unique(Final.Endoscopy),
        initial.path=sort(Pathology, decreasing=F)[1],
        highest.path=sort(Pathology, decreasing=T)[1],
        first.batch=sort(Batch.Name)[1],
        med.cellularity=median(Barretts.Cellularity,na.rm=T),
        Initial.Analysis=(sort(Batch.Name)[1] %in% levels(Batch.Name)[1:5]))
  return(arrange(sum.pt, Status, total.samples, highest.path))
}

read.patient.info<-function(file, file2=NULL, set='All', sheet=NULL) {
  if (!file.exists(file))
    stop(paste(file, "doesn't exist or is not readable."))
  
  if (grepl('\\.xlsx$', basename(file))) {
    library(readxl)
    
    patient.info = NULL;  
      if (is.null(sheet)) sheet = 'All combined'
      ws = readxl::read_xlsx(file, sheet = sheet)

      slx.cols = grep('SLX', colnames(ws), value=T)
      slx.cols = slx.cols[order(slx.cols)]
      
      if (length(slx.cols) > 1) {
        slx.rows = ws[,slx.cols]
        ws[slx.cols] = lapply(ws[slx.cols], as.null)
        ws$SLX.ID = sub('_NA$', '', strip.whitespace(apply(slx.rows, 1, function(x) paste(x, collapse='_'))))
      }
      
      ## Order columns
      ws = ws[,c(grep('Patient',colnames(ws), ignore.case=T),
                 grep('^Hospital', colnames(ws), ignore.case=T),
                 grep('Path ID',colnames(ws), ignore.case=T),
                 grep('Progressor',colnames(ws), ignore.case=T),
                 grep('Endoscopy year',colnames(ws), ignore.case=T),
                 grep('Pathology$',colnames(ws),ignore.case=T),
                 grep('Plate Index',colnames(ws),ignore.case=T),
                 grep('SLX',colnames(ws),ignore.case=T),
                 grep('cellularity',colnames(ws),ignore.case=T),
                 grep('p53 status', colnames(ws),ignore.case=F),
                 grep('Number of', colnames(ws), ignore.case=F),
                 grep('Batch', colnames(ws), ignore.case=F), 
                 grep('Set', colnames(ws), ignore.case=F),
                 grep('Block',colnames(ws),ignore.case=T)) ]
      
      colnames(ws) = c('Patient', 'Hospital.Research.ID', 'Path.ID','Status','Endoscopy.Year','Pathology','Plate.Index','SLX.ID','Barretts.Cellularity', 'p53.Status', 'Total.Reads', 'Batch.Name', 'Set', 'Block')
      
        patient.info = ws
  } else {
    patient.info = read.table(file, header=T, sep='\t', stringsAsFactors=F, quote="")

    slx.cols = grep('^SLX', colnames(patient.info), value=T)
    if (length(slx.cols) > 1) {
      slx.rows = patient.info[,slx.cols]
      patient.info[slx.cols] = lapply(patient.info[slx.cols], as.null)
      patient.info$SLX.ID = strip.whitespace(apply(slx.rows, 1, function(x) paste(x, collapse='_')))
      patient.info$SLX.ID = sub('_$', '', patient.info$SLX.ID)
    }
    patient.info = patient.info[, c(grep('Patient',colnames(patient.info), ignore.case=T),
                                    grep('Hospital.Research.ID', colnames(patient.info), ignore.case=T),
                                    grep('Path\\.ID',colnames(patient.info), ignore.case=T),
                                    grep('^(Progressor|Status)',colnames(patient.info), ignore.case=T),
                                    grep('Endoscopy.year',colnames(patient.info), ignore.case=T),
                                    grep('^Pathology',colnames(patient.info),ignore.case=T),
                                    grep('Plate.Index',colnames(patient.info),ignore.case=T),
                                    grep('SLX',colnames(patient.info),ignore.case=T),
                                    grep('cellularity',colnames(patient.info),ignore.case=T),
                                    grep('p53', colnames(patient.info),ignore.case=F),
                                    grep('Number.of|Total.Reads', colnames(patient.info), ignore.case=F), 
                                    grep('Batch', colnames(patient.info), ignore.case=F),
                                    grep('Set', colnames(patient.info), ignore.case=F),
                                    grep('Block', colnames(patient.info), ignore.case=F)) ]
    colnames(patient.info) = c('Patient','Hospital.Research.ID', 'Path.ID','Status','Endoscopy.Year','Pathology','Plate.Index','SLX.ID','Barretts.Cellularity', 'p53.Status', 'Total.Reads', 'Batch.Name', 'Set', 'Block')
  }
  
  patient.info$SLX.ID = gsub('SLX-', '', strip.whitespace( patient.info$SLX.ID ) )
  patient.info$Plate.Index = strip.whitespace(patient.info$Plate.Index)

  patient.info[] = lapply(patient.info[], strip.whitespace)
  
  patient.info[c('Status','Pathology','p53.Status')] = 
    lapply(patient.info[c('Status','Pathology','p53.Status')], function(x) factor(strip.whitespace(x)))
  
  patient.info[c('Endoscopy.Year','Barretts.Cellularity','Total.Reads','Block')] = 
    lapply(patient.info[c('Endoscopy.Year','Barretts.Cellularity','Total.Reads','Block')], as.numeric)
  
  patient.info$Samplename = gsub('-', '_', paste(patient.info$Plate.Index,strip.whitespace(patient.info$SLX.ID),sep="_"))
  ## TODO Mistake in earlier version of the patient file caused this
  snames = grep('_1072(5|9)$', patient.info$Samplename)
  if ( length(snames > 0) ) 
    patient.info$Samplename[snames] =  paste( sub('-','_', patient.info$Plate.Index[snames]), '10725_10729', sep='_' )

  # Remove 'normal' samples
  removed = subset(patient.info, grepl('D2|Gastriccardia|normal', Pathology, ignore.case=T))
  
  patient.info = subset(patient.info, !grepl('D2|Gastriccardia|normal', Pathology, ignore.case=T))

  patient.info$Pathology = droplevels(patient.info$Pathology)
  patient.info$Pathology = ordered( patient.info$Pathology, levels=c("BE","ID","LGD","HGD","IMC" ))
  
  patient.info$Batch.Name = as.factor(patient.info$Batch.Name)
  
  patient.info = arrange(patient.info, Status, Patient, Endoscopy.Year, Pathology)
  
  if (!grepl('all', set, ignore.case = T)) {
    patient.info = subset(patient.info, Set == set)
    message(paste("Returning only the", set, "set", sep=" "))
  } else {
    message(paste("Returning all patient data."))
  }

  patient.info$Hospital.Research.ID = gsub('/', '_',patient.info$Hospital.Research.ID)
  removed$Hospital.Research.ID = gsub('/', '_',removed$Hospital.Research.ID)
  
  message(paste(length(unique(patient.info$Hospital.Research.ID)), 'unique patient IDs'))
  
  if (!is.null(file2)) {
    patient.info = add.demographics(file2, patient.info)
  }
  
  patient.info$PID = sub('_$', '', unlist(lapply(patient.info$Path.ID, function(x) unlist(strsplit(x, 'B'))[1])))
  
  patient.info$Patient = as.integer(patient.info$Patient)
  
  return(list('info'=as_tibble(patient.info), 'normal'=as_tibble(removed)))
}

add.demographics<-function(file, patient.info) {
  bmicalc<-function(x) {
    ht = x[[1]]; wt = x[[2]]
    
    if (is.na(ht) || is.na(wt) || ht < 0 || wt < 0)
      return(NA)
    return( (wt/(ht/100))/(ht/100) ) 
  }
  
  if (grepl('\\.xlsx$', basename(file))) {
    require(readxl)
    
    demog = read_excel(file, sheet=1)
    #demog = read.xlsx2(file, sheetIndex = 1, colIndex = c(1:21), endRow=91, stringsAsFactors=F)
    colnames(demog) = gsub(' ', '.',colnames(demog))
  } else {
    demog = read.table(file, header=T, sep='\t', quote="", stringsAsFactors=F)
    demog = demog[1:91,1:21]
    colnames(demog) = sub('\\.\\.', '\\.', colnames(demog))
    colnames(demog) = gsub(' ', '.',colnames(demog))
    demog[demog == ""] = NA
  }

  demog$Study.Number = gsub('/', '_', demog$Study.Number)
  demog$Hospital.Research.ID =  strip.whitespace( demog$Study.Number )
  demog$Alternate.Study.Number = strip.whitespace( gsub('/', '_', demog$Alternate.Study.Number) )
  
  rows = with(demog, which(Alternate.Study.Number != ""))
  demog$Hospital.Research.ID[rows] = demog$Alternate.Study.Number[rows]
  
  # good...
  #length(which(unique(patient.info$Patient) %in% demog$Patient.ID))
  
  dates = grep('^date', colnames(demog), ignore.case = T, value=T)
  demog[dates] = lapply(demog[dates], function(x) as.Date(gsub('/', '-', x), '%d-%M-%Y'))
  
  demog$Year.Birth = format(as.Date(demog$Date.of.birth), "%Y")
  demog$Date.of.birth = NULL
  
  cols = grep('Year', colnames(demog), value=T)
  demog[cols] = lapply(demog[cols], as.numeric)
  
  cols = grep('Weight|Height|Age|Circumference|Maximal', colnames(demog), value=T)
  demog[cols] = lapply(demog[cols], as.numeric)
  
  demog = demog[,c('Hospital.Research.ID', grep('Sex|Circumference|Maximal|Date|Year|Age|Weight|Height|Smoking', colnames(demog), value=T))]
  colnames(demog) = gsub("'", '', colnames(demog))

  colnames(demog)[grep('Smoking', colnames(demog))] = 'Smoking'
  #head(demog)
  demog$BMI = apply(demog[c('Height.(cm)', 'Weight.(kg)')], 1, bmicalc)

  if (!is.null(patient.info)) {
    demog = merge(patient.info, demog, by='Hospital.Research.ID')
    demog = demog %>% group_by(Hospital.Research.ID) %>% dplyr::mutate( 
                  'Initial.Endoscopy'=min(Endoscopy.Year, Year.of.1st.Endoscopy), 'Final.Endoscopy'=max(Endoscopy.Year, Year.of.Endpoint) )
  } else {
    colnames(demog)[which(colnames(demog) == 'Hospital.Research.ID')] = 'Hospital.Research.ID'
  }

  demog[c('Sex','Smoking')] = lapply(demog[c('Sex','Smoking')], as.factor)
  
  return(demog)
}


