
strip.whitespace <- function (x) gsub("\\s+|\\s+", "", x)

summarise.patient.info<-function(df) {
  require(plyr)
  
  sum.pt = df %>% dplyr::group_by(Patient, Hospital.Research.ID, Status) %>% 
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

    slx.cols = grep('^SLX', colnames(ws), value=T)
    slx.cols = slx.cols[order(slx.cols)]
      
      if (length(slx.cols) > 1) {
        slx.rows = ws[,slx.cols]
        ws[slx.cols] = lapply(ws[slx.cols], as.null)
        ws$SLX.ID = sub('_NA$', '', strip.whitespace(apply(slx.rows, 1, function(x) paste(x, collapse='_'))))
      }
      
      ## Order columns
      ws = ws %>% dplyr::select(
        matches('^Patient', ignore.case=T), matches('^Hospital', ignore.case=T), matches('^Path ID', ignore.case=T), 
        matches('Progressor', ignore.case=T), matches('Endoscopy( |\\.|_)year', ignore.case=T),
        matches('Pathology$', ignore.case=T), matches('^Plate Index', ignore.case=T), 
        matches('^SLX', ignore.case=T), matches('cellularity', ignore.case=T), matches('p53 status', ignore.case=T),
        matches('Number of', ignore.case=T), matches('Batch', ignore.case=T), matches('Set', ignore.case=T),
        matches('Block', ignore.case=T), matches('Replicate')
      ) %>% 
        rlang::set_names(c('Patient', 'Hospital.Research.ID', 'Path.ID','Status','Endoscopy.Year','Pathology','Plate.Index','SLX.ID','Barretts.Cellularity', 'p53.Status', 'Total.Reads', 'Batch.Name', 'Set', 'Block', 'Technical.Replicate'))
      
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
  
  #patient.info$SLX.ID = gsub('SLX-', '', strip.whitespace( patient.info$SLX.ID ) )
  #patient.info = patient.info %>% mutate(SLX.ID = strip.whitespace( SLX.ID ) , Plate.Index = strip.whitespace(Plate.Index))
  
  patient.info = patient.info %>% mutate_all(strip.whitespace) %>% 
    mutate(Hospital.Research.ID = gsub('/', '_',Hospital.Research.ID))
  
  patient.info = patient.info %>% 
    mutate_at(vars(Status, Pathology, p53.Status), list(as.factor)) %>%
    mutate_at(vars(Endoscopy.Year, Barretts.Cellularity ,Total.Reads, Block), list(as.numeric)) %>% 
    mutate(Samplename = paste0(SLX.ID, '.', gsub('-','_',Plate.Index)) ) 
  
  patient.info$PID = sub('_$', '', unlist(lapply(patient.info$Path.ID, function(x) unlist(strsplit(x, 'B'))[1])))

  # Remove 'normal' samples
  normals = patient.info %>% filter(grepl('D2|Gastriccardia|normal', Pathology, ignore.case=T))
  patient.info = patient.info %>% filter(!grepl('D2|Gastriccardia|normal', Pathology, ignore.case=T))

  patient.info = patient.info %>% 
    mutate(Pathology = ordered(droplevels(Pathology),levels=c("BE","ID","LGD","HGD","IMC" )), Batch.Name = as.factor(Batch.Name), Patient = as.integer(Patient)) %>%
    arrange(Status, Patient, Endoscopy.Year, Pathology)
  
  
  # if (!grepl('all', set, ignore.case = T)) {
  #   patient.info = patient.info %>% filter(Set == set)
  #   message(paste("Returning only the", set, "set", sep=" "))
  # } else {
  #   message(paste("Returning all patient data."))
  # }
  message(paste(length(unique(patient.info$Hospital.Research.ID)), 'unique patient IDs'))
  
  if (!is.null(file2)) {
    patient.info = add.demographics(file2, patient.info)
  }

  # Remove replicates
  replicates = patient.info %>% filter(!is.na(`Technical.Replicate`)) %>% dplyr::mutate(Patient = as.integer(Patient))
  patient.info =  patient.info %>% filter(is.na(`Technical.Replicate`))
  
  
  return(list('info'=as_tibble(patient.info), 'normal'=as_tibble(normals), 'replicates'=as_tibble(replicates)))
}

add.demographics<-function(file, patient.info) {
  bmicalc<-function(ht,wt) {
    #ht = x[[1]]; wt = x[[2]]
    
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

  demog = demog %>% rowwise() %>% dplyr::mutate(
    Study.Number = gsub('/', '_', Study.Number),
    Hospital.Research.ID =  strip.whitespace(Study.Number),
    Alternate.Study.Number = strip.whitespace( gsub('/', '_', Alternate.Study.Number) )
  ) 
  rows = with(demog, which(Alternate.Study.Number != ""))
  demog$Hospital.Research.ID[rows] = demog$Alternate.Study.Number[rows]
  
  demog = demog %>% ungroup %>% dplyr::mutate(Year.Birth = format(as.Date(Date.of.birth), "%Y")) %>% 
    mutate_at(vars(matches('Weight|Height|Age|Circumference|Maximal|Year',ignore.case=T)), as.numeric) %>%
    dplyr::select(-Date.of.birth,-matches('further|notes')) %>% 
    dplyr::rename_all(list(~gsub("'",'',.))) %>%
    dplyr::rename_at(vars(matches('Smoking',ignore.case=T)), list(~sub('.*','Smoking',.))) %>%
    dplyr::mutate(BMI = bmicalc(`Height.(cm)`,`Weight.(kg)`)) %>%
    dplyr::mutate_at(vars(Sex,Smoking), list(~as.factor(.)))
    
      
  if (!is.null(patient.info)) {
    demog = left_join(patient.info, demog, by='Hospital.Research.ID') %>% 
      group_by(Hospital.Research.ID) %>% dplyr::mutate( 'Initial.Endoscopy'=min(Endoscopy.Year, Year.of.1st.Endoscopy), 'Final.Endoscopy'=max(Endoscopy.Year, Year.of.Endpoint) )
  } else {
    colnames(demog)[which(colnames(demog) == 'Hospital.Research.ID')] = 'Hospital.Research.ID'
  }

  return(demog)
}


