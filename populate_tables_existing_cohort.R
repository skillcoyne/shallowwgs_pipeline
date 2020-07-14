
library(tidyverse)
library(RSQLite)


conn <- dbConnect(RSQLite::SQLite(), "~/OneDrive - MRC Cancer Unit at University of Cambridge/Progressor Study/Patient Database/BE_Progression_Project.db")

delete_all<-function(conn) {
  tables = dbListTables(conn)
  for (t in tables) {
    print(t)
    dbExecute(conn, paste0("DELETE FROM ", t))
  }
}

row_counts<-function(conn, tables) {
  for (t in tables)
    message(paste0(t, ': ',dbGetQuery(conn, paste0("SELECT COUNT(*) FROM ", t)) ))
}

get_patient_id<-function(conn, id) {
  query = paste0("SELECT DISTINCT PatientID FROM AlternatePatientID WHERE AlternateID = ",shQuote(id)," AND AltIDStudy == 'excel sheet id'")
  dbGetQuery(conn, query) %>% pull
}


tables = c('SequenceMetadata','EndoscopySample', 'Endoscopy','AlternatePatientID','Patient', 'sqlite_sequence')

train = T
val = T

if (train) {
  ### Training Cohort
  all_pt = readxl::read_xlsx('~/Data/BarrettsProgressionRisk/QDNAseq/training/All_patient_info.xlsx', sheet = 'All combined') %>% arrange(`Patient ID`) %>%
    separate(`Path ID`, c('PathCaseID', 'Block', 'Notes'), sep = '(_|-)', remove = F, fill='right', extra='merge') %>% 
    dplyr::mutate_at(vars(PathCaseID, Block), list(~str_trim(., 'both'))) %>%
    dplyr::mutate_at(vars(matches('SLX|Index')), list(~gsub(' ', '', .)))  %>%
    rowwise() %>%
    dplyr::mutate( Samplename = paste( paste(na.omit(c(`SLX ID 1`, `SLX ID 2`)), collapse='_'), sub('-', '_',`Plate Index`), sep='.' ) )
  
  demo = readxl::read_xlsx('~/Data/BarrettsProgressionRisk/QDNAseq/training/Demographics_full.xlsx') %>% 
    dplyr::mutate(Sex = recode(Sex, 'M'='male','F'='female'), `P/NP` = dplyr::recode(`P/NP`, 'NP'='non-progressor', 'P'='progressor'))
  
  ## Patient table
  for (id in unique(all_pt$`Patient ID`)) {
    print(id)  
    pt = all_pt %>% filter(`Patient ID` == id) 
    did = pt %>% dplyr::select(`Hospital Research ID`) %>% distinct() %>% pull
    dem = demo %>% filter(`Study Number` == did | `Alternate Study Number` == did) %>% dplyr::select(`P/NP`, Sex, `Age at diagnosis`, `Smoking Status`,`further smoking info`, `Height (cm)`, `Weight (kg)`, `Year of 'Endpoint'`, `Date of First Endoscopy`, Circumference, Maximal)
    
    dem = dem %>% mutate_if( .predicate = ~ (is.numeric(.) | is.character(.)),  list(~replace_na(.,'null')))
    
    statement = 'INSERT INTO Patient (Gender, AgeAtDiagnosis, SmokingStatus, Height, Weight, Status, DateInitialDiagnosis, DateProgressed, Circumference, Maximal) VALUES'
    
    values = dem %>% dplyr::mutate(smoking = case_when( (dem$`Smoking Status` =='Y' & grepl('current',dem$`further smoking info`,ignore.case = T)) ~ "'current'",
                                               (dem$`Smoking Status` =='Y' & !grepl('current',dem$`further smoking info`,ignore.case = T)) ~ "'former'",
                                               dem$`Smoking Status` =='N' ~ "'never'", 
                                               T ~ dem$`Smoking Status`),
                          Initial = shQuote(min(c(`Date of First Endoscopy`, ISOdate(pt$`Endoscopy Year`, '01','01' ) ))), 
                          Final = case_when(`P/NP` == 'progressor' ~ shQuote(max(c(ISOdate(dem$`Year of 'Endpoint'`,'01','01'), ISOdate(pt$`Endoscopy Year`, '01','01' ) ))), T ~ 'null'),
      value = paste0('(',paste( c(shQuote(Sex), `Age at diagnosis`, smoking, `Height (cm)`, `Weight (kg)`, shQuote(dem$`P/NP`), Initial, Final, Circumference, Maximal), collapse = ','),')') ,
      ) %>% dplyr::select(value) %>% pull
  
    statement = paste(statement, paste(values,sep=',') )
    #print(statement)
    dbExecute(conn, statement)
    
    pid = dbGetQuery(conn, 'SELECT MAX(ID) FROM Patient') %>% pull
  
    ## Alternate patient IDs  
    did = pt %>% dplyr::select(`Hospital Research ID`) %>% distinct() %>% pull
    dem = demo %>% filter(`Study Number` == did | `Alternate Study Number` == did) %>% dplyr::select(matches('Study Number')) 
    altids = na.omit(unique(c(did, dem$`Study Number`, dem$`Alternate Study Number`)))
    
    statement = paste0("INSERT INTO AlternatePatientID (PatientID, AlternateID) VALUES",paste(paste0('(',paste(pid,shQuote(na.omit(altids)),sep=','),')'), collapse=','))
    #print(statement)
    dbExecute(conn, statement)
    
    statement = paste0("INSERT INTO AlternatePatientID (PatientID, AlternateID, AltIDStudy) VALUES",paste0('(', paste(c(pid,shQuote(id), shQuote('excel sheet id')), collapse=','), ')'))
    #print(statement)
    dbExecute(conn, statement)
  }
  
  row_counts(conn,rev(tables)[-1])
  
  ## Endoscopy & EndoscopySample
  for (id in unique(all_pt$`Patient ID`)) {
    print(id)  
    
    pid = get_patient_id(conn, id)
    
    pt = all_pt %>% filter(`Patient ID` == id) %>% mutate(DBPatientID = pid)
    
    did = pt %>% dplyr::select(`Hospital Research ID`) %>% distinct() %>% pull
    dem = demo %>% filter(`Study Number` == did | `Alternate Study Number` == did) %>% dplyr::select(matches('Study Number')) 
    
    tmp = pt %>% filter(`Patient ID` == id) %>% dplyr::select(DBPatientID, `Patient ID`, `Path ID`, `PathCaseID`, `Endoscopy Year`, Block, Notes, Pathology, `Pathology notes`, `p53 status`) %>% 
      arrange(`Endoscopy Year`) %>%
      dplyr::mutate(`Endoscopy Year` = ISOdate(`Endoscopy Year`, 01, 01)) %>%
      group_by(`Patient ID`, PathCaseID, `Endoscopy Year`) %>% 
      dplyr::mutate(Notes = paste(`Path ID`, collapse='; ')) %>%
      dplyr::mutate_at(vars(Notes, Block,), list(~replace_na(.,'null'))) %>% ungroup
    
     values = tmp %>% dplyr::select( `DBPatientID`, PathCaseID, `Endoscopy Year`, Notes) %>%     
       distinct %>% dplyr::mutate( value = paste0('(',paste(`DBPatientID`, shQuote(`PathCaseID`), shQuote(as.character(`Endoscopy Year`)),shQuote(Notes), sep=','),')')) %>% dplyr::select(value) %>% pull
  
     statement = paste("INSERT INTO Endoscopy (PatientID, PathCaseID, EndoscopyDate, Notes) VALUES", paste(unique(values), collapse = ','))
     dbExecute(conn, statement)
     
     select = paste0("SELECT * FROM Endoscopy WHERE PatientID == ", pid, " AND PathCaseID IN (", paste(shQuote(tmp$PathCaseID), collapse=','), ')')
     rows = dbGetQuery(conn, select) %>% arrange(ID) %>% dplyr::rename(EndoscopyID = ID)
     if (nrow(rows) <= 0)
       stop(paste("No rows in 'Endoscopy' for PatientID", pid, 'and PathCaseID', paste(tmp$PathCaseID, collapse=',')))
  
     tmp = tmp %>% left_join(rows %>% dplyr::select(-Notes), by='PathCaseID') %>% 
       dplyr::mutate(`p53 status` = recode(`p53 status`, '0'=shQuote('normal'), '1'=shQuote('aberrant')), Pathology = recode(Pathology, 'BE'='NDBE'),
                     `Pathology notes` = case_when(!is.na(`Pathology notes`) ~ shQuote(`Pathology notes`), T ~ `Pathology notes` )) %>% 
       dplyr::mutate_at(vars(matches('p53|notes')), list(~replace_na(.,'null')))
  
     values = tmp %>% mutate(value = paste0('(',paste(DBPatientID, EndoscopyID, shQuote(Block), `p53 status`, shQuote(Pathology), `Pathology notes`, sep=',' ),')')) %>%
       dplyr::select(value) %>% pull
     
     statement = paste0("INSERT INTO EndoscopySample (PatientID, EndoscopyID, Block, p53IHC, Pathology, PathNotes) VALUES ", paste(values, collapse=', ') )
     dbExecute(conn, statement)
     
     # sWGS info 
     statement = 'INSERT INTO SequenceMetadata (PatientID, SampleTable, SampleID,Samplename,SLX,PlateIndex,Replicate,TotalReads) VALUES'
     
     values = left_join(dplyr::select(rows, -Notes, -EndoscopyDate), dplyr::select(pt, PathCaseID, matches('DB|SLX|Index|Reads|Replicate|Samplename')), by='PathCaseID') %>% 
       dplyr::group_by(EndoscopyID, Samplename) %>%
       dplyr::mutate(SLX = paste(na.omit(c(`SLX ID 1`, `SLX ID 2`)), collapse='; '), Replicate = !is.na(`Technical Replicate`)) %>% ungroup %>%
       dplyr::mutate_all(list(~replace_na(.,'null'))) %>%
       dplyr::mutate(value = paste0('(',paste(DBPatientID, shQuote('EndoscopySample'), EndoscopyID, shQuote(Samplename), shQuote(SLX), shQuote(`Plate Index`), shQuote(Replicate), `Number of reads`, sep=', '), ')') ) %>% ungroup %>%
       dplyr::select(value) %>% pull
    
     statement = paste0(statement, paste(values,collapse=','))
     #print(statement)
     dbExecute(conn, statement)
  }
  
  rm(all_pt, demo)
}
row_counts(conn,rev(tables)[-1])

if (val) {
  ### Validation Cohort
  all_pt = readxl::read_xlsx('~/Data/BarrettsProgressionRisk/QDNAseq/validation/sWGS_validation_batches.xlsx', sheet = 'Final Validation Samples') %>% arrange(`Patient`) %>%
    separate(`Block ID`, c('PathCaseID', 'Block', 'Notes'), sep = '(_|[:blank:]|-)+', remove = F, fill='right', extra='merge') %>% 
    dplyr::mutate_at(vars(PathCaseID, Block), list(~str_trim(., 'both'))) %>%
    dplyr::mutate_at(vars(matches('SLX|Index')), list(~gsub(' ', '', .)))  %>%
    rowwise() %>%
    dplyr::mutate( Samplename = paste( `SLX-ID`, sub('-', '_',`Index Sequence`), sep='.' ) ) %>%
    dplyr::rename('Patient ID'='Patient')
  
  demo = readxl::read_xlsx('~/Data/BarrettsProgressionRisk/QDNAseq/validation/sWGS_validation_batches.xlsx', sheet = 'Demographics') %>% arrange(`Patient`) %>%
    dplyr::mutate(Sex = recode(Sex, 'M'='male','F'='female'), Status = recode(Status, 'P'='progressor','NP'='non-progressor'))
  
  ## Patient table
  for (id in unique(all_pt$`Patient ID`)) {
    print(id)  
    pt = all_pt %>% filter(`Patient ID` == id) 
    altids = t(unlist( sapply( pt %>% dplyr::select(`Hospital Research ID`) %>% distinct() %>% pull,  function(x) str_trim(str_split(x, pattern=',', simplify = T), 'both') ) ))
    
    dem = demo %>% filter(Patient == id) %>% dplyr::select(Patient, Status, Sex, `Age at diagnosis`, `Smoking`,`Smoking notes`, `Initial Diagnosis`, `Last Endoscopy`, Circumferential, Maximal) %>%
      dplyr::mutate_if( ~(is.character(.) | is.numeric(.)), list(~replace_na(.,'null'))) 
      
    statement = 'INSERT INTO Patient (Gender, AgeAtDiagnosis, SmokingStatus, Status, DateInitialDiagnosis, DateProgressed, Circumference, Maximal) VALUES'
    
    values = dem %>% dplyr::mutate(smoking = case_when( (dem$`Smoking` =='Y' & grepl('current',dem$`Smoking notes`,ignore.case = T)) ~ "'current'",
                                                        (dem$`Smoking` =='Y' & !grepl('current',dem$`Smoking notes`,ignore.case = T)) ~ "'former'",
                                                         dem$`Smoking` =='N' ~ "'never'", T ~ dem$`Smoking`),
                                   Initial = shQuote(`Initial Diagnosis`), 
                                   Final = case_when(`Status` == 'progressor' ~ shQuote(`Last Endoscopy`), T ~ 'null'),
                                   value = paste0('(',paste( c(shQuote(Sex), `Age at diagnosis`, smoking,  shQuote(dem$`Status`), Initial, Final, Circumferential, Maximal), collapse = ','),')') ,
    ) %>% dplyr::select(value) %>% pull
    
    statement = paste(statement, paste(values,collapse=',') )
    dbExecute(conn, statement)
    
    pid = dbGetQuery(conn, 'SELECT MAX(ID) FROM Patient') %>% pull
    
    ## Alternate patient IDs  
    statement = paste0("INSERT INTO AlternatePatientID (PatientID, AlternateID) VALUES",paste(paste0('(',paste(pid,shQuote(na.omit(altids)),sep=','),')'), collapse=','))
    #print(statement)
    dbExecute(conn, statement)
    
    statement = paste0("INSERT INTO AlternatePatientID (PatientID, AlternateID, AltIDStudy) VALUES",paste0('(', paste(c(pid,shQuote(id), shQuote('excel sheet id')), collapse=','), ')'))
    #print(statement)
    dbExecute(conn, statement)
  }
  
  row_counts(conn,rev(tables)[-1])
  
  ## Endoscopy & EndoscopySample
  for (id in unique(all_pt$`Patient ID`)) {
    print(id)  
    
    pid = get_patient_id(conn,id)
    
    pt = all_pt %>% filter(`Patient ID` == id) %>% mutate(DBPatientID = pid)
    dem = demo %>% filter(Patient == id) %>% mutate(DBPatientID = pid)
    
    tmp = pt %>% dplyr::select(DBPatientID, `Patient ID`, `Block ID`, `PathCaseID`, `Endoscopy`, Block, Notes, Pathology, `Path Notes`, `p53`) %>% 
      arrange(`Endoscopy`) %>%
      #dplyr::mutate(`Endoscopy` = ISOdate(`Endoscopy Year`, 01, 01)) %>%
      group_by(`DBPatientID`, PathCaseID, `Endoscopy`) %>% 
      dplyr::mutate(Notes = paste(`Block ID`, collapse='; ')) %>%
      dplyr::mutate_at(vars(Notes, Block,), list(~replace_na(.,'null'))) %>% ungroup
    
    values = tmp %>% dplyr::select( `DBPatientID`, PathCaseID, `Endoscopy`, Notes) %>%     
      distinct %>% dplyr::mutate( value = paste0('(',paste(`DBPatientID`, shQuote(`PathCaseID`), shQuote(as.character(`Endoscopy`)),shQuote(Notes), sep=','),')')) %>% dplyr::select(value) %>% pull
    
    statement = paste("INSERT INTO Endoscopy (PatientID, PathCaseID, EndoscopyDate, Notes) VALUES", paste(unique(values), collapse = ','))
    dbExecute(conn, statement)
    
    select = paste0("SELECT * FROM Endoscopy WHERE PatientID == ", pid, " AND PathCaseID IN (", paste(shQuote(tmp$PathCaseID), collapse=','), ')')
    rows = dbGetQuery(conn, select) %>% arrange(ID) %>% dplyr::rename(EndoscopyID = ID)
    if (nrow(rows) <= 0)
      stop(paste('No rows in Endoscopy table for patient',pid, "and PathCaseID ", paste(unique(tmp$PathCaseID),collapse=',')))
    
    tmp = tmp %>% left_join(rows %>% dplyr::select(-Notes), by='PathCaseID') %>% 
      dplyr::mutate( Pathology = recode(Pathology, 'BE'='NDBE'),
                     `Path Notes` = case_when(!is.na(`Path Notes`) ~ shQuote(`Path Notes`), T ~ 'null' )) %>% 
      dplyr::mutate_at(vars(matches('p53|notes')), list(~replace_na(.,'null')))
    
    values = tmp %>% mutate(value = paste0('(',paste(DBPatientID, EndoscopyID, shQuote(Block), shQuote(Pathology), `Path Notes`, sep=',' ),')')) %>%
      dplyr::select(value) %>% pull
    
    statement = paste0("INSERT INTO EndoscopySample (PatientID, EndoscopyID, Block, Pathology, PathNotes) VALUES ", paste(values, collapse=', ') )
    dbExecute(conn, statement)
    
    # sWGS 
    statement = 'INSERT INTO SequenceMetadata (PatientID, SampleTable, SampleID,Samplename,SLX,PlateIndex) VALUES'
    
    values = left_join(dplyr::select(rows, -Notes, -EndoscopyDate), dplyr::select(pt, PathCaseID, matches('DB|SLX|Index|Reads|Replicate|Samplename')), by='PathCaseID') %>% 
      dplyr::mutate_all(list(~replace_na(.,'null'))) %>%
      dplyr::group_by(EndoscopyID, Samplename) %>%
      dplyr::mutate(value = paste0('(',paste(DBPatientID, shQuote('EndoscopySample'), EndoscopyID, shQuote(Samplename), shQuote(`SLX-ID`), shQuote(`Index Sequence`), sep=', '), ')') ) %>% ungroup %>%
      dplyr::select(value) %>% pull
    
    statement = paste(statement, paste(values,collapse=','))
    #print(statement)
    dbExecute(conn, statement)
  }
  
  rm(all_pt, demo)
}

row_counts(conn,rev(tables)[-1])

