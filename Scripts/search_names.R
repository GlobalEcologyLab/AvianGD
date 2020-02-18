library(readxl)

#### import the data ####
taxonomy  <- read_excel("~/Documents/UCPH/Paper/Data/Birds/IOC_8.1_vs_other_lists.xlsx", col_names=TRUE)

Sp       <- read.delim("~/Documents/UCPH/Paper/Data/Birds/list_species", header = F)

## LOOP
index <- 1
while ( index <= nrow(Sp) ) { # nrow(Sp)
  cat (index, "\n")
  sp_name <- Sp[index, 1] # prende il nome da cercare
  # il seguente while testa colonna per colonna se il nome viene trovato
  # in questo modo esce dal while appena trova il nome 
  # eseguendo un numero inferiore di ricerche (raddoppia la velocit?? della precedente che testa sempre tutte le colonne)
  colindex <- 1
  while ( length( rownum <- which(taxonomy[,colindex] == sp_name, arr.ind = TRUE) ) == 0 && colindex < ncol(taxonomy)) 
  { colindex <- colindex + 1 }
  
  #  rownum <- which(taxonomy[,c(1:ncol(taxonomy))] == sp_name, arr.ind = TRUE) #cerca il nome in taxonomy
  #  result <- ifelse(nrow(rownum) > 0, taxonomy[rownum[1,1],1], "NA") # salva in result il numero di riga oppure NA
  result <- ifelse(length(rownum) > 0, taxonomy[rownum,2], "NA") # salva in result il numero di riga oppure NA
  Sp[index,2] <- result # scrive il risultato nella colonna di Sp alla riga corrispondente all'indice index
  
  # 'fintanto che' il nome successivo == al nome attuale (sp_name) esegui le seguenti istruzioni:
  while ( (next_name <- Sp[index+1, 1]) == sp_name ) {
    index <- index + 1 # incrementa l'indice
    Sp[index,2] <- result # scrive il risultato nella riga corrispondente all'indice index
    
    # ritorna in testa al 'while' per valutare se il nuovo next_name == all'attuale
    # se e' diverso esce dal while per incrementare l'indice
  }
  index <- index + 1 # incrementa l'indice
}

rm(sp_name, next_name, result, rownum, colindex, index)
write.csv(file="~/Desktop/Master_DataBase_Edit.txt", Sp)
## save(Sp, file="~/Desktop/Master_DB_Edited.RData")

new_names <- read.csv("~/Documents/UCPH/Paper/Data/Birds/Master_DataBase_Edit.txt", header = T, sep = ",", stringsAsFactors = F)
new_names <- new_names[,2:3]

NA_list <- new_names[which(is.na(new_names$V2)),]
# subset the database to take names that had no match

data_tesi <- read.csv("/Volumes/Elisabetta/UCPH/Thesis(complete)/Data/GD_calculation/cytb_all_species/grid_ids.csv", header = F, sep = "\t", stringsAsFactors = F)
data_tesi <- data_tesi[, 2:4]
data_tesi$V2 <- gsub('_', ' ', data_tesi$V2)
data_tesi$V3 <- gsub('_', ' ', data_tesi$V3)
data_tesi$V4 <- gsub('_', ' ', data_tesi$V4)
# remove "_" and replace with space " "

# find names that didn't match in table used for thesis
index <- 1
spname <- NA_list$V1
for (s in seq_along(sp_name)){
  numrow <- which(data_tesi[,1] == spname[s])
  name <- data_tesi[numrow[1],3]
  NA_list[index,2] <- name
  index <- index+1
}

no_match <- NA_list[which(is.na(NA_list$V2)),1]
matched <- NA_list[which(!is.na(NA_list$V2)),]

# re-match new found names with IOC newest version list
index <- 1
while ( index <= nrow(matched) ) { # nrow(Sp)
  cat (index, "\n")
  sp_name <- matched[index, 2] # prende il nome da cercare
  # il seguente while testa colonna per colonna se il nome viene trovato
  # in questo modo esce dal while appena trova il nome 
  # eseguendo un numero inferiore di ricerche (raddoppia la velocit?? della precedente che testa sempre tutte le colonne)
  colindex <- 1
  while ( length( rownum <- which(taxonomy[,colindex] == sp_name, arr.ind = TRUE) ) == 0 && colindex < ncol(taxonomy)) 
  { colindex <- colindex + 1 }
  
  #  rownum <- which(taxonomy[,c(1:ncol(taxonomy))] == sp_name, arr.ind = TRUE) #cerca il nome in taxonomy
  #  result <- ifelse(nrow(rownum) > 0, taxonomy[rownum[1,1],1], "NA") # salva in result il numero di riga oppure NA
  result <- ifelse(length(rownum) > 0, taxonomy[rownum,2], "NA") # salva in result il numero di riga oppure NA
  matched[index,3] <- result # scrive il risultato nella colonna di Sp alla riga corrispondente all'indice index
  
  # 'fintanto che' il nome successivo == al nome attuale (sp_name) esegui le seguenti istruzioni:
  while ( (next_name <- matched[index+1, 2]) == sp_name ) {
    index <- index + 1 # incrementa l'indice
    matched[index,3] <- result # scrive il risultato nella riga corrispondente all'indice index
    
    # ritorna in testa al 'while' per valutare se il nuovo next_name == all'attuale
    # se e' diverso esce dal while per incrementare l'indice
  }
  index <- index + 1 # incrementa l'indice
}

write.csv(file = "~/Documents/UCPH/Paper/Data/Birds/matched_list.csv", matched)

index <- 1
for (s in seq_along(spname)){
  Sp[which(Sp$CHECKLIST_NAME == spname[s]),7] <- subset(no_match, no_match$CHECKLIST_NAME == spname[s])[1,3]
  index <- index+1
}