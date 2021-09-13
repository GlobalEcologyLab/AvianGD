## This script has the purpose of identifying the number of realms covered by the species for which we have values of genetic diversity ##
## The list of bird species in each zoogegographic realm comes from the paper ().
## They use IOC Taxonomy, so we have to use the IOC names instead of the ones defined by BirdLifeInternational v2

# set up libraries
library(plyr)
library(tidyverse)
library(data.table)

# import master dataset, which contains both IOC and BirdLifeInternational v2 taxonomies.
mdb <- fread("./Data/Master_db.csv")
cols <- c("ACC_NUM", "IOC_SPECIES", "FINAL_NAME", "RED_LIST_JETZ.2017.", "NOTES")
mdb <- mdb[, ..cols]
mdb$FINAL_NAME <- ifelse(mdb$FINAL_NAME == "", NA, mdb$FINAL_NAME)

# import dataset used in the analysis
d <- read_delim(file.choose(), delim = ',') %>% .[-which(.$GD == 0),]

# remove sequences that have been discarded from the calculation of GD
data2 <- mdb[-which(mdb$NOTES == "EX" | mdb$NOTES == "DOM" | mdb$NOTES == "HYB" | 
                     mdb$NOTES == "NN" | mdb$NOTES == "Newly discovered taxon" | 
                     mdb$NOTES == "IDENTICAL" | mdb$NOTES == "REMOVE"),]
data2 <- data2[!is.na(FINAL_NAME)]

# remove duplicates (each row has one species)
data2 <- as.data.frame(cbind(gsub(" ", "_", data2$IOC_SPECIES), gsub(" ", "_", data2$FINAL_NAME), data2$RED_LIST_JETZ.2017., data2$NOTES))
colnames(data2) <- c("IOC_SPECIES", "FINAL_NAME", "RED_LIST_JETZ.2017.", "NOTES")
data2 <- data2[!duplicated(data2),]

# keep only the species analysed
final <- data2[match(d$SP, data2$FINAL_NAME),]

# import realms datasets
path <- "~/Library/Mobile Documents/com~apple~CloudDocs/PhD/Chapter1/GeneticGap/Data/bird realm spp/"
files <- list.files(path = path, pattern = ".csv", full.names = T)
realms <- lapply(files, read.csv, stringsAsFactors=F)

final[which(final$IOC_SPECIES %in% realms[[1]]$X),3] %>% length()*100/nrow(realms[[1]])
final[which(final$IOC_SPECIES %in% realms[[2]]$X),3] %>% length()*100/nrow(realms[[2]])
final[which(final$IOC_SPECIES %in% realms[[3]]$X),3] %>% length()*100/nrow(realms[[3]])
final[which(final$IOC_SPECIES %in% realms[[4]]$X),3] %>% length()*100/nrow(realms[[4]])
final[which(final$IOC_SPECIES %in% realms[[5]]$X),3] %>% length()*100/nrow(realms[[5]])
final[which(final$IOC_SPECIES %in% realms[[6]]$X),3] %>% length()*100/nrow(realms[[6]])
final[which(final$IOC_SPECIES %in% realms[[7]]$X),3] %>% length()*100/nrow(realms[[7]])
final[which(final$IOC_SPECIES %in% realms[[8]]$X),3] %>% length()*100/nrow(realms[[8]])
final[which(final$IOC_SPECIES %in% realms[[9]]$X),3] %>% length()*100/nrow(realms[[9]])
final[which(final$IOC_SPECIES %in% realms[[10]]$X),3] %>% length()*100/nrow(realms[[10]])
final[which(final$IOC_SPECIES %in% realms[[11]]$X),3] %>% length()*100/nrow(realms[[11]])




