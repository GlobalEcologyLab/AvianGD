## This script has the purpose of identifying the number of realms covered by the species for which we have values of genetic diversity ##
## The list of bird species in each zoogegographic realm comes from the paper ().
## They use IOC Taxonomy, so we have to use the IOC names instead of the ones defined by BirdLifeInternational v2

# set up libraries
library(plyr)
library(tidyverse)

# import master dataset, which contains both IOC and BirdLifeInternational v2 taxonomies.
mdb <- read.csv("~/Library/Mobile Documents/com~apple~CloudDocs/PhD/Chapter1/PaperCh1-master/Dataset/Master_db.csv", sep = ";", stringsAsFactors = F)
mdb <- mdb[,c(3,5,7,9,10)]

# import dataset used in the analysis
d <- read_delim(file.choose(), delim = ',') %>% .[-which(.$GD == 0),]

# remove sequences that have been discarded from the calculation of GD
data2 <- mdb[-which(mdb$NOTES == "EX" | mdb$NOTES == "DOM" | mdb$NOTES == "HYB" | 
                     mdb$NOTES == "NN" | mdb$NOTES == "Newly discovered taxon" | 
                     mdb$NOTES == "IDENTICAL"),]

# remove duplicates (each row has one species)
data2 <- as.data.frame(cbind(gsub(" ", "_", data2$GB_SPECIES), gsub(" ", "_", data2$IOC_SPECIES), gsub(" ", "_", data2$JETZ_NEW), data2$RED_LIST_JETZ.2017., data2$NOTES))
colnames(data2) <- c("GB_SPECIES", "IOC_SPECIES", "JETZ_NEW", "RED_LIST_JETZ.2017.", "NOTES")
data2 <- data2[!duplicated(data2),]

# keep only the species analysed
final <- data2[match(d$SP, data2$JETZ_NEW),]

# import realms datasets
path <- "~/Library/Mobile Documents/com~apple~CloudDocs/PhD/Chapter1/PaperCh1-master/Dataset/bird realm spp/"
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




