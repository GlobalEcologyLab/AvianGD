## CREATION OF DATASET FOR ANALYSES ##
library(tidyverse)

#load data tables
matlab <- read_delim(file.choose(), delim = ',') #dataset_for_matlab.csv
output <- read_delim(file.choose(), delim=',') #matlab_output.csv

# Column names were added manually to the output table from Matlab. The columns are:
# SP = species name
# GD = nucleotide diversity value
# NUM_MUT = number of mutations
# DROPPED_PAIRS = number of paired sequences not included in the calculation
# SEQS = total number of sequences for the species

# assign IUCN categories to matlab output table
iucn <- unique(matlab[,2:3])
iucn$JETZ_NEW <- gsub(' ', '_', iucn$JETZ_NEW)
d <- cbind(output,iucn[match(output$SP, iucn$JETZ_NEW),2])
colnames(d) <- c(colnames(d)[1:5], "IUCN")

# remove species with categories EW, NR and DD
table(d$IUCN)
d <- d[-which(d$IUCN == "EW" | d$IUCN == "NR" | d$IUCN == "DD"),]
d <- d[-which(is.na(d$IUCN)),]

# group into CR category
d[which(d$IUCN == "CR (PE)" | d$IUCN == "CR (PEW)"), 6] <- "CR"
table(d$IUCN)

# remove NA values of GD
d <- d[-which(is.nan(d$GD)),]

# keep only species with >5 sequences
d <- d[which(d$SEQS > 5),]

# write file
write_delim(d, "~/Documents/PhD/Chapter1/22May2020/Data/dataset_analysis.csv", delim = ',')
