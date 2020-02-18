## CREATION OF DATASET TO INPUT IN MATLAB ##
library(tidyverse)

# import the data
data_birds <- read.csv("../Dataset/Master_db.csv", header = T, sep = ";")
table(data_birds$NOTES)
data <- data_birds[-which(data_birds$NOTES == "EX" | data_birds$NOTES == "DOM" | data_birds$NOTES == "HYB" | 
                            data_birds$NOTES == "NN" | data_birds$NOTES == "Newly discovered taxon" | 
                            data_birds$NOTES == "IDENTICAL"),]

table(data$NOTES)
data <- data[-which(is.na(data$JETZ_NEW)),]
data <- data[,c(1,7,9,10,11,12)]

# write file
write_csv(data, "..//GD/matlab_input.csv")