library(phytools)
library(diversitree)
library(geiger)
library(ggtree)

# import the data
data_birds <- read.csv("../Dataset/Master_db.csv", header = T, sep = ";")
data <- data_birds[-which(data_birds$NOTES == "EX" | data_birds$NOTES == "DOM" | data_birds$NOTES == "HYB" | 
                            data_birds$NOTES == "NN" | data_birds$NOTES == "Newly discovered taxon"),]
data <- data[-which(is.na(data$JETZ_TAXONOMY)),]
data <- data_birds[,c(5,6,7,8,9)]
# remove exctinct, domesticated, hybrid species, species with uncertain nomenclature, newly discovered species 
# that are not included in the phylogeny

data$tiplabel <- NA
data$tiplabel <- gsub(' ', '_', data$JETZ_TAXONOMY) # create tip label --> species names as "Genus_species"
data <- data[,c(1,2,6,3,4,5)] # re-order

length(which(data$RED_LIST_IOC.2017. != data$RED_LIST_JETZ.2017.))
diff_IUCN <- data[(which(data$RED_LIST_IOC.2017. != data$RED_LIST_JETZ.2017.)),]
data <- data[-which(data$RED_LIST_IOC.2017. != data$RED_LIST_JETZ.2017.),]
# remove species with differences in categories due to differences in taxonomy

db <- data[!duplicated(data$tiplabel),] 
# gives a subset of the table that is without duplicates of the variables in the column specified

seqs <- as.data.frame(table(data$tiplabel)) # create database with number of sequences per species
seqs <- seqs[order(seqs$Var1),]
db <- db[order(db$tiplabel),]

db$seqsnum <- seqs$Freq # add number of sequences per species to the dataset

subdb <- subset(db, db$seqsnum > 1)
subdb <- subdb[,c(2,3,5,7)]
# write.csv(subdb, "./Birds/Dataset/subdb_phylogeny_R.csv", row.names = F, quote = F)
# take only species that have more than 1 sequence

# import phylogeny
ladtree <- ladderize(read.nexus("../Dataset/MCC_Hackett2_Full.tre"))


# look if species names are present in tip.labels of tree
which(subdb$tiplabel %in% ladtree$tip.label == FALSE)
length(which(ladtree$tip.label%in%subdb[,2]))##all match

clcolr <- rep("grey90", dim(ladtree$edge)[1]) 

for(i in 1:length(subdb[,1])){
  correct.tip<-which(ladtree$tip.label%in%subdb[i,2]) ##matches the tip to the column name of the output table
  # get number of tip
  x<-which.edge(ladtree, correct.tip) # get node number from the number of tip
  clcolr[x] <- "red"
}

jetz.tip.df<-matrix(NA,ncol=1,nrow=length(ladtree$tip.label))
rownames(jetz.tip.df)<-ladtree$tip.label

for(i in 1:length(ladtree$tip.label)){
  
  if(rownames(jetz.tip.df)[i]%in%subdb[,2]){
    jetz.tip.df[i,1]<-"red"
  }
  
  else{
    jetz.tip.df[i,1]<-"white"
  }
  
}

#pdf("Jetz_GeneticDiv_Samples_EdgeTip.pdf",height=7,width=7)
plot(ladtree,edge.color=clcolr,show.tip.label=FALSE)
segments(121, 1:nrow(jetz.tip.df), 123 + 1, 1:nrow(jetz.tip.df),col=jetz.tip.df[,1], lwd = 0.1)
#dev.off()

subdb <- subdb[-which(subdb$RED_LIST_JETZ.2017. == "EW"),]
subdb <- subdb[-which(subdb$RED_LIST_JETZ.2017. == "NR"),]
subdb <- subdb[-which(is.na(subdb$RED_LIST_JETZ.2017.)),]

subdb$coltip <- NA
subdb[which(subdb$RED_LIST_JETZ.2017. == "LC"),5] <- "#C1FFC1"
subdb[which(subdb$RED_LIST_JETZ.2017. == "NT"),5] <- "#C1FFC1"
subdb[which(subdb$RED_LIST_JETZ.2017. == "VU"),5] <- "#FFD700"
subdb[which(subdb$RED_LIST_JETZ.2017. == "EN"),5] <- "#EE7600"
subdb[which(subdb$RED_LIST_JETZ.2017. == "CR"),5] <- "#CD2626"
subdb[which(subdb$RED_LIST_JETZ.2017. == "DD"),5] <- "grey90"
subdb[which(subdb$RED_LIST_JETZ.2017. == "CR (PE)"),5] <- "#CD2626"
subdb[which(subdb$RED_LIST_JETZ.2017. == "CR (PEW)"),5] <- "#CD2626"

jetz.tip.df2<-matrix(NA,ncol=1,nrow=length(ladtree$tip.label))
rownames(jetz.tip.df2)<-ladtree$tip.label

for(i in 1:length(ladtree$tip.label)){
  
  if(rownames(jetz.tip.df2)[i]%in%subdb[,2]){
    sub <- subset(subdb, subdb$tiplabel == rownames(jetz.tip.df2)[i])
    jetz.tip.df2[i,1]<-sub$coltip
  }
  
  else{
    jetz.tip.df2[i,1]<-"white"
  }
  
}

plot(ladtree,edge.color=clcolr,show.tip.label=FALSE)
segments(121, 1:nrow(jetz.tip.df2), 123 + 1, 1:nrow(jetz.tip.df2),col=jetz.tip.df2[,1], lwd = 0.1)


# END #
#####################################################################################################


pruned.tree <- drop.tip(ladtree, ladtree$tip.label[-match(subdb$tiplabel, ladtree$tip.label)])
# tree with only the species in my data


# colour branches of the tree that end in the species present in my data
wh <- which.edge(ladtree, c(match(subdb$tiplabel, ladtree$tip.label)))
# match function gives tips numbers --> which.edge gives edge numbers linked to tips numbers
colo <- rep("grey90", dim(ladtree$edge)[1]) # colour all the branches in black
colo[wh] <- "red" # the edge numbers relative to my data are coloured in red

p <- plot(ladtree, type= 'fan', edge.color = colo, show.tip.label = F)
# difficult to understand when visualised

## other way to plot it
myd <- ladtree$tip.label[c(match(subdb$tiplabel, ladtree$tip.label))]
# get tiplabel associated to my data
group <- groupOTU(ladtree, myd)
# create a phylogeny with a group representing my data
pl <- ggtree(group, aes(color=group), layout = 'circular')
# plot tree with branches for my data 



##### TRYING TO ANNOTATE FAMILIES TO PLOT

fam <- split(subdb$tiplabel, as.factor(subdb$FAMILY))
# split species into their families --> creates a list of families with correspondent species
fam <- fam[2:183] # remove #N/A

MRCA(ladtree, tip = c(fam[[1]])) #gives the node of the species in one family

# obtain the nodes for all the families
nodes <- c()
for(s in seq_along(names(fam))){
  v <- MRCA(ladtree, tip = c(fam[[s]]))
  if (is.null(v) == T) {
    nodes <- c(nodes, "NULL") # I paste NULL so the number of nodes is the same as the number of species
  } else{
    nodes <- c(nodes, v)
  }
  cat(s, names(fam[s]), v, "\n")
}

nodes

df <- as.data.frame(matrix(nrow = length(nodes))) # both names(fam) and nodes are of the same length
df$fam <- names(fam)
df$nodes <- nodes
df <- df[-which(df$nodes == "NULL"),] # remove NULLs
df <- df[,2:3]

# pl + geom_cladelabel2(node = c(df$nodes), label = c(df$fam), offset = .8)
