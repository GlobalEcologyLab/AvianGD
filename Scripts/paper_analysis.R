## Main analysis in paper ##
library(ape)
library(phytools)
library(caper)
library(rcompanion)
library(vioplot)
library(tidyverse)

### --- Load data --- ###
# output table from Matlab (dataset_paper_analysis.R)
data <- read_delim(file.choose(), delim = ',')
# phylogenetic tree (Jetz et al. 2012)
tree <- ladderize(read.nexus(file.choose()))

# remove zero values
# zeros are just an artefact of having sequences with different lengths and 
# of the methodology used to calculate genetic diversity
d <- data[-which(data$GD == 0),]
d <- data.frame(d)
rownames(d) <- d$SP
phylo <- drop.tip(tree, tree$tip.label[-match(d$SP, tree$tip.label)])

# calculate phylogenetic signal
GD <- d$GD
names(GD) <- d$SP
sig <- phylosig(phylo, GD, method="lambda", test=TRUE)
sig$lambda
sig$P

# transform values to reach normality
d$GDt = transformTukey(d$GD,plotit=TRUE)
#d$GDt <- d$GD^0.175

# add Threatened and Non-threatened conservation status
d$status <- ifelse(d$IUCN == "LC" | d$IUCN == "NT", "Non-threatened", "Threatened")
table(d$status)

### --- Model, resampling with modified phylANOVA --- ###
source("~/Library/Mobile Documents/com~apple~CloudDocs/PhD/Chapter1/GeneticGap/Scripts/resampling_function.R")

TR <- d %>% filter(status=="Threatened") %>% .$GDt
names(TR) <- d %>% filter(status=="Threatened") %>% .$SP

NT <- d %>% filter(status=="Non-threatened") %>% .$GDt
names(NT) <- d %>% filter(status=="Non-threatened") %>% .$SP

m <- resamp_phylanova(nonthreatened = NT, threatened = TR, nsamp = 50, nreps = 1000)

P_val <- vector()
for(i in 1:1000){
  pv <- m %>% .[[3]] %>% .[[i]] %>% .[[4]] %>% .$Pf
  P_val[i] <- pv 
}

length(which(P_val < 0.05))
length(which(P_val >= 0.05))
mean(P_val)
sd(P_val)
mean(P_val) + sd(P_val)
median(P_val)
mad(P_val, constant = 1)

F.obs_vector <- vector()
for(i in 1:1000){
  fval <- m %>% .[[3]] %>% .[[i]] %>% .[[4]] %>% .$F
  F.obs_vector[i] <- fval 
}

F.null_list <- list()
for (l in 1:1000) {
  fn <- m %>% .[[3]] %>% .[[l]] %>% .[[3]]
  F.null_list[[l]] <- fn
}

F_df <- data.frame(F.v = unlist(F.null_list))
F_df$col <- 'Fnull'
F_obs_df <- data.frame(F.v = F.obs_vector)
F_obs_df$col <- 'Fobs'

F_hist <- rbind(F_df, F_obs_df)

## Effect size ##
SS <- list()
for(i in 1:1000){
  ss <- m %>% .[[3]] %>% .[[i]] %>% .[[4]] %>% .$`Sum Sq`
  SS[[i]] <- ss 
}

eta.sq <- lapply(SS, function(x){(x[1]/sum(x))})
hist(unlist(eta.sq), breaks=100)
mean(unlist(eta.sq))
sd(unlist(eta.sq))
max(unlist(eta.sq))
median(unlist(eta.sq))
mad(unlist(eta.sq), constant=1)

# distributions plot
p <- ggplot(F_hist, aes(F.v, fill=col)) + 
  geom_density(alpha=0.6) + xlim(0,55) +
  ggtitle("Observed vs. Null F-distribution") + 
  xlab("F-value") + guides(fill=guide_legend(title = "Legend"))

p + theme(axis.text = element_text(size=20), axis.title = element_text(size=24),
          axis.title.x = element_text(margin = margin(30,0,0,0)),
          axis.title.y = element_text(margin = margin(0,30,0,0)),
          plot.title = element_text(size=30, hjust = 0.5, margin = margin(0,0,30,0)),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20, face = 'bold'),
          legend.box.background = element_rect(colour = "black"),
          panel.border = element_blank(),
          panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))

# violin plot with the whole dataset
ggplot(d, aes(factor(status), GDt)) + geom_violin(aes(fill=factor(status))) + 
  scale_fill_manual(values=c("#009E73", "#D55E00")) + geom_boxplot(width=0.2, fill="grey90") + 
  xlab("Status") + ggtitle("Genetic diversity between groups") + theme_bw() + 
  # theme(axis.text = element_blank(),
  #       axis.title = element_blank(),
  #       plot.title = element_blank(),
  #       legend.position = "none") +
  theme(axis.text = element_text(size=20),
    axis.title.x = element_text(margin = margin(30,0,0,0)),
    axis.title.y = element_text(margin = margin(0,30,0,0)),
    axis.title = element_text(size=24, margin = margin(0,0,30,0)),
    plot.title = element_text(size=30, hjust = 0.5, margin = margin(0,0,30,0)),
    legend.position = "none") +
  scale_x_discrete(labels=c("Non-threatened", "Threatened"))

### --- Calculating percentiles --- ###
gd_pct <- d %>% .$GD %>% quantile(., prob = seq(0, 1, length = 21), type = 1)

# number of species in the <5th percentile
d %>% filter(status == "Threatened") %>% 
  filter(SEQS > 10) %>% 
  filter(GD <= gd_pct['5%']) %>% 
  nrow()
d %>% filter(status == "Non-threatened") %>% 
  filter(SEQS > 10) %>% 
  filter(GD <= gd_pct['5%']) %>% 
  nrow()
# number of species in the >95th percentile
d %>% filter(status == "Threatened") %>% 
  filter(SEQS > 10) %>% 
  filter(GD >= gd_pct['95%']) %>% 
  nrow()
d %>% filter(status == "Non-threatened") %>% 
  filter(SEQS > 10) %>% 
  filter(GD >= gd_pct['95%']) %>% 
  nrow()

# number of species in the <10th percentile
d %>% filter(status == "Threatened") %>% 
  filter(SEQS > 10) %>% 
  filter(GD <= gd_pct['10%']) %>% 
  nrow()
d %>% filter(status == "Non-threatened") %>% 
  filter(SEQS > 10) %>% 
  filter(GD <= gd_pct['10%']) %>% 
  nrow()
# number of species in the >90th percentile
d %>% filter(status == "Threatened") %>% 
  filter(SEQS > 10) %>% 
  filter(GD >= gd_pct['90%']) %>% 
  nrow()
d %>% filter(status == "Non-threatened") %>% 
  filter(SEQS > 10) %>% 
  filter(GD >= gd_pct['90%']) %>% 
  nrow()

# number of species in the <20th percentile
d %>% filter(status == "Threatened") %>% 
  filter(SEQS > 10) %>% 
  filter(GD <= gd_pct['20%']) %>% 
  nrow()
d %>% filter(status == "Non-threatened") %>% 
  filter(SEQS > 10) %>% 
  filter(GD <= gd_pct['20%']) %>% 
  nrow()
# number of species in the >80th percentile
d %>% filter(status == "Threatened") %>% 
  filter(SEQS > 10) %>% 
  filter(GD >= gd_pct['80%']) %>% 
  nrow()
d %>% filter(status == "Non-threatened") %>% 
  filter(SEQS > 10) %>% 
  filter(GD >= gd_pct['80%']) %>% 
  nrow()

# species in each percentile (<10th and >90th)
d %>% filter(status == "Threatened") %>% 
  filter(SEQS > 10) %>% 
  arrange(GD) %>% 
  filter(GD <= gd_pct['10%']) %>% 
  .$SP %>% 
  gsub("_", " ", .)
d %>% filter(status == "Non-threatened") %>% 
  filter(SEQS > 10) %>% 
  arrange(GD) %>% 
  filter(GD <= gd_pct['10%']) %>% 
  .$SP %>% 
  gsub("_", " ", .)
d %>% filter(status == "Threatened") %>% 
  filter(SEQS > 10) %>% 
  arrange(GD) %>% 
  filter(GD >= gd_pct['90%']) %>% 
  .$SP %>% 
  gsub("_", " ", .)
d %>% filter(status == "Non-threatened") %>% 
  filter(SEQS > 10) %>% 
  arrange(GD) %>% 
  filter(GD >= gd_pct['90%']) %>% 
  .$SP %>% 
  gsub("_", " ", .)

### --- correlation tests number of sequences and sequence lengths --- ###
library(ggpubr)
#mdb <- read.csv("~/Documents/PhD/Chapter1/22May2020/Data/Master_db.csv", stringsAsFactors = F)
seq_lengths <- read.csv("~/Documents/PhD/Chapter1/22May2020/Data/SpeciesSequenceLengths.csv", stringsAsFactors = F)
seq_lengths$SP <- seq_lengths$SP %>% str_split(., "_") %>% purrr::map(~.[1:2]) %>% purrr::map(lift(paste), sep="_") %>% unlist()

d <- d %>% merge(., seq_lengths, by="SP")
all.equal(d$SEQS.x, d$SEQS.y)
d <- d[, -ncol(d)]
d <- rename(d, "SEQS"=SEQS.x)
d$AVG_LENGTH <- apply(d[,9:10], 1, mean) %>% round()

round(mean(d$AVG_LENGTH)); round(mean(d$MIN_LENGTH)); round(mean(d$MAX_LENGTH)); round(sd(d$AVG_LENGTH))

#plots
p1 <- ggscatter(data = d, y="GDt", x="SEQS", ylab="Genetic diversity", xlab="Number of Sequences", cor.coef = T, cor.method = "kendall",
                cor.coef.coord = c(380, 0.1), add = 'reg.line', conf.int = T, add.params = list(color = 'red'))
p2 <- ggplot(data=d, aes(x=SEQS)) + 
  geom_histogram(fill='grey80', col='black', bins = 100) + 
  xlab("Number of Sequences") +
  ylab("Count") +
  theme_pubr()


p4 <- ggplot(data = d, aes(x=AVG_LENGTH)) +
  geom_histogram(fill='grey80', col='black', bins = 100) + 
  xlab("Average Sequence Length") +
  ylab("Count") +
  theme_pubr()

p3 <- ggscatter(data = d, y="GDt", x="AVG_LENGTH", ylab = "Genetic diversity", xlab="Average Sequence Length", cor.coef = T, 
                cor.method = "kendall", cor.coef.coord = c(850, 0.1), add = 'reg.line', conf.int = T, add.params = list(color = 'red'))

ggarrange(p1, NULL,  p2, NULL, NULL, NULL, p3, NULL, p4, ncol = 3, nrow = 3, widths = c(1, 0.2, 1), heights = c(1,0.2, 1), 
          labels = c('A', '', 'B', '', '', '', 'C', '', 'D'))

