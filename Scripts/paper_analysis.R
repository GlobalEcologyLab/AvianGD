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
d <- d[-which(is.na(d$GD)),]
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
T_tuk = transformTukey(d$GD,plotit=TRUE)
d$GDt <- d$GD^0.125

# add Threatened and Non-threatened conservation status
d$status <- ifelse(d$IUCN == "LC" | d$IUCN == "NT", "Non-threatened", "Threatened")
table(d$status)

### --- Model, resampling with modified phylANOVA --- ###
source("~/Library/Mobile Documents/com~apple~CloudDocs/PhD/Chapter1/PaperCh1-master/Scripts/resampling_function.R")

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
d %>% .$GD %>% quantile(., prob = seq(0, 1, length = 21), type = 1)

# number of species in the <5th percentile
d %>% filter(status == "Threatened") %>% filter(SEQS > 10) %>% filter(GD <= 0.0009685561) %>% nrow()
d %>% filter(status == "Non-threatened") %>% filter(SEQS > 10) %>% filter(GD <= 0.0009685561) %>% nrow()
# number of species in the >95th percentile
d %>% filter(status == "Threatened") %>% filter(SEQS > 10) %>% filter(GD >= 0.0424227955) %>% nrow()
d %>% filter(status == "Non-threatened") %>% filter(SEQS > 10) %>% filter(GD >= 0.0424227955) %>% nrow()


# number of species in the <10th percentile
d %>% filter(status == "Threatened") %>% filter(SEQS > 10) %>% filter(GD <= 0.0015759963) %>% nrow()
d %>% filter(status == "Non-threatened") %>% filter(SEQS > 10) %>% filter(GD <= 0.0015759963) %>% nrow()
# number of species in the >90th percentile
d %>% filter(status == "Threatened") %>% filter(SEQS > 10) %>% filter(GD >= 0.0327944040) %>% nrow()
d %>% filter(status == "Non-threatened") %>% filter(SEQS > 10) %>% filter(GD >= 0.0327944040) %>% nrow()


# number of species in the <20th percentile
d %>% filter(status == "Threatened") %>% filter(SEQS > 10) %>% filter(GD <= 0.0025488466) %>% nrow()
d %>% filter(status == "Non-threatened") %>% filter(SEQS > 10) %>% filter(GD <= 0.0025488466) %>% nrow()
# number of species in the >80th percentile
d %>% filter(status == "Threatened") %>% filter(SEQS > 10) %>% filter(GD >= 0.0221316500) %>% nrow()
d %>% filter(status == "Non-threatened") %>% filter(SEQS > 10) %>% filter(GD >= 0.0221316500) %>% nrow()

# species in each percentile (<10th and >90th)
d %>% filter(status == "Threatened") %>% filter(SEQS > 10) %>% arrange(GD) %>% filter(GD <= 0.0015759963) %>% .$SP %>% gsub("_", " ", .)
d %>% filter(status == "Non-threatened") %>% filter(SEQS > 10) %>% arrange(GD) %>% filter(GD <= 0.0015759963) %>% .$SP %>% gsub("_", " ", .)
d %>% filter(status == "Threatened") %>% filter(SEQS > 10) %>% arrange(GD) %>% filter(GD >= 0.0327944040) %>% .$SP %>% gsub("_", " ", .)
d %>% filter(status == "Non-threatened") %>% filter(SEQS > 10) %>% arrange(GD) %>% filter(GD >= 0.0327944040) %>% .$SP %>% gsub("_", " ", .)


