## PhylANOVA function ##

source("~/Library/Mobile Documents/com~apple~CloudDocs/PhD/Chapter1/PaperCh1-master/Scripts/phylANOVA_modified.R")

# the function here integrates resampling and phylANOVA. 

resamp_phylanova <- function(nonthreatened, threatened, nsamp, nreps) {
  samples_nt <- rerun(nreps, sample(nonthreatened, nsamp))
  samples_t <- rerun(nreps, sample(threatened, nsamp))
  species_names <- map2(samples_nt, samples_t, ~c(names(.x), names(.y)))
  phylos <- purrr::map(species_names, ~drop.tip(tree, tree$tip.label[-match(.x, tree$tip.label)]))
  group_vectors <- purrr::map(species_names, 
                       ~set_names(c(rep("Non-threatened", nsamp), rep("Threatened", nsamp)),
                                  .x))
  values <- list(samples_nt, samples_t) %>% transpose() %>% purrr::map(flatten_dbl)
  outputs <- pmap(list(tree = phylos, x = group_vectors, y = values), phylANOVA_modified)
  return(list(group_vectors, values, outputs))
}

# nonthreatened = named vector with values of genetic diversity for non-threatened species
# threatened = named vector with values of genetic diversity for threatened species
# nsamp = number of samples
# nreps = number of repetitions