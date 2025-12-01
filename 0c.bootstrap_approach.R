# packages
source("0b.functions.R")
library(lme4)
library(glmmTMB)
library(tidyverse)
library(gridExtra)
# run on parallel cores
library(future)
library(future.apply)

########################################################
#### Approach to obtain parametric bootstrap values ####
########################################################

bblist_D0 <- list() # example for 0D, the species/taxon richness, but idem for all other hypotheses

df <- readRDS("diversitylist_D0.rds")[["nematodes_q0"]] %>% # example with nematode 0D (taxon richness)
  filter(!is.na(D0))
fitted_model <- glmmTMB(D0 ~ landuse + (1|location),
                        family = Gamma(link="log"),
                        data = df)

plan(multisession, workers = parallel::detectCores() - 1)

# parametric bootstrapping
bb <- bootMer(fitted_model,
              FUN = function(x) {
                ff <- fixef(x)[["cond"]] 
                ff 
              },
              nsim = 5000,
              parallel = "future",
              ncpus = parallel::detectCores() - 1)


bblist_D0[["nematodes"]] <- bb 


########################################################
##### Approach to check parametric bootstrap values ####
########################################################

# Look at the histograms of the bootstrap values to see if they are approximately normally distributed
# They should be (approximately) normally distributed.
# In some cases there may be some outlying values, 
# so we also show the histograms for the data without these outliers (trimming at 1.5x the interquartile range above and below the first and third quartile).

for(m in names(bblist_D0)){
  bb_t <- as.data.frame(bblist_D0[[m]]$t) %>%
    pivot_longer(cols = everything())
  
  for(n in unique(bb_t$name)){
    # overall histograms
    print(bb_t %>% 
            filter(name == n) %>%
            ggplot(aes(x = value)) +
            geom_histogram() +
            labs(title = paste0(m, "_", n)))
    # histograms without outliers (using trim function)
    print(bb_t %>%
            filter(name == n) %>%
            trim(value) %>%
            ggplot(aes(x = value)) +
            geom_histogram() +
            labs(title = paste0(m, "_", n, "_trimmed")))
  }
}  


########################################################
############ Save parametric bootstrap values ##########
########################################################

saveRDS(bblist_D0, "bblist_D0.rds")
