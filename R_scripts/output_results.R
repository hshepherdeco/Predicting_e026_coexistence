rm(list=ls())

library(tidyverse)
library(tibble)

setwd("~/SimulateCoexistence")

devtools::load_all()
### load datasets
load("data/total_biomass_data.rda")
load("data/trait_data.rda")
load("switches.rda")


##### Function to collate results from simulated outputs  
collate_results <- function(files, vari = "replicates", year = 6) { #

  sp_names <- unlist(dimnames(biomass_data_total)[3])

  bstar0 <- get_trait_vector(trait_data, sp_names, "Bm") %>% as_tibble(rownames = "species")

  fnames <- files$full_name

  # year, subplot, species, plot
  results <- map_dfr(fnames, function(x) {
    dat <- readRDS(x)
    dat <- dat[year+1,,,]
    dat <- as.data.frame.table(dat) %>%
      mutate(subplot = as.numeric(Var1),
             plot = as.numeric(Var3),
             species = factor(Var2, labels = sp_names),
             replicate = x) %>%
      select(subplot, plot, species, pred = Freq, replicate)
    return(dat)}
  )

  # divide prediction by bstar0 for each species to make them comparable
  results <- results %>% inner_join(bstar0) %>%
    mutate(pred = pred/value)

  obs <- biomass_data_total[year,,,]
  obs <- as.data.frame.table(obs) %>%
    mutate(subplot = as.numeric(Var1),
           plot = as.numeric(Var3),
           species = factor(Var2, labels = sp_names)) %>%
    # divide observed by bstar0 for each species to make them comparable
    inner_join(bstar0) %>%
    mutate(obs = Freq/value) %>%
    select(subplot, plot, species, obs) %>%
    na.omit

  planted_richness = apply(presence_data_total, c(1,3), sum)
  planted_richness <- as.data.frame.table(planted_richness) %>%
    mutate(subplot = as.numeric(Var1),
           plot = as.numeric(Var2)) %>%
    select(subplot, plot, sp_rich = Freq)



  df <- results %>%
    inner_join(obs) %>%
    inner_join(planted_richness) %>%
    mutate(diff = obs - pred,
           rmse = sqrt(diff^2))

  return(df)
}

# se takes n as the number of plots as the replicates inflate the SE
se <- function(x, num_plots) sd(x, na.rm = TRUE) / sqrt(num_plots)


#### Note - the next section of code requires simulation outputs that are not provided due to their large data size
#### Code to produce the simulations is provided in file: xxxxx.R

args <- as.data.frame(list.files("results"))  %>% separate(1, into = c("x", "y"), sep="_") %>%
  distinct() %>%
  as.matrix()

sp_names <- unlist(dimnames(biomass_data_total)[3])
dat <- as.data.frame.table(seeding_data_total) %>%
  mutate(subplot = as.numeric(Var1),
         plot = as.numeric(Var3),
         species = factor(Var2, labels = sp_names)) %>%
  select(subplot, plot, species, seed_mass = Freq) %>% 
  filter(seed_mass != 0)

for(i in 0:11) {

f <- list.files("results_full", pattern = paste0("^", i, "_"))
f_full <- list.files("results_full", pattern = paste0("^", i, "_"), full.names = TRUE)
f_params <- list.files("params_full", pattern = paste0("^", i, "_"), full.names = TRUE)

files <- tibble(fname = f, full_name = f_full, params = f_params) %>%
  separate(fname, into = c("mod_name", "replicate"), sep = "_(?=[^_]+$)") %>%
  mutate(switches = map(params, function(x) {
    sw <- readRDS(x)
    sw <- sort(sw)
    sw <- paste(sw, collapse = "_")
    return(sw)
  }
  )
  ) %>%
  unnest(cols = c(switches))

non_zero_replicates <- files %>%
  group_by(switches) %>%
  nest() %>%
  mutate(output = map(data, function(x) {
    collate_results(x, vari = "replicates")
    })) %>%
  unnest(output) %>%
  select(-data) %>%
  mutate(nswitch = i) %>%
  inner_join(dat)

save(out_replicates, file = paste0("Data/non_zero_replicates", i, "year6.rda"))

}

### create combined files document

non_zero_replicates <- do.call(rbind, data_list)
remove(data_list)


for(i in 0:11) {
  
  f <- list.files("results_full", pattern = paste0("^", i, "_"))
  f_full <- list.files("results_full", pattern = paste0("^", i, "_"), full.names = TRUE)
  f_params <- list.files("params_full", pattern = paste0("^", i, "_"), full.names = TRUE)
  
  files <- tibble(fname = f, full_name = f_full, params = f_params) %>%
    separate(fname, into = c("mod_name", "replicate"), sep = "_(?=[^_]+$)") %>%
    mutate(switches = map(params, function(x) {
      sw <- readRDS(x)
      sw <- sort(sw)
      sw <- paste(sw, collapse = "_")
      return(sw)
    }
    )
    ) %>%
    unnest(cols = c(switches))
 
  data_list[[i + 1]] <- files
   
} 


files <- do.call(rbind, data_list)
save(files, file = "Data/files_combined.rda")
