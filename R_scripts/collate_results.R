devtools::load_all("../simulateCoexistence/")
library(tidyverse)

args = commandArgs(TRUE)

sw = args[1]

get_mae_year <- function(fname) {
  d = readRDS(fname)  

  ## compare goodness of fit
  # take mean across years
  mae = apply(abs(d[-1,,,]-biomass_data_total), 2:4, function(x) mean(x, na.rm=T))
  
  mae_df <- as.data.frame.table(mae)
  
  names(mae_df) <- c("subplot", "species", "plot", "mae")
  
  ## get rid of methodological zeros (meaningless)
  dat <- as.data.frame.table(seeding_data_total)
  names(dat) <- c("subplot", "species", "plot", "seed_mass")
  dat <- filter(dat, seed_mass != 0)
  
  mae_df <- inner_join(mae_df, dat)
  
  planted_richness <- apply(presence_data_total, c(1,3), sum)
  planted_richness <- as.data.frame.table(planted_richness)
  names(planted_richness) <- c("subplot", "plot", "sp_rich")
  
  mae_df <- inner_join(mae_df, planted_richness)
  
}

f <- list.files("../simulateCoexistence/results", pattern = paste0("^", sw, "_"))
f_full <- list.files("../simulateCoexistence/results", pattern = paste0("^", sw, "_"), full.names = TRUE)
f_params <- list.files("../simulateCoexistence/params", pattern = paste0("^", sw, "_"), full.names = TRUE)

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
  unnest(switches)

out_replicates <- files %>%
  mutate(output = map(full_name, function(x) {
    get_mae_year(x)
  })) %>% 
  unnest(output) %>% 
  select(mod_name, switches, subplot, species, plot, sp_rich, mae)

save(out_replicates, file = paste0("results_replicates_", sw, ".rda"))

# non-zero results (we are good at predicting zeros, what about the rest)
dat <- as.data.frame.table(seeding_data_total)
names(dat) <- c("subplot", "species", "plot", "seed_mass")


non_zero_replicates <- out_replicates %>% inner_join(dat %>% filter(seed_mass != 0))

save(non_zero_replicates, file = paste0("non_zero_replicates_", sw, ".rda"))
