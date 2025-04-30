rm(list=ls())
setwd("simulateCoexistence/")

devtools::load_all()
library(tidyverse)
library(tibble)
### load datasets
load("data/total_biomass_data.rda")
load("data/trait_data.rda")
load("switches.rda")
### Results directory

setwd("simulateCoexistence/")


collate_results <- function(files, vari = "replicates", year = 6) {

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

  # if(vari == "replicates") {
  #      df <- df %>%
  #    group_by(species, sp_rich, plot, subplot) %>%
  #    summarise(mean_diff = mean(diff, na.rm = TRUE),
  #             sd_diff = sd(diff, na.rm = TRUE),
  #              mean_rmse = mean(rmse, na.rm = TRUE),
  #               sd_rmse = sd(rmse, na.rm = TRUE),
  #             # r2 = 1-mean((obs - pred)^2, na.rm=TRUE)/mean((obs - mean(obs, na.rm=T))^2, na.rm=TRUE)
  #     )
  # }
  
  # if(vari == "plots") {
  #   df <- df %>%
  #     group_by(species, sp_rich, replicate) %>%
  #     summarise(mean_diff = mean(diff, na.rm = TRUE),
  #               sd_diff = sd(diff, na.rm = TRUE),
  #               mean_rmse = mean(rmse, na.rm = TRUE),
  #               sd_rmse = sd(rmse, na.rm = TRUE),
  #               #r2 = 1-mean((obs - pred)^2, na.rm=TRUE)/mean((obs - mean(obs, na.rm=T))^2, na.rm=TRUE)
  #     )
  # }

  return(df)
}

# se takes n as the number of plots as the replicates inflate the SE
se <- function(x, num_plots) sd(x, na.rm = TRUE) / sqrt(num_plots)
setwd("C:/Users/k2256287/OneDrive - King's College London/Documents/WP5/Original_manuscript")

#f <- list.files("results")
#f_full <- list.files("results", full.names = TRUE)
#f_params <- list.files("params", full.names = TRUE)

args <- as.data.frame(list.files("results"))  %>% separate(1, into = c("x", "y"), sep="_") %>%
  distinct() %>%
  as.matrix()

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

out_replicates <- files %>%
  group_by(switches) %>%
  nest() %>%
  mutate(output = map(data, function(x) {
    collate_results(x, vari = "replicates")
    })) %>%
  unnest(output) %>%
  select(-data)

save(out_replicates, file = paste0("results_replicates_131224", i, ".rda"))

}

### Combine to make one full results_replicates document

load("results_replicates_0_year6.rda")
R0 <- out_replicates %>% mutate(nswitch = 0)

load("results_replicates_1_year6.rda")
R1 <- out_replicates %>% mutate(nswitch = 1)

load("results_replicates_2_year6.rda")
R2 <- out_replicates %>% mutate(nswitch = 2)

load("results_replicates_3_year6.rda")
R3 <- out_replicates %>% mutate(nswitch = 3)

load("results_replicates_4_year6.rda")
R4 <- out_replicates %>% mutate(nswitch = 4)

load("results_replicates_5_year6.rda")
R5 <- out_replicates %>% mutate(nswitch = 5)

load("results_replicates_6_year6.rda")
R6 <- out_replicates %>% mutate(nswitch = 6)

load("results_replicates_7_year6.rda")
R7 <- out_replicates %>% mutate(nswitch = 7)

load("results_replicates_8_year6.rda")
R8 <- out_replicates %>% mutate(nswitch = 8)

load("results_replicates_9_year6.rda")
R9 <- out_replicates %>% mutate(nswitch = 9)

load("results_replicates_10.rda")
R10 <- out_replicates %>% mutate(nswitch = 10)

load("results_replicates_11_year6.rda")
R11 <- out_replicates %>% mutate(nswitch = 11)

out_replicates <- rbind(R0, R1, R2, R3, R4, R5,
                        R6, R7, R8, R9, R10, R11)  

remove(R0, R1, R2, R3, R4, R5,
       R6, R7, R8, R9, R10, R11)

save(out_replicates, file = "results_replicates_all_year6.rda")

#### Produce non-zero replicates of each mod number

# non-zero results (we are good at predicting zeros, what about the rest)
sp_names <- unlist(dimnames(biomass_data_total)[3])
dat <- as.data.frame.table(seeding_data_total) %>%
  mutate(subplot = as.numeric(Var1),
         plot = as.numeric(Var3),
         species = factor(Var2, labels = sp_names)) %>%
  select(subplot, plot, species, seed_mass = Freq) %>% 
  filter(seed_mass != 0)

# 0 switches
load("results_replicates_0_year6.rda")
non_zero_replicates <- out_replicates %>% inner_join(dat)
save(non_zero_replicates, file = "non_zero_replicates_0_year6.rda")

# 1 switch
load("results_replicates_1_year6.rda")
non_zero_replicates <- out_replicates %>% inner_join(dat)
save(non_zero_replicates, file = "non_zero_replicates_1_year6.rda")

# 2 switch
load("results_replicates_2_year6.rda")
non_zero_replicates <- out_replicates %>% inner_join(dat)
save(non_zero_replicates, file = "non_zero_replicates_2_year6.rda")

# 3 switch
load("results_replicates_3_year6.rda")
non_zero_replicates <- out_replicates%>% inner_join(dat) 
save(non_zero_replicates, file = "non_zero_replicates_3_year6.rda")

# 4 switch
load("results_replicates_4_year6.rda")
non_zero_replicates <- out_replicates %>% inner_join(dat)
save(non_zero_replicates, file = "non_zero_replicates_4.rda")

# 5 switch
load("results_replicates_5_year6.rda")
non_zero_replicates <- out_replicates %>% inner_join(dat) 
save(non_zero_replicates, file = "non_zero_replicates_5_year6.rda")

# 6 switch
load("results_replicates_6_year6.rda")
non_zero_replicates <- out_replicates %>% inner_join(dat) 
save(non_zero_replicates, file = "non_zero_replicates_6_year6.rda")

# 7 switch
load("results_replicates_7_year6.rda")
non_zero_replicates <- out_replicates %>% inner_join(dat) 
save(non_zero_replicates, file = "non_zero_replicates_7_year6.rda")

# 8 switch
load("results_replicates_8_year6.rda")
non_zero_replicates <- out_replicates %>% inner_join(dat) 
save(non_zero_replicates, file = "non_zero_replicates_8_year6.rda")

# 9 switch
load("results_replicates_9_year6.rda")
non_zero_replicates <- out_replicates %>% inner_join(dat) 
save(non_zero_replicates, file = "non_zero_replicates_9_year6.rda")

# 10 switch
load("results_replicates_10_year6.rda")
non_zero_replicates <- out_replicates %>% inner_join(dat)
save(non_zero_replicates, file = "non_zero_replicates_10_year6.rda")

# 11 switch
load("results_replicates_11_year6.rda")
non_zero_replicates <- out_replicates %>% inner_join(dat) 
save(non_zero_replicates, file = "non_zero_replicates_11_year6.rda")

#### Produce combined 

load("non_zero_replicates_0_year6.rda")
NZR0 <- non_zero_replicates 

load("non_zero_replicates_1_year6.rda")
NZR1 <- non_zero_replicates 

load("non_zero_replicates_2_year6.rda")
NZR2 <- non_zero_replicates 

load("non_zero_replicates_3_year6.rda")
NZR3 <- non_zero_replicates 

load("non_zero_replicates_4_year6.rda")
NZR4 <- non_zero_replicates 

load("non_zero_replicates_5_year6.rda")
NZR5 <- non_zero_replicates 

load("non_zero_replicates_6_year6.rda")
NZR6 <- non_zero_replicates 

load("non_zero_replicates_7_year6.rda")
NZR7 <- non_zero_replicates 

load("non_zero_replicates_8_year6.rda")
NZR8 <- non_zero_replicates 

load("non_zero_replicates_9_year6.rda")
NZR9 <- non_zero_replicates 

load("non_zero_replicates_10_year6.rda")
NZR10 <- non_zero_replicates 

load("non_zero_replicates_11_year6.rda")
NZR11 <- non_zero_replicates

non_zero_replicates <- rbind(NZR0, NZR1, NZR2, NZR3, NZR4, NZR5,
      NZR6, NZR7, NZR8, NZR9, NZR10, NZR11)

remove(NZR0, NZR1, NZR2, NZR3, NZR4, NZR5,
       NZR6, NZR7, NZR8, NZR9, NZR10, NZR11)

save(non_zero_replicates, file = "Non_zero_replicates_all_year6.rda")


# which switches are in the top 5 percent of outputs
load("Non_zero_replicates_all.rda")
top5_pc <- non_zero_replicates %>%
  #filter(sp_rich > 1) %>%
  group_by(switches) %>%
  mutate(ae = abs(diff)) %>%
  mutate(pc95 = quantile(ae, 0.05)) %>%
  filter(ae <= pc95)

save(top5_pc, file = paste0("top5pc_all.rda"))

top5_pc_nomono <- non_zero_replicates %>%
  filter(sp_rich > 1) %>%
  group_by(switches) %>%
  mutate(ae = abs(diff)) %>%
  mutate(pc95 = quantile(ae, 0.05)) %>%
  filter(ae <= pc95)

save(top5_pc_nomono, file = paste0("top5pc_nomono.rda"))

# which mechanisms are in the top 5 percent (for overall)
top_5pc_mech <- top5_pc %>%
  select(switches) %>%
  count(switches) %>%
  separate(switches, into = paste0("sw", 1:11)) %>%
  gather(delete, mech, -n) %>%
  select(-delete) %>%
  na.omit() %>%
  group_by(mech) %>%
  summarise(n = sum(n))

save(top_5pc_mech, file = paste0("top5pc_mech_all.rda"))

### no monoculture
top_5pc_nomono_mech <- top5_pc_nomono %>%
  select(switches) %>%
  count(switches) %>%
  separate(switches, into = paste0("sw", 1:11)) %>%
  gather(delete, mech, -n) %>%
  select(-delete) %>%
  na.omit() %>%
  group_by(mech) %>%
  summarise(n = sum(n))

save(top_5pc_nomono_mech, file = paste0("top5pc_mech_nomono.rda"))


### create combined files document

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
 
  save(files, file = paste0("files_", i, ".rda"))
   
} 


load("files_0.rda")
F0 <- files

load("files_1.rda")
F1 <- files

load("files_2.rda")
F2 <- files

load("files_3.rda")
F3 <- files

load("files_4.rda")
F4 <- files

load("files_5.rda")
F5 <- files

load("files_6.rda")
F6 <- files

load("files_7.rda")
F7 <- files

load("files_8.rda")
F8 <- files

load("files_9.rda")
F9 <- files

load("files_10.rda")
F10 <- files

load("files_11.rda")
F11 <- files

files <- rbind(F0, F1, F2, F3, F4, F5,
                        F6, F7, F8, F9, F10, F11)

remove(F0, F1, F2, F3, F4, F5,
       F6, F7, F8, F9, F10, F11)

save(files, file = "files_combined.rda")
