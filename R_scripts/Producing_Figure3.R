###### Figure 3: Catford et al., Mechanistic model

rm(list=ls())

setwd("~/Predicting_e026_coexistence")

library(tidyverse)
library(rcartocolor)
library(ggplot2)

se <- function(x, num_plots) sd(x, na.rm = TRUE) / sqrt(num_plots)

by_switch_number <- tibble(mod_name = numeric(),
                           mean_rmse = numeric(),
                           se_rmse = numeric(),
                           upr = numeric(),
                           lwr = numeric())

### i = number of switches
for(i in 0:11) { ## this just uses replicates where species are planted
  
  load(paste0("Data/non_zero_replicates_", i, "_year6.rda"))
  
  non_zero2 <- non_zero_replicates %>% filter(sp_rich > 1) 
  
  num_plots <- nrow(distinct(non_zero2, subplot, plot))
  
  df <- non_zero2 %>%
    group_by(species, sp_rich, plot, subplot, switches) %>%
    summarise(rmse = mean(rmse)) %>%
    ungroup() %>%
    group_by(sp_rich, plot, subplot, switches) %>%
    summarise(rmse = mean(rmse)) %>%
    ungroup() %>%
    summarise(mean_rmse = mean(rmse, na.rm = TRUE),
              se_rmse = se(rmse, num_plots)) %>%
    mutate(upr = mean_rmse + se_rmse,
           lwr = mean_rmse - se_rmse,
           mod_name = i)
  
  by_switch_number <- bind_rows(by_switch_number, df)
  
}



(Panel3A <- ggplot(by_switch_number, aes(x = mod_name, y = mean_rmse, ymin = lwr, ymax = upr)) +
    geom_ribbon(fill = "grey70", alpha = 0.5) +
    geom_line(linewidth = 1) +
    ylab("RMSE (±1 SE)") + 
    theme_bw() + 
    xlab("Number of attributes") + # 500 x 350
    theme(legend.title = element_text(size=13),
          legend.text=element_text(size=13),
          axis.text = element_text(size = 22),
          axis.title = element_text(size=22))) 



species_by_switch_number <- tibble(species = character(),
                                   mod_name = numeric(),
                                   mean_rmse = numeric(),
                                   se_rmse = numeric(),
                                   upr = numeric(),
                                   lwr = numeric())

for(i in 0:11) { ## this just uses replicates where species are planted
  
  load(paste0("Data/non_zero_replicates_", i, "_year6.rda"))
  
  non_zero2 <- non_zero_replicates %>% filter(sp_rich > 1) 
  
  num_plots <- nrow(distinct( non_zero2, subplot, plot))
  
  df <-  non_zero2 %>%
    filter(sp_rich > 1) %>%
    group_by(species, sp_rich, plot, subplot, switches) %>%
    summarise(rmse = mean(rmse)) %>%
    ungroup() %>%
    group_by(species) %>%
    summarise(mean_rmse = mean(rmse, na.rm = TRUE),
              ncov = n(),
              se_rmse = se(rmse, ncov)) %>%
    mutate(upr = mean_rmse + se_rmse,
           lwr = mean_rmse - se_rmse,
           mod_name = i)
  
  species_by_switch_number <- bind_rows(species_by_switch_number, df)
  
}

species_by_switch_number$species <- gsub("Agropyron repens", "Elymus repens", species_by_switch_number$species)

(Panel3B <- ggplot(species_by_switch_number, aes(x = mod_name, y = mean_rmse, ymin = lwr,
                                                 ymax = upr, group = species)) +
    geom_ribbon(alpha = 0.3,  aes(fill = species)) +
    geom_line(linewidth = 1, aes(colour = species)) +
    scale_fill_manual(values = c("Agrostis scabra" = "#FEC47D" , "Andropogon gerardi" = "#65A688",
                                 "Elymus repens" = "#B54445", "Poa pratensis" = "#3192C4", 
                                 "Schizachyrium scoparium" = "#633184")) +
    scale_colour_manual(values = c("Agrostis scabra" = "#FEC47D" , "Andropogon gerardi" = "#65A688",
                                   "Elymus repens" = "#B54445", "Poa pratensis" = "#3192C4", 
                                 "Schizachyrium scoparium" = "#633184")) +
    ylab("RMSE (±1 SE)") + 
    labs(fill = "Species", colour = "Species") +
    xlab("Number of attributes") + # 500 x 350
    theme_bw() +
    theme(legend.position = c(0.77, 0.8),
          legend.key.size = unit(0.7,"line"),
          legend.title = element_text(size=13),
          legend.text=element_text(size=13),
          axis.text = element_text(size = 22),
          axis.title = element_text(size=22)))


richness_by_switch_number <- tibble(sp_rich = character(),
                                    mod_name = numeric(),
                                    mean_rmse = numeric(),
                                    se_rmse = numeric(),
                                    upr = numeric(),
                                    lwr = numeric())

for(i in 0:11) { ## this just uses replicates where species are planted
  
  load(paste0("Data/non_zero_replicates_", i, "_year6.rda"))
  
  non_zero2 <- non_zero_replicates
  
  num_plots <- nrow(distinct(non_zero2, subplot, plot))
  
  df <- non_zero2 %>%
    group_by(species, sp_rich, plot, subplot, switches) %>%
    summarise(rmse = mean(rmse)) %>%
    ungroup() %>%
    group_by(sp_rich, plot, subplot, switches) %>%
    summarise(rmse = mean(rmse)) %>%
    mutate(sp_rich = as.character(sp_rich)) %>%
    group_by(sp_rich) %>%
    summarise(mean_rmse = mean(rmse, na.rm = TRUE),
              ncov = n(),
              se_rmse = se(rmse, ncov)) %>%
    mutate(upr = mean_rmse + se_rmse,
           lwr = mean_rmse - se_rmse,
           mod_name = i)
  
  richness_by_switch_number <- bind_rows(richness_by_switch_number, df)
  
}

(Panel3C <- ggplot(richness_by_switch_number, aes(x = mod_name, y = mean_rmse, ymin = lwr,
                                                  ymax = upr, group = sp_rich)) +
    geom_ribbon(alpha = 0.3, aes(fill = sp_rich)) +
    geom_line(linewidth = 1, aes(colour = sp_rich)) +
    ylab("RMSE (±1 SE)") + 
    scale_fill_manual(values = c("1" = "#FEC47D", "2" = "#65A688",
                                 "3" = "#B54445", "5" = "#3192C4")) +
    scale_colour_manual(values = c("1" = "#FEC47D", "2" = "#65A688",
                                   "3" = "#B54445", "5" = "#3192C4")) +
    xlab("Number of attributes") + # 500 x 350
    labs(fill = "Richness", colour = "Richness") +
    theme_bw() +
    theme(legend.position = c(0.85, 0.8),
          legend.key.size = unit(0.7,"line"),
          legend.title = element_text(size=13),
          legend.text=element_text(size=13),
          axis.text = element_text(size = 22),
          axis.title = element_text(size=22)))

#### Read in all files together
data_list <- list()

for (i in 0:11) {
  file_name <- paste0("Data/non_zero_replicates_", i, "_year6.rda")
  load(file_name)
  data_list[[i + 1]] <- non_zero_replicates  # Make sure `your_object` is the correct variable in each .rda file
}

non_zero_replicates <- do.call(rbind, data_list)

### switches column refers to which switches have been turned OFF - which in turn switches that component OFF
## switches are ordered alphabetically when they occur, with full 11 switches on coded as: binit_bstar_dispersal_height_lot_mor_pheno_rep_rgr_root_rstar

load("Data/switches.rda")
switches <- df

##### Creating a switch dataframe that can determine whether a group is included in a model or not
### if a group is included, that means the mechanism, or some of the mechanisms in that group, are functioning

nonedf <- switches %>% 
  filter(str_detect(switch, "none")) %>% select(switch) ### all mechanisms are switched off
nonedf$none = 0

binitdf <- switches %>% 
  filter(str_detect(switch, "binit"))%>% select(switch)
binitdf$binit = 1

bstardf <- switches %>% 
  filter(str_detect(switch, "bstar"))%>% select(switch)
bstardf$bstar = 1

dispersaldf <- switches %>% 
  filter(str_detect(switch, "dispersal"))%>% select(switch)  
dispersaldf$dispersal = 1

heightdf <- switches %>% 
  filter(str_detect(switch, "height")) %>% select(switch)
heightdf$height = 1

lotdf <- switches %>% 
  filter(str_detect(switch, "lot")) %>% select(switch)
lotdf$lot = 1

mordf <- switches %>% 
  filter(str_detect(switch, "mor")) %>% select(switch)
mordf$mor = 1

repdf <- switches %>% 
  filter(str_detect(switch, "rep")) %>% select(switch)
repdf$rep = 1

rgrdf <- switches %>% 
  filter(str_detect(switch, "rgr"))%>% select(switch)
rgrdf$rgr = 1

rootdf <- switches %>% 
  filter(str_detect(switch, "root"))%>% select(switch)
rootdf$root = 1

rstardf <- switches %>% 
  filter(str_detect(switch, "rstar")) %>% select(switch)
rstardf$rstar = 1

phenodf <- switches %>% 
  filter(str_detect(switch, "pheno")) %>% select(switch)
phenodf$pheno = 1

switches2 <- switches %>%
  full_join(nonedf,  by="switch") %>%
  full_join(binitdf, by="switch") %>%
  full_join(bstardf, by="switch") %>%
  full_join(dispersaldf, by="switch") %>%
  full_join(heightdf, by="switch") %>%
  full_join(lotdf, by="switch") %>%
  full_join(mordf, by="switch") %>%
  full_join(repdf, by="switch") %>%
  full_join(rgrdf, by="switch") %>%
  full_join(rootdf, by="switch") %>%
  full_join(rstardf, by="switch") %>%
  full_join(phenodf, by="switch") %>%
  replace(is.na(.), 0) %>%
  mutate(ngrowth = rgr + mor) %>%
  mutate(noverlap = pheno + root + height) %>%
  mutate(ncomp = bstar + rstar) %>%
  mutate(ncolo = rep + binit + lot + dispersal)

#### If one is on then a group is on

grn1 <- switches2 %>% filter(ngrowth > 0) %>% mutate(growth = 1) %>% select(switch, growth) %>% mutate(N1 = "Growth")
ovr1 <- switches2 %>% filter(noverlap > 0) %>% mutate(overlap = 1) %>% select(switch, overlap) %>% mutate(N2 = "Overlap")
com1 <- switches2 %>% filter(ncomp > 0) %>% mutate(comp = 1) %>% select(switch, comp) %>% mutate(N3 = "Competition")
col1 <- switches2 %>% filter(ncolo > 0) %>% mutate(colo = 1) %>% select(switch, colo) %>% mutate(N4 = "Colonisation")
non1 <- switches2 %>% filter(ngrowth == 0 & noverlap == 0 & ncomp == 0 & ncolo ==0) %>% select(switch) %>% mutate(N5 = "None")


singleon <- grn1 %>% full_join(ovr1, by = "switch") %>% 
  full_join(com1, by = "switch") %>% 
  full_join(col1, by = "switch") %>% 
  full_join(non1, by = "switch") %>%
  unite(GroupNm, c("N3", "N4", "N1", "N2"), remove = TRUE) %>%
  mutate(GroupNm = str_remove_all(GroupNm, "_NA")) %>%
  mutate(GroupNm = str_remove_all(GroupNm, "NA_")) %>%
  mutate(GroupNm = replace(GroupNm, GroupNm == "NA", "None")) %>%
  select(-N5) %>%
  replace(is.na(.), 0) %>%
  mutate(group_if1 = growth + overlap + comp + colo) %>%
  select(switch, GroupNm1 = GroupNm, group_if1)

#### If two groups being on is needed

grn2 <- switches2 %>% filter(ngrowth > 1) %>% mutate(growth = 1) %>% select(switch, growth) %>% mutate(N1 = "Growth")
ovr2 <- switches2 %>% filter(noverlap > 1) %>% mutate(overlap = 1) %>% select(switch, overlap) %>% mutate(N2 = "Overlap")
com2 <- switches2 %>% filter(ncomp > 1) %>% mutate(comp = 1) %>% select(switch, comp) %>% mutate(N3 = "Competition")
col2 <- switches2 %>% filter(ncolo > 1) %>% mutate(colo = 1) %>% select(switch, colo) %>% mutate(N4 = "Colonisation")
non2 <- switches2 %>% filter(ngrowth < 2 & noverlap < 2 & ncomp < 2 & ncolo < 2) %>% select(switch) %>% mutate(N5 = "None")

twoon <- grn2 %>% full_join(ovr2, by = "switch") %>% 
  full_join(com2, by = "switch") %>% 
  full_join(col2, by = "switch") %>%  
  full_join(non2, by = "switch") %>%
  unite(GroupNm, c("N3", "N4", "N1", "N2"), remove = TRUE) %>%
  mutate(GroupNm = str_remove_all(GroupNm, "_NA")) %>%
  mutate(GroupNm = str_remove_all(GroupNm, "NA_")) %>%
  mutate(GroupNm = replace(GroupNm, GroupNm == "NA", "None")) %>%
  select(-N5) %>% replace(is.na(.), 0) %>% 
  mutate(group_if2 = growth + overlap + comp + colo) %>%
  select(switch, GroupNm2 =GroupNm, group_if2)

#### If three groups are on (or two when only 2 available)

grn3 <- switches2 %>% filter(ngrowth > 1) %>% mutate(growth = 1) %>% select(switch, growth) %>% mutate(N1 = "Growth")
ovr3 <- switches2 %>% filter(noverlap > 2) %>% mutate(overlap = 1) %>% select(switch, overlap) %>% mutate(N2 = "Overlap")
com3 <- switches2 %>% filter(ncomp > 1) %>% mutate(comp = 1) %>% select(switch, comp)%>% mutate(N3 = "Competition")
col3 <- switches2 %>% filter(ncolo > 2) %>% mutate(colo = 1) %>% select(switch, colo)%>% mutate(N4 = "Colonisation")
non3 <- switches2 %>% filter(ngrowth < 2 & noverlap < 3 & ncomp < 2 & ncolo < 3) %>% select(switch) %>% mutate(N5 = "None")

threeon <- grn3 %>% full_join(ovr3, by = "switch") %>% 
  full_join(com3, by = "switch") %>% 
  full_join(col3, by = "switch") %>% 
  full_join(non3, by = "switch") %>%
  unite(GroupNm, c("N3", "N4", "N1", "N2"), remove = TRUE) %>%
  mutate(GroupNm = str_remove_all(GroupNm, "_NA")) %>%
  mutate(GroupNm = str_remove_all(GroupNm, "NA_")) %>%
  mutate(GroupNm = replace(GroupNm, GroupNm == "NA", "None")) %>%
  select(-N5) %>% replace(is.na(.), 0) %>% 
  replace(is.na(.), 1) %>%
  mutate(group_if3 = growth + overlap + comp + colo) %>%
  select(switch, GroupNm3 =GroupNm, group_if3)

#### If four groups are on (or two when only 2 available)

grn4 <- switches2 %>% filter(ngrowth > 1) %>% mutate(growth = 1) %>% select(switch, growth)%>% mutate(N1 = "Growth")
ovr4 <- switches2 %>% filter(noverlap > 2) %>% mutate(overlap = 1) %>% select(switch, overlap)%>% mutate(N2 = "Overlap")
com4 <- switches2 %>% filter(ncomp > 1) %>% mutate(comp = 1) %>% select(switch, comp)%>% mutate(N3 = "Competition")
col4 <- switches2 %>% filter(ncolo > 3) %>% mutate(colo = 1) %>% select(switch, colo)%>% mutate(N4 = "Colonisation")
non4 <- switches2 %>% filter(ngrowth < 2 & noverlap < 3 & ncomp < 2 & ncolo < 4) %>% select(switch)%>% mutate(N5 = "None")

fouron <- grn4 %>% full_join(ovr3, by = "switch") %>% 
  full_join(com4, by = "switch") %>% 
  full_join(col4, by = "switch") %>% 
  full_join(non4, by = "switch") %>% 
  unite(GroupNm, c("N3", "N4", "N1", "N2"), remove = TRUE) %>%
  mutate(GroupNm = str_remove_all(GroupNm, "_NA")) %>%
  mutate(GroupNm = str_remove_all(GroupNm, "NA_")) %>%
  mutate(GroupNm = replace(GroupNm, GroupNm == "NA", "None")) %>%
  select(-N5) %>% replace(is.na(.), 0) %>%
  mutate(group_if4 = growth + overlap + comp + colo) %>%
  select(switch, GroupNm4 =GroupNm, group_if4)

#### Combining into 1 switch df

switches_wcats <- switches %>% full_join(singleon, by="switch") %>%
  full_join(twoon, by="switch") %>%
  full_join(threeon, by="switch") %>%
  full_join(fouron, by="switch") %>%  replace(is.na(.), 1)%>%
  rename(switches = switch)

### Names of switches dont match as switches does not use alphabetical orders for names

switchesNZ <- non_zero_replicates %>% select(switches) %>% distinct() ### codes are always alphabetically arranged

switchesWC <- switches_wcats %>% select(switches) %>% distinct() ## code orders are grouped somewhat 

match(switchesNZ$switches, switchesWC$switches)

## rearrange words in switches column to alphabetical

sort_words <- function(x) {
  sorted_words <- str_split(x, "_", simplify = TRUE) %>%
    sort() %>%
    paste(collapse = "_")
  return(sorted_words)
}

switches_wcats2 <- switches_wcats %>%
  mutate(switches = sapply(switches, sort_words))

#### Getting more intricate switch details attached

switches3 <- switches2 %>%
  rename(switches = switch) %>%
  mutate(switches = sapply(switches, sort_words)) 

### attaching together

non_zero_wswitches  <- non_zero_replicates %>%
  left_join(switches_wcats2, by="switches") %>%
  left_join(switches3, by="switches")

#### Only requiring 1 variable for a group to be included

non_zero_G1 <- non_zero_wswitches %>%
  filter(ngrowth < 2 & noverlap < 2 & ncolo < 2 & ncomp < 2) %>%
  select(switches:seed_mass, GroupNm1, group_if1)

se <- function(x, num_plots) sd(x, na.rm = TRUE) / sqrt(num_plots)

by_group_numberG1 <- tibble(mod_name = numeric(),
                            mean_rmse = numeric(),
                            se_rmse = numeric(),
                            upr = numeric(),
                            lwr = numeric())

for(i in 0:4){
  
  dfb2 <- non_zero_G1 %>% filter(group_if1 == i) %>%  filter(sp_rich > 1)
  num_plotsb <- nrow(distinct(dfb2, subplot, plot))
  
  
  dfb3 <- dfb2 %>%
    group_by(species, sp_rich, plot, subplot, switches) %>%
    summarise(rmse = mean(rmse)) %>%
    ungroup()  %>%
    group_by(sp_rich, plot, subplot, switches) %>%
    summarise(rmse = mean(rmse)) %>%
    ungroup()  %>%
    summarise(mean_rmse = mean(rmse, na.rm = TRUE),
              se_rmse = se(rmse, num_plotsb)) %>%
    mutate(upr = mean_rmse + se_rmse,
           lwr = mean_rmse - se_rmse, 
           mod_name = i)
  
  by_group_numberG1 <- bind_rows(by_group_numberG1, dfb3)
}



non_zero_G2 <- non_zero_wswitches %>%
  filter(ngrowth < 3 & noverlap < 3 & ncolo < 3 & ncomp < 3 &
           ngrowth != 1 & noverlap != 1 & ncolo != 1 & ncomp != 1)

by_group_numberG2 <- tibble(mod_name = numeric(),
                            mean_rmse = numeric(),
                            se_rmse = numeric(),
                            upr = numeric(),
                            lwr = numeric())

for(i in 0:4){
  
  dfb2 <- non_zero_G2 %>% filter(group_if2 == i) %>%  filter(sp_rich > 1)
  num_plotsb <- nrow(distinct(dfb2, subplot, plot))
  
  
  dfb3 <- dfb2 %>%
    group_by(species, sp_rich, plot, subplot, switches) %>%
    summarise(rmse = mean(rmse)) %>%
    ungroup()  %>%
    group_by(sp_rich, plot, subplot, switches) %>%
    summarise(rmse = mean(rmse)) %>%
    ungroup()  %>%
    summarise(mean_rmse = mean(rmse, na.rm = TRUE),
              se_rmse = se(rmse, num_plotsb)) %>%
    mutate(upr = mean_rmse + se_rmse,
           lwr = mean_rmse - se_rmse, 
           mod_name = i)
  
  by_group_numberG2 <- bind_rows(by_group_numberG2, dfb3)
}


non_zero_G3 <- non_zero_wswitches %>%
  filter(ngrowth < 3 & noverlap < 4 & ncolo < 4 & ncomp < 3 &
           ngrowth != 1 & !noverlap %in% c(1,2) & !ncolo %in% c(1,2) & ncomp != 1)

by_group_numberG3 <- tibble(mod_name = numeric(),
                            mean_rmse = numeric(),
                            se_rmse = numeric(),
                            upr = numeric(),
                            lwr = numeric())

for(i in 0:4){
  
  dfb2 <- non_zero_G3 %>% filter(group_if3 == i) %>%  filter(sp_rich > 1)
  num_plotsb <- nrow(distinct(dfb2, subplot, plot))
  
  
  dfb3 <- dfb2 %>%
    group_by(species, sp_rich, plot, subplot, switches) %>%
    summarise(rmse = mean(rmse)) %>%
    ungroup()  %>%
    group_by(sp_rich, plot, subplot, switches) %>%
    summarise(rmse = mean(rmse)) %>%
    ungroup()  %>%
    summarise(mean_rmse = mean(rmse, na.rm = TRUE),
              se_rmse = se(rmse, num_plotsb)) %>%
    mutate(upr = mean_rmse + se_rmse,
           lwr = mean_rmse - se_rmse, 
           mod_name = i)
  
  by_group_numberG3 <- bind_rows(by_group_numberG3, dfb3)
}



non_zero_G4 <- non_zero_wswitches %>%
  filter(ngrowth != 1 & !noverlap %in% c(1,2) & !ncolo %in% c(1,2,3) & ncomp != 1)



by_group_numberG4 <- tibble(mod_name = numeric(),
                            mean_rmse = numeric(),
                            se_rmse = numeric(),
                            upr = numeric(),
                            lwr = numeric())

for(i in 0:4){
  
  dfb2 <- non_zero_G4 %>% filter(group_if4 == i) %>%  filter(sp_rich > 1)
  num_plotsb <- nrow(distinct(dfb2, subplot, plot))
  
  dfb3 <- dfb2 %>%
    group_by(species, sp_rich, plot, subplot, switches) %>%
    summarise(rmse = mean(rmse)) %>%
    ungroup()  %>%
    group_by(sp_rich, plot, subplot, switches) %>%
    summarise(rmse = mean(rmse)) %>%
    ungroup()  %>%
    summarise(mean_rmse = mean(rmse, na.rm = TRUE),
              se_rmse = se(rmse, num_plotsb)) %>%
    mutate(upr = mean_rmse + se_rmse,
           lwr = mean_rmse - se_rmse, 
           mod_name = i)
  
  by_group_numberG4 <- bind_rows(by_group_numberG4, dfb3)
}


comb_avgs <- bind_rows((by_group_numberG1 %>% mutate(factors = "1")),
                       (by_group_numberG2 %>% mutate(factors = "2")),
                       (by_group_numberG3 %>% mutate(factors = "3")),
                       (by_group_numberG4 %>% mutate(factors = "4")))


Panel3d <- ggplot(comb_avgs, aes(x = mod_name, y = mean_rmse, ymin = lwr, ymax = upr, fill = factors)) +
  geom_ribbon(alpha = 0.3, outline.type = "both") +
  geom_line(linewidth = 1, aes(colour = factors)) +
  scale_fill_manual(values = c("1" = "#FEC47D", "2" = "#65A688",
                               "3" = "#B54445", "4" = "#3192C4")) +
  scale_colour_manual(values = c("1" = "#FEC47D", "2" = "#65A688",
                                 "3" = "#B54445", "4" = "#3192C4")) +
  ylab("RMSE (±1 SE)") +
  xlab("Number of mechanisms") + # 500 x 350
  theme_bw() +
  ylim(0.65, 1.15) +
  guides(fill=guide_legend("Attributes required", order = 1),
         colour=guide_legend("Attributes required", order = 1))  +
  theme(legend.position = c(0.8, 0.8),
        legend.key.size = unit(0.7,"line"),
        legend.title = element_text(size=13),
        legend.text=element_text(size=13),
        axis.text = element_text(size = 22),
        axis.title = element_text(size=22))



#### group plots
species_by_group_numberG2 <- tibble(species = character(),
                                    mod_name = numeric(),
                                    mean_rmse = numeric(),
                                    se_rmse = numeric(),
                                    upr = numeric(),
                                    lwr = numeric())

for(i in 0:4) { ## this just uses replicates where species are planted
  
  df2 <- non_zero_G2 %>% filter(group_if1 == i) %>%  filter(sp_rich > 1)
  num_plots <- nrow(distinct(df2, subplot, plot))
  
  
  df3 <- df2 %>%
    group_by(species, sp_rich, plot, subplot, switches) %>%
    summarise(rmse = mean(rmse)) %>%
    ungroup()  %>%
    group_by(species) %>%
    summarise(mean_rmse = mean(rmse, na.rm = TRUE),
              se_rmse = se(rmse, num_plots)) %>%
    mutate(upr = mean_rmse + se_rmse,
           lwr = mean_rmse - se_rmse, 
           mod_name = i)
  
  species_by_group_numberG2 <- bind_rows(species_by_group_numberG2, df3)
  
}

species_by_group_numberG2 <- species_by_group_numberG2 %>%
  mutate(species = replace(species, species == "Agropyron repens", "Elymus repens"))

(Panel3e <- ggplot(species_by_group_numberG2, aes(x = mod_name, y = mean_rmse, ymin = lwr,
                                                   ymax = upr, group = species)) +
    geom_ribbon(aes(fill = species), alpha = 0.3) +
    geom_line(linewidth = 1, aes(colour = species)) +
    scale_fill_manual(values = c("Agrostis scabra" = "#FEC47D" , "Andropogon gerardi" = "#65A688",
                                 "Elymus repens" = "#B54445", "Poa pratensis" = "#3192C4", 
                                 "Schizachyrium scoparium" = "#633184")) +
    scale_colour_manual(values = c("Agrostis scabra" = "#FEC47D" , "Andropogon gerardi" = "#65A688",
                                   "Elymus repens" = "#B54445", "Poa pratensis" = "#3192C4", 
                                   "Schizachyrium scoparium" = "#633184")) +
    ylab("RMSE (±1 SE)") +
    labs(fill = "Species", colour = "Species") +
    xlab("Number of mechanisms") + # 500 x 350
    theme_bw() +
    ylim(0.25, 1.65) +
    theme(legend.position = c(0.77, 0.8),
          legend.key.size = unit(0.7,"line"),
          legend.title = element_text(size=13),
          legend.text=element_text(size=13),
          axis.text = element_text(size = 22),
          axis.title = element_text(size=22))) ## weird response by Poa - others seem consistent with above

richness_by_group_numberG2 <- tibble(sp_rich = character(),
                                     mod_name = numeric(),
                                     mean_rmse = numeric(),
                                     se_rmse = numeric(),
                                     upr = numeric(),
                                     lwr = numeric())

for(i in 0:4) { 
  
  df2 <- non_zero_G2 %>% filter(group_if1 == i)
  num_plots <- nrow(distinct(df2, subplot, plot))
  
  
  df3 <- df2 %>%
    group_by(species, sp_rich, plot, subplot, switches) %>%
    summarise(rmse = mean(rmse)) %>%
    ungroup()  %>%
    group_by(sp_rich, plot, subplot, switches) %>%
    summarise(rmse = mean(rmse)) %>%
    mutate(sp_rich = as.character(sp_rich)) %>%
    group_by(sp_rich) %>%
    summarise(mean_rmse = mean(rmse, na.rm = TRUE),
              se_rmse = se(rmse, num_plots)) %>%
    mutate(upr = mean_rmse + se_rmse,
           lwr = mean_rmse - se_rmse, 
           mod_name = i)
  
  richness_by_group_numberG2 <- bind_rows(richness_by_group_numberG2, df3)
  
}

(Panel3f <- ggplot(richness_by_group_numberG2, aes(x = mod_name, y = mean_rmse, ymin = lwr,
                                                    ymax = upr, group = sp_rich)) +
    geom_ribbon(aes(fill = sp_rich), alpha = 0.3) +
    geom_line(linewidth = 1, aes(colour = sp_rich)) +
    scale_fill_manual(values = c("1" = "#FEC47D", "2" = "#65A688",
                                 "3" = "#B54445", "5" = "#3192C4")) +
    scale_colour_manual(values = c("1" = "#FEC47D", "2" = "#65A688",
                                   "3" = "#B54445", "5" = "#3192C4")) +
    ylab("RMSE (±1 SE)") +
    labs(fill = "Richness", colour = "Richness") +
    xlab("Number of mechanisms") + # 500 x 350
    theme_bw() +
    theme(legend.position = c(0.85, 0.8),
          legend.key.size = unit(0.7,"line"),
          legend.title = element_text(size=13),
          legend.text=element_text(size=13),
          axis.text = element_text(size = 22),
          axis.title = element_text(size=22)))


library(ggpubr)
ggarrange(Panel3A, Panel3B, Panel3C, Panel3d, Panel3e, Panel3f,
          labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 3, nrow = 2, align = "hv",
          font.label = list(size = 22)) ## 1800 x 800


##### Alternative version - shortening species names - matching terms in manuscript

species_by_group_numberG2b <- species_by_group_numberG2 %>%
  mutate(species = replace(species, species == "Elymus repens", "Elymus"))%>%
  mutate(species = replace(species, species == "Andropogon gerardi", "Andropogon"))%>%
  mutate(species = replace(species, species == "Agrostis scabra", "Agrostis"))%>%
  mutate(species = replace(species, species == "Poa pratensis", "Poa"))%>%
  mutate(species = replace(species, species == "Schizachyrium scoparium", "Schizachyrium"))

(Panel3eb <- ggplot(species_by_group_numberG2b, aes(x = mod_name, y = mean_rmse, ymin = lwr,
                                                  ymax = upr, group = species)) +
    geom_ribbon(aes(fill = species), alpha = 0.3) +
    geom_line(linewidth = 1, aes(colour = species)) +
    scale_fill_manual(values = c("Agrostis" = "#FEC47D" , "Andropogon" = "#65A688",
                                 "Elymus" = "#B54445", "Poa" = "#3192C4", 
                                 "Schizachyrium" = "#633184")) +
    scale_colour_manual(values = c("Agrostis" = "#FEC47D" , "Andropogon" = "#65A688",
                                   "Elymus" = "#B54445", "Poa" = "#3192C4", 
                                   "Schizachyrium" = "#633184")) +
    ylab("RMSE (±1 SE)") +
    labs(fill = "Species", colour = "Species") +
    xlab("Number of mechanisms") + # 500 x 350
    theme_bw() +
    ylim(0.25, 1.65) +
    theme(legend.position = c(0.82, 0.77),
          legend.key.size = unit(0.7,"line"),
          legend.title = element_text(size=15),
          legend.text=element_text(size=15),
          axis.text = element_text(size = 22),
          axis.title = element_text(size=22))) 

species_by_switch_numberb <- species_by_switch_number  %>%
  mutate(species = replace(species, species == "Elymus repens", "Elymus"))%>%
  mutate(species = replace(species, species == "Andropogon gerardi", "Andropogon"))%>%
  mutate(species = replace(species, species == "Agrostis scabra", "Agrostis"))%>%
  mutate(species = replace(species, species == "Poa pratensis", "Poa"))%>%
  mutate(species = replace(species, species == "Schizachyrium scoparium", "Schizachyrium"))

(Panel3Bb <- ggplot(species_by_switch_numberb, aes(x = mod_name, y = mean_rmse, ymin = lwr,
                                                 ymax = upr, group = species)) +
    geom_ribbon(alpha = 0.3,  aes(fill = species)) +
    geom_line(linewidth = 1, aes(colour = species)) +
    scale_fill_manual(values = c("Agrostis" = "#FEC47D" , "Andropogon" = "#65A688",
                                 "Elymus" = "#B54445", "Poa" = "#3192C4", 
                                 "Schizachyrium" = "#633184")) +
    scale_colour_manual(values = c("Agrostis" = "#FEC47D" , "Andropogon" = "#65A688",
                                   "Elymus" = "#B54445", "Poa" = "#3192C4", 
                                   "Schizachyrium" = "#633184")) +
    ylab("RMSE (±1 SE)") + 
    labs(fill = "Species", colour = "Species") +
    xlab("Number of attributes") + # 500 x 350
    theme_bw() +
    theme(legend.position = c(0.82, 0.77),
          legend.key.size = unit(0.7,"line"),
          legend.title = element_text(size=15),
          legend.text=element_text(size=15),
          axis.text = element_text(size = 22),
          axis.title = element_text(size=22)))

(Panel3Cb <- ggplot(richness_by_switch_number, aes(x = mod_name, y = mean_rmse, ymin = lwr,
                                                  ymax = upr, group = sp_rich)) +
    geom_ribbon(alpha = 0.3, aes(fill = sp_rich)) +
    geom_line(linewidth = 1, aes(colour = sp_rich)) +
    ylab("RMSE (±1 SE)") + 
    scale_fill_manual(values = c("1" = "#FEC47D", "2" = "#65A688",
                                 "3" = "#B54445", "5" = "#3192C4")) +
    scale_colour_manual(values = c("1" = "#FEC47D", "2" = "#65A688",
                                   "3" = "#B54445", "5" = "#3192C4")) +
    xlab("Number of attributes") + # 500 x 350
    labs(fill = "Richness", colour = "Richness") +
    theme_bw() +
    theme(legend.position = c(0.85, 0.8),
          legend.key.size = unit(0.7,"line"),
          legend.title = element_text(size=15),
          legend.text=element_text(size=15),
          axis.text = element_text(size = 22),
          axis.title = element_text(size=22)))



(Panel3fb <- ggplot(richness_by_group_numberG2, aes(x = mod_name, y = mean_rmse, ymin = lwr,
                                                   ymax = upr, group = sp_rich)) +
    geom_ribbon(aes(fill = sp_rich), alpha = 0.3) +
    geom_line(linewidth = 1, aes(colour = sp_rich)) +
    scale_fill_manual(values = c("1" = "#FEC47D", "2" = "#65A688",
                                 "3" = "#B54445", "5" = "#3192C4")) +
    scale_colour_manual(values = c("1" = "#FEC47D", "2" = "#65A688",
                                   "3" = "#B54445", "5" = "#3192C4")) +
    ylab("RMSE (±1 SE)") +
    labs(fill = "Richness", colour = "Richness") +
    xlab("Number of mechanisms") + # 500 x 350
    theme_bw() +
    theme(legend.position = c(0.85, 0.8),
          legend.key.size = unit(0.7,"line"),
          legend.title = element_text(size=15),
          legend.text=element_text(size=15),
          axis.text = element_text(size = 22),
          axis.title = element_text(size=22)))

Panel3db <- ggplot(comb_avgs, aes(x = mod_name, y = mean_rmse, ymin = lwr, ymax = upr, fill = factors)) +
  geom_ribbon(alpha = 0.3, outline.type = "both") +
  geom_line(linewidth = 1, aes(colour = factors)) +
  scale_fill_manual(values = c("1" = "#FEC47D", "2" = "#65A688",
                               "3" = "#B54445", "4" = "#3192C4")) +
  scale_colour_manual(values = c("1" = "#FEC47D", "2" = "#65A688",
                                 "3" = "#B54445", "4" = "#3192C4")) +
  ylab("RMSE (±1 SE)") +
  xlab("Number of mechanisms") + # 500 x 350
  theme_bw() +
  ylim(0.65, 1.2) +
  guides(fill=guide_legend("Attributes required", order = 1),
         colour=guide_legend("Attributes required", order = 1))  +
  theme(legend.position = c(0.8, 0.8),
        legend.key.size = unit(0.7,"line"),
        legend.title = element_text(size=15),
        legend.text=element_text(size=15),
        axis.text = element_text(size = 22),
        axis.title = element_text(size=22))


library(ggpubr)
ggarrange(Panel3A, Panel3Bb, Panel3Cb, Panel3db, Panel3eb, Panel3fb,
          labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 3, nrow = 2, align = "hv",
          font.label = list(size = 22)) ## 1700 x 900 or 18.21 x 9.75


