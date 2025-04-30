###### Producing supplementary figures: Catford et al., Mechanistic model paper


rm(list=ls())

library(tidyverse)
library(rcartocolor)


##### Table S4 - correlations amongst traits 

### Correlations between overlap variables

setwd("simulateCoexistence/")
load("data/trait_data.rda")
load("data/biomass_data.rda")
source("R/competition_functions.R")


sp_list <- dimnames(biomass_data)[[3]]

root_overlap <- get_overlap_spat(trait_data$rootingdepth, switch_off_overlap = FALSE)
dimnames(root_overlap) <- list(sp_list, sp_list)
height_overlap <- get_overlap_spat(trait_data$h, switch_off_overlap = FALSE)
dimnames(height_overlap) <- list(sp_list, sp_list)
pheno_overlap <- get_overlap_pheno(trait_data[,c("PhenStart", "PhenEnd")], switch_off_overlap = FALSE)
dimnames(pheno_overlap) <- list(sp_list, sp_list)

library(reshape2)
root_df <- melt(root_overlap, varnames = c("Species1", "Species2")) %>% dplyr::rename("Root" = value)
height_df <- melt(height_overlap, varnames = c("Species1", "Species2")) %>% dplyr::rename("Height" = value)
pheno_df <- melt(pheno_overlap, varnames = c("Species1", "Species2")) %>% dplyr::rename("Temporal" = value)


comb <- cbind(root_df, height_df$Height, pheno_df$Temporal) %>% filter(!Species1 == Species2)
colnames(comb) <- c("Sp1", "Sp2", "Root", "Height", "Pheno")
cor(comb[,3:5])

### averages for each species
overlap_avgs <- comb %>% group_by(Sp1) %>%
  summarise(Root_mean = mean(Root, na.rm = TRUE),
            Height_mean = mean(Height, na.rm = TRUE),
            Pheno_mean = mean(Pheno, na.rm = TRUE)) %>% rename(species = Sp1)


init_df <- read.csv("e026_seed_weights.csv") %>% select(plot = Plot.number, subplot = Subplot.number, 7:11)
colnames(init_df) <- c("plot", "subplot", "Agropyron repens", "Agrostis scabra", "Andropogon gerardi", "Schizachyrium scoparium", "Poa pratensis")

init_long <- init_df %>% tidyr::gather("species", "binit", 3:7) %>% filter(binit > 0 & binit != "NA") %>% select(species, plot, subplot, binit)

binit_avg <- read.csv("E26_modelled_results_year6.csv") %>%
  filter(sown == 1) %>% select(species, plot, subplot) %>% inner_join(init_long) %>% group_by(species) %>% summarise(binit = mean(binit))


attributes <- trait_data %>% select(species, Bm, mor, R, rf, rgr) %>% left_join(overlap_avgs, by="species") %>% left_join(binit_avg, by="species") %>%
  select(bstar = Bm, rstar = R, mor, rgr, root = Root_mean, height = Height_mean, pheno = Pheno_mean, binit, rep = rf)


as.data.frame(cor(attributes)) %>% round(2)

###### Figure S1: Model error with differing number of attributes in groups
load("Non_zero_replicates_all_year6.rda")

### switches column refers to which switches have been turned OFF - which in turn switches that component OFF
## switches are ordered alphabetically when they occur, with full 11 switches on coded as: binit_bstar_dispersal_height_lot_mor_pheno_rep_rgr_root_rstar

load("switches.rda")
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


PanelS1A <- ggplot(comb_avgs[comb_avgs$factors == "1",], aes(x = mod_name, y = mean_rmse, ymin = lwr, ymax = upr)) +
  geom_ribbon(alpha = 0.3, outline.type = "both") +
  geom_line(linewidth = 1) +
  ylab("RMSE (±1 SE)") +
  xlab("Number of mechanisms") + # 500 x 350
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size=12))

#### group plots

species_by_group_numberG1 <- tibble(species = character(),
                                    mod_name = numeric(),
                                    mean_rmse = numeric(),
                                    se_rmse = numeric(),
                                    upr = numeric(),
                                    lwr = numeric())

for(i in 0:4) { ## this just uses replicates where species are planted

  df2 <- non_zero_G1 %>% filter(group_if1 == i) %>%  filter(sp_rich > 1)
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

  species_by_group_numberG1 <- bind_rows(species_by_group_numberG1, df3)

}

species_by_group_numberG1 <- species_by_group_numberG1 %>%
  mutate(species = replace(species, species == "Agropyron repens", "Elymus repens"))


(PanelS1B <- ggplot(species_by_group_numberG1, aes(x = mod_name, y = mean_rmse, ymin = lwr,
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
    ylim(0.25, 1.85) +
    theme(legend.position = c(0.7, 0.8),
          legend.key.size = unit(0.7,"line"),
          axis.text = element_text(size = 10),
          axis.title = element_text(size=12))) ## weird response by Poa - others seem consistent with above

richness_by_group_numberG1 <- tibble(sp_rich = character(),
                                     mod_name = numeric(),
                                     mean_rmse = numeric(),
                                     se_rmse = numeric(),
                                     upr = numeric(),
                                     lwr = numeric())

for(i in 0:4) {

  df2 <- non_zero_G1 %>% filter(group_if1 == i)
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

  richness_by_group_numberG1 <- bind_rows(richness_by_group_numberG1, df3)

}

(PanelS1C <- ggplot(richness_by_group_numberG1, aes(x = mod_name, y = mean_rmse, ymin = lwr,
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
          axis.text = element_text(size = 10),
          axis.title = element_text(size=12)))


### Three attributes

PanelS1D <- ggplot(comb_avgs[comb_avgs$factors == "3",], aes(x = mod_name, y = mean_rmse, ymin = lwr, ymax = upr)) +
  geom_ribbon(alpha = 0.3, outline.type = "both") +
  geom_line(linewidth = 1) +
  ylab("RMSE (±1 SE)") +
  xlab("Number of mechanisms") + # 500 x 350
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size=12))

#### group plots

species_by_group_numberG3 <- tibble(species = character(),
                                    mod_name = numeric(),
                                    mean_rmse = numeric(),
                                    se_rmse = numeric(),
                                    upr = numeric(),
                                    lwr = numeric())

for(i in 0:4) { ## this just uses replicates where species are planted

  df2 <- non_zero_G3 %>% filter(group_if1 == i) %>%  filter(sp_rich > 1)
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

  species_by_group_numberG3 <- bind_rows(species_by_group_numberG3, df3)

}

species_by_group_numberG3 <- species_by_group_numberG3 %>%
  mutate(species = replace(species, species == "Agropyron repens", "Elymus repens"))

(PanelS1E <- ggplot(species_by_group_numberG3, aes(x = mod_name, y = mean_rmse, ymin = lwr,
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
          axis.text = element_text(size = 10),
          axis.title = element_text(size=12))) ## weird response by Poa - others seem consistent with above


richness_by_group_numberG3 <- tibble(sp_rich = character(),
                                     mod_name = numeric(),
                                     mean_rmse = numeric(),
                                     se_rmse = numeric(),
                                     upr = numeric(),
                                     lwr = numeric())

for(i in 0:4) {

  df2 <- non_zero_G3 %>% filter(group_if1 == i)
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

  richness_by_group_numberG3 <- bind_rows(richness_by_group_numberG3, df3)

}

(PanelS1F <- ggplot(richness_by_group_numberG3, aes(x = mod_name, y = mean_rmse, ymin = lwr,
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
          axis.text = element_text(size = 10),
          axis.title = element_text(size=12)))

### four attributes per plot

PanelS1G <- ggplot(comb_avgs[comb_avgs$factors == "4",], aes(x = mod_name, y = mean_rmse, ymin = lwr, ymax = upr)) +
  geom_ribbon(alpha = 0.3, outline.type = "both") +
  geom_line(linewidth = 1) +
  ylab("RMSE (±1 SE)") +
  xlab("Number of mechanisms") + # 500 x 350
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size=12))

#### group plots

species_by_group_numberG4 <- tibble(species = character(),
                                    mod_name = numeric(),
                                    mean_rmse = numeric(),
                                    se_rmse = numeric(),
                                    upr = numeric(),
                                    lwr = numeric())

for(i in 0:4) { ## this just uses replicates where species are planted

  df2 <- non_zero_G4 %>% filter(group_if1 == i) %>%  filter(sp_rich > 1)
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

  species_by_group_numberG4 <- bind_rows(species_by_group_numberG4, df3)

}

species_by_group_numberG4 <- species_by_group_numberG4 %>%
  mutate(species = replace(species, species == "Agropyron repens", "Elymus repens"))

(PanelS1H <- ggplot(species_by_group_numberG4, aes(x = mod_name, y = mean_rmse, ymin = lwr,
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
          axis.text = element_text(size = 10),
          axis.title = element_text(size=12))) ## weird response by Poa - others seem consistent with above


richness_by_group_numberG4 <- tibble(sp_rich = character(),
                                     mod_name = numeric(),
                                     mean_rmse = numeric(),
                                     se_rmse = numeric(),
                                     upr = numeric(),
                                     lwr = numeric())

for(i in 0:4) {

  df2 <- non_zero_G4 %>% filter(group_if1 == i)
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

  richness_by_group_numberG4 <- bind_rows(richness_by_group_numberG4, df3)

}

(PanelS1I <- ggplot(richness_by_group_numberG4, aes(x = mod_name, y = mean_rmse, ymin = lwr,
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
          axis.text = element_text(size = 10),
          axis.title = element_text(size=12)))



###### Figure S1
ggarrange(PanelS1A, PanelS1B, PanelS1C, PanelS1D,
          PanelS1E, PanelS1F, PanelS1G, PanelS1H, PanelS1I,
          ncol = 3, nrow = 3, align = "hv",
          labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I")) ## 1400 x 1100


##### Figure S2: Frequency top mechanism occurs in each top 5pc of best model sets
load("data/total_biomass_data.rda")
sp_names <- unlist(dimnames(biomass_data_total)[3])

seed_plots <- as.data.frame.table(seeding_data_total) %>%
  mutate(subplot = as.numeric(Var1),
         plot = as.numeric(Var3),
         species = factor(Var2, labels = sp_names)) %>%
  select(subplot, plot, species, seed_mass = Freq) %>%
  filter(seed_mass != 0) %>% spread(species, seed_mass) %>%
  unite(plot_code, c("plot", "subplot"))

load(paste0("non_zero_replicates_", 0, "_year6.rda"))

plot_codes <- non_zero_replicates %>% filter(sp_rich == 2) %>%
  ungroup() %>%
  select(plot, subplot) %>% ##105 plots
  distinct() %>% unite(plot_code, c("plot", "subplot"), remove = FALSE) %>%
  left_join(seed_plots, by="plot_code") %>%
  mutate(Agrre = if_else(`Agropyron repens` > 1, "Agrre", NA_character_)) %>%
  mutate(Agrsc = if_else(`Agrostis scabra` > 1, "Agrsc", NA_character_)) %>%
  mutate(Andge = if_else(`Andropogon gerardi` > 1, "Andge", NA_character_)) %>%
  mutate(Poapr = if_else(`Poa pratensis` > 1, "Poapr", NA_character_)) %>%
  mutate(Schsc = if_else(`Schizachyrium scoparium` > 1, "Schsc", NA_character_)) %>%
  rowwise() %>%
  mutate(code = paste(na.omit(c(Agrre, Agrsc, Andge, Poapr, Schsc)), collapse = "_")) %>%
  select(plot_code, code)

load("Non_zero_replicates_all_year6.rda")
load("files_combined.rda")

df_Split <- non_zero_replicates  %>%
  inner_join(files %>% select(full_name, mod_name), by = c("replicate" = "full_name")) %>%
  filter(sp_rich == 2) %>%
  unite(plot_code, c("plot", "subplot"), remove = FALSE) %>%
  left_join(plot_codes, by="plot_code") %>%
  mutate(mod_name = as.numeric(mod_name)) %>%
  group_by(species, sp_rich, plot, subplot, switches, mod_name, code) %>%
  summarise(rmse = mean(rmse)) %>%
  ungroup() %>%
  group_by(sp_rich, plot, subplot, switches, mod_name, code) %>%
  summarise(rmse = mean(rmse)) %>%
  ungroup()## gives dataframe of output replicates with model names


top_5pc_Split <- df_Split %>%
  group_by(mod_name, code) %>%
  mutate(pc95 = quantile(rmse, 0.05)) %>%
  filter(rmse <= pc95) ### top 5% quantile of RMSE: i.e. best performing models, from each amount of switches

top_5pc_sw_Split <- top_5pc_Split %>%
  filter(mod_name == 1) %>% ##filters out single switch models
  group_by(code) %>%
  count(switches) %>%
  mutate(switches = as.factor(switches),
         switches = fct_reorder(switches, n, .desc = TRUE)) %>%
  arrange(-n) %>%
  slice(1:11)

top5pc_mech_Split <- top_5pc_Split %>%
  select(switches, mod_name, code) %>%
  filter(mod_name %in% 1:11) %>% ## gets rid of all switches on (full model)
  separate(switches, into = paste0("sw", 1:11)) %>% ## separates out models
  gather(delete, mech, -mod_name, -code) %>%
  select(-delete) %>%
  na.omit()

top_5pc_sw_Split_switches <- top_5pc_sw_Split %>% group_by(code) %>%
  mutate(switch_total= sum(n)) %>% ungroup() %>%
  mutate(perc = (n/switch_total)*100)

dfZG2_code <- non_zero_G2 %>%
  unite(plot_code, c("plot", "subplot"), remove = FALSE) %>%
  left_join(plot_codes, by="plot_code") %>%
  rename(groups_ch = GroupNm2, Groups_on = group_if2) %>%
  mutate(Groups_on = as.numeric(Groups_on)) %>%
  filter(sp_rich ==2) %>%
  group_by(species, sp_rich, plot, subplot, groups_ch, Groups_on, code) %>%
  summarise(rmse = mean(rmse)) %>%
  ungroup()  %>%
  group_by(sp_rich, plot, subplot, groups_ch, Groups_on, code) %>%
  summarise(rmse = mean(rmse)) %>%
  ungroup() ## gives dataframe of output replicates with model names

top_5pcG2_code <- dfZG2_code %>%
  group_by(Groups_on, code) %>%
  mutate(pc95 = quantile(rmse, 0.05)) %>%
  filter(rmse <= pc95) ### top 5% quantile of RMSE: i.e. best performing models, from each amount of switches

top_5pc_swG2_code <- top_5pcG2_code %>%
  filter(Groups_on == 1) %>% ##filters out single group models
  count(groups_ch) %>%
  mutate(groups_ch = as.factor(groups_ch),
         groups_ch = fct_reorder(groups_ch, n, .desc = TRUE)) %>%
  arrange(-n)

top_5pc_mechG2_code <- top_5pcG2_code %>%
  select(code, groups_ch, Groups_on) %>%
  filter(Groups_on %in% 1:4) %>% ## gets rid of all switches on (full model)
  separate(groups_ch, into = paste0("sw", 1:4)) %>% ## separates out models
  gather(delete, mech, -Groups_on, -code) %>%
  select(-delete) %>%
  na.omit()

perc_top5pc_mechs <- top_5pc_mechG2_code %>% mutate(count = 1) %>% group_by(code, Groups_on, mech) %>%
  summarise(grp_total= sum(count)) %>% ungroup() %>%
  group_by(code, Groups_on) %>% mutate(total_count = sum(grp_total)) %>%
  mutate(perc = (grp_total/total_count)*100)

top5pc_mech_Split$mech <- factor(top5pc_mech_Split$mech, levels = c("binit", "lot", "rep", "dispersal", "bstar", "rstar",
                                                                    "mor", "rgr", "root", "height", "pheno"))

top_5pc_sw_Split_switches$switches <- factor(top_5pc_sw_Split_switches$switches, levels = c("binit", "lot", "rep", "dispersal", "bstar", "rstar",
                                                                                            "mor", "rgr", "root", "height", "pheno"))

df_row1 <- bind_rows((top_5pc_swG2_code[top_5pc_swG2_code$code == "Agrre_Agrsc",] %>% mutate(Groups_on = as.character(Groups_on)) %>%
                        mutate(groups_ch = as.character(groups_ch))),
                     bind_cols(Groups_on = c("1"), code = c("Agrre_Agrsc"),
                               groups_ch = c("Growth"), n = c(0)))
df_row1$groups_ch <- factor(df_row1$groups_ch, levels = c("Colonisation", "Competition", "Growth", "Overlap"))

df_row1_4 <- bind_rows((top_5pc_sw_Split_switches[top_5pc_sw_Split_switches$code == "Agrre_Agrsc",] %>% mutate(switches = as.character(switches)) %>%
                          mutate(switches = as.character(switches))),
                       bind_cols(code = c("Agrre_Agrsc"), switches = c("bstar"),
                                 n = c(0), mech = c("bstar"), switch_total = c(10), perc = c(0)))
df_row1_4$switches <- factor(df_row1_4$switches, levels = c("binit", "lot", "rep", "dispersal", "bstar", "rstar", "mor", "rgr", "root", "height", "pheno"))


df_row2 <- bind_rows((top_5pc_swG2_code[top_5pc_swG2_code$code == "Agrsc_Andge",] %>% mutate(Groups_on = as.character(Groups_on)) %>%
                        mutate(groups_ch = as.character(groups_ch))),
                     bind_cols(Groups_on = c("1", "1", "1"), code = c("Agrsc_Andge", "Agrsc_Andge", "Agrsc_Andge"),
                               groups_ch = c("Competition", "Growth", "Overlap"), n = c(0, 0, 0)))
df_row2$groups_ch <- factor(df_row2$groups_ch, levels = c("Colonisation", "Competition", "Growth", "Overlap"))

df_row2_4 <- bind_rows((top_5pc_sw_Split_switches[top_5pc_sw_Split_switches$code == "Agrre_Poapr",] %>% mutate(switches = as.character(switches)) %>%
                          mutate(switches = as.character(switches))),
                       bind_cols(code = c("Agrre_Poapr"), switches = c("pheno", "rstar"),
                                 n = c(0, 0), mech = c("pheno", "rstar"), switch_total = c(10, 10), perc = c(0, 0)))
df_row2_4$switches <- factor(df_row2_4$switches, levels = c("binit", "lot", "rep", "dispersal", "bstar", "rstar", "mor", "rgr", "root", "height", "pheno"))

df_row3 <- bind_rows((top_5pc_swG2_code[top_5pc_swG2_code$code == "Agrre_Poapr",] %>% mutate(Groups_on = as.character(Groups_on)) %>%
                        mutate(groups_ch = as.character(groups_ch))),
                     bind_cols(Groups_on = c("1"), code = c("Agrre_Poapr"),
                               groups_ch = c("Overlap"), n = c(0)))
df_row3$groups_ch <- factor(df_row3$groups_ch, levels = c("Colonisation", "Competition", "Growth", "Overlap"))

df_row3_4 <- bind_rows((top_5pc_sw_Split_switches[top_5pc_sw_Split_switches$code == "Agrre_Schsc",] %>% mutate(switches = as.character(switches)) %>%
                          mutate(switches = as.character(switches))),
                       bind_cols(code = c("Agrre_Schsc"), switches = c("binit", "dispersal", "height", "lot", "mor",
                                                                       "pheno", "rep", "rgr", "rstar"),
                                 n = c(0,0,0,0,0,0,0,0,0), mech = c("binit", "dispersal", "height", "lot", "mor",
                                                                    "pheno", "rep", "rgr", "rstar"), switch_total = c(0,0,0,0,0,0,0,0,0), perc = c(0,0,0,0,0,0,0,0,0)))
df_row3_4$switches <- factor(df_row3_4$switches, levels = c("binit", "lot", "rep", "dispersal", "bstar", "rstar", "mor", "rgr", "root", "height", "pheno"))

df_row4 <- bind_rows((top_5pc_swG2_code[top_5pc_swG2_code$code == "Agrre_Schsc",] %>% mutate(Groups_on = as.character(Groups_on)) %>%
                        mutate(groups_ch = as.character(groups_ch))),
                     bind_cols(Groups_on = c("1", "1"), code = c("Agrre_Schsc","Agrre_Schsc"),
                               groups_ch = c("Competition", "Colonisation"), n = c(0, 0)))
df_row4$groups_ch <- factor(df_row4$groups_ch, levels = c("Colonisation", "Competition", "Growth", "Overlap"))

df_row4_4 <- bind_rows((top_5pc_sw_Split_switches[top_5pc_sw_Split_switches$code == "Agrsc_Andge",] %>% mutate(switches = as.character(switches)) %>%
                          mutate(switches = as.character(switches))),
                       bind_cols(code = c("Agrsc_Andge"), switches = c("dispersal", "lot", "mor",
                                                                       "pheno", "rep", "rgr", "rstar"),
                                 n = c(0,0,0,0,0,0,0), mech = c("dispersal", "lot", "mor","pheno", "rep", "rgr", "rstar"), switch_total = c(0,0,0,0,0,0,0), perc = c(0,0,0,0,0,0,0)))
df_row4_4$switches <- factor(df_row4_4$switches, levels = c("binit", "lot", "rep", "dispersal", "bstar", "rstar", "mor", "rgr", "root", "height", "pheno"))

df_row5 <- bind_rows((top_5pc_swG2_code[top_5pc_swG2_code$code == "Poapr_Schsc",] %>% mutate(Groups_on = as.character(Groups_on)) %>%
                        mutate(groups_ch = as.character(groups_ch))),
                     bind_cols(Groups_on = c("1", "1"), code = c("Poapr_Schsc", "Poapr_Schsc"),
                               groups_ch = c("Competition", "Overlap"), n = c(0, 0)))
df_row5$groups_ch <- factor(df_row5$groups_ch, levels = c("Colonisation", "Competition", "Growth", "Overlap"))

df_row5_4 <- bind_rows((top_5pc_sw_Split_switches[top_5pc_sw_Split_switches$code == "Poapr_Schsc",] %>% mutate(switches = as.character(switches)) %>%
                          mutate(switches = as.character(switches))),
                       bind_cols(code = c("Poapr_Schsc"), switches = c("dispersal", "lot", "mor",
                                                                       "pheno", "rep", "rgr", "rstar"),
                                 n = c(0,0,0,0,0,0,0), mech = c("dispersal", "lot", "mor","pheno", "rep", "rgr", "rstar"), switch_total = c(0,0,0,0,0,0,0), perc = c(0,0,0,0,0,0,0)))
df_row5_4$switches <- factor(df_row5_4$switches, levels = c("binit", "lot", "rep", "dispersal", "bstar", "rstar", "mor", "rgr", "root", "height", "pheno"))

df_row6 <- bind_rows((top_5pc_swG2_code[top_5pc_swG2_code$code == "Agrsc_Schsc",] %>% mutate(Groups_on = as.character(Groups_on)) %>%
                        mutate(groups_ch = as.character(groups_ch))),
                     bind_cols(Groups_on = c("1", "1"), code = c("Agrsc_Schsc", "Agrsc_Schsc"),
                               groups_ch = c("Competition", "Growth"), n = c(0, 0)))
df_row6$groups_ch <- factor(df_row6$groups_ch, levels = c("Colonisation", "Competition", "Growth", "Overlap"))

df_row6_4 <- bind_rows((top_5pc_sw_Split_switches[top_5pc_sw_Split_switches$code == "Agrsc_Schsc",] %>% mutate(switches = as.character(switches)) %>%
                          mutate(switches = as.character(switches))),
                       bind_cols(code = c("Agrsc_Schsc"), switches = c("lot", "rep", "dispersal", "rstar", "mor", "rgr"),
                                 n = c(0,0,0,0,0,0), mech = c("lot", "rep", "dispersal", "rstar", "mor", "rgr"), switch_total = c(0,0,0,0,0,0), perc = c(0,0,0,0,0,0)))
df_row6_4$switches <- factor(df_row5_4$switches, levels = c("binit", "lot", "rep", "dispersal", "bstar", "rstar", "mor", "rgr", "root", "height", "pheno"))


attribute_colours <- c("Binit" =  "#ADD1EC", "lottery" = "#6699CC", "fecun" = "#0072B2", "dispersal" = "#608B8B",
                       "B*" = "#CC6666", "R*" = "#DF4A6B",
                       "mortality" = "#EFCAB1", "RGR" = "#D55E00",
                       "root" = "#DFCBDE", "height" = "#CC79A7" , "pheno"= "#AA4499")

mech_colours <- c("Colonisation" = "#ADD1EC", "Growth" = "#EFCAB1",
                  "Niche \ndifferentiation" = "#DFCBDE", "Competition" = "#CC6666")

#### Agrostis scabra and Elymus repens
df_row1_4 <- df_row1_4 %>% dplyr::mutate(switches = as.character(switches)) %>%
  dplyr::mutate(switches = replace(switches, switches == "binit", "Binit"))%>%
  dplyr::mutate(switches = replace(switches, switches == "lot", "lottery"))%>%
  dplyr::mutate(switches = replace(switches, switches == "rep", "fecun"))%>%
  dplyr::mutate(switches = replace(switches, switches == "bstar", "B*"))%>%
  dplyr::mutate(switches = replace(switches, switches == "rstar", "R*"))%>%
  dplyr::mutate(switches = replace(switches, switches == "mor", "mortality"))%>%
  dplyr::mutate(switches = replace(switches, switches == "rgr", "RGR"))

df_row1_4$switches <- factor(df_row1_4$switches, levels = c("Binit", "lottery", "fecun", "dispersal", "B*", "R*", "mortality", "RGR", "root", "height", "pheno"))

FigS2a <- ggplot(df_row1_4, aes(x = switches, y = perc, fill = switches)) +
  geom_bar(stat = "identity") + scale_fill_manual(values = attribute_colours)  +
  theme_bw() + theme(legend.position = "none") + xlab("Attributes") + ylab("n") +
  theme(legend.position = "none", axis.text = element_text(size = 10),
        axis.title = element_text(size=12)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

top5pc_mech_Split <- top5pc_mech_Split %>% dplyr::mutate(mech = as.character(mech)) %>%
  dplyr::mutate(mech = replace(mech, mech == "binit", "Binit"))%>%
  dplyr::mutate(mech = replace(mech, mech == "lot", "lottery"))%>%
  dplyr::mutate(mech = replace(mech, mech == "rep", "fecun"))%>%
  dplyr::mutate(mech = replace(mech, mech == "bstar", "B*"))%>%
  dplyr::mutate(mech = replace(mech, mech == "rstar", "R*"))%>%
  dplyr::mutate(mech = replace(mech, mech == "mor", "mortality"))%>%
  dplyr::mutate(mech = replace(mech, mech == "rgr", "RGR"))

top5pc_mech_Split$mech <- factor(top5pc_mech_Split$mech, levels = c("Binit", "lottery", "fecun", "dispersal", "B*", "R*", "mortality", "RGR", "root", "height", "pheno"))

FigS2b <- ggplot(top5pc_mech_Split[top5pc_mech_Split$code == "Agrre_Agrsc",], aes(x = as.factor(mod_name), fill = mech)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = attribute_colours) +
  xlab("Model complexity (number of attributes)") + ylab("% in top models") +
  theme_bw() + theme(legend.position = "NULL")  +
  theme(legend.position = "none", axis.text = element_text(size = 10),
        axis.title = element_text(size=12)) +
  scale_y_continuous(breaks = c(0, 0.25,0.5,0.75, 1),
                     labels = c(0,25, 50,75, 100))

df_row1 <- df_row1 %>% dplyr::mutate(groups_ch = as.character(groups_ch)) %>%
  dplyr::mutate(groups_ch = replace(groups_ch, groups_ch == "Overlap", "Niche \ndifferentiation")) %>%
  mutate(groups_ch = as.factor(groups_ch))

FigS2c <- ggplot(df_row1, aes(x = groups_ch, y = n, fill = groups_ch)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = mech_colours) +
  theme_bw() + xlab("Mechanism") + theme(legend.position = "none") +
  theme(legend.position = "none", axis.text = element_text(size = 10),
        axis.title = element_text(size=12))

perc_top5pc_mechs <- perc_top5pc_mechs %>% dplyr::mutate(mech = as.character(mech)) %>%
  dplyr::mutate(mech = replace(mech, mech == "Overlap", "Niche \ndifferentiation"))

FigS2d <- ggplot(perc_top5pc_mechs[perc_top5pc_mechs$code == "Agrre_Agrsc",], aes(x = as.factor(Groups_on), y = perc, fill = mech)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = mech_colours) +
  xlab("Number of mechanisms") + ylab("% in top models") + theme_bw() +
  theme(legend.position = "none") +
  theme(legend.position = "none", axis.text = element_text(size = 10),
        axis.title = element_text(size=12))


##### Elymus repens and poa pratensis
df_row2_4 <- df_row2_4 %>% dplyr::mutate(switches = as.character(switches)) %>%
  dplyr::mutate(switches = replace(switches, switches == "binit", "Binit"))%>%
  dplyr::mutate(switches = replace(switches, switches == "lot", "lottery"))%>%
  dplyr::mutate(switches = replace(switches, switches == "rep", "fecun"))%>%
  dplyr::mutate(switches = replace(switches, switches == "bstar", "B*"))%>%
  dplyr::mutate(switches = replace(switches, switches == "rstar", "R*"))%>%
  dplyr::mutate(switches = replace(switches, switches == "mor", "mortality"))%>%
  dplyr::mutate(switches = replace(switches, switches == "rgr", "RGR"))

df_row2_4$switches <- factor(df_row2_4$switches, levels = c("Binit", "lottery", "fecun", "dispersal", "B*", "R*", "mortality", "RGR", "root", "height", "pheno"))


FigS2e <- ggplot(df_row2_4, aes(x = switches, y = perc, fill = switches)) +
  geom_bar(stat = "identity") + scale_fill_manual(values = attribute_colours)  +
  theme_bw() + theme(legend.position = "none") + xlab("Attributes") + ylab("n")  +
  theme(legend.position = "none", axis.text = element_text(size = 10),
        axis.title = element_text(size=12)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

FigS2f <- ggplot(top5pc_mech_Split[top5pc_mech_Split$code == "Agrre_Poapr",], aes(x = as.factor(mod_name), fill = mech)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = attribute_colours) +
  xlab("Model complexity (number of attributes)") + ylab("% in top models") +
  theme_bw() + theme(legend.position = "NULL")  +
  theme(legend.position = "none", axis.text = element_text(size = 10),
        axis.title = element_text(size=12)) +
  scale_y_continuous(breaks = c(0, 0.25,0.5,0.75, 1),
                     labels = c(0,25, 50,75, 100))

df_row3 <- bind_rows((top_5pc_swG2_code[top_5pc_swG2_code$code == "Agrre_Poapr",] %>% mutate(Groups_on = as.character(Groups_on)) %>%
                        mutate(groups_ch = as.character(groups_ch))),
                     bind_cols(Groups_on = c("1"), code = c("Agrre_Poapr"),
                               groups_ch = c("Overlap"), n = c(0)))

df_row3 <- df_row3 %>% dplyr::mutate(groups_ch = as.character(groups_ch)) %>%
  dplyr::mutate(groups_ch = replace(groups_ch, groups_ch == "Overlap", "Niche \ndifferentiation")) %>%
  mutate(groups_ch = as.factor(groups_ch))

FigS2g <- ggplot(df_row3, aes(x = groups_ch, y = n, fill = groups_ch)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = mech_colours) +
  theme_bw() + xlab("Mechanism") + theme(legend.position = "none") +
  theme(legend.position = "none", axis.text = element_text(size = 10),
        axis.title = element_text(size=12))

FigS2h <- ggplot(perc_top5pc_mechs[perc_top5pc_mechs$code == "Agrre_Poapr",], aes(x = as.factor(Groups_on), y = perc, fill = mech)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = mech_colours) +
  xlab("Number of mechanisms") + ylab("% in top models") + theme_bw() +
  theme(legend.position = "none") +
  theme(legend.position = "none", axis.text = element_text(size = 10),
        axis.title = element_text(size=12))



##### Elymus repens and Schizachyrium scoparium

df_row3_4 <- df_row3_4 %>% dplyr::mutate(switches = as.character(switches)) %>%
  dplyr::mutate(switches = replace(switches, switches == "binit", "Binit"))%>%
  dplyr::mutate(switches = replace(switches, switches == "lot", "lottery"))%>%
  dplyr::mutate(switches = replace(switches, switches == "rep", "fecun"))%>%
  dplyr::mutate(switches = replace(switches, switches == "bstar", "B*"))%>%
  dplyr::mutate(switches = replace(switches, switches == "rstar", "R*"))%>%
  dplyr::mutate(switches = replace(switches, switches == "mor", "mortality"))%>%
  dplyr::mutate(switches = replace(switches, switches == "rgr", "RGR"))

df_row3_4$switches <- factor(df_row3_4$switches, levels = c("Binit", "lottery", "fecun", "dispersal", "B*", "R*", "mortality", "RGR", "root", "height", "pheno"))

FigS2i <- ggplot(df_row3_4, aes(x = switches, y = perc, fill = switches)) +
  geom_bar(stat = "identity") + scale_fill_manual(values = attribute_colours)  +
  theme_bw() + theme(legend.position = "none") + xlab("Attributes") + ylab("n")  +
  theme(legend.position = "none", axis.text = element_text(size = 10),
        axis.title = element_text(size=12)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

FigS2j <- ggplot(top5pc_mech_Split[top5pc_mech_Split$code == "Agrre_Schsc",], aes(x = as.factor(mod_name), fill = mech)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = attribute_colours) +
  xlab("Model complexity (number of attributes)") + ylab("% in top models") +
  theme_bw() + theme(legend.position = "NULL")  +
  theme(legend.position = "none", axis.text = element_text(size = 10),
        axis.title = element_text(size=12)) +
  scale_y_continuous(breaks = c(0, 0.25,0.5,0.75, 1),
                     labels = c(0,25, 50,75, 100))

df_row4 <- bind_rows((top_5pc_swG2_code[top_5pc_swG2_code$code == "Agrre_Schsc",] %>% mutate(Groups_on = as.character(Groups_on)) %>%
                        mutate(groups_ch = as.character(groups_ch))),
                     bind_cols(Groups_on = c("1", "1"), code = c("Agrre_Schsc","Agrre_Schsc"),
                               groups_ch = c("Competition", "Colonisation"), n = c(0, 0)))
df_row4$groups_ch <- factor(df_row4$groups_ch, levels = c("Colonisation", "Competition", "Growth", "Overlap"))

df_row4 <- df_row4 %>% dplyr::mutate(groups_ch = as.character(groups_ch)) %>%
  dplyr::mutate(groups_ch = replace(groups_ch, groups_ch == "Overlap", "Niche \ndifferentiation")) %>%
  mutate(groups_ch = as.factor(groups_ch))

FigS2k <- ggplot(df_row4, aes(x = groups_ch, y = n, fill = groups_ch)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = mech_colours) +
  theme_bw() + xlab("Mechanism") + theme(legend.position = "none") +
  theme(legend.position = "none", axis.text = element_text(size = 10),
        axis.title = element_text(size=12))

FigS2l <- ggplot(perc_top5pc_mechs[perc_top5pc_mechs$code == "Agrre_Schsc",], aes(x = as.factor(Groups_on), y = perc, fill = mech)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = mech_colours) +
  xlab("Number of mechanisms") + ylab("% in top models") + theme_bw() +
  theme(legend.position = "none") +
  theme(legend.position = "none", axis.text = element_text(size = 10),
        axis.title = element_text(size=12))



##### Agrostis & Andropogon
df_row4_4 <- df_row4_4 %>% dplyr::mutate(switches = as.character(switches)) %>%
  dplyr::mutate(switches = replace(switches, switches == "binit", "Binit"))%>%
  dplyr::mutate(switches = replace(switches, switches == "lot", "lottery"))%>%
  dplyr::mutate(switches = replace(switches, switches == "rep", "fecun"))%>%
  dplyr::mutate(switches = replace(switches, switches == "bstar", "B*"))%>%
  dplyr::mutate(switches = replace(switches, switches == "rstar", "R*"))%>%
  dplyr::mutate(switches = replace(switches, switches == "mor", "mortality"))%>%
  dplyr::mutate(switches = replace(switches, switches == "rgr", "RGR"))

df_row4_4$switches <- factor(df_row4_4$switches, levels = c("Binit", "lottery", "fecun", "dispersal", "B*", "R*", "mortality", "RGR", "root", "height", "pheno"))

FigS2m <- ggplot(df_row4_4, aes(x = switches, y = perc, fill = switches)) +
  geom_bar(stat = "identity") + scale_fill_manual(values = attribute_colours)  +
  theme_bw() + theme(legend.position = "none") + xlab("Attributes") + ylab("n")  +
  theme(legend.position = "none", axis.text = element_text(size = 10),
        axis.title = element_text(size=12)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

FigS2n <- ggplot(top5pc_mech_Split[top5pc_mech_Split$code == "Agrsc_Andge",], aes(x = as.factor(mod_name), fill = mech)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = attribute_colours) +
  xlab("Model complexity (number of attributes)") + ylab("% in top models") +
  theme_bw() + theme(legend.position = "NULL")  +
  theme(legend.position = "none", axis.text = element_text(size = 10),
        axis.title = element_text(size=12)) +
  scale_y_continuous(breaks = c(0, 0.25,0.5,0.75, 1),
                     labels = c(0,25, 50,75, 100))

df_row2 <- df_row2 %>% dplyr::mutate(groups_ch = as.character(groups_ch)) %>%
  dplyr::mutate(groups_ch = replace(groups_ch, groups_ch == "Overlap", "Niche \ndifferentiation")) %>%
  mutate(groups_ch = as.factor(groups_ch))

FigS2o <- ggplot(df_row2, aes(x = groups_ch, y = n, fill = groups_ch)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = mech_colours) +
  theme_bw() + xlab("Mechanism") + theme(legend.position = "none") +
  theme(legend.position = "none", axis.text = element_text(size = 10),
        axis.title = element_text(size=12))

FigS2p <- ggplot(perc_top5pc_mechs[perc_top5pc_mechs$code == "Agrsc_Andge",], aes(x = as.factor(Groups_on), y = perc, fill = mech)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = mech_colours) +
  xlab("Number of mechanisms") + ylab("% in top models") + theme_bw() +
  theme(legend.position = "none") +
  theme(legend.position = "none", axis.text = element_text(size = 10),
        axis.title = element_text(size=12))

### Poa & Schizachyrium
df_row5_4 <- df_row5_4 %>% dplyr::mutate(switches = as.character(switches)) %>%
  dplyr::mutate(switches = replace(switches, switches == "binit", "Binit"))%>%
  dplyr::mutate(switches = replace(switches, switches == "lot", "lottery"))%>%
  dplyr::mutate(switches = replace(switches, switches == "rep", "fecun"))%>%
  dplyr::mutate(switches = replace(switches, switches == "bstar", "B*"))%>%
  dplyr::mutate(switches = replace(switches, switches == "rstar", "R*"))%>%
  dplyr::mutate(switches = replace(switches, switches == "mor", "mortality"))%>%
  dplyr::mutate(switches = replace(switches, switches == "rgr", "RGR"))

df_row5_4$switches <- factor(df_row5_4$switches, levels = c("Binit", "lottery", "fecun", "dispersal", "B*", "R*", "mortality", "RGR", "root", "height", "pheno"))

FigS2q <- ggplot(df_row5_4, aes(x = switches, y = perc, fill = switches)) +
  geom_bar(stat = "identity") + scale_fill_manual(values = attribute_colours)  +
  theme_bw() + theme(legend.position = "none") + xlab("Attributes") + ylab("n")   +
  theme(legend.position = "none", axis.text = element_text(size = 10),
        axis.title = element_text(size=12))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

FigS2r <- ggplot(top5pc_mech_Split[top5pc_mech_Split$code == "Poapr_Schsc",], aes(x = as.factor(mod_name), fill = mech)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = attribute_colours) +
  xlab("Model complexity (number of attributes)") + ylab("% in top models") +
  theme_bw() + theme(legend.position = "NULL")  +
  theme(legend.position = "none", axis.text = element_text(size = 10),
        axis.title = element_text(size=12)) +
  scale_y_continuous(breaks = c(0, 0.25,0.5,0.75, 1),
                     labels = c(0,25, 50,75, 100))

df_row5 <- bind_rows((top_5pc_swG2_code[top_5pc_swG2_code$code == "Poapr_Schsc",] %>% mutate(Groups_on = as.character(Groups_on)) %>%
                        mutate(groups_ch = as.character(groups_ch))),
                     bind_cols(Groups_on = c("1", "1"), code = c("Poapr_Schsc", "Poapr_Schsc"),
                               groups_ch = c("Competition", "Overlap"), n = c(0, 0)))
df_row5$groups_ch <- factor(df_row5$groups_ch, levels = c("Colonisation", "Competition", "Growth", "Overlap"))

df_row5 <- df_row5 %>% dplyr::mutate(groups_ch = as.character(groups_ch)) %>%
  dplyr::mutate(groups_ch = replace(groups_ch, groups_ch == "Overlap", "Niche \ndifferentiation")) %>%
  mutate(groups_ch = as.factor(groups_ch))

FigS2s <- ggplot(df_row5, aes(x = groups_ch, y = n, fill = groups_ch)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = mech_colours) +
  theme_bw() + xlab("Mechanism") + theme(legend.position = "none") +
  theme(legend.position = "none", axis.text = element_text(size = 10),
        axis.title = element_text(size=12))

FigS2t <- ggplot(perc_top5pc_mechs[perc_top5pc_mechs$code == "Poapr_Schsc",], aes(x = as.factor(Groups_on), y = perc, fill = mech)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = mech_colours) +
  xlab("Number of mechanisms") + ylab("% in top models") + theme_bw() +
  theme(legend.position = "none") +
  theme(legend.position = "none", axis.text = element_text(size = 10),
        axis.title = element_text(size=12))

### Agrostis scabra and Schizachyrium scoparium
df_row6_4 <- df_row6_4 %>% dplyr::mutate(switches = as.character(switches)) %>%
  dplyr::mutate(switches = replace(switches, switches == "binit", "Binit"))%>%
  dplyr::mutate(switches = replace(switches, switches == "lot", "lottery"))%>%
  dplyr::mutate(switches = replace(switches, switches == "rep", "fecun"))%>%
  dplyr::mutate(switches = replace(switches, switches == "bstar", "B*"))%>%
  dplyr::mutate(switches = replace(switches, switches == "rstar", "R*"))%>%
  dplyr::mutate(switches = replace(switches, switches == "mor", "mortality"))%>%
  dplyr::mutate(switches = replace(switches, switches == "rgr", "RGR"))

df_row6_4$switches <- factor(df_row6_4$switches, levels = c("Binit", "lottery", "fecun", "dispersal", "B*", "R*", "mortality", "RGR", "root", "height", "pheno"))


FigS2u <- ggplot(df_row6_4, aes(x = switches, y = perc, fill = switches)) +
  geom_bar(stat = "identity") + scale_fill_manual(values = attribute_colours)  +
  theme_bw() + theme(legend.position = "none") + xlab("Attributes") + ylab("n")  +
  theme(legend.position = "none", axis.text = element_text(size = 10),
        axis.title = element_text(size=12)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

FigS2v <- ggplot(top5pc_mech_Split[top5pc_mech_Split$code == "Agrsc_Schsc",], aes(x = as.factor(mod_name), fill = mech)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = attribute_colours) +
  xlab("Model complexity (number of attributes)") + ylab("% in top models") +
  theme_bw() + theme(legend.position = "NULL")  +
  theme(legend.position = "none", axis.text = element_text(size = 10),
        axis.title = element_text(size=12)) +
  scale_y_continuous(breaks = c(0, 0.25,0.5,0.75, 1),
                     labels = c(0,25, 50,75, 100))

df_row6 <- bind_rows((top_5pc_swG2_code[top_5pc_swG2_code$code == "Agrsc_Schsc",] %>% mutate(Groups_on = as.character(Groups_on)) %>%
                        mutate(groups_ch = as.character(groups_ch))),
                     bind_cols(Groups_on = c("1", "1"), code = c("Agrsc_Schsc", "Agrsc_Schsc"),
                               groups_ch = c("Competition", "Growth"), n = c(0, 0)))
df_row6$groups_ch <- factor(df_row6$groups_ch, levels = c("Colonisation", "Competition", "Growth", "Overlap"))



df_row6 <- df_row6 %>% dplyr::mutate(groups_ch = as.character(groups_ch)) %>%
  dplyr::mutate(groups_ch = replace(groups_ch, groups_ch == "Overlap", "Niche \ndifferentiation")) %>%
  mutate(groups_ch = as.factor(groups_ch))

FigS2w <- ggplot(df_row6, aes(x = groups_ch, y = n, fill = groups_ch)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = mech_colours) +
  theme_bw() + xlab("Mechanism") + theme(legend.position = "none") +
  theme(legend.position = "none", axis.text = element_text(size = 10),
        axis.title = element_text(size=12))

FigS2x <- ggplot(perc_top5pc_mechs[perc_top5pc_mechs$code == "Agrsc_Schsc",], aes(x = as.factor(Groups_on), y = perc, fill = mech)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = mech_colours) +
  xlab("Number of mechanisms") + ylab("% in top models") + theme_bw() +
  theme(legend.position = "none") +
  theme(legend.position = "none", axis.text = element_text(size = 10),
        axis.title = element_text(size=12))



ggarrange(((annotate_figure((ggarrange(FigS2a, FigS2b, FigS2c, FigS2d, ncol = 4, align = "hv", labels = c("A", "B", "C", "D"))),
                            top = text_grob(expression(italic("Elymus repens") ~ "and" ~ italic("Agrostis scabra")), color = "black", size = 16)))),
          ((annotate_figure((ggarrange(FigS2e, FigS2f, FigS2g, FigS2h, ncol = 4, align = "hv", labels = c("E", "F", "G", "H"))),
                            top = text_grob(expression(italic("Elymus repens") ~ "and" ~ italic("Poa pratensis")), color = "black", size = 16)))),
          ((annotate_figure((ggarrange(FigS2i, FigS2j, FigS2k, FigS2l, ncol = 4, align = "hv", labels = c("I", "J", "K", "L"))),
                            top = text_grob(expression(italic("Elymus repens") ~ "and" ~ italic("Schizachyrium scoparium")), color = "black", size = 16)))),
          ((annotate_figure((ggarrange(FigS2m, FigS2n, FigS2o, FigS2p, ncol = 4, align = "hv", labels = c("M", "N", "O", "P"))),
                            top = text_grob(expression(italic("Agrostis scabra") ~ "and" ~ italic("Andropogon gerardi")), color = "black", size = 16)))),
          ((annotate_figure((ggarrange(FigS2q, FigS2r, FigS2s, FigS2t, ncol = 4, align = "hv", labels = c("M", "N", "O", "P"))),
                            top = text_grob(expression(italic("Poa pratensis") ~ "and" ~ italic("Schizachyrium scoparium")), color = "black", size = 16)))),
          ((annotate_figure((ggarrange(FigS2u, FigS2v, FigS2w, FigS2x, ncol = 4, align = "hv", labels = c("U", "V", "W", "X"))),
                            top = text_grob(expression(italic("Agrostis scabra") ~ "and" ~ italic("Schizachyrium scoparium")), color = "black", size = 16)))),
          align = "hv", ncol = 1, nrow = 6)


###### Figure S3 - frequency of each attribute for each species
library(ggradar)
library(scales)
library(showtext)

results_species <- results_main <- results_sp_rich <-
  tibble(mod_name = numeric(),
         binit = numeric(),
         bstar = numeric(),
         dispersal = numeric(),
         height = numeric(),
         lot = numeric(),
         mor = numeric(),
         pheno = numeric(),
         rep = numeric(),
         rgr = numeric(),
         root = numeric(),
         rstar = numeric())

for(i in 1:11) {
  
  load(paste0("non_zero_replicates_", i, "_year6.rda"))
  
  non_zero2 <- non_zero_replicates
  
  radar_main <- non_zero2 %>%
    group_by(species, sp_rich, plot, subplot, switches) %>%
    summarise(rmse = mean(rmse, na.rm = TRUE)) %>% ## addition: gives average RMSE for each combination (average of 10 simulations)
    ungroup() %>%
    group_by(sp_rich, plot, subplot, switches) %>%
    summarise(mean_rmse = mean(rmse, na.rm = TRUE)) %>% ## addition: gives average RMSE for each combination (average of 10 simulations)
    ungroup() %>%
    filter(sp_rich > 1) %>%
    mutate(mae_top5 = quantile(mean_rmse, 0.05)) %>%
    filter(mean_rmse <= mae_top5) ## then finds top 5% performing models
  
  
  mech_main <- radar_main %>%
    select(switches) %>%
    count(switches) %>% ## how many times does each appear
    separate(switches, into = paste0("sw", 1:11)) %>%
    gather(delete, mech, -n) %>%
    select(-delete) %>%
    na.omit() %>%
    group_by(mech) %>%
    summarise(n = sum(n)) %>%
    mutate(n = n / nrow(radar_main)) %>%
    spread(mech, n) %>%
    mutate(mod_name = i)
  
  results_main <- bind_rows(results_main, mech_main)
  
  radar_species <- non_zero2 %>%
    group_by(species, sp_rich, plot, subplot, switches) %>%
    summarise(mean_rmse = mean(rmse, na.rm = TRUE)) %>% ## addition: gives average RMSE for each combination (average of 10 simulations) %>%
    ungroup() %>%
    filter(sp_rich > 1) %>%
    group_by(species) %>%
    mutate(mae_top5 = quantile(mean_rmse, 0.05)) %>%
    filter(mean_rmse <= mae_top5)
  
  mod_per_species <- radar_species %>% group_by(species) %>% summarise(count = n())
  
  mech_species <- radar_species %>%
    select(switches) %>%
    count(switches) %>%
    separate(switches, into = paste0("sw", 1:11)) %>%
    gather(delete, mech, -n, -species) %>%
    select(-delete) %>%
    na.omit() %>%
    group_by(mech, species) %>%
    summarise(n = sum(n)) %>%
    inner_join(mod_per_species) %>%
    mutate(n = n /count) %>%
    spread(mech, n) %>%
    mutate(mod_name = i) %>%
    select(-count)
  
  results_species <- bind_rows(results_species, mech_species)
  
  radar_sp_rich <- non_zero2 %>%
    group_by(species, sp_rich, plot, subplot, switches) %>%
    summarise(rmse = mean(rmse, na.rm = TRUE)) %>% ## addition: gives average RMSE for each combination (average of 10 simulations) %>%
    ungroup() %>%
    group_by(sp_rich, plot, subplot, switches) %>%
    summarise(mean_rmse = mean(rmse, na.rm = TRUE)) %>% ## addition: gives average RMSE for each combination (average of 10 simulations)
    ungroup() %>%
    group_by(sp_rich) %>%
    mutate(mae_top5 = quantile(mean_rmse, 0.05)) %>%
    filter(mean_rmse <= mae_top5)
  
  mod_per_sp_rich <- radar_sp_rich %>% group_by(sp_rich) %>% summarise(count = n())
  
  mech_sp_rich <- radar_sp_rich %>%
    select(switches) %>%
    count(switches) %>%
    separate(switches, into = paste0("sw", 1:11)) %>%
    gather(delete, mech, -n, -sp_rich) %>%
    select(-delete) %>%
    na.omit() %>%
    group_by(mech, sp_rich) %>%
    summarise(n = sum(n)) %>%
    inner_join(mod_per_sp_rich) %>%
    mutate(n = n / count) %>%
    spread(mech, n) %>%
    mutate(mod_name = i) %>%
    select(-count)
  
  results_sp_rich <- bind_rows(results_sp_rich, mech_sp_rich)
  
} ## see how to deal with NAs


##### Figure S4 ### species level radar plots for groups

results_species <- results_species %>%
  mutate(species = replace(species, species == "Agropyron repens", "Elymus repens"))

results_species2 <- results_species %>%
  rename(Binit = binit, lottery = lot, fecun = rep, `B*` = bstar, `R*` = rstar, mortality = mor, RGR = rgr)

results_species2$mod_name <- factor(results_species2$mod_name, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"))

plt_species <- results_species2 %>%
  replace_na(list(Binit = 0,
                  `B*` = 0,
                  dispersal = 0,
                  height = 0,
                  lottery = 0,
                  mortality = 0,
                  pheno = 0,
                  fecun = 0,
                  RGR = 0,
                  root = 0,
                  `R*` = 0)) %>%
  group_by(species) %>%
  nest() %>%
  mutate(plt = map2(data, species, function(x, y) {
    ggradar(x, plot.title = paste0(y),
            group.point.size = 0.7,
            grid.label.size = 4,
            axis.label.size = 4,
            group.line.width = 1.3) +
      theme(plot.title = element_text(size = 15,
                                      face = "bold"))}))


ggarrange(plt_species$plt[[2]], plt_species$plt[[3]], plt_species$plt[[1]], plt_species$plt[[4]], plt_species$plt[[5]], common.legend = TRUE, legend = "right")



###### Figure S4 - top 5% of best performing models


## 1groups

non_zero_G1 <- non_zero_wswitches %>%
  filter(ngrowth < 2 & noverlap < 2 & ncolo < 2 & ncomp < 2) %>%
  select(switches:seed_mass, GroupNm1, group_if1)


dfZG1 <- non_zero_G1 %>%
  rename(groups_ch = GroupNm1, Groups_on = group_if1) %>%
  mutate(Groups_on = as.numeric(Groups_on)) %>%
  group_by(species, sp_rich, plot, subplot, groups_ch, Groups_on) %>%
  filter(sp_rich > 1) %>%
  summarise(rmse = mean(rmse)) %>%
  ungroup() %>%
  group_by(sp_rich, plot, subplot, groups_ch, Groups_on) %>%
  summarise(rmse = mean(rmse)) %>%
  ungroup() ## gives dataframe of output replicates with model names

top_5pcG1 <- dfZG1 %>%
  group_by(Groups_on) %>%
  mutate(pc95 = quantile(rmse, 0.05)) %>%
  filter(rmse <= pc95) ### top 5% quantile of RMSE: i.e. best performing models, from each amount of switches

top_5pc_swG1 <- top_5pcG1 %>%
  filter(Groups_on == 1) %>% ##filters out single switch models
  count(groups_ch) %>%
  mutate(totaln = sum(n)) %>%
  mutate(perc = (n/totaln)*100) %>%
  mutate(groups_ch = as.factor(groups_ch),
         groups_ch = fct_reorder(groups_ch, perc, .desc = TRUE))

top_5pc_swG1 <- top_5pc_swG1 %>%
  mutate(groups_ch = as.character(groups_ch)) %>%
  dplyr::mutate(groups_ch = replace(groups_ch, groups_ch == "Overlap", "Niche differentiation")) %>%
  mutate(groups_ch = as.factor(groups_ch))

top_5pc_swG1$groups_ch <- factor(top_5pc_swG1$groups_ch, levels = c("Competition", "Niche differentiation",
                                                                    "Colonisation", "Growth"))

(PanelS4B <- ggplot(top_5pc_swG1, aes(x = groups_ch, y = perc, colour = groups_ch, fill = groups_ch)) +
    geom_bar(stat = "identity") +
    scale_colour_manual(values = c("Colonisation" = "#ADD1EC", "Growth" = "#EFCAB1",
                                   "Niche differentiation" = "#DFCBDE", "Competition" = "#CC6666")) +
    scale_fill_manual(values = c("Colonisation" = "#ADD1EC", "Growth" = "#EFCAB1",
                                 "Niche differentiation" = "#DFCBDE", "Competition" = "#CC6666")) +
    theme_bw() +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(size = 9),
          axis.text.y = element_text(size = 11),
          axis.title.x = element_text(size=12),
          axis.title.y = element_text(size = 12)) +
    xlab("Mechanism") + ylab("% in top models")) # 750 x 300

top_5pc_mechG1 <- top_5pcG1 %>%
  select(groups_ch, Groups_on) %>%
  filter(Groups_on %in% 1:4) %>% ## gets rid of all switches on (full model)
  separate(groups_ch, into = paste0("sw", 1:4)) %>% ## separates out models
  gather(delete, mech, -Groups_on) %>%
  select(-delete) %>%
  na.omit() %>%
  mutate(count = 1, Groups_on = as.factor(Groups_on)) %>%
  group_by(Groups_on, mech) %>%
  summarise(total_count = sum(count)) %>% ungroup() %>%
  na.omit() %>% group_by(Groups_on) %>% mutate(sum = sum(total_count)) %>%
  mutate(perc = (total_count/sum)*100)

top_5pc_mechG1 <- top_5pc_mechG1 %>%
  mutate(mech = as.character(mech)) %>%
  dplyr::mutate(mech = replace(mech, mech == "Overlap", "Niche differentiation")) %>%
  mutate(mech = as.factor(mech))

top_5pc_mechG1$mech <- factor(top_5pc_mechG1$mech, levels = c("Competition", "Niche differentiation",
                                                              "Colonisation", "Growth"))

(PanelS4c <- ggplot(top_5pc_mechG1, aes(x = as.factor(Groups_on), y = perc, fill = mech)) +
    geom_bar(stat = "identity") +
    scale_colour_manual(values = c("Colonisation" = "#ADD1EC", "Growth" = "#EFCAB1",
                                   "Niche differentiation" = "#DFCBDE", "Competition" = "#CC6666")) +
    scale_fill_manual(values = c("Colonisation" = "#ADD1EC", "Growth" = "#EFCAB1",
                                 "Niche differentiation" = "#DFCBDE", "Competition" = "#CC6666")) +
    xlab("Number of mechanisms") + # 750 x 350
    ylab("") + theme_bw() +
    ylab("% in top models") +
    guides(fill=guide_legend("Mechanism"), colour = guide_legend("Mechanism"))+
    theme(legend.position = "right",
          axis.text = element_text(size = 9),
          axis.title = element_text(size=11)) +
    guides(fill=guide_legend("Mechanism", order = 1, ncol = 1))) ## switches that appear in the top 5% of each number of switches

(PanelS4c_noleg <- ggplot(top_5pc_mechG1, aes(x = as.factor(Groups_on), y = perc, fill = mech)) +
    geom_bar(stat = "identity") +
    scale_colour_manual(values = c("Colonisation" = "#ADD1EC", "Growth" = "#EFCAB1",
                                   "Niche differentiation" = "#DFCBDE", "Competition" = "#CC6666")) +
    scale_fill_manual(values = c("Colonisation" = "#ADD1EC", "Growth" = "#EFCAB1",
                                 "Niche differentiation" = "#DFCBDE", "Competition" = "#CC6666")) +
    xlab("Number of mechanisms") + # 750 x 350
    ylab("") + theme_bw() +
    ylab("% in top models") +
    guides(fill=guide_legend("Mechanism"), colour = guide_legend("Mechanism"))+
    theme(legend.position = "none",
          axis.text = element_text(size = 9),
          axis.title = element_text(size=11)))

results_mainG1 <- tibble(mod_name = numeric(),
                         Growth = as.numeric(),
                         Colonisation = as.numeric(),
                         Competition = as.numeric(),
                         Overlap = as.numeric())

for(i in 1:4) {

  dfZZG1 <- dfZG1 %>% filter(Groups_on == i) %>% dplyr::filter(groups_ch != "None")

  radar_mainG1 <- dfZZG1 %>%
    filter(sp_rich > 1) %>%
    mutate(mae_top5 = quantile(rmse, 0.05)) %>%
    filter(rmse <= mae_top5) ## then finds top 5% performing models


  mech_mainG1 <- radar_mainG1 %>%
    select(groups_ch) %>%
    count(groups_ch) %>% ## how many times does each appear
    separate(groups_ch, into = paste0("sw", 1:4)) %>%
    gather(delete, mech, -n) %>%
    select(-delete) %>%
    na.omit() %>%
    group_by(mech) %>%
    summarise(n = sum(n)) %>%
    mutate(n = n / nrow(radar_mainG1)) %>%
    spread(mech, n) %>%
    mutate(mod_name = i)

  results_mainG1 <- bind_rows(results_mainG1, mech_mainG1)

} ## see how to deal with NAs

results_mainG1$mod_name <- factor(results_mainG1$mod_name, levels = c("1", "2", "3", "4"))
results_mainG1 <- results_mainG1 %>% select(mod_name, Colonisation, Growth, `Niche \ndifferentiation`= Overlap, Comp = Competition)

PanelS4A <- ggradar(results_mainG1,
                    axis.label.size = 6, # Afftects the names of the variables
                    grid.label.size = 4,
                    group.point.size = 3,# Simply the size of the point
                    legend.position = "bottom")  +
  theme_minimal() +
  theme(axis.line = element_blank(),   # Remove axis lines
        axis.text = element_blank(),   # Remove axis text (numbers)
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),   # Remove major grid lines
        panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  #  scale_colour_manual(values = c("1" = "#FEC47D", "2" = "#65A688",
  #                                 "3" = "#B54445", "4" = "#3192C4")) +
  theme(panel.grid=element_line(size=0.1)) +
  labs(color = "Number of \nmechanisms") # 700 x 450


## switches that appear in the top 5% of each number of switches


### 3 groups

non_zero_G3 <- non_zero_wswitches %>%
  filter(ngrowth < 3 & noverlap < 4 & ncolo < 4 & ncomp < 3 &
           ngrowth != 1 & !noverlap %in% c(1,2) & !ncolo %in% c(1,2) & ncomp != 1)

dfZG3 <- non_zero_G3 %>%
  rename(groups_ch = GroupNm1, Groups_on = group_if1) %>%
  mutate(Groups_on = as.numeric(Groups_on)) %>%
  group_by(species, sp_rich, plot, subplot, groups_ch, Groups_on) %>%
  filter(sp_rich > 1) %>%
  summarise(rmse = mean(rmse)) %>%
  ungroup() %>%
  group_by(sp_rich, plot, subplot, groups_ch, Groups_on) %>%
  summarise(rmse = mean(rmse)) %>%
  ungroup() ## gives dataframe of output replicates with model names

top_5pcG3 <- dfZG3 %>%
  group_by(Groups_on) %>%
  mutate(pc95 = quantile(rmse, 0.05)) %>%
  filter(rmse <= pc95) ### top 5% quantile of RMSE: i.e. best performing models, from each amount of switches

top_5pc_swG3 <- top_5pcG3 %>%
  filter(Groups_on == 1) %>% ##filters out single switch models
  count(groups_ch) %>%
  mutate(totaln = sum(n)) %>%
  mutate(perc = (n/totaln)*100) %>%
  mutate(groups_ch = as.factor(groups_ch),
         groups_ch = fct_reorder(groups_ch, perc, .desc = TRUE))

top_5pc_swG3 <- top_5pc_swG3 %>%
  mutate(groups_ch = as.character(groups_ch)) %>%
  dplyr::mutate(groups_ch = replace(groups_ch, groups_ch == "Overlap", "Niche differentiation")) %>%
  mutate(groups_ch = as.factor(groups_ch))

top_5pc_swG3$groups_ch <- factor(top_5pc_swG3$groups_ch, levels = c("Competition", "Niche differentiation",
                                                                    "Colonisation", "Growth"))

(PanelS4E <- ggplot(top_5pc_swG3, aes(x = groups_ch, y = perc, colour = groups_ch, fill = groups_ch)) +
    geom_bar(stat = "identity") +
    scale_colour_manual(values = c("Colonisation" = "#ADD1EC", "Growth" = "#EFCAB1",
                                   "Niche differentiation" = "#DFCBDE", "Competition" = "#CC6666")) +
    scale_fill_manual(values = c("Colonisation" = "#ADD1EC", "Growth" = "#EFCAB1",
                                 "Niche differentiation" = "#DFCBDE", "Competition" = "#CC6666")) +
    theme_bw() +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(size = 9),
          axis.text.y = element_text(size = 11),
          axis.title.x = element_text(size=12),
          axis.title.y = element_text(size = 12)) +
    xlab("Mechanism") + ylab("% in top models")) # 750 x 300

top_5pc_mechG3 <- top_5pcG3 %>%
  select(groups_ch, Groups_on) %>%
  filter(Groups_on %in% 1:4) %>% ## gets rid of all switches on (full model)
  separate(groups_ch, into = paste0("sw", 1:4)) %>% ## separates out models
  gather(delete, mech, -Groups_on) %>%
  select(-delete) %>%
  na.omit() %>%
  mutate(count = 1, Groups_on = as.factor(Groups_on)) %>%
  group_by(Groups_on, mech) %>%
  summarise(total_count = sum(count)) %>% ungroup() %>%
  na.omit() %>% group_by(Groups_on) %>% mutate(sum = sum(total_count)) %>%
  mutate(perc = (total_count/sum)*100)

top_5pc_mechG3 <- top_5pc_mechG3 %>%
  mutate(mech = as.character(mech)) %>%
  dplyr::mutate(mech = replace(mech, mech == "Overlap", "Niche differentiation")) %>%
  mutate(mech = as.factor(mech))

top_5pc_mechG3$mech <- factor(top_5pc_mechG3$mech, levels = c("Competition", "Niche differentiation",
                                                              "Colonisation", "Growth"))

(PanelS4F <- ggplot(top_5pc_mechG3, aes(x = as.factor(Groups_on), y = perc, fill = mech)) +
    geom_bar(stat = "identity") +
    scale_colour_manual(values = c("Colonisation" = "#ADD1EC", "Growth" = "#EFCAB1",
                                   "Niche differentiation" = "#DFCBDE", "Competition" = "#CC6666")) +
    scale_fill_manual(values = c("Colonisation" = "#ADD1EC", "Growth" = "#EFCAB1",
                                 "Niche differentiation" = "#DFCBDE", "Competition" = "#CC6666")) +
    xlab("Number of mechanisms") + # 750 x 350
    ylab("") + theme_bw() +
    ylab("% in top models") +
    guides(fill=guide_legend("Mechanism"), colour = guide_legend("Mechanism"))+
    theme(legend.position = "right",
          axis.text = element_text(size = 9),
          axis.title = element_text(size=11)) +
    guides(fill=guide_legend("Mechanism", order = 1, ncol = 1))) ## switches that appear in the top 5% of each number of switches


results_mainG3 <- tibble(mod_name = numeric(),
                         Growth = as.numeric(),
                         Colonisation = as.numeric(),
                         Competition = as.numeric(),
                         Overlap = as.numeric())

for(i in 1:4) {

  dfZZG3 <- dfZG3 %>% filter(Groups_on == i) %>% dplyr::filter(groups_ch != "None")

  radar_mainG3 <- dfZZG3 %>%
    filter(sp_rich > 1) %>%
    mutate(mae_top5 = quantile(rmse, 0.05)) %>%
    filter(rmse <= mae_top5) ## then finds top 5% performing models


  mech_mainG3 <- radar_mainG3 %>%
    select(groups_ch) %>%
    count(groups_ch) %>% ## how many times does each appear
    separate(groups_ch, into = paste0("sw", 1:4)) %>%
    gather(delete, mech, -n) %>%
    select(-delete) %>%
    na.omit() %>%
    group_by(mech) %>%
    summarise(n = sum(n)) %>%
    mutate(n = n / nrow(radar_mainG3)) %>%
    spread(mech, n) %>%
    mutate(mod_name = i)

  results_mainG3 <- bind_rows(results_mainG3, mech_mainG3)

} ## see how to deal with NAs

results_mainG3$mod_name <- factor(results_mainG3$mod_name, levels = c("1", "2", "3", "4"))
results_mainG3 <- results_mainG3 %>% select(mod_name, Colonisation, Growth, `Niche \ndifferentiation`= Overlap, Comp = Competition)

PanelS4D <- ggradar(results_mainG3,
                    axis.label.size = 6, # Afftects the names of the variables
                    grid.label.size = 4,
                    group.point.size = 3,# Simply the size of the point
                    legend.position = "bottom")  +
  theme_minimal() +
  theme(axis.line = element_blank(),   # Remove axis lines
        axis.text = element_blank(),   # Remove axis text (numbers)
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),   # Remove major grid lines
        panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  #  scale_colour_manual(values = c("1" = "#FEC47D", "2" = "#65A688",
  #                                 "3" = "#B54445", "4" = "#3192C4")) +
  theme(panel.grid=element_line(size=0.1)) +
  labs(color = "Number of \nmechanisms") # 700 x 450


(PanelS4F_noleg <- ggplot(top_5pc_mechG3, aes(x = as.factor(Groups_on), y = perc, fill = mech)) +
    geom_bar(stat = "identity") +
    scale_colour_manual(values = c("Colonisation" = "#ADD1EC", "Growth" = "#EFCAB1",
                                   "Niche differentiation" = "#DFCBDE", "Competition" = "#CC6666")) +
    scale_fill_manual(values = c("Colonisation" = "#ADD1EC", "Growth" = "#EFCAB1",
                                 "Niche differentiation" = "#DFCBDE", "Competition" = "#CC6666")) +
    xlab("Number of mechanisms") + # 750 x 350
    ylab("") + theme_bw() +
    ylab("% in top models") +
    guides(fill=guide_legend("Mechanism"), colour = guide_legend("Mechanism"))+
    theme(legend.position = "none",
          axis.text = element_text(size = 11),
          axis.title = element_text(size=12)) +
    guides(fill=guide_legend("Mechanism", order = 1, ncol = 1)))



#### 4 attributes

non_zero_G4 <- non_zero_wswitches %>%
  filter(ngrowth != 1 & !noverlap %in% c(1,2) & !ncolo %in% c(1,2,3) & ncomp != 1)

dfZG4 <- non_zero_G4 %>%
  rename(groups_ch = GroupNm1, Groups_on = group_if1) %>%
  mutate(Groups_on = as.numeric(Groups_on)) %>%
  group_by(species, sp_rich, plot, subplot, groups_ch, Groups_on) %>%
  filter(sp_rich > 1) %>%
  summarise(rmse = mean(rmse)) %>%
  ungroup() %>%
  group_by(sp_rich, plot, subplot, groups_ch, Groups_on) %>%
  summarise(rmse = mean(rmse)) %>%
  ungroup() ## gives dataframe of output replicates with model names

top_5pcG4 <- dfZG4 %>%
  group_by(Groups_on) %>%
  mutate(pc95 = quantile(rmse, 0.05)) %>%
  filter(rmse <= pc95) ### top 5% quantile of RMSE: i.e. best performing models, from each amount of switches

top_5pc_swG4 <- top_5pcG4 %>%
  filter(Groups_on == 1) %>% ##filters out single switch models
  count(groups_ch) %>%
  mutate(totaln = sum(n)) %>%
  mutate(perc = (n/totaln)*100) %>%
  mutate(groups_ch = as.factor(groups_ch),
         groups_ch = fct_reorder(groups_ch, perc, .desc = TRUE))

top_5pc_swG4 <- top_5pc_swG4 %>%
  mutate(groups_ch = as.character(groups_ch)) %>%
  dplyr::mutate(groups_ch = replace(groups_ch, groups_ch == "Overlap", "Niche differentiation")) %>%
  mutate(groups_ch = as.factor(groups_ch))

top_5pc_swG4$groups_ch <- factor(top_5pc_swG4$groups_ch, levels = c("Competition", "Niche differentiation",
                                                                    "Colonisation", "Growth"))

(PanelS4H <- ggplot(top_5pc_swG4, aes(x = groups_ch, y = perc, colour = groups_ch, fill = groups_ch)) +
    geom_bar(stat = "identity") +
    scale_colour_manual(values = c("Colonisation" = "#ADD1EC", "Growth" = "#EFCAB1",
                                   "Niche differentiation" = "#DFCBDE", "Competition" = "#CC6666")) +
    scale_fill_manual(values = c("Colonisation" = "#ADD1EC", "Growth" = "#EFCAB1",
                                 "Niche differentiation" = "#DFCBDE", "Competition" = "#CC6666")) +
    theme_bw() +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(size = 9),
          axis.text.y = element_text(size = 11),
          axis.title.x = element_text(size=12),
          axis.title.y = element_text(size = 12)) +
    xlab("Mechanism") + ylab("% in top models")) # 750 x 300

top_5pc_mechG4 <- top_5pcG4 %>%
  select(groups_ch, Groups_on) %>%
  filter(Groups_on %in% 1:4) %>% ## gets rid of all switches on (full model)
  separate(groups_ch, into = paste0("sw", 1:4)) %>% ## separates out models
  gather(delete, mech, -Groups_on) %>%
  select(-delete) %>%
  na.omit() %>%
  mutate(count = 1, Groups_on = as.factor(Groups_on)) %>%
  group_by(Groups_on, mech) %>%
  summarise(total_count = sum(count)) %>% ungroup() %>%
  na.omit() %>% group_by(Groups_on) %>% mutate(sum = sum(total_count)) %>%
  mutate(perc = (total_count/sum)*100)

top_5pc_mechG4 <- top_5pc_mechG4 %>%
  mutate(mech = as.character(mech)) %>%
  dplyr::mutate(mech = replace(mech, mech == "Overlap", "Niche differentiation")) %>%
  mutate(mech = as.factor(mech))

top_5pc_mechG4$mech <- factor(top_5pc_mechG4$mech, levels = c("Competition", "Niche differentiation",
                                                              "Colonisation", "Growth"))

(PanelS4I <- ggplot(top_5pc_mechG4, aes(x = as.factor(Groups_on), y = perc, fill = mech)) +
    geom_bar(stat = "identity") +
    scale_colour_manual(values = c("Colonisation" = "#ADD1EC", "Growth" = "#EFCAB1",
                                   "Niche differentiation" = "#DFCBDE", "Competition" = "#CC6666")) +
    scale_fill_manual(values = c("Colonisation" = "#ADD1EC", "Growth" = "#EFCAB1",
                                 "Niche differentiation" = "#DFCBDE", "Competition" = "#CC6666")) +
    xlab("Number of mechanisms") + # 750 x 350
    ylab("") + theme_bw() +
    ylab("% in top models") +
    guides(fill=guide_legend("Mechanism"), colour = guide_legend("Mechanism"))+
    theme(legend.position = "right",
          axis.text = element_text(size = 9),
          axis.title = element_text(size=11)) +
    guides(fill=guide_legend("Mechanism", order = 1, ncol = 1))) ## switches that appear in the top 5% of each number of switches


results_mainG4 <- tibble(mod_name = numeric(),
                         Growth = as.numeric(),
                         Colonisation = as.numeric(),
                         Competition = as.numeric(),
                         Overlap = as.numeric())

for(i in 1:4) {

  dfZZG4 <- dfZG4 %>% filter(Groups_on == i) %>% dplyr::filter(groups_ch != "None")

  radar_mainG4 <- dfZZG4 %>%
    filter(sp_rich > 1) %>%
    mutate(mae_top5 = quantile(rmse, 0.05)) %>%
    filter(rmse <= mae_top5) ## then finds top 5% performing models


  mech_mainG4 <- radar_mainG4 %>%
    select(groups_ch) %>%
    count(groups_ch) %>% ## how many times does each appear
    separate(groups_ch, into = paste0("sw", 1:4)) %>%
    gather(delete, mech, -n) %>%
    select(-delete) %>%
    na.omit() %>%
    group_by(mech) %>%
    summarise(n = sum(n)) %>%
    mutate(n = n / nrow(radar_mainG4)) %>%
    spread(mech, n) %>%
    mutate(mod_name = i)

  results_mainG4 <- bind_rows(results_mainG4, mech_mainG4)

} ## see how to deal with NAs

results_mainG4$mod_name <- factor(results_mainG4$mod_name, levels = c("1", "2", "3", "4"))
results_mainG4 <- results_mainG4 %>% select(mod_name, Colonisation, Growth, `Niche \ndifferentiation` = Overlap, Comp = Competition)

PanelS4G <- ggradar(results_mainG4,
                    axis.label.size = 6, # Afftects the names of the variables
                    grid.label.size = 4,
                    group.point.size = 3,# Simply the size of the point
                    legend.position = "bottom")  +
  theme_minimal() +
  theme(axis.line = element_blank(),   # Remove axis lines
        axis.text = element_blank(),   # Remove axis text (numbers)
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),   # Remove major grid lines
        panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  #  scale_colour_manual(values = c("1" = "#FEC47D", "2" = "#65A688",
  #                                 "3" = "#B54445", "4" = "#3192C4")) +
  theme(panel.grid=element_line(size=0.1)) +
  labs(color = "Number of \nmechanisms") ## 700 x 400


(PanelS4I_noleg <- ggplot(top_5pc_mechG4, aes(x = as.factor(Groups_on), y = perc, fill = mech)) +
    geom_bar(stat = "identity") +
    scale_colour_manual(values = c("Colonisation" = "#ADD1EC", "Growth" = "#EFCAB1",
                                   "Niche differentiation" = "#DFCBDE", "Competition" = "#CC6666")) +
    scale_fill_manual(values = c("Colonisation" = "#ADD1EC", "Growth" = "#EFCAB1",
                                 "Niche differentiation" = "#DFCBDE", "Competition" = "#CC6666")) +
    xlab("Number of mechanisms") + # 750 x 350
    ylab("") + theme_bw() +
    ylab("% in top models") +
    guides(fill=guide_legend("Mechanism"), colour = guide_legend("Mechanism"))+
    theme(legend.position = "none",
          axis.text = element_text(size = 11),
          axis.title = element_text(size=12)) +
    guides(fill=guide_legend("Mechanism", order = 1, ncol = 1)))


PanelS4A # 850 x 625
PanelS4B # 450 x 300
PanelS4c_noleg # 450 x 300
PanelS4D # 850 x 625
PanelS4E # 450 x 300
PanelS4F_noleg # 450 x 300
PanelS4G # 850 x 625
PanelS4H # 450 x 300
PanelS4I_noleg # 450 x 300



##### Figure S5
results_sp_rich2 <- results_sp_rich %>%
  rename(Binit = binit, lottery = lot, fecun = rep, `B*` = bstar, `R*` = rstar, mortality = mor, RGR = rgr)

results_sp_rich2$mod_name <- factor(results_sp_rich2$mod_name, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"))

plt_sp_rich <- results_sp_rich2 %>%
  replace_na(list(Binit = 0,
                  `B*` = 0,
                  dispersal = 0,
                  height = 0,
                  lottery = 0,
                  mortality = 0,
                  pheno = 0,
                  fecun = 0,
                  RGR = 0,
                  root = 0,
                  `R*` = 0)) %>%
  group_by(sp_rich) %>%
  nest() %>%
  mutate(plt = map2(data, sp_rich, function(x, y) {
    ggradar(x, plot.title = paste0("Richness = ", y),
            group.point.size = 0.7,
            grid.label.size = 4,
            axis.label.size = 4,
            group.line.width = 1.3) +
      theme(plot.title = element_text(size = 15,
                                      face = "bold"))}))

## plt_sp_rich$plt[[1]] + plt_sp_rich$plt[[2]] + plt_sp_rich$plt[[3]] + plt_sp_rich$plt[[4]] + plot_layout(guides = "collect")
ggarrange(plt_sp_rich$plt[[1]], plt_sp_rich$plt[[2]], plt_sp_rich$plt[[3]], plt_sp_rich$plt[[4]], common.legend = TRUE, legend = "right")




#### Figure S6 - radar plots of groups split by species richness

results_sp_richG1 <-
  tibble(mod_name = numeric(),
         Growth = as.numeric(),
         Colonisation = as.numeric(),
         Competition = as.numeric(),
         Overlap = as.numeric())

for(i in 1:4) {

  dfZZG1 <- dfZG1 %>% filter(Groups_on == i) %>% dplyr::filter(groups_ch != "None")

  radar_sp_richG1 <- dfZZG1 %>%
    group_by(sp_rich, plot, subplot, groups_ch) %>%
    summarise(mean_rmse = mean(rmse, na.rm = TRUE)) %>%
    group_by(sp_rich) %>%
    mutate(mae_top5 = quantile(mean_rmse, 0.05)) %>%
    filter(mean_rmse <= mae_top5)

  mod_per_sp_richG1 <- radar_sp_richG1 %>% group_by(sp_rich) %>% summarise(count = n())

  mech_sp_richG1 <- radar_sp_richG1 %>%
    select(groups_ch) %>%
    count(groups_ch) %>%
    separate(groups_ch, into = paste0("sw", 1:4)) %>%
    gather(delete, mech, -n, -sp_rich) %>%
    select(-delete) %>%
    na.omit() %>%
    group_by(mech, sp_rich) %>%
    summarise(n = sum(n)) %>%
    inner_join(mod_per_sp_richG1) %>%
    mutate(n = n / count) %>%
    spread(mech, n) %>%
    mutate(mod_name = i) %>%
    select(-count)  %>%
    replace(is.na(.), 0)

  results_sp_richG1 <- bind_rows(results_sp_richG1, mech_sp_richG1)

}

results_sp_richG1 <- results_sp_richG1 %>% select(mod_name, Colonisation, Growth, `Niche \ndifferentiation` = Overlap, Comp = Competition, sp_rich)
results_sp_richG1$mod_name <- factor(results_sp_richG1$mod_name, levels = c("1", "2", "3", "4", "5", "6",
                                                                            "7", "8", "9", "10", "11"))
plt_sp_richG1 <- results_sp_richG1 %>%
  replace_na(list(binit = 0,
                  bstar = 0,
                  dispersal = 0,
                  height = 0,
                  lot = 0,
                  mor = 0,
                  pheno = 0,
                  rep = 0,
                  rgr = 0,
                  root = 0,
                  rstar = 0)) %>%
  group_by(sp_rich) %>%
  nest() %>%
  mutate(plt = map2(data, sp_rich, function(x, y) {
    ggradar(x, plot.title = paste0("Richness = ", y),
            group.point.size = 0.7,
            grid.label.size = 4,
            axis.label.size = 5,
            group.line.width = 1.3) +
      theme(plot.title = element_text(size = 15,
                                      face = "bold",
      ),
      legend.position = "none"
      )}))

# plt_sp_richG1$plt[[1]] +plt_sp_richG1$plt[[2]] + plt_sp_richG1$plt[[3]]  + plot_layout(guides = "collect")
ggarrange(plt_sp_richG1$plt[[1]], plt_sp_richG1$plt[[2]], plt_sp_richG1$plt[[3]],  ncol  = 3)
### 1250 x 400

#### G2

dfZG2 <- non_zero_G2 %>%
  rename(groups_ch = GroupNm2, Groups_on = group_if2) %>%
  mutate(Groups_on = as.numeric(Groups_on)) %>%
  group_by(species, sp_rich, plot, subplot, groups_ch, Groups_on) %>%
  filter(sp_rich > 1) %>%
  summarise(rmse = mean(rmse)) %>%
  ungroup() %>%
  group_by(sp_rich, plot, subplot, groups_ch, Groups_on) %>%
  summarise(rmse = mean(rmse)) %>%
  ungroup() 

results_sp_richG2 <-
  tibble(mod_name = numeric(),
         Growth = as.numeric(),
         Colonisation = as.numeric(),
         Competition = as.numeric(),
         Overlap = as.numeric())

for(i in 1:4) {

  dfZZG2 <- dfZG2 %>% filter(Groups_on == i) %>% dplyr::filter(groups_ch != "None")

  radar_sp_richG2 <- dfZZG2 %>%
    group_by(sp_rich, plot, subplot, groups_ch) %>%
    summarise(mean_rmse = mean(rmse, na.rm = TRUE)) %>%
    group_by(sp_rich) %>%
    mutate(mae_top5 = quantile(mean_rmse, 0.05)) %>%
    filter(mean_rmse <= mae_top5)

  mod_per_sp_richG2 <- radar_sp_richG2 %>% group_by(sp_rich) %>% summarise(count = n())

  mech_sp_richG2 <- radar_sp_richG2 %>%
    select(groups_ch) %>%
    count(groups_ch) %>%
    separate(groups_ch, into = paste0("sw", 1:4)) %>%
    gather(delete, mech, -n, -sp_rich) %>%
    select(-delete) %>%
    na.omit() %>%
    group_by(mech, sp_rich) %>%
    summarise(n = sum(n)) %>%
    inner_join(mod_per_sp_richG2) %>%
    mutate(n = n / count) %>%
    spread(mech, n) %>%
    mutate(mod_name = i) %>%
    select(-count)  %>%
    replace(is.na(.), 0)

  results_sp_richG2 <- bind_rows(results_sp_richG2, mech_sp_richG2)

}

results_sp_richG2 <- results_sp_richG2 %>% select(mod_name, Colonisation, Growth, `Niche \ndifferentiation` = Overlap, Comp = Competition, sp_rich)
results_sp_richG2$mod_name <- factor(results_sp_richG2$mod_name, levels = c("1", "2", "3", "4", "5", "6",
                                                                            "7", "8", "9", "10", "11"))
plt_sp_richG2 <- results_sp_richG2 %>%
  replace_na(list(binit = 0,
                  bstar = 0,
                  dispersal = 0,
                  height = 0,
                  lot = 0,
                  mor = 0,
                  pheno = 0,
                  rep = 0,
                  rgr = 0,
                  root = 0,
                  rstar = 0)) %>%
  group_by(sp_rich) %>%
  nest() %>%
  mutate(plt = map2(data, sp_rich, function(x, y) {
    ggradar(x, plot.title = paste0("Richness = ", y),
            group.point.size = 0.7,
            grid.label.size = 4,
            axis.label.size = 5,
            group.line.width = 1.3) +
      theme(plot.title = element_text(size = 15,
                                      face = "bold"),
            legend.position = "none"
      )}))

# plt_sp_richG2$plt[[1]] +plt_sp_richG2$plt[[2]] + plt_sp_richG2$plt[[3]]  + plot_layout(guides = "collect")
ggarrange(plt_sp_richG2$plt[[1]], plt_sp_richG2$plt[[2]], plt_sp_richG2$plt[[3]], ncol  = 3)
### 1250 x 400


results_sp_richG3 <-
  tibble(mod_name = numeric(),
         Growth = as.numeric(),
         Colonisation = as.numeric(),
         Competition = as.numeric(),
         Overlap = as.numeric())

for(i in 1:4) {

  dfZZG3 <- dfZG3 %>% filter(Groups_on == i) %>% dplyr::filter(groups_ch != "None")

  radar_sp_richG3 <- dfZZG3 %>%
    group_by(sp_rich, plot, subplot, groups_ch) %>%
    summarise(mean_rmse = mean(rmse, na.rm = TRUE)) %>%
    group_by(sp_rich) %>%
    mutate(mae_top5 = quantile(mean_rmse, 0.05)) %>%
    filter(mean_rmse <= mae_top5)

  mod_per_sp_richG3 <- radar_sp_richG3 %>% group_by(sp_rich) %>% summarise(count = n())

  mech_sp_richG3 <- radar_sp_richG3 %>%
    select(groups_ch) %>%
    count(groups_ch) %>%
    separate(groups_ch, into = paste0("sw", 1:4)) %>%
    gather(delete, mech, -n, -sp_rich) %>%
    select(-delete) %>%
    na.omit() %>%
    group_by(mech, sp_rich) %>%
    summarise(n = sum(n)) %>%
    inner_join(mod_per_sp_richG3) %>%
    mutate(n = n / count) %>%
    spread(mech, n) %>%
    mutate(mod_name = i) %>%
    select(-count)  %>%
    replace(is.na(.), 0)

  results_sp_richG3 <- bind_rows(results_sp_richG3, mech_sp_richG3)

}

results_sp_richG3 <- results_sp_richG3 %>% select(mod_name, Colonisation, Growth, `Niche \ndifferentiation` = Overlap, Comp = Competition, sp_rich)
results_sp_richG3$mod_name <- factor(results_sp_richG3$mod_name, levels = c("1", "2", "3", "4", "5", "6",
                                                                            "7", "8", "9", "10", "11"))
plt_sp_richG3 <- results_sp_richG3 %>%
  replace_na(list(binit = 0,
                  bstar = 0,
                  dispersal = 0,
                  height = 0,
                  lot = 0,
                  mor = 0,
                  pheno = 0,
                  rep = 0,
                  rgr = 0,
                  root = 0,
                  rstar = 0)) %>%
  group_by(sp_rich) %>%
  nest() %>%
  mutate(plt = map2(data, sp_rich, function(x, y) {
    ggradar(x, plot.title = paste0("Richness = ", y),
            group.point.size = 0.7,
            grid.label.size = 4,
            axis.label.size = 5,
            group.line.width = 1.3) +
      theme(plot.title = element_text(size = 15,
                                      face = "bold"),
            legend.position = "none"
      )}))

# plt_sp_richG3$plt[[1]] +plt_sp_richG3$plt[[2]] + plt_sp_richG3$plt[[3]]  + plot_layout(guides = "collect")
ggarrange(plt_sp_richG3$plt[[1]], plt_sp_richG3$plt[[2]], plt_sp_richG3$plt[[3]], ncol  = 3)



results_sp_richG4 <-
  tibble(mod_name = numeric(),
         Growth = as.numeric(),
         Colonisation = as.numeric(),
         Competition = as.numeric(),
         Overlap = as.numeric())

for(i in 1:4) {

  dfZZG4 <- dfZG4 %>% filter(Groups_on == i) %>% dplyr::filter(groups_ch != "None")

  radar_sp_richG4 <- dfZZG4 %>%
    group_by(sp_rich, plot, subplot, groups_ch) %>%
    summarise(mean_rmse = mean(rmse, na.rm = TRUE)) %>%
    group_by(sp_rich) %>%
    mutate(mae_top5 = quantile(mean_rmse, 0.05)) %>%
    filter(mean_rmse <= mae_top5)

  mod_per_sp_richG4 <- radar_sp_richG4 %>% group_by(sp_rich) %>% summarise(count = n())

  mech_sp_richG4 <- radar_sp_richG4 %>%
    select(groups_ch) %>%
    count(groups_ch) %>%
    separate(groups_ch, into = paste0("sw", 1:4)) %>%
    gather(delete, mech, -n, -sp_rich) %>%
    select(-delete) %>%
    na.omit() %>%
    group_by(mech, sp_rich) %>%
    summarise(n = sum(n)) %>%
    inner_join(mod_per_sp_richG4) %>%
    mutate(n = n / count) %>%
    spread(mech, n) %>%
    mutate(mod_name = i) %>%
    select(-count)  %>%
    replace(is.na(.), 0)

  results_sp_richG4 <- bind_rows(results_sp_richG4, mech_sp_richG4)

}

results_sp_richG4 <- results_sp_richG4 %>% select(mod_name, Colonisation, Growth, `Niche \ndifferentiation` = Overlap, Comp = Competition, sp_rich)
results_sp_richG4$mod_name <- factor(results_sp_richG4$mod_name, levels = c("1", "2", "3", "4", "5", "6",
                                                                            "7", "8", "9", "10", "11"))
plt_sp_richG4 <- results_sp_richG4 %>%
  replace_na(list(binit = 0,
                  bstar = 0,
                  dispersal = 0,
                  height = 0,
                  lot = 0,
                  mor = 0,
                  pheno = 0,
                  rep = 0,
                  rgr = 0,
                  root = 0,
                  rstar = 0)) %>%
  group_by(sp_rich) %>%
  nest() %>%
  mutate(plt = map2(data, sp_rich, function(x, y) {
    ggradar(x, plot.title = paste0("Richness = ", y),
            group.point.size = 0.7,
            grid.label.size = 4,
            axis.label.size = 5,
            group.line.width = 1.3) +
      theme(plot.title = element_text(size = 15,
                                      face = "bold"),
            legend.position = "none"
      )}))

# plt_sp_richG4$plt[[1]] +plt_sp_richG4$plt[[2]] + plt_sp_richG4$plt[[3]]  + plot_layout(guides = "collect")
ggarrange(plt_sp_richG4$plt[[1]], plt_sp_richG4$plt[[2]], plt_sp_richG4$plt[[3]], ncol  = 3)



########### Figure S7 - speices x group radar plot

## 1 group per species
results_speciesG1 <-
  tibble(mod_name = numeric(),
         Growth = as.numeric(),
         Colonisation = as.numeric(),
         Competition = as.numeric(),
         Overlap = as.numeric())

for(i in 1:4) {

  dfZZG1_sp <- non_zero_G1 %>%
    rename(groups_ch = GroupNm1, Groups_on = group_if1) %>%
    mutate(Groups_on = as.numeric(Groups_on)) %>%
    group_by(species, sp_rich, plot, subplot, groups_ch, Groups_on) %>%
    filter(sp_rich > 1) %>%
    summarise(mean_rmse = mean(rmse)) %>%
    ungroup() %>% filter(Groups_on == i) %>% dplyr::filter(groups_ch != "None")

  radar_speciesG1 <- dfZZG1_sp %>%
    group_by(species) %>%
    mutate(mae_top5 = quantile(mean_rmse, 0.05)) %>%
    filter(mean_rmse <= mae_top5)

  mod_per_speciesG1 <- radar_speciesG1 %>% group_by(species) %>% summarise(count = n())

  mech_speciesG1 <- radar_speciesG1 %>%
    select(groups_ch) %>%
    count(groups_ch) %>%
    separate(groups_ch, into = paste0("sw", 1:4)) %>%
    gather(delete, mech, -n, -species) %>%
    select(-delete) %>%
    na.omit() %>%
    group_by(mech, species) %>%
    summarise(n = sum(n)) %>%
    inner_join(mod_per_speciesG1) %>%
    mutate(n = n / count) %>%
    spread(mech, n) %>%
    mutate(mod_name = i) %>%
    select(-count)  %>%
    replace(is.na(.), 0)

  results_speciesG1 <- bind_rows(results_speciesG1, mech_speciesG1)

}


results_speciesG1 <- results_speciesG1 %>%
  mutate(species = replace(species, species == "Agropyron repens", "Elymus repens"))

results_speciesG1 <- results_speciesG1 %>% select(mod_name, Colonisation, Growth, `Niche \ndifferentiation` = Overlap, Comp = Competition, species)

results_speciesG1$mod_name <- factor(results_speciesG1$mod_name, levels = c("1", "2", "3", "4", "5", "6",
                                                                            "7", "8", "9", "10", "11"))
plt_speciesG1 <- results_speciesG1 %>%
  replace_na(list(binit = 0,
                  bstar = 0,
                  dispersal = 0,
                  height = 0,
                  lot = 0,
                  mor = 0,
                  pheno = 0,
                  rep = 0,
                  rgr = 0,
                  root = 0,
                  rstar = 0)) %>%
  group_by(species) %>%
  nest() %>%
  mutate(plt = map2(data, species, function(x, y) {
    ggradar(x, plot.title = paste0(y),
            group.point.size = 0.7,
            grid.label.size = 4,
            axis.label.size = 4,
            group.line.width = 1.3) +
      theme(plot.title = element_text(size = 15,
                                      face = "bold"),
            legend.position = "none"
      )}))

#plt_speciesG1$plt[[2]] +plt_speciesG1$plt[[3]] + plt_speciesG1$plt[[1]] +
#  plt_speciesG1$plt[[4]] + plt_speciesG1$plt[[5]]+ plot_layout(guides = "collect")

ggarrange(plt_speciesG1$plt[[2]], plt_speciesG1$plt[[3]], plt_speciesG1$plt[[1]], plt_speciesG1$plt[[4]], plt_speciesG1$plt[[5]], ncol = 5)
### 2000 x 800


### 2 groups per species
results_speciesG2 <-
  tibble(mod_name = numeric(),
         Growth = as.numeric(),
         Colonisation = as.numeric(),
         Competition = as.numeric(),
         Overlap = as.numeric())

for(i in 1:4) {

  dfZZG2_sp <- non_zero_G2 %>%
    rename(groups_ch = GroupNm2, Groups_on = group_if2) %>%
    mutate(Groups_on = as.numeric(Groups_on)) %>%
    group_by(species, sp_rich, plot, subplot, groups_ch, Groups_on) %>%
    filter(sp_rich > 1) %>%
    summarise(mean_rmse = mean(rmse)) %>%
    ungroup() %>% filter(Groups_on == i) %>% dplyr::filter(groups_ch != "None")

  radar_speciesG2 <- dfZZG2_sp %>%
    group_by(species) %>%
    mutate(mae_top5 = quantile(mean_rmse, 0.05)) %>%
    filter(mean_rmse <= mae_top5)

  mod_per_speciesG2 <- radar_speciesG2 %>% group_by(species) %>% summarise(count = n())

  mech_speciesG2 <- radar_speciesG2 %>%
    select(groups_ch) %>%
    count(groups_ch) %>%
    separate(groups_ch, into = paste0("sw", 1:4)) %>%
    gather(delete, mech, -n, -species) %>%
    select(-delete) %>%
    na.omit() %>%
    group_by(mech, species) %>%
    summarise(n = sum(n)) %>%
    inner_join(mod_per_speciesG2) %>%
    mutate(n = n / count) %>%
    spread(mech, n) %>%
    mutate(mod_name = i) %>%
    select(-count)  %>%
    replace(is.na(.), 0)

  results_speciesG2 <- bind_rows(results_speciesG2, mech_speciesG2)

}


results_speciesG2 <- results_speciesG2 %>%
  mutate(species = replace(species, species == "Agropyron repens", "Elymus repens"))

results_speciesG2 <- results_speciesG2 %>% select(mod_name, Colonisation, Growth, `Niche \ndifferentiation` = Overlap, Comp = Competition, species)

results_speciesG2$mod_name <- factor(results_speciesG2$mod_name, levels = c("1", "2", "3", "4", "5", "6",
                                                                            "7", "8", "9", "10", "11"))
plt_speciesG2 <- results_speciesG2 %>%
  replace_na(list(binit = 0,
                  bstar = 0,
                  dispersal = 0,
                  height = 0,
                  lot = 0,
                  mor = 0,
                  pheno = 0,
                  rep = 0,
                  rgr = 0,
                  root = 0,
                  rstar = 0)) %>%
  group_by(species) %>%
  nest() %>%
  mutate(plt = map2(data, species, function(x, y) {
    ggradar(x, plot.title = paste0(y),
            group.point.size = 0.7,
            grid.label.size = 4,
            axis.label.size = 4,
            group.line.width = 1.3) +
      theme(plot.title = element_text(size = 15,
                                      face = "bold"),
            legend.position = "none"
      )}))

#plt_speciesG2$plt[[2]] +plt_speciesG2$plt[[3]] + plt_speciesG2$plt[[1]] +
#  plt_speciesG2$plt[[4]] + plt_speciesG2$plt[[5]]+ plot_layout(guides = "collect")
ggarrange(plt_speciesG2$plt[[2]], plt_speciesG2$plt[[3]], plt_speciesG2$plt[[1]], plt_speciesG2$plt[[4]], plt_speciesG2$plt[[5]], ncol = 5)
### 2000 x 800

#### G3

results_speciesG3 <-
  tibble(mod_name = numeric(),
         Growth = as.numeric(),
         Colonisation = as.numeric(),
         Competition = as.numeric(),
         Overlap = as.numeric())

for(i in 1:4) {

  dfZZG3_sp <- non_zero_G3 %>%
    rename(groups_ch = GroupNm3, Groups_on = group_if3) %>%
    mutate(Groups_on = as.numeric(Groups_on)) %>%
    group_by(species, sp_rich, plot, subplot, groups_ch, Groups_on) %>%
    filter(sp_rich > 1) %>%
    summarise(mean_rmse = mean(rmse)) %>%
    ungroup() %>% filter(Groups_on == i) %>% dplyr::filter(groups_ch != "None")

  radar_speciesG3 <- dfZZG3_sp %>%
    group_by(species) %>%
    mutate(mae_top5 = quantile(mean_rmse, 0.05)) %>%
    filter(mean_rmse <= mae_top5)

  mod_per_speciesG3 <- radar_speciesG3 %>% group_by(species) %>% summarise(count = n())

  mech_speciesG3 <- radar_speciesG3 %>%
    select(groups_ch) %>%
    count(groups_ch) %>%
    separate(groups_ch, into = paste0("sw", 1:4)) %>%
    gather(delete, mech, -n, -species) %>%
    select(-delete) %>%
    na.omit() %>%
    group_by(mech, species) %>%
    summarise(n = sum(n)) %>%
    inner_join(mod_per_speciesG3) %>%
    mutate(n = n / count) %>%
    spread(mech, n) %>%
    mutate(mod_name = i) %>%
    select(-count)  %>%
    replace(is.na(.), 0)

  results_speciesG3 <- bind_rows(results_speciesG3, mech_speciesG3)

}


results_speciesG3 <- results_speciesG3 %>%
  mutate(species = replace(species, species == "Agropyron repens", "Elymus repens"))

results_speciesG3 <- results_speciesG3 %>% select(mod_name, Colonisation, Growth, `Niche \ndifferentiation` = Overlap, Comp = Competition, species)

results_speciesG3$mod_name <- factor(results_speciesG3$mod_name, levels = c("1", "2", "3", "4", "5", "6",
                                                                            "7", "8", "9", "10", "11"))
plt_speciesG3 <- results_speciesG3 %>%
  replace_na(list(binit = 0,
                  bstar = 0,
                  dispersal = 0,
                  height = 0,
                  lot = 0,
                  mor = 0,
                  pheno = 0,
                  rep = 0,
                  rgr = 0,
                  root = 0,
                  rstar = 0)) %>%
  group_by(species) %>%
  nest() %>%
  mutate(plt = map2(data, species, function(x, y) {
    ggradar(x, plot.title = paste0(y),
            group.point.size = 0.7,
            grid.label.size = 4,
            axis.label.size = 4,
            group.line.width = 1.3) +
      theme(plot.title = element_text(size = 15,
                                      face = "bold"),
            legend.position = "none"
      )}))

#plt_speciesG3$plt[[2]] +plt_speciesG3$plt[[3]] + plt_speciesG3$plt[[1]] +
#  plt_speciesG3$plt[[4]] + plt_speciesG3$plt[[5]]+ plot_layout(guides = "collect")
ggarrange(plt_speciesG3$plt[[2]], plt_speciesG3$plt[[3]], plt_speciesG3$plt[[1]], plt_speciesG3$plt[[4]], plt_speciesG3$plt[[5]], ncol = 5)
# 2000 x 800


### 1250 x 500

results_speciesG4 <-
  tibble(mod_name = numeric(),
         Growth = as.numeric(),
         Colonisation = as.numeric(),
         Competition = as.numeric(),
         Overlap = as.numeric())

for(i in 1:4) {

  dfZZG4_sp <- non_zero_G4 %>%
    rename(groups_ch = GroupNm4, Groups_on = group_if4) %>%
    mutate(Groups_on = as.numeric(Groups_on)) %>%
    group_by(species, sp_rich, plot, subplot, groups_ch, Groups_on) %>%
    filter(sp_rich > 1) %>%
    summarise(mean_rmse = mean(rmse)) %>%
    ungroup() %>% filter(Groups_on == i) %>% dplyr::filter(groups_ch != "None")

  radar_speciesG4 <- dfZZG4_sp %>%
    group_by(species) %>%
    mutate(mae_top5 = quantile(mean_rmse, 0.05)) %>%
    filter(mean_rmse <= mae_top5)

  mod_per_speciesG4 <- radar_speciesG4 %>% group_by(species) %>% summarise(count = n())

  mech_speciesG4 <- radar_speciesG4 %>%
    select(groups_ch) %>%
    count(groups_ch) %>%
    separate(groups_ch, into = paste0("sw", 1:4)) %>%
    gather(delete, mech, -n, -species) %>%
    select(-delete) %>%
    na.omit() %>%
    group_by(mech, species) %>%
    summarise(n = sum(n)) %>%
    inner_join(mod_per_speciesG4) %>%
    mutate(n = n / count) %>%
    spread(mech, n) %>%
    mutate(mod_name = i) %>%
    select(-count)  %>%
    replace(is.na(.), 0)

  results_speciesG4 <- bind_rows(results_speciesG4, mech_speciesG4)

}


results_speciesG4 <- results_speciesG4 %>%
  mutate(species = replace(species, species == "Agropyron repens", "Elymus repens"))

results_speciesG4 <- results_speciesG4 %>% select(mod_name, Colonisation, Growth, `Niche \ndifferentition` = Overlap, Comp = Competition, species)


results_speciesG4$mod_name <- factor(results_speciesG4$mod_name, levels = c("1", "2", "3", "4", "5", "6",
                                                                            "7", "8", "9", "10", "11"))
plt_speciesG4 <- results_speciesG4 %>%
  replace_na(list(binit = 0,
                  bstar = 0,
                  dispersal = 0,
                  height = 0,
                  lot = 0,
                  mor = 0,
                  pheno = 0,
                  rep = 0,
                  rgr = 0,
                  root = 0,
                  rstar = 0)) %>%
  group_by(species) %>%
  nest() %>%
  mutate(plt = map2(data, species, function(x, y) {
    ggradar(x, plot.title = paste0(y),
            group.point.size = 0.7,
            grid.label.size = 4,
            axis.label.size = 4,
            group.line.width = 1.3) +
      theme(plot.title = element_text(size = 15,
                                      face = "bold"),
            legend.position = "none"
      )}))

#plt_speciesG4$plt[[2]] +plt_speciesG4$plt[[3]] + plt_speciesG4$plt[[1]] +
#  plt_speciesG4$plt[[4]] + plt_speciesG4$plt[[5]]+ plot_layout(guides = "collect")

ggarrange(plt_speciesG4$plt[[2]], plt_speciesG4$plt[[3]], plt_speciesG4$plt[[1]], plt_speciesG4$plt[[4]], plt_speciesG4$plt[[5]], ncol = 5)


######## Figure S8 - radar plots connecting the co-occurance of mechanisms per group

load("Non_zero_replicates_all_year6.rda")
load("files_combined.rda")

non_zero2 <- non_zero_replicates %>% filter(sp_rich > 1)

df <- non_zero2 %>% inner_join(files %>% select(full_name, mod_name), by = c("replicate" = "full_name")) %>%
  mutate(mod_name = as.numeric(mod_name)) %>%
  group_by(species, sp_rich, plot, subplot, switches, mod_name) %>%
  ungroup() %>% ## get rid of replicates
  group_by(sp_rich, plot, subplot, switches, mod_name) %>%
  summarise(rmse = mean(rmse)) %>%
  ungroup() ## gives data frame of output replicates with model names

top_5pc <- df %>%
  group_by(mod_name) %>%
  mutate(pc95 = quantile(rmse, 0.05)) %>%
  filter(rmse <= pc95) 

##### Chord diagram - Fig. S11
switches <- c("binit", "lot", "rep", "dispersal", "bstar", "rstar",  "mor", "rgr", "root", "height", "pheno")

out2 <- combn(switches, 2) %>%
  t() %>%
  as.data.frame() %>%
  mutate(freq = 0)


top5_pc_pairs <- map_df(top_5pc, function(x) x)

top_5pc_pairs <- top5_pc_pairs %>%
  mutate(replicate = 1:n()) %>%
  select(switches, replicate, subplot, plot) %>%
  separate(switches, into = paste0("sw", 1:11)) %>%
  gather(key, value, -replicate, -subplot, -plot) %>%
  na.omit()

for(i in 1:nrow(out2)) {
  df <- top_5pc_pairs %>%
    group_by(replicate, subplot, plot) %>%
    mutate(test = sum(value %in% out2[i,]) == 2) %>%
    filter(test)

  out2[i,3] <- nrow(distinct(df, replicate, subplot, plot))
}


switches <- c("binit", "lot", "rep", "dispersal", "bstar", "rstar",  "mor", "rgr", "root", "height", "pheno")

out2 <- combn(switches, 2) %>%
  t() %>%
  as.data.frame() %>%
  mutate(freq = 0)


top5_pc_pairs <- map_df(top_5pc, function(x) x)

top_5pc_pairs <- top5_pc_pairs %>%
  mutate(replicate = 1:n()) %>%
  select(switches, replicate, subplot, plot) %>%
  separate(switches, into = paste0("sw", 1:11)) %>%
  gather(key, value, -replicate, -subplot, -plot) %>%
  na.omit()

for(i in 1:nrow(out2)) {
  df <- top_5pc_pairs %>%
    group_by(replicate, subplot, plot) %>%
    mutate(test = sum(value %in% out2[i,]) == 2) %>%
    filter(test)
  
  out2[i,3] <- nrow(distinct(df, replicate, subplot, plot))
}


out2 <- mutate(out2, freq = (freq - min(freq))/(max(freq) - min(freq)))
colours <- structure(c("#ADD1EC75", "#6699CC75", "#0072B275", "#00008B50", "#CC6666", "#DF4A4A",
                       "#EFCAB1", "#D55E0075","#DFCBDE","#CC79A775" , "#AA449975"), names = switches)
groups <- structure(c("Colonisation",  "Colonisation", "Colonisation",  "Colonisation",
                      "Competition", "Competition", "Growth", "Growth",
                      "Niche differentiation", "Niche differentiation", "Niche differentiation"), names = switches)
library(circlize)
library(rcartocolor)
chordDiagram(out2,
             group = groups,
             grid.col = colours,
             annotationTrack = c("grid", "axis"),
             preAllocateTracks = list(track.height = mm_h(4),track.margin = c(mm_h(4), 0)),
             transparency = 0.5,   # Set transparency for a lighter background
             annotationTrackHeight = c(0.05, 0.1),
             link.lwd = 2)

circos.track(track.index = 2, panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.8, niceFacing = TRUE)
}, bg.border = NA)

highlight.sector(c("root", "height", "pheno"), track.index = 1, col = "#DFCBDE",
                 text = "Overlap", cex = 1, text.col = "black", niceFacing = TRUE, font = 2)
highlight.sector(c("binit", "dispersal", "lot", "rep"), track.index = 1, col = "#ADD1EC",
                 text = "Colonisation", cex = 0.8, text.col = "black", niceFacing = TRUE, font = 2)
highlight.sector(c("mor", "rgr"), track.index = 1, col = "#EFCAB1",
                 text = "Growth", cex = 0.8, text.col = "black", niceFacing = TRUE, font = 2)
highlight.sector(c("bstar", "rstar"), track.index = 1, col = "#CC666699",
                 text = "Competition", cex = 0.8, text.col = "black", niceFacing = TRUE, font = 2)
### not used in MS

top_5pcG2 <- dfZG2 %>%
  group_by(Groups_on) %>%
  mutate(pc95 = quantile(rmse, 0.05)) %>%
  filter(rmse <= pc95)

groups <- c("Colonisation", "Competition", "Growth", "Overlap")

out2 <- combn(groups, 2) %>%
  t() %>%
  as.data.frame() %>%
  mutate(freq = 0)

top5pc_G2 <- top_5pcG2 %>% filter(Groups_on == 2)

top5_pc_G2_pairs <- map_df(top5pc_G2, function(x) x)

top_5pc_G2_pairs <- top5_pc_G2_pairs %>%
  mutate(replicate = 1:n()) %>%
  select(groups_ch, replicate, subplot, plot) %>%
  separate(groups_ch, into = paste0("sw", 1:4)) %>%
  gather(key, value, -replicate, -subplot, -plot) %>%
  na.omit()

for(i in 1:nrow(out2)) {
  df <- top_5pc_G2_pairs %>%
    group_by(replicate, subplot, plot) %>%
    mutate(test = sum(value %in% out2[i,]) == 2) %>%
    filter(test)

  out2[i,3] <- nrow(distinct(df, replicate, subplot, plot))
}


out2 <- out2 %>% dplyr::mutate(V1 = as.character(V1)) %>%
  dplyr::mutate(V1 = replace(V1, V1 == "Overlap", "Niche differentiation"))%>%
  dplyr::mutate(V2 = as.character(V2)) %>%
  dplyr::mutate(V2 = replace(V2, V2 == "Overlap", "Niche differentiation"))


groups2 <- c("Colonisation", "Competition", "Growth", "Niche differentiation")
#out2 <- mutate(out2, freq = (freq - min(freq))/(max(freq) - min(freq)))
colours <- structure(c("#ADD1EC", "#CC333399", "#EFCAB1", "#DFCBDE"), names = groups2)
groups <- structure(c("Colonisation", "Competition", "Growth", "Niche differentiation"), names = groups2)
library(circlize)
library(rcartocolor)
chordDiagram(out2,
             group = groups,
             grid.col = colours,
             annotationTrack = c("grid", "axis"),
             preAllocateTracks = list(track.height = mm_h(6),track.margin = c(mm_h(4), 0)),
             transparency = 0.5,   # Set transparency for a lighter background
             annotationTrackHeight = c(0.1, 0.15),
             link.lwd = 2)


circos.track(track.index = 2, panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.text(mean(xlim), mean(ylim), sector.index, cex = 1.2, niceFacing = TRUE)
}, bg.border = NA)


#### three groups
dfZG3 <- non_zero_G3 %>%
  rename(groups_ch = GroupNm2, Groups_on = group_if2) %>%
  mutate(Groups_on = as.numeric(Groups_on)) %>%
  group_by(species, sp_rich, plot, subplot, groups_ch, Groups_on) %>%
  filter(sp_rich > 1) %>%
  summarise(rmse = mean(rmse)) %>%
  ungroup() %>%
  group_by(sp_rich, plot, subplot, groups_ch, Groups_on) %>%
  summarise(rmse = mean(rmse)) %>%
  ungroup() ## gives dataframe of output replicates with model names

top_5pcG3 <- dfZG3 %>%
  group_by(Groups_on) %>%
  mutate(pc95 = quantile(rmse, 0.05)) %>%
  filter(rmse <= pc95)

groups <- c("Colonisation", "Competition", "Growth", "Overlap")

out2 <- combn(groups, 2) %>%
  t() %>%
  as.data.frame() %>%
  mutate(freq = 0)

top5pc_G3 <- top_5pcG3 %>% filter(Groups_on == 2)

top5_pc_G3_pairs <- map_df(top5pc_G3, function(x) x)

top_5pc_G3_pairs <- top5_pc_G3_pairs %>%
  mutate(replicate = 1:n()) %>%
  select(groups_ch, replicate, subplot, plot) %>%
  separate(groups_ch, into = paste0("sw", 1:4)) %>%
  gather(key, value, -replicate, -subplot, -plot) %>%
  na.omit()

for(i in 1:nrow(out2)) {
  df <- top_5pc_G3_pairs %>%
    group_by(replicate, subplot, plot) %>%
    mutate(test = sum(value %in% out2[i,]) == 2) %>%
    filter(test)

  out2[i,3] <- nrow(distinct(df, replicate, subplot, plot))
}

out2 <- out2 %>% dplyr::mutate(V1 = as.character(V1)) %>%
  dplyr::mutate(V1 = replace(V1, V1 == "Overlap", "Niche differentiation"))%>%
  dplyr::mutate(V2 = as.character(V2)) %>%
  dplyr::mutate(V2 = replace(V2, V2 == "Overlap", "Niche differentiation"))


groups2 <- c("Colonisation", "Competition", "Growth", "Niche differentiation")
#out2 <- mutate(out2, freq = (freq - min(freq))/(max(freq) - min(freq)))
colours <- structure(c("#ADD1EC", "#CC333399", "#EFCAB1", "#DFCBDE"), names = groups2)
groups <- structure(c("Colonisation", "Competition", "Growth", "Niche differentiation"), names = groups2)


library(circlize)
library(rcartocolor)
chordDiagram(out2,
             group = groups,
             grid.col = colours,
             annotationTrack = c("grid", "axis"),
             preAllocateTracks = list(track.height = mm_h(6),track.margin = c(mm_h(4), 0)),
             transparency = 0.5,   # Set transparency for a lighter background
             annotationTrackHeight = c(0.1, 0.15),
             link.lwd = 2)


circos.track(track.index = 2, panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.text(mean(xlim), mean(ylim), sector.index, cex = 1.2, niceFacing = TRUE)
}, bg.border = NA)



### 4 attributes per group
dfZG4 <- non_zero_G4 %>%
  rename(groups_ch = GroupNm2, Groups_on = group_if2) %>%
  mutate(Groups_on = as.numeric(Groups_on)) %>%
  group_by(species, sp_rich, plot, subplot, groups_ch, Groups_on) %>%
  filter(sp_rich > 1) %>%
  summarise(rmse = mean(rmse)) %>%
  ungroup() %>%
  group_by(sp_rich, plot, subplot, groups_ch, Groups_on) %>%
  summarise(rmse = mean(rmse)) %>%
  ungroup() ## gives dataframe of output replicates with model names

top_5pcG4 <- dfZG4 %>%
  group_by(Groups_on) %>%
  mutate(pc95 = quantile(rmse, 0.05)) %>%
  filter(rmse <= pc95)

groups <- c("Colonisation", "Competition", "Growth", "Overlap")

out2 <- combn(groups, 2) %>%
  t() %>%
  as.data.frame() %>%
  mutate(freq = 0)

top5pc_G4 <- top_5pcG4 %>% filter(Groups_on == 2)

top5_pc_G4_pairs <- map_df(top5pc_G4, function(x) x)

top_5pc_G4_pairs <- top5_pc_G4_pairs %>%
  mutate(replicate = 1:n()) %>%
  select(groups_ch, replicate, subplot, plot) %>%
  separate(groups_ch, into = paste0("sw", 1:4)) %>%
  gather(key, value, -replicate, -subplot, -plot) %>%
  na.omit()

for(i in 1:nrow(out2)) {
  df <- top_5pc_G4_pairs %>%
    group_by(replicate, subplot, plot) %>%
    mutate(test = sum(value %in% out2[i,]) == 2) %>%
    filter(test)

  out2[i,3] <- nrow(distinct(df, replicate, subplot, plot))
}

out2 <- out2 %>% dplyr::mutate(V1 = as.character(V1)) %>%
  dplyr::mutate(V1 = replace(V1, V1 == "Overlap", "Niche differentiation"))%>%
  dplyr::mutate(V2 = as.character(V2)) %>%
  dplyr::mutate(V2 = replace(V2, V2 == "Overlap", "Niche differentiation"))


groups2 <- c("Colonisation", "Competition", "Growth", "Niche differentiation")
#out2 <- mutate(out2, freq = (freq - min(freq))/(max(freq) - min(freq)))
colours <- structure(c("#ADD1EC", "#CC333399", "#EFCAB1", "#DFCBDE"), names = groups2)
groups <- structure(c("Colonisation", "Competition", "Growth", "Niche differentiation"), names = groups2)
library(rcartocolor)
chordDiagram(out2,
             group = groups,
             grid.col = colours,
             annotationTrack = c("grid", "axis"),
             preAllocateTracks = list(track.height = mm_h(6),track.margin = c(mm_h(4), 0)),
             transparency = 0.5,   # Set transparency for a lighter background
             annotationTrackHeight = c(0.1, 0.15),
             link.lwd = 2)


circos.track(track.index = 2, panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.text(mean(xlim), mean(ylim), sector.index, cex = 1.2, niceFacing = TRUE)
}, bg.border = NA)


##### Figure S9: Pairwise radar plots


results_mainG2_code <- results_mainG2_code %>% select(mod_name, Colonisation, Growth, `Niche \ndifferentiation` = Overlap, Comp = Competition,  code, count)

pair_rad1 <- ggradar((results_mainG2_code[results_mainG2_code$code == "Agrre_Agrsc",] %>% select(-code, -count)),
                     axis.label.size = 4, # Afftects the names of the variables
                     grid.label.size = 4,
                     group.point.size = 3,# Simply the size of the point
                     legend.position = "bottom")  +
  theme_minimal() +
  theme(axis.line = element_blank(),   # Remove axis lines
        axis.text = element_blank(),   # Remove axis text (numbers)
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),   # Remove major grid lines
        panel.grid.minor = element_blank()) +
  theme(panel.grid=element_line(size=0.1),
        legend.position = "none") +
  labs(color = "Number of \nmechanisms") +
  ggtitle("Elymus repens and Agrostis scabra")

pair_rad2 <- ggradar((results_mainG2_code[results_mainG2_code$code == "Agrre_Poapr",] %>% select(-code, -count)),
                     axis.label.size = 4, # Afftects the names of the variables
                     grid.label.size = 4,
                     group.point.size = 3,# Simply the size of the point
                     legend.position = "bottom")  +
  theme_minimal() +
  theme(axis.line = element_blank(),   # Remove axis lines
        axis.text = element_blank(),   # Remove axis text (numbers)
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),   # Remove major grid lines
        panel.grid.minor = element_blank()) +
  theme(panel.grid=element_line(size=0.1),
        legend.position = "none") +
  labs(color = "Number of \nmechanisms") +
  ggtitle("Elymus repens and Poa pratensis")

pair_rad3 <- ggradar((results_mainG2_code[results_mainG2_code$code == "Agrre_Schsc",] %>% select(-code, -count)),
                     axis.label.size = 4, # Afftects the names of the variables
                     grid.label.size = 4,
                     group.point.size = 3,# Simply the size of the point
                     legend.position = "bottom")  +
  theme_minimal() +
  theme(axis.line = element_blank(),   # Remove axis lines
        axis.text = element_blank(),   # Remove axis text (numbers)
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),   # Remove major grid lines
        panel.grid.minor = element_blank()) +
  theme(panel.grid=element_line(size=0.1),
        legend.position = "none") +
  labs(color = "Number of \nmechanisms") +
  ggtitle("Elymus repens and Schizachyrium scoparium")

pair_rad4 <- ggradar((results_mainG2_code[results_mainG2_code$code == "Agrsc_Andge",] %>% select(-code, -count)),
                     axis.label.size = 4, # Afftects the names of the variables
                     grid.label.size = 4,
                     group.point.size = 3,# Simply the size of the point
                     legend.position = "bottom")  +
  theme_minimal() +
  theme(axis.line = element_blank(),   # Remove axis lines
        axis.text = element_blank(),   # Remove axis text (numbers)
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),   # Remove major grid lines
        panel.grid.minor = element_blank()) +
  theme(panel.grid=element_line(size=0.1),
        legend.position = "none") +
  labs(color = "Number of \nmechanisms") +
  ggtitle("Agrostis scabra and Andropogon gerardi")

pair_rad5 <- ggradar((results_mainG2_code[results_mainG2_code$code == "Poapr_Schsc",] %>% select(-code, -count)),
                     axis.label.size = 4, # Afftects the names of the variables
                     grid.label.size = 4,
                     group.point.size = 3,# Simply the size of the point
                     legend.position = "bottom")  +
  theme_minimal() +
  theme(axis.line = element_blank(),   # Remove axis lines
        axis.text = element_blank(),   # Remove axis text (numbers)
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),   # Remove major grid lines
        panel.grid.minor = element_blank()) +
  theme(panel.grid=element_line(size=0.1),
        legend.position = "none") +
  labs(color = "Number of \nmechanisms") +
  ggtitle("Poa pratensis and Schizachyrium scoparium")

pair_rad6 <- ggradar((results_mainG2_code[results_mainG2_code$code == "Agrsc_Schsc",] %>% select(-code, -count)),
                     axis.label.size = 4, # Afftects the names of the variables
                     grid.label.size = 4,
                     group.point.size = 3,# Simply the size of the point
                     legend.position = "bottom")  +
  theme_minimal() +
  theme(axis.line = element_blank(),   # Remove axis lines
        axis.text = element_blank(),   # Remove axis text (numbers)
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),   # Remove major grid lines
        panel.grid.minor = element_blank()) +
  theme(panel.grid=element_line(size=0.1),
        legend.position = "none") +
  labs(color = "Number of \nmechanisms") +
  ggtitle("Agrostis scabra and Schizachyrium scoparium")

ggarrange(pair_rad1, pair_rad2, pair_rad3, pair_rad4, pair_rad5, pair_rad6,
          ncol= 3, nrow = 2, align = "hv")




##### Figure S10: Pairwise species errors


sp_names <- unlist(dimnames(biomass_data_total)[3])

seed_plots <- as.data.frame.table(seeding_data_total) %>%
  mutate(subplot = as.numeric(Var1),
         plot = as.numeric(Var3),
         species = factor(Var2, labels = sp_names)) %>%
  select(subplot, plot, species, seed_mass = Freq) %>%
  filter(seed_mass != 0) %>% spread(species, seed_mass) %>%
  unite(plot_code, c("plot", "subplot"))

load(paste0("non_zero_replicates_", 0, "_year6.rda"))

non_zero2 <- non_zero_replicates %>%
  unite(code, c("plot", "subplot"), remove = FALSE) %>%
  filter(!code %in% excl_plots$code) %>%
  select(-code)

plot_codes <- non_zero2 %>% filter(sp_rich == 2) %>%
  ungroup() %>%
  select(plot, subplot) %>% ##105 plots
  distinct() %>% unite(plot_code, c("plot", "subplot"), remove = FALSE) %>%
  left_join(seed_plots, by="plot_code") %>%
  mutate(Agrre = if_else(`Agropyron repens` > 1, "Agrre", NA_character_)) %>%
  mutate(Agrsc = if_else(`Agrostis scabra` > 1, "Agrsc", NA_character_)) %>%
  mutate(Andge = if_else(`Andropogon gerardi` > 1, "Andge", NA_character_)) %>%
  mutate(Poapr = if_else(`Poa pratensis` > 1, "Poapr", NA_character_)) %>%
  mutate(Schsc = if_else(`Schizachyrium scoparium` > 1, "Schsc", NA_character_)) %>%
  rowwise() %>%
  mutate(code = paste(na.omit(c(Agrre, Agrsc, Andge, Poapr, Schsc)), collapse = "_")) %>%
  select(plot_code, code)

plot_codes %>% select(code) %>% distinct() ### 6 different combinations

se <- function(x, num_plots) sd(x, na.rm = TRUE) / sqrt(num_plots)

by_switch_number_code <- tibble(mod_name = numeric(),
                                mean_rmse = numeric(),
                                se_rmse = numeric(),
                                upr = numeric(),
                                lwr = numeric(),
                                code = character())

### i = number of switches
for(i in 0:11) { ## this just uses replicates where species are planted
  
  load(paste0("non_zero_replicates_", i, "_year6.rda"))
  
  non_zero2 <- non_zero_replicates  %>%
    filter(sp_rich == 2) %>%
    unite(plot_code, c("plot", "subplot"), remove = FALSE) %>%
    left_join(plot_codes, by="plot_code")
  
  num_plots <- nrow(distinct(non_zero2, subplot, plot, code))
  
  df <- non_zero2 %>%
    group_by(species, plot, subplot, switches, code) %>%
    summarise(rmse = mean(rmse)) %>%
    ungroup() %>%
    group_by(plot, subplot, switches, code) %>%
    summarise(rmse = mean(rmse)) %>%
    ungroup() %>%
    group_by(code) %>%
    summarise(mean_rmse = mean(rmse, na.rm = TRUE),
              se_rmse = se(rmse, num_plots)) %>%
    mutate(upr = mean_rmse + se_rmse,
           lwr = mean_rmse - se_rmse,
           mod_name = i)
  
  by_switch_number_code <- bind_rows(by_switch_number_code, df)
  
}

(code1 <- ggplot(by_switch_number_code, aes(x = mod_name, y = mean_rmse, ymin = lwr, ymax = upr)) +
    geom_ribbon(fill = "grey70") +
    geom_line(linewidth = 1) +
    ylab("RMSE (±1 SE)") +
    xlab("Number of switches") +
    facet_wrap(~code) + # 500 x 350
    theme_bw()) ## this looks correct: increasing number of mechanisms increases model performance


by_switch_number_species_code <- tibble(species = character(),
                                        mod_name = numeric(),
                                        mean_rmse = numeric(),
                                        se_rmse = numeric(),
                                        upr = numeric(),
                                        lwr = numeric(),
                                        code = character())

### i = number of switches
for(i in 0:11) { ## this just uses replicates where species are planted
  
  load(paste0("non_zero_replicates_", i, "_year6.rda"))
  
  non_zero2 <- non_zero_replicates  %>% filter(sp_rich == 2)
  
  num_plots <- nrow(distinct(non_zero2, species, subplot, plot))
  
  df <- non_zero2 %>%
    group_by(species, sp_rich, plot, subplot, switches) %>%
    summarise(rmse = mean(rmse)) %>%
    ungroup() %>%
    group_by(species) %>%
    summarise(mean_rmse = mean(rmse, na.rm = TRUE),
              se_rmse = se(rmse, num_plots)) %>%
    mutate(upr = mean_rmse + se_rmse,
           lwr = mean_rmse - se_rmse,
           mod_name = i)
  
  by_switch_number_species_code <- bind_rows(by_switch_number_species_code, df)
  
}

(code2 <- ggplot(by_switch_number_species_code, aes(x = mod_name, y = mean_rmse, ymin = lwr, ymax = upr,
                                                    group = species, colour = species, fill = species)) +
    geom_ribbon(alpha = 0.5) +
    geom_line(linewidth = 1) +
    ylab("RMSE (±1 SE)") +
    xlab("Number of switches") +
    facet_wrap(~code, scale = "free_y") + # 500 x 350
    theme_bw())



###### number of groups - reduce to 2

non_zero_G2 <- non_zero_wswitches %>%
  filter(ngrowth < 3 & noverlap < 3 & ncolo < 3 & ncomp < 3 &
           ngrowth != 1 & noverlap != 1 & ncolo != 1 & ncomp != 1)

by_group_numberG2_code <- tibble(mod_name = numeric(),
                                 mean_rmse = numeric(),
                                 se_rmse = numeric(),
                                 upr = numeric(),
                                 lwr = numeric(),
                                 code = character())

for(i in 0:4){
  
  dfb2 <- non_zero_G2 %>% filter(group_if2 == i)%>% filter(sp_rich == 2) %>%
    unite(plot_code, c("plot", "subplot"), remove = FALSE) %>%
    left_join(plot_codes, by="plot_code")
  
  num_plotsb <- nrow(distinct(dfb2, subplot, plot, code))
  
  dfb3 <- dfb2 %>%
    group_by(species, sp_rich, plot, subplot, switches, code) %>%
    summarise(rmse = mean(rmse)) %>%
    ungroup()  %>%
    group_by(sp_rich, plot, subplot, switches, code) %>%
    summarise(rmse = mean(rmse)) %>%
    ungroup()  %>%
    group_by(code) %>%
    summarise(mean_rmse = mean(rmse, na.rm = TRUE),
              se_rmse = se(rmse, num_plotsb)) %>%
    mutate(upr = mean_rmse + se_rmse,
           lwr = mean_rmse - se_rmse,
           mod_name = i)
  
  by_group_numberG2_code <- bind_rows(by_group_numberG2_code, dfb3)
}

(G2A_code <- ggplot(by_group_numberG2_code, aes(x = mod_name, y = mean_rmse, ymin = lwr, ymax = upr)) +
    geom_ribbon(fill = "grey70") +
    geom_line(linewidth = 1) +
    ylab("RMSE (±1 SE)") +
    xlab("Number of groups") +
    facet_wrap(~code) + # 500 x 350
    theme_bw()) ## increasing number of groups continually increases performance

species_by_group_numberG2_code <- tibble(species = character(),
                                         mod_name = numeric(),
                                         mean_rmse = numeric(),
                                         se_rmse = numeric(),
                                         upr = numeric(),
                                         lwr = numeric(),
                                         code = character())

for(i in 0:4) { ## this just uses replicates where species are planted
  
  df2 <- non_zero_G2 %>% filter(group_if2 == i) %>% filter(sp_rich == 2) %>%
    unite(plot_code, c("plot", "subplot"), remove = FALSE) %>%
    left_join(plot_codes, by="plot_code")
  num_plots <- nrow(distinct(df2, subplot, plot))
  
  
  df3 <- df2 %>%
    group_by(species, sp_rich, plot, subplot, switches, code) %>%
    summarise(rmse = mean(rmse)) %>%
    ungroup()  %>%
    group_by(species, code) %>%
    summarise(mean_rmse = mean(rmse, na.rm = TRUE),
              se_rmse = se(rmse, num_plots)) %>%
    mutate(upr = mean_rmse + se_rmse,
           lwr = mean_rmse - se_rmse,
           mod_name = i)
  
  species_by_group_numberG2_code <- bind_rows(species_by_group_numberG2_code, df3)
  
}

species_by_group_numberG2_code <- species_by_group_numberG2_code %>%
  mutate(species = replace(species, species == "Agropyron repens", "Elymus repens"))

(G2B <- ggplot(species_by_group_numberG2_code, aes(x = mod_name, y = mean_rmse, ymin = lwr,
                                                   ymax = upr, group = species, fill = species, colour = species)) +
    geom_ribbon(alpha = 0.5) +
    geom_line(linewidth = 1) +
    ylab("RMSE (±1 SE)") +
    scale_fill_manual(values = c("Agrostis scabra" = "#FEC47D" , "Andropogon gerardi" = "#65A688",
                                 "Elymus repens" = "#B54445", "Poa pratensis" = "#3192C4",
                                 "Schizachyrium scoparium" = "#633184")) +
    scale_colour_manual(values = c("Agrostis scabra" = "#FEC47D" , "Andropogon gerardi" = "#65A688",
                                   "Elymus repens" = "#B54445", "Poa pratensis" = "#3192C4",
                                   "Schizachyrium scoparium" = "#633184")) +
    xlab("Number of mechanisms") + # 500 x 350
    theme_bw() +
    facet_wrap(~code)) ## weird response by Poa - others seem consistent with above

#### putting into rows

by_switch_number_species_code <- by_switch_number_species_code %>%
  mutate(species = replace(species, species == "Agropyron repens", "Elymus repens"))

(FigS10A <- ggplot(by_switch_number_species_code[by_switch_number_species_code$code == "Agrre_Agrsc",],
                   aes(x = mod_name, y = mean_rmse, ymin = lwr,
                       ymax = upr, group = species,
                       fill = species, colour = species)) +
    geom_ribbon(alpha = 0.5) +
    geom_line(linewidth = 1) +
    ylab("RMSE (±1 SE)") +
    scale_fill_manual(values = c("Agrostis scabra" = "#FEC47D" , "Andropogon gerardi" = "#65A688",
                                 "Elymus repens" = "#B54445", "Poa pratensis" = "#3192C4",
                                 "Schizachyrium scoparium" = "#633184")) +
    scale_colour_manual(values = c("Agrostis scabra" = "#FEC47D" , "Andropogon gerardi" = "#65A688",
                                   "Elymus repens" = "#B54445", "Poa pratensis" = "#3192C4",
                                   "Schizachyrium scoparium" = "#633184")) +
    xlab("Number of attributes") + # 500 x 350
    theme_bw() +
    theme(legend.position="none",
          axis.text = element_text(size = 12),
          axis.title = element_text(size=14)))
(FigS10B <- ggplot(species_by_group_numberG2_code[species_by_group_numberG2_code$code == "Agrre_Agrsc",],
                   aes(x = mod_name, y = mean_rmse, ymin = lwr,
                       ymax = upr, group = species,
                       fill = species, colour = species)) +
    geom_ribbon(alpha = 0.5) +
    geom_line(linewidth = 1) +
    ylab("RMSE (±1 SE)") +
    scale_fill_manual(values = c("Agrostis scabra" = "#FEC47D" , "Andropogon gerardi" = "#65A688",
                                 "Elymus repens" = "#B54445", "Poa pratensis" = "#3192C4",
                                 "Schizachyrium scoparium" = "#633184")) +
    scale_colour_manual(values = c("Agrostis scabra" = "#FEC47D" , "Andropogon gerardi" = "#65A688",
                                   "Elymus repens" = "#B54445", "Poa pratensis" = "#3192C4",
                                   "Schizachyrium scoparium" = "#633184")) +
    xlab("Number of mechanisms") + # 500 x 350
    theme_bw() +
    theme(legend.position="none",
          axis.text = element_text(size = 12),
          axis.title = element_text(size=14)))


(FigS10C <- ggplot(by_switch_number_species_code[by_switch_number_species_code$code == "Agrre_Poapr",],
                   aes(x = mod_name, y = mean_rmse, ymin = lwr,
                       ymax = upr, group = species,
                       fill = species, colour = species)) +
    geom_ribbon(alpha = 0.5) +
    geom_line(linewidth = 1) +
    ylab("RMSE (±1 SE)") +
    scale_fill_manual(values = c("Agrostis scabra" = "#FEC47D" , "Andropogon gerardi" = "#65A688",
                                 "Elymus repens" = "#B54445", "Poa pratensis" = "#3192C4",
                                 "Schizachyrium scoparium" = "#633184")) +
    scale_colour_manual(values = c("Agrostis scabra" = "#FEC47D" , "Andropogon gerardi" = "#65A688",
                                   "Elymus repens" = "#B54445", "Poa pratensis" = "#3192C4",
                                   "Schizachyrium scoparium" = "#633184")) +
    xlab("Number of attributes") + # 500 x 350
    theme_bw() +
    theme(legend.position="none",
          axis.text = element_text(size = 12),
          axis.title = element_text(size=14)))

(FigS10D <- ggplot(species_by_group_numberG2_code[species_by_group_numberG2_code$code == "Agrre_Poapr",],
                   aes(x = mod_name, y = mean_rmse, ymin = lwr,
                       ymax = upr, group = species,
                       fill = species, colour = species)) +
    geom_ribbon(alpha = 0.5) +
    geom_line(linewidth = 1) +
    ylab("RMSE (±1 SE)") +
    scale_fill_manual(values = c("Agrostis scabra" = "#FEC47D" , "Andropogon gerardi" = "#65A688",
                                 "Elymus repens" = "#B54445", "Poa pratensis" = "#3192C4",
                                 "Schizachyrium scoparium" = "#633184")) +
    scale_colour_manual(values = c("Agrostis scabra" = "#FEC47D" , "Andropogon gerardi" = "#65A688",
                                   "Elymus repens" = "#B54445", "Poa pratensis" = "#3192C4",
                                   "Schizachyrium scoparium" = "#633184")) +
    xlab("Number of mechanisms") + # 500 x 350
    theme_bw() +
    theme(legend.position="none",
          axis.text = element_text(size = 12),
          axis.title = element_text(size=14)))


(FigS10E <- ggplot(by_switch_number_species_code[by_switch_number_species_code$code == "Agrre_Schsc",],
                   aes(x = mod_name, y = mean_rmse, ymin = lwr,
                       ymax = upr, group = species,
                       fill = species, colour = species)) +
    geom_ribbon(alpha = 0.5) +
    geom_line(linewidth = 1) +
    ylab("RMSE (±1 SE)") +
    scale_fill_manual(values = c("Agrostis scabra" = "#FEC47D" , "Andropogon gerardi" = "#65A688",
                                 "Elymus repens" = "#B54445", "Poa pratensis" = "#3192C4",
                                 "Schizachyrium scoparium" = "#633184")) +
    scale_colour_manual(values = c("Agrostis scabra" = "#FEC47D" , "Andropogon gerardi" = "#65A688",
                                   "Elymus repens" = "#B54445", "Poa pratensis" = "#3192C4",
                                   "Schizachyrium scoparium" = "#633184")) +
    xlab("Number of attributes") + # 500 x 350
    theme_bw() +
    theme(legend.position="none",
          axis.text = element_text(size = 12),
          axis.title = element_text(size=14)))

(FigS10F <- ggplot(species_by_group_numberG2_code[species_by_group_numberG2_code$code == "Agrre_Schsc",],
                   aes(x = mod_name, y = mean_rmse, ymin = lwr,
                       ymax = upr, group = species,
                       fill = species, colour = species)) +
    geom_ribbon(alpha = 0.5) +
    geom_line(linewidth = 1) +
    ylab("RMSE (±1 SE)") +
    scale_fill_manual(values = c("Agrostis scabra" = "#FEC47D" , "Andropogon gerardi" = "#65A688",
                                 "Elymus repens" = "#B54445", "Poa pratensis" = "#3192C4",
                                 "Schizachyrium scoparium" = "#633184")) +
    scale_colour_manual(values = c("Agrostis scabra" = "#FEC47D" , "Andropogon gerardi" = "#65A688",
                                   "Elymus repens" = "#B54445", "Poa pratensis" = "#3192C4",
                                   "Schizachyrium scoparium" = "#633184")) +
    xlab("Number of mechanisms") + # 500 x 350
    theme_bw() +
    theme(legend.position="none",
          axis.text = element_text(size = 12),
          axis.title = element_text(size=14)))


(FigS10G <- ggplot(by_switch_number_species_code[by_switch_number_species_code$code == "Agrsc_Andge",],
                   aes(x = mod_name, y = mean_rmse, ymin = lwr,
                       ymax = upr, group = species,
                       fill = species, colour = species)) +
    geom_ribbon(alpha = 0.5) +
    geom_line(linewidth = 1) +
    ylab("RMSE (±1 SE)") +
    scale_fill_manual(values = c("Agrostis scabra" = "#FEC47D" , "Andropogon gerardi" = "#65A688",
                                 "Elymus repens" = "#B54445", "Poa pratensis" = "#3192C4",
                                 "Schizachyrium scoparium" = "#633184")) +
    scale_colour_manual(values = c("Agrostis scabra" = "#FEC47D" , "Andropogon gerardi" = "#65A688",
                                   "Elymus repens" = "#B54445", "Poa pratensis" = "#3192C4",
                                   "Schizachyrium scoparium" = "#633184")) +
    xlab("Number of attributes") + # 500 x 350
    theme_bw() +
    theme(legend.position="none",
          axis.text = element_text(size = 12),
          axis.title = element_text(size=14)))

(FigS10H <- ggplot(species_by_group_numberG2_code[species_by_group_numberG2_code$code == "Agrsc_Andge",],
                   aes(x = mod_name, y = mean_rmse, ymin = lwr,
                       ymax = upr, group = species,
                       fill = species, colour = species)) +
    geom_ribbon(alpha = 0.5) +
    geom_line(linewidth = 1) +
    ylab("RMSE (±1 SE)") +
    scale_fill_manual(values = c("Agrostis scabra" = "#FEC47D" , "Andropogon gerardi" = "#65A688",
                                 "Elymus repens" = "#B54445", "Poa pratensis" = "#3192C4",
                                 "Schizachyrium scoparium" = "#633184")) +
    scale_colour_manual(values = c("Agrostis scabra" = "#FEC47D" , "Andropogon gerardi" = "#65A688",
                                   "Elymus repens" = "#B54445", "Poa pratensis" = "#3192C4",
                                   "Schizachyrium scoparium" = "#633184")) +
    xlab("Number of mechanisms") + # 500 x 350
    theme_bw() +
    theme(legend.position="none",
          axis.text = element_text(size = 12),
          axis.title = element_text(size=14)))


(FigS10I <- ggplot(by_switch_number_species_code[by_switch_number_species_code$code == "Poapr_Schsc",],
                   aes(x = mod_name, y = mean_rmse, ymin = lwr,
                       ymax = upr, group = species,
                       fill = species, colour = species)) +
    geom_ribbon(alpha = 0.5) +
    geom_line(linewidth = 1) +
    ylab("RMSE (±1 SE)") +
    scale_fill_manual(values = c("Agrostis scabra" = "#FEC47D" , "Andropogon gerardi" = "#65A688",
                                 "Elymus repens" = "#B54445", "Poa pratensis" = "#3192C4",
                                 "Schizachyrium scoparium" = "#633184")) +
    scale_colour_manual(values = c("Agrostis scabra" = "#FEC47D" , "Andropogon gerardi" = "#65A688",
                                   "Elymus repens" = "#B54445", "Poa pratensis" = "#3192C4",
                                   "Schizachyrium scoparium" = "#633184")) +
    xlab("Number of attributes") + # 500 x 350
    theme_bw() +
    theme(legend.position="none",
          axis.text = element_text(size = 12),
          axis.title = element_text(size=14)))

(FigS10J <- ggplot(species_by_group_numberG2_code[species_by_group_numberG2_code$code == "Poapr_Schsc",],
                   aes(x = mod_name, y = mean_rmse, ymin = lwr,
                       ymax = upr, group = species,
                       fill = species, colour = species)) +
    geom_ribbon(alpha = 0.5) +
    geom_line(linewidth = 1) +
    ylab("RMSE (±1 SE)") +
    scale_fill_manual(values = c("Agrostis scabra" = "#FEC47D" , "Andropogon gerardi" = "#65A688",
                                 "Elymus repens" = "#B54445", "Poa pratensis" = "#3192C4",
                                 "Schizachyrium scoparium" = "#633184")) +
    scale_colour_manual(values = c("Agrostis scabra" = "#FEC47D" , "Andropogon gerardi" = "#65A688",
                                   "Elymus repens" = "#B54445", "Poa pratensis" = "#3192C4",
                                   "Schizachyrium scoparium" = "#633184")) +
    xlab("Number of mechanisms") + # 500 x 350
    theme_bw() +
    theme(legend.position="none",
          axis.text = element_text(size = 12),
          axis.title = element_text(size=14)))

(FigS10K <- ggplot(by_switch_number_species_code[by_switch_number_species_code$code == "Agrsc_Schsc",],
                   aes(x = mod_name, y = mean_rmse, ymin = lwr,
                       ymax = upr, group = species,
                       fill = species, colour = species)) +
    geom_ribbon(alpha = 0.5) +
    geom_line(linewidth = 1) +
    ylab("RMSE (±1 SE)") +
    scale_fill_manual(values = c("Agrostis scabra" = "#FEC47D" , "Andropogon gerardi" = "#65A688",
                                 "Elymus repens" = "#B54445", "Poa pratensis" = "#3192C4",
                                 "Schizachyrium scoparium" = "#633184")) +
    scale_colour_manual(values = c("Agrostis scabra" = "#FEC47D" , "Andropogon gerardi" = "#65A688",
                                   "Elymus repens" = "#B54445", "Poa pratensis" = "#3192C4",
                                   "Schizachyrium scoparium" = "#633184")) +
    xlab("Number of attributes") + # 500 x 350
    theme_bw() +
    theme(legend.position="none",
          axis.text = element_text(size = 12),
          axis.title = element_text(size=14)))

(FigS10L <- ggplot(species_by_group_numberG2_code[species_by_group_numberG2_code$code == "Agrsc_Schsc",],
                   aes(x = mod_name, y = mean_rmse, ymin = lwr,
                       ymax = upr, group = species,
                       fill = species, colour = species)) +
    geom_ribbon(alpha = 0.5) +
    geom_line(linewidth = 1) +
    ylab("RMSE (±1 SE)") +
    scale_fill_manual(values = c("Agrostis scabra" = "#FEC47D" , "Andropogon gerardi" = "#65A688",
                                 "Elymus repens" = "#B54445", "Poa pratensis" = "#3192C4",
                                 "Schizachyrium scoparium" = "#633184")) +
    scale_colour_manual(values = c("Agrostis scabra" = "#FEC47D" , "Andropogon gerardi" = "#65A688",
                                   "Elymus repens" = "#B54445", "Poa pratensis" = "#3192C4",
                                   "Schizachyrium scoparium" = "#633184")) +
    xlab("Number of mechanisms") + # 500 x 350
    theme_bw() +
    theme(legend.position="none",
          axis.text = element_text(size = 12),
          axis.title = element_text(size=14)))


### making figure

ggarrange((annotate_figure((ggarrange(FigS10A, FigS10B, ncol = 2, align = "hv", labels = c("A", "B"))),
                           top = text_grob(expression(italic("Elymus repens") ~ "and" ~ italic("Agrostis scabra")), color = "black", size = 16))),
          (annotate_figure((ggarrange(FigS10C, FigS10D, ncol = 2, align = "hv", labels = c("C", "D"))),
                           top = text_grob(expression(italic("Elymus repens") ~ "and" ~ italic("Poa pratensis")), color = "black", size = 16))),
          (annotate_figure((ggarrange(FigS10E, FigS10F, ncol = 2, align = "hv", labels = c("E", "F"))),
                           top = text_grob(expression(italic("Elymus repens") ~ "and" ~ italic("Schizachyrium scoparium")), color = "black", size = 16))),
          (annotate_figure((ggarrange(FigS10G, FigS10H, ncol = 2, align = "hv", labels = c("G", "H"))),
                           top = text_grob(expression(italic("Agrostis scabra") ~ "and" ~ italic("Andropogon gerardi")), color = "black", size = 16))),
          (annotate_figure((ggarrange(FigS10I, FigS10J, ncol = 2, align = "hv", labels = c("I", "J"))),
                           top = text_grob(expression(italic("Poa pratensis") ~ "and" ~ italic("Schizachyrium scoparium")), color = "black", size = 16))),
          (annotate_figure((ggarrange(FigS10K, FigS10L, ncol = 2, align = "hv", labels = c("K", "L"))),
                           top = text_grob(expression(italic("Agrostis scabra") ~ "and" ~ italic("Schizachyrium scoparium")), color = "black", size = 16))),
          ncol = 1, nrow = 6, align = "hv")




##### Figure S11: Model error exluding Poa plots
load("data/total_biomass_data.rda")

sp_names <- unlist(dimnames(biomass_data_total)[3])

PP_plots <- as.data.frame.table(seeding_data_total) %>%
  mutate(subplot = as.numeric(Var1),
         plot = as.numeric(Var3),
         species = factor(Var2, labels = sp_names)) %>%
  select(subplot, plot, species, seed_mass = Freq) %>%
  filter(seed_mass != 0) %>% filter(species == "Poa pratensis") %>%
  unite(plot_code, c("plot", "subplot")) %>% select(plot_code)

se <- function(x, num_plots) sd(x, na.rm = TRUE) / sqrt(num_plots)



by_switch_number <- tibble(mod_name = numeric(),
                           mean_rmse = numeric(),
                           se_rmse = numeric(),
                           upr = numeric(),
                           lwr = numeric())

### i = number of switches
for(i in 0:11) { ## this just uses replicates where species are planted
  
  load(paste0("non_zero_replicates_", i, "_year6.rda"))
  
  non_zero2 <- non_zero_replicates %>%
    filter(sp_rich > 1) %>%
    unite(code, c("plot", "subplot"), remove = FALSE) %>%
    filter(!code %in% PP_plots$plot_code) %>%
    select(-code)
  
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


(PanelS11A <- ggplot(by_switch_number, aes(x = mod_name, y = mean_rmse, ymin = lwr, ymax = upr)) +
    geom_ribbon(fill = "grey70", alpha = 0.5) +
    geom_line(linewidth = 1) +
    ylab("RMSE (±1 SE)") +
    theme_bw() +
    xlab("Number of attributes") + # 500 x 350
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size=12))) ## this looks correct: increasing number of mechanisms increases model performance


species_by_switch_number <- tibble(species = character(),
                                   mod_name = numeric(),
                                   mean_rmse = numeric(),
                                   se_rmse = numeric(),
                                   upr = numeric(),
                                   lwr = numeric())

for(i in 0:11) { ## this just uses replicates where species are planted
  
  load(paste0("non_zero_replicates_", i, "_year6.rda"))
  
  non_zero2 <- non_zero_replicates %>%
    filter(sp_rich > 1) %>%
    unite(code, c("plot", "subplot"), remove = FALSE) %>%
    filter(!code %in% PP_plots$plot_code) %>%
    select(-code)
  
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

(PanelS11B <- ggplot(species_by_switch_number, aes(x = mod_name, y = mean_rmse, ymin = lwr,
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
    theme(legend.position = c(0.7, 0.8),
          legend.key.size = unit(0.7,"line"),
          axis.text = element_text(size = 10),
          axis.title = element_text(size=12)))


richness_by_switch_number <- tibble(sp_rich = character(),
                                    mod_name = numeric(),
                                    mean_rmse = numeric(),
                                    se_rmse = numeric(),
                                    upr = numeric(),
                                    lwr = numeric())

for(i in 0:11) { ## this just uses replicates where species are planted
  
  load(paste0("non_zero_replicates_", i, "_year6.rda"))
  
  non_zero2 <- non_zero_replicates %>%
    unite(code, c("plot", "subplot"), remove = FALSE) %>%
    filter(!code %in% PP_plots$plot_code) %>%
    select(-code)
  
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

(PanelS11C <- ggplot(richness_by_switch_number, aes(x = mod_name, y = mean_rmse, ymin = lwr,
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
          axis.text = element_text(size = 10),
          axis.title = element_text(size=12)))

#### groups


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
  
  dfb2 <- non_zero_G1 %>% filter(group_if1 == i) %>%  filter(sp_rich > 1) %>%
    unite(code, c("plot", "subplot"), remove = FALSE) %>%
    filter(!code %in% PP_plots$plot_code) %>%
    select(-code)
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
  
  dfb2 <- non_zero_G2 %>% filter(group_if2 == i) %>%  filter(sp_rich > 1)  %>%
    unite(code, c("plot", "subplot"), remove = FALSE) %>%
    filter(!code %in% PP_plots$plot_code) %>%
    select(-code)
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
  
  dfb2 <- non_zero_G3 %>% filter(group_if3 == i) %>%  filter(sp_rich > 1)  %>%
    unite(code, c("plot", "subplot"), remove = FALSE) %>%
    filter(!code %in% PP_plots$plot_code) %>%
    select(-code)
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
  
  dfb2 <- non_zero_G4 %>% filter(group_if4 == i) %>%  filter(sp_rich > 1)  %>%
    unite(code, c("plot", "subplot"), remove = FALSE) %>%
    filter(!code %in% PP_plots$plot_code) %>%
    select(-code)
  
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


PanelS11D <- ggplot(comb_avgs, aes(x = mod_name, y = mean_rmse, ymin = lwr, ymax = upr, fill = factors)) +
  geom_ribbon(alpha = 0.3, outline.type = "both") +
  geom_line(linewidth = 1, aes(colour = factors)) +
  scale_fill_manual(values = c("1" = "#FEC47D", "2" = "#65A688",
                               "3" = "#B54445", "4" = "#3192C4")) +
  scale_colour_manual(values = c("1" = "#FEC47D", "2" = "#65A688",
                                 "3" = "#B54445", "4" = "#3192C4")) +
  ylab("RMSE (±1 SE)") +
  xlab("Number of mechanisms") + # 500 x 350
  theme_bw() +
  ylim(0.55, 1.3) +
  guides(fill=guide_legend("Attributes required", order = 1),
         colour=guide_legend("Attributes required", order = 1))  +
  theme(legend.position = c(0.8, 0.8),
        legend.key.size = unit(0.7,"line"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size=12))



#### group plots

species_by_group_numberG2 <- tibble(species = character(),
                                    mod_name = numeric(),
                                    mean_rmse = numeric(),
                                    se_rmse = numeric(),
                                    upr = numeric(),
                                    lwr = numeric())

for(i in 0:4) { ## this just uses replicates where species are planted
  
  df2 <- non_zero_G2 %>% filter(group_if2 == i) %>%  filter(sp_rich > 1)  %>%
    unite(code, c("plot", "subplot"), remove = FALSE) %>%
    filter(!code %in% PP_plots$plot_code) %>%
    select(-code)
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

(PanelS11e <- ggplot(species_by_group_numberG2, aes(x = mod_name, y = mean_rmse, ymin = lwr,
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
    ylim(0.25, 2.2) +
    theme(legend.position = c(0.70, 0.8),
          legend.key.size = unit(0.7,"line"),
          axis.text = element_text(size = 10),
          axis.title = element_text(size=12))) ## weird response by Poa - others seem consistent with above


richness_by_group_numberG2 <- tibble(sp_rich = character(),
                                     mod_name = numeric(),
                                     mean_rmse = numeric(),
                                     se_rmse = numeric(),
                                     upr = numeric(),
                                     lwr = numeric())

for(i in 0:4) {
  
  df2 <- non_zero_G2 %>% filter(group_if2 == i)  %>%
    unite(code, c("plot", "subplot"), remove = FALSE) %>%
    filter(!code %in% PP_plots$plot_code) %>%
    select(-code)
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

(PanelS11f <- ggplot(richness_by_group_numberG2, aes(x = mod_name, y = mean_rmse, ymin = lwr,
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
          axis.text = element_text(size = 10),
          axis.title = element_text(size=12)))

library(ggpubr)
ggarrange(PanelS11A, PanelS11B, PanelS11C, PanelS11D, PanelS11e, PanelS11f,
          labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 3, nrow = 2, align = "hv")


##### One attribute included - not in main text

load(paste0("non_zero_replicates_1_year6.rda"))
non_zero_one_on <- non_zero_replicates %>%
  mutate(species = replace(species, species == "Agropyron repens", "Elymus repens")) %>%
  mutate(sp_rich = as.factor(sp_rich),
         species = as.factor(species))

non_zero_one_on$switches <- factor(non_zero_one_on$switches,
                                   levels = c("binit", "lot", "rep", "dispersal", "bstar", "rstar",
                                              "mor", "rgr", "root", "height", "pheno"))

ggplot() + # 1200 x 750
  geom_boxplot(aes(x=switches, y=rmse, colour = switches),
               data = non_zero_one_on) +
  facet_grid(sp_rich~species, scales = "free") +
  theme_bw() +
  scale_colour_manual(values = c("binit" =  "#ADD1EC", "lot" = "#6699CC", "rep" = "#0072B2", "dispersal" = "darkblue",
                                 "bstar" = "#CC3333", "rstar" = "#BB1036",
                                 "mor" = "#EFCAB1", "rgr" = "#D55E00",
                                 "root" = "#DFCBDE", "height" = "#CC79A7" , "pheno"= "#AA4499")) +
  ylab("Mean RMSE") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())


##### Figure S13 - one attribute not included

load(paste0("non_zero_replicates_10_year6.rda"))
non_zero_ten_on <- non_zero_replicates %>%
  mutate(species = replace(species, species == "Agropyron repens", "Elymus repens")) %>%
  mutate(sp_rich = as.factor(sp_rich),
         species = as.factor(species))

non_zero_ten_on <- non_zero_ten_on %>%
  mutate(switches = replace(switches, switches == "bstar_dispersal_height_lot_mor_pheno_rep_rgr_root_rstar", "binit"))%>%
  mutate(switches = replace(switches, switches == "binit_dispersal_height_lot_mor_pheno_rep_rgr_root_rstar", "bstar"))%>%
  mutate(switches = replace(switches, switches == "binit_bstar_height_lot_mor_pheno_rep_rgr_root_rstar", "dispersal"))%>%
  mutate(switches = replace(switches, switches == "binit_bstar_dispersal_lot_mor_pheno_rep_rgr_root_rstar", "height"))%>%
  mutate(switches = replace(switches, switches == "binit_bstar_dispersal_height_mor_pheno_rep_rgr_root_rstar", "lot"))%>%
  mutate(switches = replace(switches, switches == "binit_bstar_dispersal_height_lot_pheno_rep_rgr_root_rstar", "mor"))%>%
  mutate(switches = replace(switches, switches == "binit_bstar_dispersal_height_lot_mor_rep_rgr_root_rstar", "pheno"))%>%
  mutate(switches = replace(switches, switches == "binit_bstar_dispersal_height_lot_mor_pheno_rgr_root_rstar", "rep"))%>%
  mutate(switches = replace(switches, switches == "binit_bstar_dispersal_height_lot_mor_pheno_rep_root_rstar", "rgr"))%>%
  mutate(switches = replace(switches, switches == "binit_bstar_dispersal_height_lot_mor_pheno_rep_rgr_rstar", "root"))%>%
  mutate(switches = replace(switches, switches == "binit_bstar_dispersal_height_lot_mor_pheno_rep_rgr_root", "rstar"))

non_zero_ten_on$switches <- factor(non_zero_ten_on$switches,
                                   levels = c("binit", "lot", "rep", "dispersal", "bstar", "rstar",
                                              "mor", "rgr", "root", "height", "pheno"))

ggplot() + # 1200 x 750
  geom_boxplot(aes(x=switches, y=rmse, colour = switches),
               data = non_zero_ten_on) +
  facet_grid(sp_rich~species, scales = "free") +
  theme_bw() +
  scale_colour_manual(values = c("binit" =  "#ADD1EC", "lot" = "#6699CC", "rep" = "#0072B2", "dispersal" = "darkblue",
                                 "bstar" = "#CC3333", "rstar" = "#BB1036",
                                 "mor" = "#EFCAB1", "rgr" = "#D55E00",
                                 "root" = "#DFCBDE", "height" = "#CC79A7" , "pheno"= "#AA4499")) +
  ylab("Mean RMSE") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())



######## Figure S14

se <- function(x, num_plots) sd(x, na.rm = TRUE) / sqrt(num_plots)

by_switch_number <- tibble(mod_name = numeric(),
                           mean_rmse = numeric(),
                           se_rmse = numeric(),
                           upr = numeric(),
                           lwr = numeric())


### i = number of switches
for(i in 0:11) { ## this just uses replicates where species are planted
  
  load(paste0("non_zero_replicates_", i, "_year6.rda"))
  
  non_zero2 <- non_zero_replicates %>% filter(sp_rich > 1) 
  
  num_plots <- nrow(distinct(non_zero2, subplot, plot))
  
  df <- non_zero2 %>%
    group_by(species, sp_rich, plot, subplot, switches) %>%
    summarise(rmse = mean(rmse)) %>%
    ungroup() %>%
    group_by(sp_rich, plot, subplot, switches) %>%
    summarise(rmse = mean(rmse)) %>%
    ungroup() %>%
    summarise(mean_rmse = median(rmse, na.rm = TRUE),
              se_rmse = se(rmse, num_plots)) %>%
    mutate(upr = mean_rmse + se_rmse,
           lwr = mean_rmse - se_rmse,
           mod_name = i)
  
  by_switch_number <- bind_rows(by_switch_number, df)
  
}



(PanelS14A <- ggplot(by_switch_number, aes(x = mod_name, y = mean_rmse, ymin = lwr, ymax = upr)) +
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
  
  load(paste0("non_zero_replicates_", i, "_year6.rda"))
  
  non_zero2 <- non_zero_replicates %>% filter(sp_rich > 1) 
  
  num_plots <- nrow(distinct( non_zero2, subplot, plot))
  
  df <-  non_zero2 %>%
    filter(sp_rich > 1) %>%
    group_by(species, sp_rich, plot, subplot, switches) %>%
    summarise(rmse = mean(rmse)) %>%
    ungroup() %>%
    group_by(species) %>%
    summarise(mean_rmse = median(rmse, na.rm = TRUE),
              ncov = n(),
              se_rmse = se(rmse, ncov)) %>%
    mutate(upr = mean_rmse + se_rmse,
           lwr = mean_rmse - se_rmse,
           mod_name = i)
  
  species_by_switch_number <- bind_rows(species_by_switch_number, df)
  
}

species_by_switch_number$species <- gsub("Agropyron repens", "Elymus repens", species_by_switch_number$species)

(PanelS14B <- ggplot(species_by_switch_number, aes(x = mod_name, y = mean_rmse, ymin = lwr,
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
  
  load(paste0("non_zero_replicates_", i, "_year6.rda"))
  
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
    summarise(mean_rmse = median(rmse, na.rm = TRUE),
              ncov = n(),
              se_rmse = se(rmse, ncov)) %>%
    mutate(upr = mean_rmse + se_rmse,
           lwr = mean_rmse - se_rmse,
           mod_name = i)
  
  richness_by_switch_number <- bind_rows(richness_by_switch_number, df)
  
}

(PanelS14C <- ggplot(richness_by_switch_number, aes(x = mod_name, y = mean_rmse, ymin = lwr,
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
    summarise(mean_rmse = median(rmse, na.rm = TRUE),
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
    summarise(mean_rmse = median(rmse, na.rm = TRUE),
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
    summarise(mean_rmse = median(rmse, na.rm = TRUE),
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
    summarise(mean_rmse = median(rmse, na.rm = TRUE),
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

PanelS14d <- ggplot(comb_avgs, aes(x = mod_name, y = mean_rmse, ymin = lwr, ymax = upr, fill = factors)) +
  geom_ribbon(alpha = 0.3, outline.type = "both") +
  geom_line(linewidth = 1, aes(colour = factors)) +
  scale_fill_manual(values = c("1" = "#FEC47D", "2" = "#65A688",
                               "3" = "#B54445", "4" = "#3192C4")) +
  scale_colour_manual(values = c("1" = "#FEC47D", "2" = "#65A688",
                                 "3" = "#B54445", "4" = "#3192C4")) +
  ylab("RMSE (±1 SE)") +
  xlab("Number of mechanisms") + # 500 x 350
  theme_bw() +
  ylim(0.5, 0.9) +
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
    summarise(mean_rmse = median(rmse, na.rm = TRUE),
              se_rmse = se(rmse, num_plots)) %>%
    mutate(upr = mean_rmse + se_rmse,
           lwr = mean_rmse - se_rmse, 
           mod_name = i)
  
  species_by_group_numberG2 <- bind_rows(species_by_group_numberG2, df3)
  
}

species_by_group_numberG2 <- species_by_group_numberG2 %>%
  mutate(species = replace(species, species == "Agropyron repens", "Elymus repens"))

(PanelS14E <- ggplot(species_by_group_numberG2, aes(x = mod_name, y = mean_rmse, ymin = lwr,
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
    ylim(0.2, 1.3) +
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
    summarise(mean_rmse = median(rmse, na.rm = TRUE),
              se_rmse = se(rmse, num_plots)) %>%
    mutate(upr = mean_rmse + se_rmse,
           lwr = mean_rmse - se_rmse, 
           mod_name = i)
  
  richness_by_group_numberG2 <- bind_rows(richness_by_group_numberG2, df3)
  
}

(PanelS14f <- ggplot(richness_by_group_numberG2, aes(x = mod_name, y = mean_rmse, ymin = lwr,
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
    theme(legend.position = c(0.77, 0.8),
          legend.key.size = unit(0.7,"line"),
          legend.title = element_text(size=13),
          legend.text=element_text(size=13),
          axis.text = element_text(size = 22),
          axis.title = element_text(size=22)))


library(ggpubr)
ggarrange(PanelS14A, PanelS14B, PanelS14C, PanelS14d, PanelS14E, PanelS14f,
          labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 3, nrow = 2, align = "hv",
          font.label = list(size = 22)) ## 1800 x 800


###### Figure S15


### Three attributes included

switches <- c("binit", "lot", "rep", "dispersal", "bstar", "rstar",  "mor", "rgr", "root", "height", "pheno")

out2 <- combn(switches, 2) %>%
  t() %>%
  as.data.frame() %>%
  mutate(freq = 0)

top5pc_A3 <- top_5pc %>% filter(mod_name == 3)

top5_pc_A3_pairs <- map_df(top5pc_A3, function(x) x)

top_5pc_A3_pairs <- top5_pc_A3_pairs %>%
  mutate(replicate = 1:n()) %>%
  select(switches, replicate, subplot, plot) %>%
  separate(switches, into = paste0("sw", 1:11)) %>%
  gather(key, value, -replicate, -subplot, -plot) %>%
  na.omit()

for(i in 1:nrow(out2)) {
  df <- top_5pc_A3_pairs %>%
    group_by(replicate, subplot, plot) %>%
    mutate(test = sum(value %in% out2[i,]) == 2) %>%
    filter(test)
  
  out2[i,3] <- nrow(distinct(df, replicate, subplot, plot))
}

#save(out2, file = "chord_diagram_pairs_nomonos.rda")
#load("chord_diagram_pairs_nomonos.rda")

out2 <- out2 %>% dplyr::mutate(V1 = as.character(V1)) %>%
  dplyr::mutate(V1 = replace(V1, V1 == "binit", "Binit"))%>%
  dplyr::mutate(V1 = replace(V1, V1 == "lot", "lottery"))%>%
  dplyr::mutate(V1 = replace(V1, V1 == "rep", "fecun"))%>%
  dplyr::mutate(V1 = replace(V1, V1 == "bstar", "B*"))%>%
  dplyr::mutate(V1 = replace(V1, V1 == "rstar", "R*"))%>%
  dplyr::mutate(V1 = replace(V1, V1 == "rgr", "RGR"))

out2 <- out2 %>% dplyr::mutate(V2 = as.character(V2)) %>%
  dplyr::mutate(V2 = replace(V2, V2 == "binit", "Binit"))%>%
  dplyr::mutate(V2 = replace(V2, V2 == "lot", "lottery"))%>%
  dplyr::mutate(V2 = replace(V2, V2 == "rep", "fecun"))%>%
  dplyr::mutate(V2 = replace(V2, V2 == "bstar", "B*"))%>%
  dplyr::mutate(V2 = replace(V2, V2 == "rstar", "R*"))%>%
  dplyr::mutate(V2 = replace(V2, V2 == "rgr", "RGR"))

switches2 <- c("Binit", "lottery", "fecun", "dispersal", "B*", "R*",  "mor", "RGR", "root", "height", "pheno")

# out2 <- mutate(out2, freq = (freq - min(freq))/(max(freq) - min(freq)))
colours <- structure(c("#ADD1EC", "#6699CC", "#0072B2", "#608B8B", "#CC6666", "#DF4A6B",
                       "#EFCAB1", "#D55E00","#DFCBDE","#CC79A7" , "#AA4499"), names = switches2)
groups <- structure(c("Colonisation",  "Colonisation", "Colonisation",  "Colonisation",
                      "Competition", "Competition", "Growth", "Growth",
                      "Niche differentiation", "Niche differentiation", "Niche differentiation"), names = switches2)

library(circlize)
library(rcartocolor)
chordDiagram(out2,
             group = groups,
             grid.col = colours,
             annotationTrack = c("grid", "axis"),
             preAllocateTracks = list(track.height = mm_h(6),track.margin = c(mm_h(4), 0)),
             transparency = 0.5,   # Set transparency for a lighter background
             annotationTrackHeight = c(0.1, 0.15),
             link.lwd = 2)


circos.track(track.index = 2, panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.text(mean(xlim), mean(ylim), sector.index, cex = 1.2, niceFacing = TRUE)
}, bg.border = NA)

highlight.sector(c("root", "height", "pheno"), track.index = 1, col = "#DFCBDE",
                 text = "Niche differentiation", cex = 1.2, text.col = "black", niceFacing = TRUE, font = 2, text.vjust = "-1mm")
highlight.sector(c("Binit", "dispersal", "lottery", "fecun"), track.index = 1, col = "#ADD1EC",
                 text = "Colonisation", cex = 1.2, text.col = "black", niceFacing = TRUE, font = 2, text.vjust = "-1mm")
highlight.sector(c("mor", "RGR"), track.index = 1, col = "#EFCAB1",
                 text = "Growth", cex = 1.2, text.col = "black", niceFacing = TRUE, font = 2)
highlight.sector(c("B*", "R*"), track.index = 1, col = "#CC666699",
                 text = "Competition", cex = 1.2, text.col = "black", niceFacing = TRUE, font = 2, text.vjust = "-1mm")

##### 4 attribute models

switches <- c("binit", "lot", "rep", "dispersal", "bstar", "rstar",  "mor", "rgr", "root", "height", "pheno")

out2 <- combn(switches, 2) %>%
  t() %>%
  as.data.frame() %>%
  mutate(freq = 0)

top5pc_A4 <- top_5pc %>% filter(mod_name == 4)

top5_pc_A4_pairs <- map_df(top5pc_A4, function(x) x)

top_5pc_A4_pairs <- top5_pc_A4_pairs %>%
  mutate(replicate = 1:n()) %>%
  select(switches, replicate, subplot, plot) %>%
  separate(switches, into = paste0("sw", 1:11)) %>%
  gather(key, value, -replicate, -subplot, -plot) %>%
  na.omit()

for(i in 1:nrow(out2)) {
  df <- top_5pc_A4_pairs %>%
    group_by(replicate, subplot, plot) %>%
    mutate(test = sum(value %in% out2[i,]) == 2) %>%
    filter(test)
  
  out2[i,3] <- nrow(distinct(df, replicate, subplot, plot))
}

#save(out2, file = "chord_diagram_pairs_nomonos.rda")
#load("chord_diagram_pairs_nomonos.rda")


out2 <- out2 %>% dplyr::mutate(V1 = as.character(V1)) %>%
  dplyr::mutate(V1 = replace(V1, V1 == "binit", "Binit"))%>%
  dplyr::mutate(V1 = replace(V1, V1 == "lot", "lottery"))%>%
  dplyr::mutate(V1 = replace(V1, V1 == "rep", "fecun"))%>%
  dplyr::mutate(V1 = replace(V1, V1 == "bstar", "B*"))%>%
  dplyr::mutate(V1 = replace(V1, V1 == "rstar", "R*"))%>%
  dplyr::mutate(V1 = replace(V1, V1 == "rgr", "RGR"))

out2 <- out2 %>% dplyr::mutate(V2 = as.character(V2)) %>%
  dplyr::mutate(V2 = replace(V2, V2 == "binit", "Binit"))%>%
  dplyr::mutate(V2 = replace(V2, V2 == "lot", "lottery"))%>%
  dplyr::mutate(V2 = replace(V2, V2 == "rep", "fecun"))%>%
  dplyr::mutate(V2 = replace(V2, V2 == "bstar", "B*"))%>%
  dplyr::mutate(V2 = replace(V2, V2 == "rstar", "R*"))%>%
  dplyr::mutate(V2 = replace(V2, V2 == "rgr", "RGR"))

switches2 <- c("Binit", "lottery", "fecun", "dispersal", "B*", "R*",  "mor", "RGR", "root", "height", "pheno")

# out2 <- mutate(out2, freq = (freq - min(freq))/(max(freq) - min(freq)))
colours <- structure(c("#ADD1EC", "#6699CC", "#0072B2", "#608B8B", "#CC6666", "#DF4A6B",
                       "#EFCAB1", "#D55E00","#DFCBDE","#CC79A7" , "#AA4499"), names = switches2)
groups <- structure(c("Colonisation",  "Colonisation", "Colonisation",  "Colonisation",
                      "Competition", "Competition", "Growth", "Growth",
                      "Niche differentiation", "Niche differentiation", "Niche differentiation"), names = switches2)

library(circlize)
library(rcartocolor)
chordDiagram(out2,
             group = groups,
             grid.col = colours,
             annotationTrack = c("grid", "axis"),
             preAllocateTracks = list(track.height = mm_h(6),track.margin = c(mm_h(4), 0)),
             transparency = 0.5,   # Set transparency for a lighter background
             annotationTrackHeight = c(0.1, 0.15),
             link.lwd = 2)


circos.track(track.index = 2, panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.text(mean(xlim), mean(ylim), sector.index, cex = 1.2, niceFacing = TRUE)
}, bg.border = NA)

highlight.sector(c("root", "height", "pheno"), track.index = 1, col = "#DFCBDE",
                 text = "Niche differentiation", cex = 1.2, text.col = "black", niceFacing = TRUE, font = 2, text.vjust = "-1mm")
highlight.sector(c("Binit", "dispersal", "lottery", "fecun"), track.index = 1, col = "#ADD1EC",
                 text = "Colonisation", cex = 1.2, text.col = "black", niceFacing = TRUE, font = 2, text.vjust = "-1mm")
highlight.sector(c("mor", "RGR"), track.index = 1, col = "#EFCAB1",
                 text = "Growth", cex = 1.2, text.col = "black", niceFacing = TRUE, font = 2)
highlight.sector(c("B*", "R*"), track.index = 1, col = "#CC666699",
                 text = "Competition", cex = 1.2, text.col = "black", niceFacing = TRUE, font = 2, text.vjust = "-1mm")



###### Extra figures

##### Top 5% increasing model complexity
by_switch_number_all <- tibble(mod_name = numeric(),
                               rmse = numeric())

### i = number of switches
for(i in 0:11) { ## this just uses replicates where species are planted
  
  load(paste0("non_zero_replicates_", i, "_year6.rda"))
  
  non_zero2 <- non_zero_replicates %>%
    filter(sp_rich > 1) %>%
    unite(code, c("plot", "subplot"), remove = FALSE) %>%
    select(-code)
  
  num_plots <- nrow(distinct(non_zero2, subplot, plot))
  
  df <- non_zero2 %>%
    group_by(species, sp_rich, plot, subplot, switches) %>%
    summarise(rmse = mean(rmse)) %>%
    ungroup() %>%
    group_by(sp_rich, plot, subplot, switches) %>%
    summarise(rmse = mean(rmse)) %>%
    ungroup() %>%
    group_by(sp_rich, plot, subplot) %>%
    summarise(rmse = mean(rmse)) %>%
    mutate(mod_name = i) %>% select(mod_name, rmse)
  
  by_switch_number_all  <- bind_rows(by_switch_number_all, df)
  
}


ggplot() + geom_histogram(aes(x=rmse), data = by_switch_number_all) +
  facet_wrap(~mod_name, scale = "free_y") + theme_bw() +
  geom_vline(data = by_switch_number, aes(xintercept = mean_rmse), color = "red")


by_switch_number_t5 <- tibble(mod_name = numeric(),
                              rmse = numeric())

### i = number of switches
for(i in 0:11) { ## this just uses replicates where species are planted
  
  load(paste0("non_zero_replicates_", i, "_year6.rda"))
  
  non_zero2 <- non_zero_replicates %>%
    filter(sp_rich > 1)
  
  num_plots <- nrow(distinct(non_zero2, subplot, plot))
  
  df <- non_zero2 %>%
    group_by(species, sp_rich, plot, subplot, switches) %>%
    summarise(rmse = mean(rmse)) %>%
    ungroup() %>%
    group_by(sp_rich, plot, subplot, switches) %>%
    summarise(rmse = mean(rmse)) %>%
    ungroup() %>%
    group_by(sp_rich, plot, subplot) %>%
    summarise(rmse = mean(rmse)) %>%
    mutate(pc95 = quantile(rmse, 0.05)) %>%
    filter(rmse <= pc95) %>%
    ungroup() %>%
    summarise(mean_rmse = mean(rmse, na.rm = TRUE),
              se_rmse = se(rmse, num_plots)) %>%
    mutate(upr = mean_rmse + se_rmse,
           lwr = mean_rmse - se_rmse,
           mod_name = i)
  
  by_switch_number_t5  <- bind_rows(by_switch_number_t5, df)
  
}


(Panel2A_tf <- ggplot(by_switch_number_t5, aes(x = mod_name, y = mean_rmse, ymin = lwr, ymax = upr)) +
    geom_ribbon(fill = "grey70", alpha = 0.5) +
    geom_line(linewidth = 1) +
    ylab("RMSE (±1 SE)") +
    theme_bw() +
    xlab("Number of attributes") + # 500 x 350
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size=12))) ## this looks correct: increasing number of mechanisms increases model performance



###### radar plot - full attributes


results_main <- results_main %>% select(mod_name, Binit = binit, lottery = lot, fecundity = rep,
                                        dispersal, `B*` = bstar, `R*` = rstar, mortality = mor, RGR = rgr,
                                        root, height, phenology = pheno)

results_main$mod_name <- factor(results_main$mod_name, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"))

(Panel4a <- ggradar(results_main,
                    axis.label.size = 4, # Afftects the names of the variables
                    grid.label.size = 4,
                    group.point.size = 3,# Simply the size of the point
                    legend.position = "bottom")  +
    theme_minimal() +
    theme(axis.line = element_blank(),   # Remove axis lines
          axis.text = element_blank(),   # Remove axis text (numbers)
          axis.ticks = element_blank(),
          panel.grid.major = element_blank(),   # Remove major grid lines
          panel.grid.minor = element_blank()) +
    theme(legend.title = element_text(size = 12, color = "black"),
          legend.position = "right") +
    theme(panel.grid=element_line(size=0.1)) +
    labs(color = "Number of \nattributes"))



## two groups
dfZG2 <- non_zero_G2 %>%
  rename(groups_ch = GroupNm2, Groups_on = group_if2) %>%
  mutate(Groups_on = as.numeric(Groups_on)) %>%
  group_by(species, sp_rich, plot, subplot, groups_ch, Groups_on) %>%
  filter(sp_rich > 1) %>%
  summarise(rmse = mean(rmse)) %>%
  ungroup() %>%
  group_by(sp_rich, plot, subplot, groups_ch, Groups_on) %>%
  summarise(rmse = mean(rmse)) %>%
  ungroup() ## gives dataframe of output replicates with model names


results_mainG2 <- tibble(mod_name = numeric(),
                         Growth = as.numeric(),
                         Colonisation = as.numeric(),
                         Competition = as.numeric(),
                         Overlap = as.numeric())

for(i in 1:4) {
  
  dfZZG2 <- dfZG2 %>% filter(Groups_on == i) %>% dplyr::filter(groups_ch != "None")
  
  radar_mainG2 <- dfZZG2 %>%
    filter(sp_rich > 1) %>%
    mutate(mae_top5 = quantile(rmse, 0.05)) %>%
    filter(rmse <= mae_top5) ## then finds top 5% performing models
  
  
  mech_mainG2 <- radar_mainG2 %>%
    select(groups_ch) %>%
    count(groups_ch) %>% ## how many times does each appear
    separate(groups_ch, into = paste0("sw", 1:4)) %>%
    gather(delete, mech, -n) %>%
    select(-delete) %>%
    na.omit() %>%
    group_by(mech) %>%
    summarise(n = sum(n)) %>%
    mutate(n = n / nrow(radar_mainG2)) %>%
    spread(mech, n) %>%
    mutate(mod_name = i)
  
  results_mainG2 <- bind_rows(results_mainG2, mech_mainG2)
  
} ## see how to deal with NAs

results_mainG2$mod_name <- factor(results_mainG2$mod_name, levels = c("1", "2", "3", "4"))
results_mainG2 <- results_mainG2 %>% select(mod_name, Colonisation, Growth, Competition, Overlap)

Panel4d <- ggradar(results_mainG2,
                   axis.label.size = 4, # Afftects the names of the variables
                   grid.label.size = 4,
                   group.point.size = 3,# Simply the size of the point
                   legend.position = "bottom")  +
  theme_minimal() +
  theme(axis.line = element_blank(),   # Remove axis lines
        axis.text = element_blank(),   # Remove axis text (numbers)
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),   # Remove major grid lines
        panel.grid.minor = element_blank()) +
  theme(legend.title = element_text(size = 12, color = "black"),
        legend.position = "right") +
  #  scale_colour_manual(values = c("1" = "#FEC47D", "2" = "#65A688",
  #                                 "3" = "#B54445", "4" = "#3192C4")) +
  theme(panel.grid=element_line(size=0.1)) +
  labs(color = "Number of \nmechanisms")


(Panel4f_noleg <- ggplot(top_5pc_mechG2, aes(x = as.factor(Groups_on), y = perc, fill = mech)) +
    geom_bar(stat = "identity") +
    scale_colour_manual(values = c("Colonisation" = "#ADD1EC", "Growth" = "#EFCAB1",
                                   "Niche differentiation" = "#DFCBDE", "Competition" = "#CC6666")) +
    scale_fill_manual(values = c("Colonisation" = "#ADD1EC", "Growth" = "#EFCAB1",
                                 "Niche differentiation" = "#DFCBDE", "Competition" = "#CC6666")) +
    xlab("Number of mechanisms") + # 750 x 350
    ylab("") + theme_bw() +
    ylab("% in top models") +
    guides(fill=guide_legend("Mechanism"), colour = guide_legend("Mechanism"))+
    theme(legend.position = "none",
          plot.tag = element_text(face = "bold", size = 24),
          axis.text = element_text(size = 22),
          axis.title = element_text(size=24)) +
    guides(fill=guide_legend("Mechanism", order = 1, ncol = 1)) +
    labs(tag = "E"))





