##### Figure 4: Catford et al., Mechanistic model

rm(list=ls())

setwd("~/Predicting_e026_coexistence")

library(tidyverse)
library(rcartocolor)
library(circlize)

#### Read in all files together
data_list <- list()

for (i in 0:11) {
  file_name <- paste0("Data/non_zero_replicates_", i, "_year6.rda")
  load(file_name)
  data_list[[i + 1]] <- non_zero_replicates  # Make sure `your_object` is the correct variable in each .rda file
}

non_zero_replicates <- do.call(rbind, data_list)
remove(data_list)


load("Data/files_combined.rda") ##links files to switches 


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
  filter(rmse <= pc95) ### top 5% quantile of RMSE: i.e. best performing models, from each amount of switches

#### figure 4a


switches <- c("binit", "lot", "rep", "dispersal", "bstar", "rstar",  "mor", "rgr", "root", "height", "pheno")

out2 <- combn(switches, 2) %>%
  t() %>%
  as.data.frame() %>%
  mutate(freq = 0)

top5pc_A2 <- top_5pc %>% filter(mod_name == 2)

top5_pc_A2_pairs <- map_df(top5pc_A2, function(x) x)

top_5pc_A2_pairs <- top5_pc_A2_pairs %>%
  mutate(replicate = 1:n()) %>%
  select(switches, replicate, subplot, plot) %>%
  separate(switches, into = paste0("sw", 1:11)) %>%
  gather(key, value, -replicate, -subplot, -plot) %>%
  na.omit()

for(i in 1:nrow(out2)) {
  df <- top_5pc_A2_pairs %>%
    group_by(replicate, subplot, plot) %>%
    mutate(test = sum(value %in% out2[i,]) == 2) %>%
    filter(test)
  
  out2[i,3] <- nrow(distinct(df, replicate, subplot, plot))
}


##rename
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

colours <- structure(c("#ADD1EC", "#6699CC", "#0072B2", "#608B8B", "#CC6666", "#DF4A6B",
                       "#EFCAB1", "#D55E00","#DFCBDE","#CC79A7" , "#AA4499"), names = switches2)
groups <- structure(c("Colonisation",  "Colonisation", "Colonisation",  "Colonisation", 
                      "Competition", "Competition", "Growth", "Growth",
                      "Niche differentiation", "Niche differentiation", "Niche differentiation"), names = switches2)

sum(out2$freq) ## total number


###### Figure 4a
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
### 1100 x 700
### Table of values

out_table <- out2 %>% rename("Attribute_1" = "V1", "Attribute_2" = "V2") ## Table of pairwise values


#### making remaining panels
switches <- c("bstar", "root", "rstar", "height", "binit", "pheno", "lot", "mor", "rep", "dispersal", "rgr") 
group <- c("Competition", "Overlap", "Competition", "Overlap", "Colonisation",
           "Overlap", "Colonisation", "Growth", "Colonisation", "Colonisation", "Growth")
switch_cats <- as.data.frame(cbind(switches, group))

top_5pc_sw <- top_5pc %>%
  filter(mod_name == 1) %>% ##filters out single switch models
  left_join(switch_cats, by="switches") %>%
  count(switches) %>%
  mutate(switches = as.factor(switches),
         switches = fct_reorder(switches, n, .desc = TRUE)) %>%
  mutate(total_n = sum(n)) %>%
  mutate(perc = (n/total_n)*100) %>%
  arrange(-perc) 

### Rename switches with desired names
top_5pc_sw <- top_5pc_sw %>% dplyr::mutate(switches = as.character(switches)) %>%
  dplyr::mutate(switches = replace(switches, switches == "binit", "Binit"))%>%
  dplyr::mutate(switches = replace(switches, switches == "lot", "lottery"))%>%
  dplyr::mutate(switches = replace(switches, switches == "rep", "fecun"))%>%
  dplyr::mutate(switches = replace(switches, switches == "bstar", "B*"))%>%
  dplyr::mutate(switches = replace(switches, switches == "rstar", "R*"))%>%
  dplyr::mutate(switches = replace(switches, switches == "mor", "mortality"))%>%
  dplyr::mutate(switches = replace(switches, switches == "rgr", "RGR"))

top_5pc_sw$switches <- factor(top_5pc_sw$switches, levels = c("Binit", "lottery", "fecun", "dispersal", "B*", "R*", 
                                                          "mortality", "RGR", "root", "height", "pheno"))

(Panel4b <- ggplot(top_5pc_sw, aes(x = switches, y = perc, fill = switches)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  xlab("Attribute") +
  ylab("% in top models") +
  theme(axis.text.x = element_text(vjust = 0.5))  +
  scale_fill_manual(values = c("Binit" =  "#ADD1EC", "lottery" = "#6699CC", "fecun" = "#0072B2", "dispersal" = "#608B8B",
                               "B*" = "#CC6666", "R*" = "#DF4A6B",
                               "mortality" = "#EFCAB1", "RGR" = "#D55E00",
                               "root" = "#DFCBDE", "height" = "#CC79A7" , "pheno"= "#AA4499")) +
  theme(legend.position = "none",
        axis.text = element_text(size = 26),
        plot.tag = element_text(face = "bold", size = 28),
        axis.title = element_text(size=28)) +
    labs(tag = "B")) # 750 x 300

top_5pc_mech <- top_5pc %>%
  select(switches, mod_name) %>%
  filter(mod_name %in% 1:11) %>% ## gets rid of all switches on (full model)
  separate(switches, into = paste0("sw", 1:11)) %>% ## separates out models
  gather(delete, mech, -mod_name) %>%
  select(-delete) %>%
  mutate(count = 1, mod_name = as.factor(mod_name)) %>%
  group_by(mod_name, mech) %>%
  summarise(total_count = sum(count)) %>% ungroup() %>%
  na.omit() %>% group_by(mod_name) %>% mutate(sum = sum(total_count)) %>%
  mutate(perc = (total_count/sum)*100)

top_5pc_mech <- top_5pc_mech %>% dplyr::mutate(mech = as.character(mech)) %>%
  dplyr::mutate(mech = replace(mech, mech == "binit", "Binit"))%>%
  dplyr::mutate(mech = replace(mech, mech == "lot", "lottery"))%>%
  dplyr::mutate(mech = replace(mech, mech == "rep", "fecun"))%>%
  dplyr::mutate(mech = replace(mech, mech == "bstar", "B*"))%>%
  dplyr::mutate(mech = replace(mech, mech == "rstar", "R*"))%>%
  dplyr::mutate(mech = replace(mech, mech == "mor", "mortality"))%>%
  dplyr::mutate(mech = replace(mech, mech == "rgr", "RGR"))

top_5pc_mech$mech <- factor(top_5pc_mech$mech, levels = c("Binit", "lottery", "fecun", "dispersal", "B*", "R*", 
                                                          "mortality", "RGR", "root", "height", "pheno"))

Panel4c <- ggplot(top_5pc_mech, aes(x = as.factor(mod_name), y = perc, fill = mech)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Binit" =  "#ADD1EC", "lottery" = "#6699CC", "fecun" = "#0072B2", "dispersal" = "#608B8B",
                               "B*" = "#CC6666", "R*" = "#DF4A6B",
                               "mortality" = "#EFCAB1", "RGR" = "#D55E00",
                               "root" = "#DFCBDE", "height" = "#CC79A7" , "pheno"= "#AA4499")) +
  xlab("Model complexity (number of attributes)") + # 750 x 350
  ylab("% in top models") + theme_bw() +
  theme(axis.text = element_text(size = 26),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        plot.tag = element_text(face = "bold", size = 28),
        axis.title = element_text(size=28),
        legend.position = "right") +
  guides(fill=guide_legend("Attribute", order = 1, ncol = 1)) ## switches that appear in the top 5% of each number of switches

### no legend - for publication
Panel4c_noleg <- ggplot(top_5pc_mech, aes(x = as.factor(mod_name), y = perc, fill = mech)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Binit" =  "#ADD1EC", "lottery" = "#6699CC", "fecundity" = "#0072B2", "dispersal" = "#608B8B",
                               "B*" = "#CC6666", "R*" = "#DF4A4A",
                               "mortality" = "#EFCAB1", "RGR" = "#D55E00",
                               "root" = "#DFCBDE", "height" = "#CC79A7" , "pheno"= "#AA4499")) +
  xlab("Model complexity (number of attributes)") + # 750 x 350
  ylab("% in top models") + theme_bw() +
  theme(axis.text = element_text(size = 26),
        plot.tag = element_text(face = "bold", size = 28),
        axis.title = element_text(size=28),
        legend.position = "none") +
  guides(fill=guide_legend("Attribute", order = 1, ncol = 1)) +
           labs(tag = "C")  ## switches that appear in the top 5% of each number of switches


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

switchesNZ <- non_zero2 %>% select(switches) %>% distinct() ### codes are always alphabetically arranged

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

non_zero_wswitches  <- non_zero2 %>%
  left_join(switches_wcats2, by="switches") %>%
  left_join(switches3, by="switches")

#### Only requiring 2 variable for a group to be included

non_zero_G2 <- non_zero_wswitches %>%
  filter(ngrowth < 3 & noverlap < 3 & ncolo < 3 & ncomp < 3 &
           ngrowth != 1 & noverlap != 1 & ncolo != 1 & ncomp != 1)

dfZG2 <- non_zero_G2 %>%
  rename(groups_ch = GroupNm1, Groups_on = group_if1) %>%
  mutate(Groups_on = as.numeric(Groups_on)) %>%
  group_by(species, sp_rich, plot, subplot, groups_ch, Groups_on) %>%
  filter(sp_rich > 1) %>%
  summarise(rmse = mean(rmse)) %>%
  ungroup() %>%
  group_by(sp_rich, plot, subplot, groups_ch, Groups_on) %>%
  summarise(rmse = mean(rmse)) %>%
  ungroup() ## gives dataframe of output replicates with model names

top_5pcG2 <- dfZG2 %>%
  group_by(Groups_on) %>%
  mutate(pc95 = quantile(rmse, 0.05)) %>%
  filter(rmse <= pc95) ### top 5% quantile of RMSE: i.e. best performing models, from each amount of switches

top_5pc_swG2 <- top_5pcG2 %>%
  filter(Groups_on == 1) %>% ##filters out single switch models
  count(groups_ch) %>%
  mutate(totaln = sum(n)) %>%
  mutate(perc = (n/totaln)*100) %>%
  mutate(groups_ch = as.factor(groups_ch),
         groups_ch = fct_reorder(groups_ch, perc, .desc = TRUE)) 

top_5pc_swG2 <- top_5pc_swG2 %>%
  mutate(groups_ch = as.character(groups_ch)) %>%
  dplyr::mutate(groups_ch = replace(groups_ch, groups_ch == "Overlap", "Niche \ndifferentiation")) %>%
  mutate(groups_ch = as.factor(groups_ch))

top_5pc_swG2$groups_ch <- factor(top_5pc_swG2$groups_ch, levels = c("Competition", "Niche \ndifferentiation",
                                                                    "Colonisation", "Growth"))

(Panel4D <- ggplot(top_5pc_swG2, aes(x = groups_ch, y = perc, colour = groups_ch, fill = groups_ch)) +
    geom_bar(stat = "identity") +
    scale_colour_manual(values = c("Colonisation" = "#ADD1EC", "Growth" = "#EFCAB1",
                                   "Niche \ndifferentiation" = "#DFCBDE", "Competition" = "#CC6666")) +
    scale_fill_manual(values = c("Colonisation" = "#ADD1EC", "Growth" = "#EFCAB1",
                                 "Niche \ndifferentiation" = "#DFCBDE", "Competition" = "#CC6666")) +
    theme_bw() +
    theme(legend.position = "none") + 
    theme(axis.text = element_text(size = 26),
          plot.tag = element_text(face = "bold", size = 28),
          axis.title = element_text(size=28)) +
    xlab("Mechanism") + ylab("% in top models") +
    labs(tag = "D")) # 750 x 300

top_5pc_mechG2 <- top_5pcG2 %>%
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

top_5pc_mechG2 <- top_5pc_mechG2 %>%
  mutate(mech = as.character(mech)) %>%
  dplyr::mutate(mech = replace(mech, mech == "Overlap", "Niche differentiation")) %>%
  mutate(mech = as.factor(mech))

top_5pc_mechG2$mech <- factor(top_5pc_mechG2$mech, levels = c("Competition", "Niche differentiation",
                                                                    "Colonisation", "Growth"))

(Panel4E <- ggplot(top_5pc_mechG2, aes(x = as.factor(Groups_on), y = perc, fill = mech)) +
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
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          axis.text = element_text(size = 26),
          plot.tag = element_text(face = "bold", size = 28),
          axis.title = element_text(size=28)) +
    guides(fill=guide_legend("Mechanism", order = 1, ncol = 1))+
  labs(tag = "E")) ## switches that appear in the top 5% of each number of switches


### no legend - for publication
(Panel4E_noleg <- ggplot(top_5pc_mechG2, aes(x = as.factor(Groups_on), y = perc, fill = mech)) +
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
          axis.text = element_text(size = 26),
          plot.tag = element_text(face = "bold", size = 28),
          axis.title = element_text(size=28)) +
    labs(tag = "E"))

library(grid) ## export as 900 x 1100 pixels
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 115, ncol = 100)))
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
}
# Arrange the plots
print(Panel4b, vp = define_region(row = 1:40, col = 1:100))   # Span over two columns
print(Panel4c_noleg, vp = define_region(row = 41:80, col = 1:100))
print(Panel4D, vp = define_region(row = 81:115, col = 2:51))
print(Panel4E_noleg, vp = define_region(row = 81:113, col = 52:100))

