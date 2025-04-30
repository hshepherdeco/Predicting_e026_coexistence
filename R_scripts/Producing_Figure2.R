### Figure 2 - Catford et al., Mechanistic model
### Comparing modelled results to observed values

rm(list=ls())

setwd("~/SimulateCoexistence") ### requires data and functions stored in from the main function 
devtools::load_all()
library(tidyverse); library(rcartocolor)

rstar_switch <- FALSE
height_switch <- FALSE
root_switch <- FALSE
pheno_switch <- FALSE
dispersal_switch <- FALSE
lot_switch <- FALSE
bstar_switch = FALSE
rgr_switch = FALSE
b_init_switch = FALSE
rep_local_diff_switch = FALSE
mor_diff_switch = FALSE
soiln_switch = FALSE

# array for saving data
# year, subplot, species, plot
sim_data_total = array(dim=dim(biomass_data_total)+c(1,0,0,0))

# which plot? (1,2,3,4,6,7,8,9,10)
#plotnum = 1
for(plotnum in c(1,2,3,4,6,7,8,9,10)) {
  
  # load data from plot
  tmp = data.frame(subplot=1:64, totaln=soiln_data_total[,which(unique_list$uplot==plotnum)])
  
  #planted_richness = apply(presence_data_total, c(1,3), sum)
  #uplanted_richness = sort(unique(c(planted_richness)))
  #siteuse = which(planted_richness[,match(plotnum, unique_list$uplot)]==1 & presence_data_total[,4,match(plotnum, unique_list$uplot)]==1)
  
  siteuse = which(is.finite(tmp$totaln)) # identifies whether plots are large or small
  tmp = tmp[siteuse,,drop=FALSE]
  
  if(soiln_switch) {
    tmp[,2] = mean(soiln_data_total)
  }
  
  totaln_data = as_tibble(tmp)
  biomass_data = biomass_data_total[,,,which(unique_list$uplot==plotnum)]
  biomass_data = biomass_data[,siteuse,,drop=FALSE]
  site_data = site_data[siteuse,,drop=FALSE]
  
  # replace missing values in first year
  #missing_rows = which(is.na(rowSums(biomass_data[1,,])))
  #if(sum(missing_rows)>0) {
  #  for(ii in 1:length(missing_rows)) {
  #    pstmp = which(is.na(biomass_data[1,missing_rows[ii],]))
  #    biomass_data[1,missing_rows[ii],pstmp] =
  #      seeding_data_total[siteuse[missing_rows[ii]],pstmp,which(unique_list$uplot==plotnum)]*
  #      trait_data[pstmp,"s"]
  #  }
  #}
  #biomass_data[1,,][is.na(biomass_data[1,,])] = 0
  
  # get initial seeded biomass
  b_init = matrix(nrow = dim(biomass_data)[2], ncol = dim(biomass_data)[3])
  for(ii in 1:nrow(b_init)) {
    b_init[ii,] =
      seeding_data_total[siteuse[ii],,which(unique_list$uplot==plotnum)]*
      trait_data[,"s"]
  }
  row.names(b_init) = dimnames(biomass_data)[[2]]
  colnames(b_init) = dimnames(biomass_data)[[3]]
  
  # make dispersal matrix
  if(dim(site_data)[1] == 64 & dim(trait_data)[1] == 5) {
    # use pre-compiled data to save times
    pd_array = pd_array_precompiled
  } else {
    pd_array = make_disp_matrix(trait_data, site_data)
  }
  
  # add in mask to dispersal function to deal with effects of weeding.
  for(i in 1:dim(pd_array)[3]) {
    notplantedsp = unname(which(b_init[i,] == 0))
    if(sum(notplantedsp)>0) {
      pd_array[notplantedsp,,i] = 0
    }
  }
  
  print(paste("plot: ", plotnum))
  #simulate_from = "simulated";switch_off_rstar = FALSE;switch_off_overlap = FALSE;switch_off_height = FALSE;switch_off_dispersal = FALSE;mod_bm = NA;mod_rgr = NA;sp_list = NA;run_in_e26 = TRUE
  community_E26param <- simulate_coexistence(biomass_data=biomass_data,
                                             totaln_data=totaln_data,
                                             site_data=site_data,
                                             trait_data=trait_data,
                                             monoculture_biomass=monoculture_biomass,
                                             rgr_data=rgr_data,
                                             simulate_from = "simulated",
                                             pd_array = pd_array,
                                             b_init = b_init,
                                             run_in_e26 = TRUE,
                                             switch_off_rstar = rstar_switch,
                                             switch_off_height = height_switch,
                                             switch_off_root = root_switch,
                                             switch_off_dispersal = dispersal_switch,
                                             switch_off_pheno = pheno_switch,
                                             switch_off_lot = lot_switch,
                                             switch_off_bstar = bstar_switch,
                                             switch_off_rgr = rgr_switch,
                                             switch_off_b_init = b_init_switch,
                                             switch_off_fecundity_diff = rep_local_diff_switch,
                                             switch_off_mor_diff = mor_diff_switch)
  sim_data_total[,,,match(plotnum, unique_list$uplot)] = community_E26param
}

#saveRDS(sim_data_total, paste0("results/", args[1], "_", jobid, ".rds"))

if(FALSE) {
  # plot time-series for individual plots
  plotnum = 3
  community_E26param = sim_data_total[,,,match(plotnum, unique_list$uplot)]
  
  par(mar=c(4,4,2,2))
  for(i in 1:length(siteuse)) {
    observeddynamics_site_i=biomass_data[,i,1:5]
    predicteddynamics_site_i_E26p=community_E26param[,i,1:5]
    
    matplot(0:9, predicteddynamics_site_i_E26p+1, type = "l",
            col=1:5, lty="solid", lwd=2,
            xlab = "time", ylab = "biomass", ylim=c(1,500),
            main =paste("plot ", plotnum, "; site ", siteuse[i], sep = ""), log = "y", axes = FALSE)
    matpoints(jitter(1:9),
              observeddynamics_site_i+1,
              type = "p", pch=1:5,col=1:5)
    
    axis(1)
    axis(2, at=c(1,c(1,2,5,10,20,50,100,200,500)+1),
         c(0,c(1,2,5,10,20,50,100,200,500)), las = 2)
    box()
    
    legend("topleft",
           #       legend = c("Species 1", "Species 2", "Species 3", "Species 4", "Species 5"),
           legend = c("Agropyron repens", "Agrostis scabra", "Andropogon g.", "Poa pratensis", "Schizachyrium s."),
           col=1:5,
           pch = 1:5,
           bty = "n",
           pt.cex = 1,
           cex = 0.8,
           text.col = "black",
           horiz = F,
           lty = 1, lwd=2,
           inset = c(0.01, 0.01))
  }
}




# plot overall fits
# get richness per treatment
planted_richness = apply(presence_data_total, c(1,3), sum)
uplanted_richness = sort(unique(c(planted_richness)))

# coefficient of efficiency
e2fun = function(x,y, dolog=FALSE) {
  if(dolog) {
    x = log(x); y = log(y)
  }
  
  1-mean((x-y)^2)/mean((y-mean(y))^2)
}

# room mean square error
rmse = function(x,y, dolog=FALSE) {
  if(dolog) {
    x = log(x); y = log(y)
  }
  
  sqrt(mean((x-y)^2))
}


#### Extracting data for figure 2 - year 6 biomass estimates and observed values


dim(sim_data_total) ## 10 year, 64 subplots, 5 species, 9 plots

fill_df <-   tibble(species = character(),
                    observed = numeric(),
                    predicted = numeric(),
                    subplot = numeric(),
                    plot = numeric(),
                    year = numeric(),
                    planted_richness = numeric(),
                    sown = numeric())

for(i in 1:length(uplanted_richness)) {
  subplots_use = row(planted_richness)[planted_richness==uplanted_richness[i]]
  plots_use = col(planted_richness)[planted_richness==uplanted_richness[i]]
  
  for(j in 1:length(subplots_use)){
    comb_bio <-  data.frame(t(rbind(biomass_data_total[6, subplots_use[j],,plots_use[j]], sim_data_total[6,subplots_use[j],,plots_use[j]])))
    colnames(comb_bio) <- c("observed", "predicted")
    comb_bio$subplot <- subplots_use[j]
    comb_bio$plot <- plots_use[j]
    comb_bio$year <- 6
    comb_bio <- comb_bio %>% rownames_to_column("species")
    comb_bio$planted_richness <- uplanted_richness[i]
    
    seeded <- data.frame((seeding_data_total[subplots_use[j],,plots_use[j]]))
    colnames(seeded) <- "init"
    seeded <- seeded %>% rownames_to_column("species")
    seeded <- seeded[seeded$init > 0,]
    seeded$sown <- 1
    seeded <- seeded[,-2]
    
    comb_bio <- comb_bio %>% left_join(seeded, by = "species") %>%
      replace(is.na(.), 0)
    
    fill_df <- bind_rows(fill_df, comb_bio)
  }}

remove(b_init, comb_bio, planted_richness, seeded, site_data, tmp, totaln_data)



#### Standardize based on monoculture biomass values ####
load("data/trait_data.rda")
sp_names <- unique(fill_df$species)

get_trait_vector <- function(trait_data, sp_list, trait_name) {
  out <- trait_data[trait_data$species %in% sp_list, trait_name]
  names(out) <- sp_list
  return(out)
}

bstar0 <- get_trait_vector(trait_data, sp_names, "Bm") %>% as_tibble(rownames = "species")

### Standardise by Bm values
df <- fill_df %>%
  filter(sown == 1) %>%
  inner_join(bstar0) %>%
  mutate(obs = observed/value) %>%
  mutate(pred = predicted/value) %>%
  dplyr::select(-observed, -predicted) %>%
  rename(observed = obs, predicted = pred)


df_sums <- fill_df %>%
  filter(sown == 1) %>%
  inner_join(bstar0) %>%
  group_by(plot, subplot, planted_richness) %>%
  summarise(observed = sum(observed), predicted = sum(predicted), value = sum(value)) %>%
  mutate(observed = observed/value, predicted = predicted/value) ### Sum together for total biomass estimate

df_zeros <- df %>% filter(observed == 0 & sown == 1) %>% unite(code, c("plot", "subplot"), remove = FALSE)

total_R1 <- df_sums[df_sums$observed > 0 & df_sums$planted_richness == 1,]
summary(lm(log(observed+1) ~ log(predicted+1), data = total_R1))$r.squared

total_R2 <- df_sums[df_sums$observed > 0 & df_sums$planted_richness == 2,]
summary(lm(log(observed+1) ~ log(predicted+1), data = total_R2))$r.squared

total_R3 <- df_sums[df_sums$observed > 0 & df_sums$planted_richness == 3,]
summary(lm(log(observed+1) ~ log(predicted+1), data = total_R3))$r.squared

total_R5 <- df_sums[df_sums$observed > 0 & df_sums$planted_richness == 5,]
summary(lm(log(observed+1) ~ log(predicted+1), data = total_R5))$r.squared

### summary R2 value for all polycultures
total_poly <- df_sums[df_sums$observed > 0 & df_sums$planted_richness != 1,]
summary(lm(log(observed+1) ~ log(predicted+1), data = total_poly))$r.squared


#### Total biomass
(Fig2a <- ggplot() +
  geom_abline(slope=1, linetype = "dashed") +
  geom_point(aes(y = observed+1, x= predicted+1,
                 colour = as.factor(planted_richness), shape = as.factor(planted_richness)),
             data = df_sums[df_sums$observed > 0,],
             alpha = 1, size = 1) + theme_bw() +
  geom_smooth(method = "lm", formula = y~x, se = FALSE,
              aes(y = observed+1, x= predicted+1,colour = as.factor(planted_richness)),
              df_sums[df_sums$observed > 0,]) +
  scale_colour_manual(labels = c("1 (0.60)", "2 (0.65)", "3 (0.67)", "5 (0.72)"),
                      values = c("1" = "#FEC47D", "2" = "#65A688", "3" = "#B54445", "5" = "#3192C4")) +
  scale_shape_manual(labels = c("1 (0.60)", "2 (0.65)", "3 (0.67)", "5 (0.72)"),
                     values = c("1" = 1, "2" = 2,"3" = 3, "5" = 4)) +
    scale_x_log10(breaks = c(1, 1 + c(1, 2, 5, 10)),
                  labels = c(0, c(1, 2, 5, 10)),
                  limits = c(1, 15)) +
    scale_y_log10(breaks = c(1, 1 + c(1, 2, 5, 10)),
                  labels = c(0, c(1, 2, 5, 10)),
                  limits = c(1, 15)) +
  ylab("Observed") + xlab("") +
 labs(shape = expression("Richness (R"^2*")"), colour = expression("Richness (R"^2*")")) +
  theme(legend.position = c(0.15, 0.8),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        panel.grid.minor = element_blank(),
        plot.tag = element_text(face = "bold", size = 18),
        plot.title = element_text(hjust = 0.5, size = 18),
        axis.text = element_text(size = 14),
        axis.title = element_text(size=16)) +
    labs(title = "Total community biomass", tag = "A"))

##### Agrostis scabra

Aggsc_R1 <- df[df$observed > 0 & df$predicted > 0 & df$species == "Agrostis scabra" &
                 df$planted_richness == 1,]
summary(lm(log(observed+1) ~ log(predicted+1), data = Aggsc_R1))$r.squared

Aggsc_R2 <- df[df$observed > 0 & df$predicted > 0 & df$species == "Agrostis scabra" &
                 df$planted_richness == 2,]
summary(lm(log(observed+1) ~ log(predicted+1), data = Aggsc_R2))$r.squared

#Aggsc_R3 <- df[df$observed > 0 & df$predicted > 0 & df$species == "Agrostis scabra" &
#                 df$planted_richness == 3,]
#summary(lm(log(observed) ~ log(predicted), data = Aggsc_R3))$r.squared

Aggsc_R5 <- df[df$observed > 0 & df$predicted > 0 & df$species == "Agrostis scabra" &
                 df$planted_richness == 5,]
summary(lm(log(observed+1) ~ log(predicted+1), data = Aggsc_R5))$r.squared


(Fig2b <- ggplot() +
  geom_abline(slope=1, linetype = "dashed") +
  geom_point(aes(y = observed+1, x= predicted+1,
                 colour = as.factor(planted_richness), shape = as.factor(planted_richness)),
             data = df[df$observed > 0 & df$predicted > 0 & df$species == "Agrostis scabra",],
             alpha = 1, size = 1) +
  geom_point(aes(y = (observed+1), x= (predicted+1),
                 colour = as.factor(planted_richness), shape = as.factor(planted_richness)),
             data = df_zeros[df_zeros$species == "Agrostis scabra",],
             alpha = 1, size = 1) +
  geom_smooth(method = "lm", formula = y~x, se = FALSE,
              aes(y = observed+1, x= predicted+1,colour = as.factor(planted_richness)),
              df[df$observed > 0 & df$predicted > 0 & df$species == "Agrostis scabra",]) +
  theme_bw() +
  scale_colour_manual(labels = c("1 (0.22)", "2 (0.04)", "5 (1)"),
                      values = c("1" = "#FEC47D", "2" = "#65A688", "5" = "#3192C4")) +
  scale_shape_manual(labels = c("1 (0.22)", "2 (0.04)", "5 (1)"),
                     values = c("1" = 1, "2" = 2, "3" = 3, "5" = 4)) +
    scale_x_log10(breaks = c(1, 1 + c(1, 2, 5, 10)),
                  labels = c(0, c(1, 2, 5, 10)),
                  limits = c(1, 15)) +
    scale_y_log10(breaks = c(1, 1 + c(1, 2, 5, 10)),
                  labels = c(0, c(1, 2, 5, 10)),
                  limits = c(1, 15)) +
    labs(shape = expression("Richness (R"^2*")"), colour = expression("Richness (R"^2*")")) +
  ylab("") + xlab("") +
    theme(legend.position = c(0.82, 0.28),
          legend.title = element_text(size = 12),
          legend.text=element_text(size=12),
          plot.tag = element_text(face = "bold", size = 18),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 18, face = "italic"),
          axis.text = element_text(size = 14),
          axis.title = element_text(size=16)) +
    labs(title = "Agrostis scabra", tag = "B"))


##### Andropogon gerardi

Andge_R1 <- df[df$observed > 0 & df$predicted > 0 & df$species == "Andropogon gerardi" &
                 df$planted_richness == 1,]
summary(lm(log(observed+1) ~ log(predicted+1), data = Andge_R1))$r.squared

Andge_R2 <- df[df$observed > 0 & df$predicted > 0 & df$species == "Andropogon gerardi" &
                 df$planted_richness == 2,]
summary(lm(log(observed+1) ~ log(predicted+1), data = Andge_R2))$r.squared

#Andge_R3 <- df[df$observed > 0 & df$predicted > 0 & df$species == "Andropogon gerardi" &
#df$planted_richness == 3,]
#summary(lm(log(observed) ~ log(predicted), data = Andge_R3))$r.squared

Andge_R5 <- df[df$observed > 0 & df$predicted > 0 & df$species == "Andropogon gerardi" &
                 df$planted_richness == 5,]
summary(lm(log(observed+1) ~ log(predicted+1), data = Andge_R5))$r.squared


(Fig2c <- ggplot() +
    geom_abline(slope=1, linetype = "dashed") +
    geom_point(aes(y = observed+1, x= predicted+1,
                   colour = as.factor(planted_richness), shape = as.factor(planted_richness)),
               data = df[df$observed > 0 & df$predicted > 0 & df$species == "Andropogon gerardi",],
               alpha = 1, size = 1) +
    geom_point(aes(y = (observed+1), x= (predicted+1),
                   colour = as.factor(planted_richness), shape = as.factor(planted_richness)),
               data = df_zeros[df_zeros$species == "Andropogon gerardi",],
               alpha = 1, size = 1) +
    geom_smooth(method = "lm", formula = y~x, se = FALSE,
                aes(y = observed+1, x= predicted+1,colour = as.factor(planted_richness)),
                df[df$observed > 0 & df$predicted > 0 & df$species == "Andropogon gerardi",]) +
    theme_bw() +
    scale_colour_manual(labels = c("1 (0.91)", "2 (0.78)", "5 (0.66)"),
                        values = c("1" = "#FEC47D", "2" = "#65A688", "5" = "#3192C4")) +
    scale_shape_manual(labels = c("1 (0.91)", "2 (0.78)", "5 (0.66)"),
                       values = c("1" = 1, "2" = 2,"3" = 3, "5" = 4)) +
    scale_x_log10(breaks = c(1, 1 + c(1, 2, 5, 10)),
                  labels = c(0, c(1, 2, 5, 10)),
                  limits = c(1, 15)) +
    scale_y_log10(breaks = c(1, 1 + c(1, 2, 5, 10)),
                  labels = c(0, c(1, 2, 5, 10)),
                  limits = c(1, 15)) +
    labs(shape = expression("Richness (R"^2*")"), colour = expression("Richness (R"^2*")")) +
    ylab("") + xlab("") +
    theme(legend.position = c(0.82, 0.28),
          legend.text=element_text(size=12),
          legend.title = element_text(size = 12),
          plot.tag = element_text(face = "bold", size = 18),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 18, face = "italic"),
          axis.text = element_text(size = 14),
          axis.title = element_text(size=14)) +
    labs(title = "Andropogon gerardi", tag = "C"))


##### Elymus repens

Elyre_R1 <- df[df$observed > 0 & df$predicted > 0 & df$species == "Agropyron repens" &
                 df$planted_richness == 1,]
summary(lm(log(observed+1) ~ log(predicted+1), data = Elyre_R1))$r.squared

Elyre_R2 <- df[df$observed > 0 & df$predicted > 0 & df$species == "Agropyron repens" &
                 df$planted_richness == 2,]
summary(lm(log(observed+1) ~ log(predicted+1), data = Elyre_R2))$r.squared

Elyre_R3 <- df[df$observed > 0 & df$predicted > 0 & df$species == "Agropyron repens" &
                 df$planted_richness == 3,]
summary(lm(log(observed+1) ~ log(predicted+1), data = Elyre_R3))$r.squared

#Elyre_R5 <- df[df$observed > 0 & df$predicted > 0 & df$species == "Agropyron repens" & df$planted_richness == 5,]
#summary(lm(log(observed+1) ~ log(predicted+1), data = Elyre_R5))$r.squared


(Fig2d <- ggplot() +
    geom_abline(slope=1, linetype = "dashed") +
    geom_point(aes(y = observed+1, x= predicted+1,
                   colour = as.factor(planted_richness), shape = as.factor(planted_richness)),
               data = df[df$observed > 0 & df$predicted > 0 & df$species == "Agropyron repens",],
               alpha = 1, size = 1) +
    geom_point(aes(y = (observed+1), x= (predicted+1),
                   colour = as.factor(planted_richness), shape = as.factor(planted_richness)),
               data = df_zeros[df_zeros$species == "Agropyron repens",],
               alpha = 1, size = 1) +
    geom_smooth(method = "lm", formula = y~x, se = FALSE,
                aes(y = observed+1, x= predicted+1,colour = as.factor(planted_richness)),
                df[df$observed > 0 & df$predicted > 0 & df$species == "Agropyron repens",]) +
    theme_bw() +
    scale_colour_manual(labels = c("1 (0.85)", "2 (0.73)", "3 (0.97)", "5 (-)"),
                        values = c("1" = "#FEC47D", "2" = "#65A688", "3" = "#B54445", "5" = "#3192C4")) +
    scale_shape_manual(labels = c("1 (0.85)", "2 (0.73)", "3 (0.97)", "5 (-)"),values = c("1" = 1, "2" = 2, "3" = 3, "5" = 4)) +
    scale_x_log10(breaks = c(1, 1 + c(1, 2, 5, 10)),
                  labels = c(0, c(1, 2, 5, 10)),
                  limits = c(1, 15)) +
    scale_y_log10(breaks = c(1, 1 + c(1, 2, 5, 10)),
                  labels = c(0, c(1, 2, 5, 10)),
                  limits = c(1, 15)) +
    labs(shape = expression("Richness (R"^2*")"), colour = expression("Richness (R"^2*")")) +
    ylab("Observed") + xlab("Predicted") +
    theme(legend.position = c(0.82, 0.33),
          legend.title = element_text(size=12),
          legend.text=element_text(size=12),
          plot.tag = element_text(face = "bold", size = 18),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 18, face = "italic"),
          axis.text = element_text(size = 14),
          axis.title = element_text(size=16)) +
    labs(title = "Elymus repens", tag = "D"))


##### Poa pratensis

Poapr_R1 <- df[df$observed > 0 & df$predicted > 0 & df$species == "Poa pratensis" &
                 df$planted_richness == 1,]
summary(lm(log(observed+1) ~ log(predicted+1), data = Poapr_R1))$r.squared

Poapr_R2 <- df[df$observed > 0 & df$predicted > 0 & df$species == "Poa pratensis" &
                 df$planted_richness == 2,]
summary(lm(log(observed+1) ~ log(predicted+1), data = Poapr_R2))$r.squared

Poapr_R3 <- df[df$observed > 0 & df$predicted > 0 & df$species == "Poa pratensis" &
                 df$planted_richness == 3,]
summary(lm(log(observed+1) ~ log(predicted+1), data = Poapr_R3))$r.squared

Poapr_R5 <- df[df$observed > 0 & df$predicted > 0 & df$species == "Poa pratensis" &
                 df$planted_richness == 5,]
summary(lm(log(observed+1) ~ log(predicted+1), data = Poapr_R5))$r.squared

(Fig2e <- ggplot() +
    geom_abline(slope=1, linetype = "dashed") +
    geom_point(aes(y = observed+1, x= predicted+1,
                   colour = as.factor(planted_richness), shape = as.factor(planted_richness)),
               data = df[df$observed > 0 & df$predicted > 0 & df$species == "Poa pratensis",],
               alpha = 1, size = 1) +
    geom_point(aes(y = (observed+1), x= (predicted+1),
                   colour = as.factor(planted_richness), shape = as.factor(planted_richness)),
               data = df_zeros[df_zeros$species == "Poa pratensis",],
               alpha = 1, size = 1) +
    geom_smooth(method = "lm", formula = y~x, se = FALSE,
                aes(y = observed+1, x= predicted+1,colour = as.factor(planted_richness)),
                df[df$observed > 0 & df$predicted > 0 & df$species == "Poa pratensis",]) +
    theme_bw() +
    scale_colour_manual(labels = c("1 (0.73)", "2 (0.36)", "3 (0.07)", "5 (0.29)"),
                        values = c("1" = "#FEC47D", "2" = "#65A688", "3" = "#B54445", "5" = "#3192C4")) +
    scale_shape_manual(labels = c("1 (0.73)", "2 (0.36)", "3 (0.07)", "5 (0.29)"),
                       values = c("1" = 1, "2" = 2, "3" = 3, "5" = 4)) +
    scale_x_log10(breaks = c(1, 1 + c(1, 2, 5, 10)),
                  labels = c(0, c(1, 2, 5, 10)),
                  limits = c(1, 15)) +
    scale_y_log10(breaks = c(1, 1 + c( 1, 2, 5, 10)),
                  labels = c(0, c( 1, 2, 5, 10)),
                  limits = c(1, 15)) +
    labs(shape = expression("Richness (R"^2*")"), colour = expression("Richness (R"^2*")")) +
    ylab("") + xlab("Predicted") +
    theme(legend.position = c(0.82, 0.33),
          legend.title = element_text(size=12),
          legend.text=element_text(size=12),
          panel.grid.minor = element_blank(),
          plot.tag = element_text(face = "bold", size = 18),
          plot.title = element_text(hjust = 0.5, size = 18, face = "italic"),
          axis.text = element_text(size = 14),
          axis.title = element_text(size=16)) +
    labs(title = "Poa pratensis", tag = "E"))


#### Schizachyrium scoparium

Schsc_R1 <- df[df$observed > 0 & df$predicted > 0 & df$species == "Schizachyrium scoparium" &
                 df$planted_richness == 1,]
summary(lm(log(observed+1) ~ log(predicted+1), data = Schsc_R1))$r.squared

Schsc_R2 <- df[df$observed > 0 & df$predicted > 0 & df$species == "Schizachyrium scoparium" &
                 df$planted_richness == 2,]
summary(lm(log(observed+1) ~ log(predicted+1), data = Schsc_R2))$r.squared

Schsc_R3 <- df[df$observed > 0 & df$predicted > 0 & df$species == "Schizachyrium scoparium" &
                 df$planted_richness == 3,]
summary(lm(log(observed+1) ~ log(predicted+1), data = Schsc_R3))$r.squared

Schsc_R5 <- df[df$observed > 0 & df$predicted > 0 & df$species == "Schizachyrium scoparium" &
                 df$planted_richness == 5,]
summary(lm(log(observed+1) ~ log(predicted+1), data = Schsc_R5))$r.squared


(Fig2f <- ggplot() +
    geom_abline(slope=1, linetype = "dashed") +
    geom_point(aes(y = observed+1, x= predicted+1,
                   colour = as.factor(planted_richness), shape = as.factor(planted_richness)),
               data = df[df$observed > 0 & df$predicted > 0 & df$species == "Schizachyrium scoparium",],
               alpha = 1, size = 1) +
    geom_point(aes(y = (observed+1), x= (predicted+1),
                   colour = as.factor(planted_richness), shape = as.factor(planted_richness)),
               data = df_zeros[df_zeros$species == "Schizachyrium scoparium",],
               alpha = 1, size = 1) +
    geom_smooth(method = "lm", formula = y~x, se = FALSE,
                aes(y = observed+1, x= predicted+1,colour = as.factor(planted_richness)),
                df[df$observed > 0 & df$predicted > 0 & df$species == "Schizachyrium scoparium",]) +
    theme_bw() +
    scale_colour_manual(labels = c("1 (0.77)", "2 (0.68)", "3 (0.65)", "5 (0.21)"),
                        values = c("1" = "#FEC47D", "2" = "#65A688", "3" = "#B54445", "5" = "#3192C4")) +
    scale_shape_manual(labels = c("1 (0.77)", "2 (0.68)", "3 (0.65)", "5 (0.21)"),
                       values = c("1" = 1, "2" = 2, "3" = 3, "5" = 4)) +
    scale_x_log10(breaks = c(1, 1 + c(1, 2, 5, 10)),
                  labels = c(0, c(1, 2, 5, 10)),
                  limits = c(1, 15)) +
    scale_y_log10(breaks = c(1, 1 + c(1, 2, 5, 10)),
                  labels = c(0, c(1, 2, 5, 10)),
                  limits = c(1, 15)) +
    ylab("") + xlab("Predicted") +
    labs(shape = expression("Richness (R"^2*")"), colour = expression("Richness (R"^2*")")) +
    theme(legend.position = c(0.82, 0.33),
          legend.title = element_text(size=12),
          legend.text=element_text(size=12),
          panel.grid.minor = element_blank(),
          plot.tag = element_text(face = "bold", size = 18),
          plot.title = element_text(hjust = 0.5, size = 18, face = "italic"),
          axis.text = element_text(size = 14),
          axis.title = element_text(size=16)) +
    labs(title = "Schizachyrium scoparium", tag = "F"))


library(grid) ## export as 1200 x 1000 pixels
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 3, ncol = 3)))
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
}
# Arrange the plots
print(Fig2a, vp = define_region(row = 1:2, col = 1:2))   # Span over two columns
print(Fig2b, vp = define_region(row = 1, col = 3))
print(Fig2c, vp = define_region(row = 2, col = 3))
print(Fig2d, vp = define_region(row = 3, col = 1))
print(Fig2e, vp = define_region(row = 3, col = 2))
print(Fig2f, vp = define_region(row = 3, col = 3))

## Export as 1250 x 1000 or 13.21 x 10.5 pdf




