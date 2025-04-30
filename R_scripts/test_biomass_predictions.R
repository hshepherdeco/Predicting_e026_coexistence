#### Catford et al., Mechanistic model
#### R script required for simulations with differing number of attributes switched on

setwd("~/SimulateCoexistence")

devtools::load_all()
library(tidyverse)
args = commandArgs(TRUE)

switches <- c("dispersal", "root", "height", "pheno", "lot", "rstar", "bstar", "rgr", "binit", "rep", "mor")
if(args == "0") {
  args2 = "none"
} else {
  args2 <- sample(switches, args)
}

jobid <- Sys.getenv("SLURM_ARRAY_TASK_ID")

rstar_switch <- TRUE
height_switch <- TRUE
root_switch <- TRUE
pheno_switch <- TRUE
dispersal_switch <- TRUE
lot_switch <- TRUE
bstar_switch <- TRUE
rgr_switch <- TRUE
binit_switch <- TRUE
rep_switch <- TRUE
mor_switch <- TRUE

for(which_switch in args2) {
  assign(paste0(which_switch, "_switch"), FALSE)
}

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
                                             switch_off_b_init = binit_switch,
                                             switch_off_fecun = rep_switch,
                                             switch_off_mor_diff = mor_switch)
  sim_data_total[,,,match(plotnum, unique_list$uplot)] = community_E26param
}

saveRDS(sim_data_total, paste0("results/", args, "_", jobid, ".rds"))
saveRDS(args2, paste0("params/", args, "_", jobid, ".rds"))
