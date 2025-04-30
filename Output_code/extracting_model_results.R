rm(list=ls())
setwd("~/simulateCoexistence")
devtools::load_all()
library(tidyverse)
#args = commandArgs(TRUE)
#args2 = str_split(args, pattern = "_")
#args2 = args2[[1]]
#jobid <- Sys.getenv("SLURM_ARRAY_TASK_ID")

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
                                             switch_off_temp = pheno_switch,
                                             switch_off_lot = lot_switch,
                                             switch_off_bstar = bstar_switch,
                                             switch_off_rgr = rgr_switch,
                                             switch_off_b_init = b_init_switch,
                                             switch_off_rep_local_diff = rep_local_diff_switch,
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

pdf("figures/goodness_of_fit_allon.pdf", width = 12, height=30)
fitdat = NULL
par(mar=c(4,4,2,2), mfcol=c(9, 4))
for(i in 1:length(uplanted_richness)) {
  subplots_use = row(planted_richness)[planted_richness==uplanted_richness[i]]
  plots_use = col(planted_richness)[planted_richness==uplanted_richness[i]]
  for(j in 1:length(unique_list$uyear)) {
    plot(c(1,501), c(1,501), type = "n",
         xlab = "estimated", ylab = "observed",
         axes = F, log="xy")
    title(paste("richness: ", uplanted_richness[i], "; year: ", unique_list$uyear[j]-1984, sep = ""))

    axis(1, at=c(1,c(1,2,5,10,20,50,100,200,500)+1),
         c(0,c(1,2,5,10,20,50,100,200,500)))
    axis(2, at=c(1,c(1,2,5,10,20,50,100,200,500)+1),
         c(0,c(1,2,5,10,20,50,100,200,500)), las = 2)
    box()

    if(i==1 & j ==1) {
      legend("topleft", unique_list$uspecies, col = 1:5, lty = 1, pch = 1:5, bty = "n")
    }

    for(k in 1:length(unique_list$uspecies)) {
      # year, subplot, species, plot
      est = sim_data_total[cbind(j,subplots_use,k,plots_use)]+1
      obs = biomass_data_total[cbind(j,subplots_use,k,plots_use)]+1

      points(est, obs, pch = k, col = k)

      kp = which(est>1 & obs >1 & !is.na(est) & !is.na(obs))
      if(length(kp)>0) {
        est = est[kp]
        obs = obs[kp]

        mod = lm(log(obs)~log(est))
        est_plot = seq(min(est, na.rm=T), max(est, na.rm=T), length=100)

        fitdat = rbind(fitdat,
                       data.frame(richness = uplanted_richness[i],
                                  year = unique_list$uyear[j],
                                  species = unique_list$uspecies[k],
                                  r2 = summary(mod)$r.squared,
                                  e2 = e2fun(est, obs),
                                  rmse = rmse(est, obs)))
        lines(est_plot, exp(suppressWarnings(predict(mod, newdata=data.frame(est=est_plot)))), col = k)
      }
    }
    abline(a=0,b=1,lwd=2,col="darkgrey", lty=2)
  }
}
dev.off()

write.csv(fitdat, "figures/goodness_of_fit_allon.csv", row.names=FALSE)







# separate by plot
# get order by totaln
plottotaln = colMeans(soiln_data_total)
plotorder = order(plottotaln)


pdf("figures/goodness_of_fit_byplot_allon.pdf", width = 12, height=30)
fitdat_plot = NULL
for(kk in plotorder) {
  par(mar=c(4,4,2,2), mfcol=c(9, 4), oma = c(0,0,4,0))
  for(i in 1:length(uplanted_richness)) {
    subplots_use = which(planted_richness[,kk]==uplanted_richness[i])
    plots_use = kk
    for(j in 1:length(unique_list$uyear)) {
      plot(c(1,501), c(1,501), type = "n",
           xlab = "estimated", ylab = "observed",
           axes = F, log="xy")
      title(paste("richness: ", uplanted_richness[i], "; year: ", unique_list$uyear[j]-1984, sep = ""))

      axis(1, at=c(1,c(1,2,5,10,20,50,100,200,500)+1),
           c(0,c(1,2,5,10,20,50,100,200,500)))
      axis(2, at=c(1,c(1,2,5,10,20,50,100,200,500)+1),
           c(0,c(1,2,5,10,20,50,100,200,500)), las = 2)
      box()

      if(i==1 & j ==1) {
        title(paste("plot = ", unique_list$uplot[kk], "; totaln = ", round(plottotaln[kk],4), sep = ""), outer=TRUE)
        legend("topleft", unique_list$uspecies, col = 1:5, lty = 1, pch = 1:5, bty = "n")
      }

      for(k in 1:length(unique_list$uspecies)) {
        est = sim_data_total[cbind(j,subplots_use,k,plots_use)]+1
        obs = biomass_data_total[cbind(j,subplots_use,k,plots_use)]+1

        points(est, obs, pch = k, col = k)

        kp = which(est>1 & obs >1 & !is.na(est) & !is.na(obs))
        if(length(kp)>0) {
          est = est[kp]
          obs = obs[kp]

          mod = lm(log(obs)~log(est))
          est_plot = seq(min(est, na.rm=T), max(est, na.rm=T), length=100)

          fitdat_plot = rbind(fitdat_plot,
                              data.frame(richness = uplanted_richness[i],
                                         year = unique_list$uyear[j],
                                         species = unique_list$uspecies[k],
                                         plot = unique_list$uplot[kk],
                                         totaln = plottotaln[kk],
                                         r2 = summary(mod)$r.squared,
                                         e2 = e2fun(est, obs),
                                         rmse = rmse(est, obs)))
          lines(est_plot, exp(suppressWarnings(predict(mod, newdata=data.frame(est=est_plot)))), col = k)
        }
      }
      abline(a=0,b=1,lwd=2,col="darkgrey", lty=2)
    }
  }
}
dev.off()

write.csv(fitdat_plot, "figures/goodness_of_fit_plot_allon.csv", row.names=FALSE)




############### Jane figures for proposal
pdf("figures/goodness_of_fit_allon_summary.pdf", width = 14, height=7)
par(mar=c(4,4,1,1), mfcol=c(1, 1))
m = cbind(c(1,1,1), c(1,1,1), c(1,1,1), c(2,4,6), c(3,5,7))
layout(m)
for(k in 0:length(unique_list$uspecies)) {
  plot(c(1,501), c(1,501), type = "n",
       xlab = "Estimated", ylab = "Observed",
       axes = F, log="xy")
  title(c("Total Community", unique_list$uspecies)[k+1])

  axis(1, at=c(1,c(1,2,5,10,20,50,100,200,500)+1),
       c(0,c(1,2,5,10,20,50,100,200,500)))
  axis(2, at=c(1,c(1,2,5,10,20,50,100,200,500)+1),
       c(0,c(1,2,5,10,20,50,100,200,500)), las = 2)
  box()
  r2lst = NULL

  for(i in 1:length(uplanted_richness)) {
    subplots_use = row(planted_richness)[planted_richness==uplanted_richness[i]]
    plots_use = col(planted_richness)[planted_richness==uplanted_richness[i]]
    for(j in (length(unique_list$uyear)-1)) {
      r2tmp = "-"
      if(k == 0) {
        # year, subplot, species, plot
        tmp = apply(sim_data_total, c(1,2,4), sum)
        est = tmp[cbind(j,subplots_use,plots_use)]+1
        tmp = apply(biomass_data_total, c(1,2,4), sum)
        obs = tmp[cbind(j,subplots_use,plots_use)]+1

        points(est, obs, pch = i, col = i)

        kp = which(est>1 & obs >1 & !is.na(est) & !is.na(obs))
        if(length(kp)>0) {
          est = est[kp]
          obs = obs[kp]

          mod = lm(log(obs)~log(est))
          est_plot = seq(min(est, na.rm=T), max(est, na.rm=T), length=100)

          lines(est_plot, exp(suppressWarnings(predict(mod, newdata=data.frame(est=est_plot)))), col = i, lwd = 1.5)
          r2tmp = round(summary(mod)$r.sq,2)
        }
      } else {
        # year, subplot, species, plot
        est = sim_data_total[cbind(j,subplots_use,k,plots_use)]+1
        obs = biomass_data_total[cbind(j,subplots_use,k,plots_use)]+1

        points(est, obs, pch = i, col = i)

        kp = which(est>1 & obs >1 & !is.na(est) & !is.na(obs))
        if(length(kp)>0) {
          est = est[kp]
          obs = obs[kp]

          mod = lm(log(obs)~log(est))
          est_plot = seq(min(est, na.rm=T), max(est, na.rm=T), length=100)

          lines(est_plot, exp(suppressWarnings(predict(mod, newdata=data.frame(est=est_plot)))), col = i, lwd = 1.5)
          r2tmp = round(summary(mod)$r.sq,2)
        }
      }

      r2lst = c(r2lst, r2tmp)
    }
  }
  abline(a=0,b=1,lwd=2,col="darkgrey", lty=2)
  legend("topleft", legend = paste(uplanted_richness, " (", r2lst, ")", sep = ""), col = 1:5, lty = 1, pch = 1:5, bty = "n", title = expression(paste("Species Richness (", "R"^2, "):")), cex = c(2,1,1,1,1,1)[k+1])
}
dev.off()

#### Plotting figure 2 - for paper
#### Just using 1990 data (year 6)

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

#write.csv(fill_df, "E26_modelled_results_year6.csv")

fit1992 <- fitdat %>% filter(year == "1992")
#write.csv(fitdat, "E26_modelled_fits_year8.csv")
