### Figure 2 - Catford et al., Mechanistic model
### Comparing modelled results to observed values

rm(list=ls())

library(tidyverse)
library(rcartocolor)

setwd("simulateCoexistence/")

#### Standardize based on monoculture biomass values ####

load("data/trait_data.rda")
sp_names <- unique((read.csv("E26_modelled_results_year6.csv"))$species)

get_trait_vector <- function(trait_data, sp_list, trait_name) {
  out <- trait_data[trait_data$species %in% sp_list, trait_name]
  names(out) <- sp_list
  return(out)
}

bstar0 <- get_trait_vector(trait_data, sp_names, "Bm") %>% as_tibble(rownames = "species")

### Standardise by Bm values
df <- read.csv("E26_modelled_results_year6.csv") %>%
  filter(sown == 1) %>%
  inner_join(bstar0) %>%
  mutate(obs = observed/value) %>%
  mutate(pred = predicted/value) %>%
  dplyr::select(-observed, -predicted) %>%
  rename(observed = obs, predicted = pred)


df_sums <- read.csv("E26_modelled_results_year6.csv") %>%
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




