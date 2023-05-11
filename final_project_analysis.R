# CSCI 5423 Final Project
# Agent-based model of food web evolution
# Analysis of simulation results
# Henry Li

#-----------------------------------------------------------------------------#
# load packages and set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(tidyverse)
source("./ABM_simulator.R")

#-----------------------------------------------------------------------------#
# run simulations and store output in a dataframe
experiment_df <- run_experiment(1, c(25, 50, 100), c(0.1, 0.2, 0.3))

save(experiment_df, file = "./experiment_data.rds")

#-----------------------------------------------------------------------------#
# plot results

# plot species abundance timeseries
species_timeseries <- experiment_df[[7]][[9]]
species_timeseries <- species_timeseries[, c(1:(ncol(species_timeseries)-1))]
species_timeseries <- pivot_longer(species_timeseries, !t, names_to = "species", values_to = "abundance")
species_timeseries$species <- as.factor(species_timeseries$species)
species_timeseries_plot <- ggplot(data = species_timeseries, aes(x = t, y = abundance, color = species)) +
  geom_line() +
  scale_x_continuous(name = "Time") +
  scale_y_continuous(name = "Species Abundance") +
  labs(title = "Species Abundances over Time",
       subtitle = "Global Community S = 100, C = 0.3") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")
species_timeseries_plot

# plot community total abundance timeseries
community_timeseries <- experiment_df[[7]][[9]]
community_timeseries <- community_timeseries[, c("t", "total")]
community_timeseries_plot <- ggplot(data = community_timeseries, aes(x = t, y = total)) +
  geom_line() +
  scale_x_continuous(name = "Time") +
  scale_y_continuous(name = "Total Community Occupancy", limits = c(0, 10000), breaks = seq(0, 10000, 1000)) +
  labs(title = "Total Community Occupancy over Time",
       subtitle = "Global Community S = 100, C = 0.3") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
community_timeseries_plot

# plot community network structure metric timeseries
metric_timeseries <- experiment_df[[6]][[9]]
SC_timeseries <- metric_timeseries[, c("t", "species", "connectance")]
TIB_timeseries <- metric_timeseries[, c("t", "top", "intermediate", "basal")]
TIB_timeseries <- pivot_longer(TIB_timeseries, !t, names_to = "Metric", values_to = "Value")
scale <- 150
SC_timeseries_plot <- ggplot(data = SC_timeseries, aes(x = t, y = connectance)) +
  geom_line(aes(color = "Connectance")) +
  geom_line(aes(y = species/scale, color = "Species")) +
  scale_x_continuous(name = "Time") +
  scale_y_continuous(name = "Network Connectance", sec.axis = sec_axis(~.*scale, name = "Network Size")) +
  labs(title = "Network Size and Connectance over Time",
       subtitle = "Global Community S = 100, C = 0.3",
       color = "Metric") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
TIB_timeseries_plot <- ggplot(data = TIB_timeseries, aes(x = t, y = Value, color = Metric)) +
  geom_line() +
  scale_x_continuous(name = "Time") +
  scale_y_continuous(name = "Proportion") +
  labs(title = "Proportion of Species Types over Time",
       subtitle = "Global Community S = 100, C = 0.3",
       color = "Species Type") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
SC_timeseries_plot
TIB_timeseries_plot