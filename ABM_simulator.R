# CSCI 5423 Final Project
# Agent-based model of food web evolution
# Henry Li

#-----------------------------------------------------------------------------#
# load packages and set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(tidyverse)
library(igraph)
library(EBImage)
source("./niche.R")
source("./niche_validation.R")

#-----------------------------------------------------------------------------#
# define helper functions

# initialize_species_set(): function to create the set of global species in the simulation and their attributes
# input: Number of species S, number of basal species B, connectance C for positive integers S and B, and C in [0,1]
# output: list(S x 6 dataframe, S x S matrix)
initialize_species_set <- function(S, B, C) {
  master_web <- make_web(S, B, C)
  species_adjacency_matrix <- master_web[[2]]
  
  # create a dataframe to store the species in the global food web
  species_df <- as.data.frame(master_web[[3]], row.names = FALSE)
  species_df$species_id <- seq(1, S, 1)
  # get the trophic level of each species
  species_df$trophic_level <- trophic_levels(species_adjacency_matrix)
  # define the metabolism of a species according to the CDF of a beta distribution, with the lowest trophic level with m = 0 and the highest trophic level with m = 1
  species_df$m_i <- pbeta((species_df$trophic_level - 1) / (max(species_df$trophic_level) - 1), 1, 2)
  # set the x-position of each species for visualization
  species_df$x_position <- runif(nrow(species_df))
  
  # set the color of each species for visualization
  species_df <- species_df[order(species_df$trophic_level), ]
  species_df$color <- NA
  number_of_plants <- sum(species_df$trophic_level == 1)
  number_of_tl <- max(floor(species_df$trophic_level))
  plant_palette <- terrain.colors(number_of_plants)
  animal_palette <- hsv(1, seq(0,1,length.out = number_of_tl + 1) , 1)
  for (i in  1:nrow(species_df)) { # assign a color to each species
    if (i < number_of_plants + 1) { # if species is a plant
      species_df$color[i] <- plant_palette[i]
    }
    else { # if the species is an animal
      species_df$color[i] <- animal_palette[floor(species_df$trophic_level[i] + 1)]
    }
  }
  species_df <- species_df[order(species_df$species_id), ]
  
  # reorder columns of species_df
  species_df <- species_df[, c(4, 1, 2, 3, 6, 5, 7, 8)]
  
  return(list(species_df, species_adjacency_matrix))
}

# initialize_environment(): function to initialize the environment
# input: list(S x 6 dataframe, S x S matrix) created by initialize_species_set()
# output: 1 x 4 dataframe
initialize_environment <- function(species_list) {
  species_df <- species_list[[1]]
  
  # create a dataframe to store the species present in the environment
  environment_df <- as.data.frame(matrix(NA, 1, 4))
  colnames(environment_df) <- c("species_id", "x", "y", "energy")
  # add a random basal species to a random location in the environment
  basal_index <- sample(which(species_df$trophic_level == 1), 1)
  environment_df$species_id[1] <- species_df$species_id[basal_index]
  environment_df$x[1] <- sample(1:100, 1)
  environment_df$y[1] <- sample(1:100, 1)
  environment_df$energy[1] <- 1
  
  return(environment_df)
}

# get_invader_location(): function to get an occupied cell for an invader
# input: integer i for species_id of invader, n x 4 dataframe for non-negative integer n, and list(S x 6 dataframe, S x S matrix) created by initialize_species_set()
# output: 2 x 1 vector [x, y] for x, y integers in 1:100
get_invader_location <- function(i, environment_df, species_list) {
  # subset environment_df to plants in the environment and animals in the environment
  species_df <- species_list[[1]]
  plant_species <- unique(species_df$species_id[species_df$trophic_level == 1])
  plant_environment_df <- environment_df[environment_df$species_id %in% plant_species, ]
  animal_environment_df <- environment_df[!(environment_df$species_id %in% plant_species), ]
  
  proposed_x <- sample(1:100, 1)
  proposed_y <- sample(1:100, 1)
  feasible <- 0
  while (feasible == 0) { # check if the proposed cell is unoccupied given the current environment_df
    if (i %in% plant_species) { # if the invader is a plant species
      if (sum(environment_df$x == proposed_x & environment_df$y == proposed_y) == 0) { # stop if the proposed cell is unoccupied
        feasible <- 1
      }
      else { # otherwise redraw a new proposed cell
        proposed_x <- sample(1:100, 1)
        proposed_y <- sample(1:100, 1)
      }
    }
    else { # if the invader is an animal species
      if ((sum(environment_df$x == proposed_x & environment_df$y == proposed_y) == 0) || (sum(environment_df$x == proposed_x & environment_df$y == proposed_y) == 1 & environment_df$species_id[which(environment_df$x == proposed_x & environment_df$y == proposed_y)] %in% plant_species)) { # stop if the proposed cell is unoccupied or occupied by a plant only
        feasible <- 1
      }
      else { # otherwise draw a new proposed cell
        proposed_x <- sample(1:100, 1)
        proposed_y <- sample(1:100, 1)
      }
    }
  }
  
  return(c(proposed_x, proposed_y))
}

# get_plant_neighbors(): function to get a dataframe of a cell's plant neighbors and attributes, including empty neighboring cells
# input: integers x and y between 1 and 100, n x 4 dataframe  for non-negative integer n, and list(S x 6 dataframe, S x S matrix) created by initialize_species_set()
# output: 8 x 4 dataframe
get_plant_neighbors <- function(x, y, environment_df, species_list) {
  # subset environment_df to just plants in the environment
  species_df <- species_list[[1]]
  plant_species <- unique(species_df$species_id[species_df$trophic_level == 1])
  environment_df <- environment_df[environment_df$species_id %in% plant_species, ]
  
  # create dataframe to store neighbors and their attributes
  neighbors_df <- as.data.frame(matrix(NA, 8, 4))
  colnames(neighbors_df) <- c("species_id", "x", "y", "energy")
  
  # populate neighbors_df
  neighbors_df$species_id <- NA
  # lord forgive me idk how to wrap the environment
  if (x == 1 & y == 1) {
    neighbors_df$x <- c(x-2, x-2, x-2, x, x, x+1, x+1, x+1) %% 101
    neighbors_df$y <- c(y-2, y, y+1, y-2, y+1, y-2, y, y+1) %% 101
  }
  else if (x == 100 & y == 1) {
    neighbors_df$x <- c(x-1, x-1, x-1, x, x, x+2, x+2, x+2) %% 101
    neighbors_df$y <- c(y-2, y, y+1, y-2, y+1, y-2, y, y+1) %% 101 
  }
  else if (x == 1 & y == 100) {
    neighbors_df$x <- c(x-2, x-2, x-2, x, x, x+1, x+1, x+1) %% 101
    neighbors_df$y <- c(y-1, y, y+2, y-1, y+2, y-1, y, y+2) %% 101
  }
  else if (x == 100 & y == 100) {
    neighbors_df$x <- c(x-1, x-1, x-1, x, x, x+2, x+2, x+2) %% 101
    neighbors_df$y <- c(y-1, y, y+2, y-1, y+2, y-1, y, y+2) %% 101
  }
  else if (x == 1 & y < 100 & y > 1) {
    neighbors_df$x <- c(x-2, x-2, x-2, x, x, x+1, x+1, x+1) %% 101
    neighbors_df$y <- c(y-1, y, y+1, y-1, y+1, y-1, y, y+1) %% 101
  }
  else if (x == 100 & y < 100 & y > 1) {
    neighbors_df$x <- c(x-1, x-1, x-1, x, x, x+2, x+2, x+2) %% 101
    neighbors_df$y <- c(y-1, y, y+1, y-1, y+1, y-1, y, y+1) %% 101
  }
  else if (y == 1 & x < 100 & x > 1) {
    neighbors_df$x <- c(x-1, x-1, x-1, x, x, x+1, x+1, x+1) %% 101
    neighbors_df$y <- c(y-2, y, y+1, y-2, y+1, y-2, y, y+1) %% 101
  }
  else if (y == 100 & x < 100 & x > 1) {
    neighbors_df$x <- c(x-1, x-1, x-1, x, x, x+1, x+1, x+1) %% 101
    neighbors_df$y <- c(y-1, y, y+2, y-1, y+2, y-1, y, y+2) %% 101
  }
  else {
    neighbors_df$x <- c(x-1, x-1, x-1, x, x, x+1, x+1, x+1) %% 101
    neighbors_df$y <- c(y-1, y, y+1, y-1, y+1, y-1, y, y+1) %% 101 
  }
  neighbors_df$energy <- 0
  for (i in 1:nrow(neighbors_df)) {
    # find the current occupant of each neighboring cell
    matching_neighbor_index <- which(environment_df$x == neighbors_df$x[i] & environment_df$y == neighbors_df$y[i] & environment_df$energy > 0)
    if (length(matching_neighbor_index) > 0) {
      # fill in with neighbor's attributes if neighboring cell is occupied
      neighbors_df$species_id[i] <- environment_df$species_id[matching_neighbor_index]
      neighbors_df$energy[i] <- environment_df$energy[matching_neighbor_index]
    }
  }
  
  return(neighbors_df)
}

# get_animal_neighbors(): function to get a dataframe of a cell's animal neighbors and attributes, including empty neighboring cells
# input: integers x and y between 1 and 100, n x 4 dataframe  for non-negative integer n, and list(S x 6 dataframe, S x S matrix) created by initialize_species_set()
# output: 8 x 4 dataframe
get_animal_neighbors <- function(x, y, environment_df, species_list) {
  # subset environment_df to just animals in the environment
  species_df <- species_list[[1]]
  animal_species <- unique(species_df$species_id[species_df$trophic_level > 1])
  environment_df <- environment_df[environment_df$species_id %in% animal_species, ]
  # create dataframe to store neighbors and their attributes
  neighbors_df <- as.data.frame(matrix(NA, 8, 4))
  colnames(neighbors_df) <- c("species_id", "x", "y", "energy")
  
  # populate neighbors_df
  neighbors_df$species_id <- NA
  # lord forgive me idk how to wrap the environment
  if (x == 1 & y == 1) {
    neighbors_df$x <- c(x-2, x-2, x-2, x, x, x+1, x+1, x+1) %% 101
    neighbors_df$y <- c(y-2, y, y+1, y-2, y+1, y-2, y, y+1) %% 101
  }
  else if (x == 100 & y == 1) {
    neighbors_df$x <- c(x-1, x-1, x-1, x, x, x+2, x+2, x+2) %% 101
    neighbors_df$y <- c(y-2, y, y+1, y-2, y+1, y-2, y, y+1) %% 101 
  }
  else if (x == 1 & y == 100) {
    neighbors_df$x <- c(x-2, x-2, x-2, x, x, x+1, x+1, x+1) %% 101
    neighbors_df$y <- c(y-1, y, y+2, y-1, y+2, y-1, y, y+2) %% 101
  }
  else if (x == 100 & y == 100) {
    neighbors_df$x <- c(x-1, x-1, x-1, x, x, x+2, x+2, x+2) %% 101
    neighbors_df$y <- c(y-1, y, y+2, y-1, y+2, y-1, y, y+2) %% 101
  }
  else if (x == 1 & y < 100 & y > 1) {
    neighbors_df$x <- c(x-2, x-2, x-2, x, x, x+1, x+1, x+1) %% 101
    neighbors_df$y <- c(y-1, y, y+1, y-1, y+1, y-1, y, y+1) %% 101
  }
  else if (x == 100 & y < 100 & y > 1) {
    neighbors_df$x <- c(x-1, x-1, x-1, x, x, x+2, x+2, x+2) %% 101
    neighbors_df$y <- c(y-1, y, y+1, y-1, y+1, y-1, y, y+1) %% 101
  }
  else if (y == 1 & x < 100 & x > 1) {
    neighbors_df$x <- c(x-1, x-1, x-1, x, x, x+1, x+1, x+1) %% 101
    neighbors_df$y <- c(y-2, y, y+1, y-2, y+1, y-2, y, y+1) %% 101
  }
  else if (y == 100 & x < 100 & x > 1) {
    neighbors_df$x <- c(x-1, x-1, x-1, x, x, x+1, x+1, x+1) %% 101
    neighbors_df$y <- c(y-1, y, y+2, y-1, y+2, y-1, y, y+2) %% 101
  }
  else {
    neighbors_df$x <- c(x-1, x-1, x-1, x, x, x+1, x+1, x+1) %% 101
    neighbors_df$y <- c(y-1, y, y+1, y-1, y+1, y-1, y, y+1) %% 101 
  }
  neighbors_df$energy <- 0
  for (i in 1:nrow(neighbors_df)) {
    # find the current occupant of each neighboring cell
    matching_neighbor_index <- which(environment_df$x == neighbors_df$x[i] & environment_df$y == neighbors_df$y[i] & environment_df$energy > 0)
    if (length(matching_neighbor_index) > 0) {
      # fill in with neighbor's attributes if neighboring cell is occupied
      neighbors_df$species_id[i] <- environment_df$species_id[matching_neighbor_index]
      neighbors_df$energy[i] <- environment_df$energy[matching_neighbor_index]
    }
  }
  
  return(neighbors_df)
}

# get_plant_availables(): function to get a dataframe of available cells to add offspring to
# input: 8 x 4 dataframe created by get_plant_neighbors() and 8 x 4 dataframe created by get_animal_neighbors()
# output: n x 4 dataframe with integer n between 0 and 8
get_plant_availables <- function(plant_neighbors, animal_neighbors) {
  available_cells <- plant_neighbors
  for (i in 1:nrow(available_cells)) {
    if (!is.na(animal_neighbors$species_id[i])) {
      available_cells$species_id[i] <- animal_neighbors$species_id[i]
    }
  }
  available_cells <- available_cells[is.na(available_cells$species_id), c("x", "y")]
  
  return(available_cells)
}

# get_animal_availables(): function to get a dataframe of available cells to move to
# input: 8 x 4 dataframe created by get_animal_neighbors()
# output: n x 4 dataframe with integer n between 0 and 8
get_animal_availables <- function(animal_neighbors) {
  available_cells <- animal_neighbors[is.na(animal_neighbors$species_id), c("x", "y")]
  
  return(available_cells)
}

# increment_timestep(): function to simulate one timestep
# input: n x 4 dataframe for non-negative integer n and list(S x 6 dataframe, S x S matrix) created by initialize_species_set()
# output: m x 4 dataframe for non-negative integer m 
increment_timestep <- function(environment_df, species_list) {
  species_df <- species_list[[1]]
  species_adjacency_matrix <- species_list[[2]]
  # get the list of plant species in the global environment
  plant_species <- unique(species_df$species_id[species_df$trophic_level == 1])
  
  # add a random invader to a random unoccupied cell in the environment
  invader_index <- sample(nrow(species_df), 1)
  invader_location <- get_invader_location(invader_index, environment_df, species_list)
  invader <- c(species_df$species_id[invader_index],
               invader_location[1],
               invader_location[2],
               1)
  environment_df <- rbind(environment_df, invader)
  
  # create a matrix indicating which plant organisms are completely surrounded and should be skipped in the timestep
  plant_environment_df <- environment_df[environment_df$species_id %in% plant_species, ]
  plant_matrix <- matrix(0, nrow = 100, ncol = 100)
  plant_matrix[as.matrix(plant_environment_df[, c("x", "y")])] <- 1
  kernel <- matrix(rep(1, 9), nrow = 3, ncol = 3)
  plant_matrix_convolution <- filter2(plant_matrix, kernel, boundary = "circular")
  plant_matrix_convolution <- round(plant_matrix_convolution)
  plant_matrix_convolution[plant_matrix_convolution < 9] <- 0
  plant_matrix_convolution[plant_matrix_convolution == 9] <- 1
  surrounded_plants <- plant_matrix * plant_matrix_convolution
  
  # create a random order of organisms in the environment
  order <- sample(1:nrow(environment_df), replace = FALSE)
  # iterate over all animals, and plants that are not surrounded in the environment
  order <- order[!(environment_df$species_id[order] %in% plant_species) |
                 ((environment_df$species_id[order] %in% plant_species) & (surrounded_plants[cbind(environment_df$x[order], environment_df$y[order])] != 1))]
  for (i in order) {
    if (environment_df$energy[i] == 0) { # if the organism is dead, skip
      next
    }
    else if (environment_df$species_id[i] %in% plant_species) { # otherwise, if the alive organism is a plant
      # get neighbors and available cells
      plant_neighbors <- get_plant_neighbors(environment_df$x[i], environment_df$y[i], environment_df, species_list)
      animal_neighbors <- get_animal_neighbors(environment_df$x[i], environment_df$y[i], environment_df, species_list)
      plant_availables <- get_plant_availables(plant_neighbors, animal_neighbors)
      
      # "feed"
      environment_df$energy[i] <- environment_df$energy[i] + 0.1
      # reproduce
      if (environment_df$energy[i] > 1 & nrow(plant_availables) > 0) {
        # randomly choose an unoccupied neighboring cell
        new_cell_index <- sample(nrow(plant_availables), 1)
        
        # create new offspring to be added
        offspring <- c(environment_df$species_id[i],
                       plant_availables$x[new_cell_index],
                       plant_availables$y[new_cell_index],
                       environment_df$energy[i] - 1)
        environment_df <- rbind(environment_df, offspring)
        
        environment_df$energy[i] <- 1
      }
    }
    else { # otherwise, the alive organism is an animal
      # get neighbors and available cells
      plant_neighbors <- get_plant_neighbors(environment_df$x[i], environment_df$y[i], environment_df, species_list)
      animal_neighbors <- get_animal_neighbors(environment_df$x[i], environment_df$y[i], environment_df, species_list)
      animal_availables <- get_animal_availables(animal_neighbors)
      # feed
      # see if the organism has any neighboring prey
      potential_local_prey <- plant_neighbors
      for (j in 1:nrow(potential_local_prey)) {
        if (!is.na(animal_neighbors$species_id[j])) {
          potential_local_prey[j, ] <- animal_neighbors[j, ]
        }
      }
      potential_local_prey <- potential_local_prey[!is.na(potential_local_prey$species_id) & potential_local_prey$energy > 0, ]
      # see if any neighboring prey are within the organism's feeding range
      potential_global_prey <- which(species_adjacency_matrix[environment_df$species_id[i], ] == 1)
      potential_local_prey <- potential_local_prey[potential_local_prey$species_id %in% potential_global_prey, ]
      
      # if there are no neighboring prey
      if (nrow(potential_local_prey) == 0) {
        if (nrow(animal_availables) > 0) { # if there are unoccupied adjacent cells
          # move to a random adjacent cell
          new_cell_index <- sample(nrow(animal_availables), 1)
          environment_df$x[i] <- animal_availables$x[new_cell_index]
          environment_df$y[i] <- animal_availables$y[new_cell_index]
        }
        else { # if all adjacent cells are occupied
          next # do nothing
        }
      }
      else { # otherwise, try to consume a random prey
        prey_index <- sample(nrow(potential_local_prey), 1)
        prey_id <- potential_local_prey$species_id[prey_index]
        predator_id <- environment_df$species_id[i]
        success_probability <- exp(-((species_df$n_i[prey_id]-species_df$c_i[predator_id])/ (0.5 * species_df$r_i[predator_id]))^2)
        if (runif(1, 0, 1) < success_probability) { # if predation is successful
          # update predator location and energy
          environment_df$x[i] <- potential_local_prey$x[prey_index]
          environment_df$y[i] <- potential_local_prey$y[prey_index]
          environment_df$energy[i] <- environment_df$energy[i] + 0.25 * potential_local_prey$energy[prey_index]
          # update prey energy
          prey_environment_index <- which(environment_df$species_id == potential_local_prey$species_id[prey_index] & environment_df$x == potential_local_prey$x[prey_index] & environment_df$y == potential_local_prey$y[prey_index])
          environment_df$energy[prey_environment_index] <- 0
        }
        else { # if predation is unsuccessful
          next # do nothing
        }
      }
      # reproduce
      animal_neighbors <- get_animal_neighbors(environment_df$x[i], environment_df$y[i], environment_df, species_list)
      animal_availables <- get_animal_availables(animal_neighbors)
      if (environment_df$energy[i] > 1 & nrow(animal_availables) > 0) {
        # randomly choose an unoccupied neighboring cell
        new_cell_index <- sample(nrow(animal_availables), 1)
        
        # create new offspring to be added
        offspring <- c(environment_df$species_id[i],
                       animal_availables$x[new_cell_index],
                       animal_availables$y[new_cell_index],
                       environment_df$energy[i] - 1)
        environment_df <- rbind(environment_df, offspring)
        
        environment_df$energy[i] <- 1
      }
      # expend energy
      metabolism <- species_df$m_i[which(species_df$species_id == environment_df$species_id[i])]
      environment_df$energy[i] <- environment_df$energy[i] - 0.1 * metabolism
    }
  }
  
  # delete organisms with 0 energy and add new organisms to environment_df
  environment_df <- environment_df[environment_df$energy > 0, ]
  
  return(environment_df)
}

# summarize_community(): function to summarize community at a given timestep
# input: n x 4 dataframe for non-negative integer n and list(S x 6 dataframe, S x S matrix) created by initialize_species_set()
# output: list(a x 2 matrix, S x S matrix)
summarize_community <- function(environment_df, species_list) {
  species_df <- species_list[[1]]
  species_adjacency_matrix <- species_list[[2]]
  # get the list of species present in the environment
  present_species <- unique(environment_df$species_id)
  
  # create a dataframe to store the abundance of each species present in the environment
  community_df <- as.data.frame(matrix(NA, length(present_species), 2))
  colnames(community_df) <- c("species_id", "abundance")
  for (i in 1:nrow(community_df)) {
    community_df$species_id[i] <- present_species[i]
    community_df$abundance[i] <- sum(environment_df$species_id == present_species[i])
  }
  # subset the global adjacency matrix to only the species present in the environment
  present_adjacency_matrix <- species_adjacency_matrix[present_species, present_species]
  
  return(list(community_df, present_adjacency_matrix))
}

# visualize_community(): # function to visualize community at a given timestep
# input: n x 4 dataframe for non-negative integer n, list(S x 6 dataframe, S x S matrix) created by initialize_species_set(), and positive integers S, C, replicate, and t
# output: .png file of visualized food web and .png file of environment, saved to working directory
visualize_community <- function(environment_df, species_list, S, C, replicate, t) {
  species_df <- species_list[[1]]
  # add a row in species_df for visualizing unoccupied cells
  species_adjacency_matrix <- species_list[[2]]
  species_edge_list <- as.data.frame(make_Edg(species_adjacency_matrix))
  
  # create the current food web for visualization
  # get the list of species present in the environment
  present_species <- c(unique(environment_df$species_id))
  present_species_df <- species_df[present_species, ]
  present_adjacency_matrix <- species_adjacency_matrix[present_species, present_species]
  present_edge_list <- as.matrix(species_edge_list[species_edge_list$consumer %in% present_species & species_edge_list$resource %in% present_species, ])
  present_node_list <- present_species
  present_graph <- graph.data.frame(present_edge_list,
                                    directed = TRUE,
                                    vertices = present_node_list)
  # organize species' height in the food web by trophic level
  graph_layout <- present_species_df[, c("x_position", "trophic_level")]
  graph_layout$trophic_level <- graph_layout$trophic_level - 1
  graph_layout <- as.matrix(graph_layout)
  # set aesthetic parameters in the food web
  E(present_graph)$color <- "grey"
  V(present_graph)$color <- "goldenrod"
  V(present_graph)$label.cex <- 1
  # export the visualized food web
  food_web_snapshot_dir <- paste("./food_web_snapshots/", S, "_", C, "_", replicate, sep = "")
  png(paste(food_web_snapshot_dir, "/food_web_", t, ".png", sep = ""), 1000, 1000)
  plot(present_graph,
       layout = graph_layout,
       vertex.label = present_node_list,
       vertex.size = 3,
       edge.arrow.size = 1,
       edge.arrow.mode = T,
       edge.width = 0.3,
       xlim = c(0, 1),
       ylim = c(0, ceiling(max(species_df$trophic_level) - 1)),
       rescale = FALSE,
       asp = 0)
  dev.off()

  # add a row in species_df for visualizing unoccupied cells
  species_df <- rbind(species_df,
                      list(0, 0, 0, 0, 0, -1, 0, "#808080"))
  species_df <- species_df[order(species_df$species_id), ]
  # get the list of plant species in the global environment
  plant_species <- unique(species_df$species_id[species_df$trophic_level == 1])
  # split environment_df into a dataframe for plant species and a dataframe for animal species
  plant_environment_df <- environment_df[environment_df$species_id %in% plant_species, ]
  animal_environment_df <- environment_df[!(environment_df$species_id %in% plant_species), ]
  # create the environment matrix for visualization
  environment_matrix <- matrix(0, 100, 100)
  for (i in 1:nrow(plant_environment_df)) {
    x <- plant_environment_df$x[i]
    y <- plant_environment_df$y[i]
    environment_matrix[x, y] <- plant_environment_df$species_id[i]
  }
  for (i in 1:nrow(animal_environment_df)) {
    x <- animal_environment_df$x[i]
    y <- animal_environment_df$y[i]
    environment_matrix[x, y] <- animal_environment_df$species_id[i]
  }
  # add 0.01 to each value in environment_matrix to make it work with how image() maps values to colors
  environment_matrix <- environment_matrix + 0.01
  # export the visualized environment
  environment_snapshot_dir <- paste("./environment_snapshots/", S, "_", C, "_", replicate, sep = "")
  png(paste(environment_snapshot_dir, "/environment_", t, ".png", sep = ""), 1000, 1000)
  image(environment_matrix,
        col = species_df$color,
        breaks = 0:length(species_df$color),
        xaxt = "n",
        yaxt = "n")
  dev.off()
}
#-----------------------------------------------------------------------------#
# master function to run the simulation, store results, and export visualizations

# run_simulation(): master function that generates the simulation
# input: Number of species S, number of basal species B, connectance C for positive integers S and B, and C in [0,1], positive integer t_limit for number of timesteps, and positive integer snapshot_increment for snapshot frequency
# output: list(list((t_limit/snapshot_increment + 1) n x 2 dataframes), list((t_limit/snapshot_increment + 1) m x m dataframes)) for non-zero integers n and m
run_simulation <- function(S, B, C, replicate, t_limit, snapshot_increment) {
  start_time <- Sys.time()
  food_web_snapshot_dir <- paste("./food_web_snapshots/", S, "_", C, "_", replicate, sep = "")
  dir.create(food_web_snapshot_dir)
  environment_snapshot_dir <- paste("./environment_snapshots/", S, "_", C, "_", replicate, sep = "")
  dir.create(environment_snapshot_dir)
  
  # initialize the set of species in the global environment
  global_species_set <- initialize_species_set(S, B, C)
  # initialize the local environment
  local_environment <- initialize_environment(global_species_set)
  # initialize lists to save snapshot community_dfs and adjacency matrices
  community_df_list <- vector(mode = "list", length = t_limit + 1)
  adjacency_matrix_list <- vector(mode = "list", length = t_limit + 1)
  community_df_list[[1]] <- summarize_community(local_environment, global_species_set)[[1]]
  adjacency_matrix_list[[1]] <- summarize_community(local_environment, global_species_set)[[2]]
  visualize_community(local_environment, global_species_set, S, C, replicate, 1)
  
  for (i in 1:t_limit) { # increment t_limit timesteps
    # increment the local environment by 1 timestep
    local_environment <- increment_timestep(local_environment, global_species_set)
    
    if (i %% snapshot_increment == 0) { # if a snapshot is required for this timestep
      print(paste("iterations completed: ", i, sep = ""))
      community_summary <- summarize_community(local_environment, global_species_set)
      # add community_df snapshot to community_df list
      community_df_list[[i+1]] <- community_summary[[1]]
      # add adjacency matrix to adjacency matrix list
      adjacency_matrix_list[[i+1]] <- community_summary[[2]]
      # generate visual snapshots
      visualize_community(local_environment, global_species_set, S, C, replicate, i+1)
    }
  }
  community_df_list <- community_df_list[lengths(community_df_list) != 0]
  adjacency_matrix_list <- adjacency_matrix_list[lengths(adjacency_matrix_list) != 0]
  
  end_time <- Sys.time()
  print(end_time - start_time)
  
  return(list(community_df_list, adjacency_matrix_list))
}

#-----------------------------------------------------------------------------#
# function to store and analyze results from simulations

# create_dfs(): function to summarize results of a simulation and organize them into two dataframes (one for network metrics over time, one for community abundances over time)
# input: Number of species S, snapshot_increment, list(list((t_limit/snapshot_increment + 1) n x 2 dataframes), list((t_limit/snapshot_increment + 1) m x m dataframes)) created by run_simulation()
# output: list((t_limit/snapshot_increment + 1) x 7 dataframe, (t_limit/snapshot_increment + 1) x (S + 1) dataframe)
create_dfs <- function(S, snapshot_increment, simulation_run) {
  community_df_list <- simulation_run[[1]]
  adjacency_matrix_list <- simulation_run[[2]]
  
  # create a dataframe tracking network metrics over time
  metrics_df <- as.data.frame(matrix(NA, nrow = length(adjacency_matrix_list), ncol = 6))
  colnames(metrics_df) <- c("t", "species", "connectance", "top", "intermediate", "basal")
  for (i in 1:nrow(metrics_df)) {
    metrics_df$t[i] <- (i - 1) * snapshot_increment + 1
    adjacency_matrix <- adjacency_matrix_list[[i]]
    if (is.null(nrow(adjacency_matrix))) {
      metrics_df[i, c(2:6)] <- c(1, 0, 0, 0, 0)
    }
    else {
      metrics_df[i, c(2:6)] <- c(species(adjacency_matrix),
                                 connectance(adjacency_matrix),
                                 top(adjacency_matrix),
                                 intermediate(adjacency_matrix),
                                 basal(adjacency_matrix))
    }
  }
  
  # create a dataframe tracking species abundances over time
  community_df <- as.data.frame(matrix(0, nrow = length(community_df_list), ncol = S + 2))
  colnames(community_df) <- c("t", c(1:S), "total")
  for (i in 1:nrow(community_df)) {
    community_df$t[i] <- (i - 1) * snapshot_increment + 1
    community_abundance <- community_df_list[[i]]
    for (j in 1:nrow(community_abundance)) {
      species_id <- community_abundance$species_id[j]
      species_abundance <- community_abundance$abundance[j]
      community_df[i, species_id + 1] <- species_abundance
    }
    community_df$total[i] <- sum(community_df[i, c(2:S+1)])
  }
  
  return(list(metrics_df, community_df))
}

# run_experiment(): function to create a factorial experimental design of simulations, run specified simulations, and return all results in a dataframe
# input: number_of_replicates, vector containing values for all numbers of species S in web, and vector containing valuse for all expected connectances C of web
# output: dataframe containing results for all combinations of specified inputs
run_experiment <- function(number_of_replicates, S, C) {
  replicate <- 1:number_of_replicates
  # create factorial dataframe of different web parameter combinations
  df <- expand_grid(S, C, replicate)
  
  # add additional columns to df
  df$community_df_list <- NA
  df$adjacency_matrix_list <- NA
  df$metrics_df <- NA
  df$community_df <- NA
  
  # run simulations specified in df
  for (i in 1:nrow(df)) {
    if (i %% 10 == 0) {print(paste("simulations completed: ", i, "/", nrow(df), sep = ""))}
    simulation <- run_simulation(df$S[i],
                                 ceiling(0.1 * df$S[i]),
                                 df$C[i],
                                 df$replicate[i],
                                 1000,
                                 10)
    df$community_df_list[i] <- list(simulation[[1]])
    df$adjacency_matrix_list[i] <- list(simulation[[2]])
    simulation_dfs <- create_dfs(df$S[i],
                                 10,
                                 simulation)
    df$metrics_df[i] <- list(simulation_dfs[[1]])
    df$community_df[i] <- list(simulation_dfs[[2]])
  }
  
  return(df)
}