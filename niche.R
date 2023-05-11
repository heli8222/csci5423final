# Niche model
# returns a model-generated food web with species S and connectance C
# see Williams and Martinez (2000) for model
# Williams and Martinez (2000): 10.1038/35004572

# load packages
library(igraph)
library(doBy)
library(tidyverse)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#-----------------------------------------------------------------------------#
# define helper functions

# create_species(): function to create and parameterize a species by providing beta for Gamma distribution
# input: positive real number beta
# output: vector of length 3
create_species <- function(beta) {
  new_species <- numeric(length = 3) # initialize a vector to describe a new species
  new_species[1] <- runif(1, min = 0, max = 1) # draw niche value n_i from Unif[0,1] and store in vector
  y <- runif(1, min = 0, max = 1) # set up to draw fundamental generality range r_i
  x <- 1 - (1 - y)^(1 / beta) # set up to draw fundamental generality range r_i
  new_species[2] <- x * new_species[1] # draw r_i, see Williams and Martinez (2000) for derivation
  new_species[3] <- runif(1, min = new_species[2] / 2, max = new_species[1]) # draw range center c_i from Unif[r_i/2, n_i] and store in vector
  return(new_species)
}

# make_Par(): function to construct a parameter matrix of n species
# input: positive integer S, positive real number beta
# output: n x 3 matrix
make_Par <- function(S, beta) {
  parameter_matrix <- matrix(NA, nrow = S, ncol = 3) # initialize a matrix to store species parameters
  colnames(parameter_matrix) <- c("n_i", "r_i", "c_i") # name matrix columns
  for (i in 1:S) { # populate matrix by looping through each species
    parameter_matrix[i,] <- t(create_species(beta)) # for each row create a new species
  }
  return(parameter_matrix)
}

# assign_basal(): function to assign B basal species to parameter matrix
# input: S x 3 matrix, B for positive integers S and B
# output: S x 3 matrix
assign_basal <- function(parameter_matrix, B) {
  basal_species <- which.minn(parameter_matrix[,1], B) # initialize a vector to store the B species with lowest n_i values
  for (i in 1:B) { # set r_i to 0 for all species in this vector
    parameter_matrix[basal_species[i],2] <- 0
  }
  return(parameter_matrix)
}

# make_Adj(): function to construct an adjacency matrix from a parameter matrix
# input: S x 3 matrix for positive integer S
# output: S x S matrix
make_Adj <- function(parameter_matrix) {
  S <- nrow(parameter_matrix) # store the number of species
  adjacency_matrix <- matrix(NA, nrow = S, ncol = S) # initialize a matrix to store network edges
  for (i in 1:S) { # calculate feeding ranges for each species
    i_lower <- parameter_matrix[i,3] - parameter_matrix[i,2] / 2 # define lower boundary of a species' feeding range
    i_upper <- parameter_matrix[i,3] + parameter_matrix[i,2] / 2 # define upper boundary of a species' feeding range
    i_eats <- integer(length = S) # initialize a vector to store network edges for a species
    for (j in 1:S) { # populate vector by looping through each species
      if (between(parameter_matrix[j,1], i_lower, i_upper)) { # determine whether a species' n_i falls within the feeding range
        i_eats[j] <- 1
      }
    }
    adjacency_matrix[i,] <- i_eats # populate matrix by inserting vector for each species
  }
  return(adjacency_matrix)
}

# make_Edg(): function to construct an edge list from an adjacency matrix
# input: S x S matrix for positive integer S
# output: a x 2 matrix for positive integer a
make_Edg <- function(adjacency_matrix) {
  n <- sum(adjacency_matrix) # store the number of edges in an adjacency matrix
  edge_list <- matrix(NA, nrow = n, ncol = 2) # initialize a matrix to store network edges
  colnames(edge_list) <- c("consumer", "resource") # name matrix columns
  row_index <- 1 # initialize a counter to store the ith interaction to the ith row of the edge list
  for (i in 1:nrow(adjacency_matrix)) { # loop through all rows of adjacency matrix
    for (j in 1:ncol(adjacency_matrix)) { # loop through all columns of adjacency matrix
      if (adjacency_matrix[i,j] == 1) { # add i,j to the edge list if species i eats species j
        edge_list[row_index, 1] <- i
        edge_list[row_index, 2] <- j
        row_index <- row_index + 1 # increment row_index to the next row of edge list
      }
    }
  }
  return(edge_list)
}

# remove_species(): function to remove a species from a parameter matrix
# input:  S x 3 matrix, i for positive integers S and i
# output: S-1 x 3 matrix
remove_species <- function(parameter_matrix, i) {
  parameter_matrix <- parameter_matrix[-i,] # remove the ith species from parameter matrix
  return(parameter_matrix)
}

# add_species(): function to add a species to a parameter matrix
# input: S x 3 matrix, beta for positive integer S and positive real number beta
# output: S+1 x 3 matrix
add_species <- function(parameter_matrix, beta) {
  new_species <- c(0,0,0) # initialize vector to describe new species
  basal_species <- which(parameter_matrix[,2] == 0)
  max_basal_ni <- max(parameter_matrix[basal_species,1])
  while (new_species[1] < max_basal_ni) { # generate a new species that has a higher n_i than all existing basal species
    new_species <- t(create_species(beta))
  }
  parameter_matrix <- rbind(parameter_matrix, new_species) # add new species to parameter matrix
  row.names(parameter_matrix) <- rep("", nrow(parameter_matrix)) # get rid of weird row name that shows up after adding the new species
  return(parameter_matrix)
}

# isolates(): function to detect identity of single isolated species in adjacency matrix
# input: S x S matrix for positive integer S
# output: vector of length a for non-negative integer a
isolates <- function(adjacency_matrix) {
  isolated_species <- numeric(length = 0) # initialize a vector to store indices of isolated species
  for (i in 1:nrow(adjacency_matrix)) { # loop through all species to identify isolated species
    if (sum(adjacency_matrix[i,]) == 0 & sum(adjacency_matrix[,i]) == 0) { # add species to vector if it eats no one and no one eats it
      isolated_species <- append(isolated_species, i)
    }
  }
  return(isolated_species)
}

# identicals(): function to detect identity of single trophically identical species in adjacency matrix
# input: S x S matrix for positive integer S
# output: vector of length a for non-negative integer a
identicals <- function(adjacency_matrix) {
  identical_species <- numeric(length = 0) # initialize a vector to store indices of identical species
  for (i in 1:nrow(adjacency_matrix)) { # loop through all species-pairs to identify identical species
    i_eats <- which(adjacency_matrix[i,] == 1) # initialize a vector to store all species that species i eats
    i_eaten <- which(adjacency_matrix[,i] == 1) # initialize a vector to store all species that eat species i
    for (j in 1:nrow(adjacency_matrix)) {
      if (j != i){
        j_eats <- which(adjacency_matrix[j,] == 1) # initialize a vector to store all species that species j eats
        j_eaten <- which(adjacency_matrix[,j] == 1) # initialize a vector to store all species that eats species j
        if (identical(i_eats, j_eats) & identical(i_eaten, j_eaten)) { # add species j to vector if it has the same eats and eaten vectors as species i
          identical_species <- append(identical_species, j)
        }
      }
    }
  }
  return(identical_species)
}

# lower_triangular(): function to make an acyclic adjacency matrix lower triangular
# input: S x S matrix for positive integer S
# output: S x S matrix
lower_triangular <- function(adjacency_matrix) {
  graph <- graph_from_adjacency_matrix(adjacency_matrix) # convert adjacency matrix into igraph object
  graph <- topo_sort(graph, mode = "in") # sort adjacency matrix based on incoming edges
  return(adjacency_matrix[graph, graph])
}

# DFS(): function to decycle an adjacency matrix
# input: S x S matrix for positive integer S
# output: S x S matrix
# taken from threat_functions.R by Gyorgy Barabas
DFS <- function(adjacency_matrix) {
  DFSCOLOR <- numeric(0)
  DFSBACKEDGE <- numeric(0)
  ORDERVERTICES <- numeric(0)
  DFSVisit <- function(adjacency_matrix, i) {
    DFSCOLOR[i] <<- 1
    for (j in 1:nrow(adjacency_matrix)) {
      if(adjacency_matrix[i,j] != 0) {
        if (DFSCOLOR[j] == 0) {
          DFSVisit(adjacency_matrix,j)
        } else {
          if(DFSCOLOR[j] == 1) {
            DFSBACKEDGE[i,j] <<- 1 # It's a back edge: list for removal
          }
        }
      }
    }
    DFSCOLOR[i] <<- 2
    ORDERVERTICES <<- c(i, ORDERVERTICES)
  }
  run_DFS <- function(adjacency_matrix) {
    S <- nrow(adjacency_matrix)
    DFSCOLOR <<- rep(0, S)
    DFSBACKEDGE <<- matrix(0, S, S)
    ORDERVERTICES <<- numeric(0)
    for (i in 1:S) if (DFSCOLOR[i] == 0) DFSVisit(adjacency_matrix, i)
    return(adjacency_matrix - DFSBACKEDGE)
  }
  return(run_DFS(adjacency_matrix))
}

# trophic_levels(): function to calculate the prey averaged trophic level for each species from an adjacency matrix
# input: S x S matrix for positive integer S
# output: S x 1 vector
# adapted from threat_functions.R by Gyorgy Barabas
trophic_levels <- function(adjacency_matrix) {
  A <- adjacency_matrix
  A <- A / rowSums(A)
  A[is.nan(A)] <- 0
  S <- nrow(A)
  tl <- solve(diag(S) - A, rep(1, S))
  return(tl)
}

#-----------------------------------------------------------------------------#
# make_web(): master function that generates niche model network
# input: Number of species S, number of basal species B, connectance C for positive integers S and B, and C in [0,1]
# output: list(a x 2 matrix, S x S matrix, S x 3 matrix) for positive integer a
make_web <- function(S, B, C) {
  beta <- (1 - 2*C) / (2*C) # store beta parameter of Gamma distribution
  parameter_matrix <- make_Par(S, beta) # make a parameter matrix with S species
  parameter_matrix <- assign_basal(parameter_matrix, B) # assign B basal species to the parameter matrix
  adjacency_matrix <- make_Adj(parameter_matrix) # make an adjacency matrix using the parameter matrix
  problem_species <- c(isolates(adjacency_matrix), identicals(adjacency_matrix)) # initialize a vector of all isolated or identical species
  
  while (length(problem_species) > 0) { # while there are isolated or identical species, adjust the parameter matrix
    i <- problem_species[1] # identify one such isolated or identical species
    parameter_matrix <- remove_species(parameter_matrix, i) # remove this species from the parameter matrix
    parameter_matrix <- add_species(parameter_matrix, beta) # add a new species to the parameter matrix
    parameter_matrix <- assign_basal(parameter_matrix, B) # reassign basal species to ensure that there are at least B basal species
    adjacency_matrix <- make_Adj(parameter_matrix) # rebuild the adjacency matrix based on the adjusted parameter matrix
    problem_species <- c(isolates(adjacency_matrix), identicals(adjacency_matrix)) # re-evaluate the vector of all isolated or identical species
  }
  
  adjacency_matrix <- DFS(adjacency_matrix) # remove all cycles from the adjacency matrix
  edge_list <- make_Edg(adjacency_matrix) # make the edge list based on the adjacency matrix
  return(list(edge_list, adjacency_matrix, parameter_matrix))
}

#-----------------------------------------------------------------------------#
# Output
# S <- 25 # input number of species
# B <- 1 # input number of basal species
# C <- .32 # input connectance of web
# 
# web <- make_web(S, B, C)
# write.csv(web[[1]], "edge_list.csv", row.names = FALSE)
# write.csv(web[[2]], "adjacency_matrix.csv", row.names = FALSE)
# write.csv(web[[3]], "parameter_matrix.csv", row.names = FALSE)
# 
# # visualize adjacency matrix
# image(t(apply(web[[2]],2,rev)))