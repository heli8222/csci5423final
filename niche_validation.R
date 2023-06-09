# Niche model validation
# compare the properties of food webs generated by niche.R to the properties of food webs generated by the niche model in Williams and Martinez (2000)
# Williams and Martinez (2000): 10.1038/35004572

source("./niche.R")
#-----------------------------------------------------------------------------#
# define functions to calculate web properties

species <- function(adjacency_matrix) {
  return(nrow(adjacency_matrix))
}

connectance <- function(adjacency_matrix) {
  return(sum(adjacency_matrix) / (nrow(adjacency_matrix) - 1)^2)
}

# T(): function to calculate the fraction of top species within a web
# input: adjacency matrix
# output: non-negative real number
top <- function(adjacency_matrix) {
  S <- nrow(adjacency_matrix) # store the number of species
  counter <- 0 # initialize counter for number of top species
  for (i in 1:S) { # loop through all species in web
    if (sum(adjacency_matrix[,i]) == 0) { # if no one eats species i, add it to counter
      counter <- counter + 1
    }
  }
  return(counter / S)
}

# I(): function to calculate the fraction of intermediate species within a web
# input: adjacency matrix
# output: non-negative real number
intermediate <- function(adjacency_matrix) {
  S <- nrow(adjacency_matrix) # store the number of species
  counter <- 0 # initialize counter for number of intermediate species
  for (i in 1:S) { # loop through all species in web
    if (sum(adjacency_matrix[i,]) > 0 & sum(adjacency_matrix[,i]) > 0) { # if species i eats someone and someone eats species i, add it to counter
      counter <- counter + 1
    }
  }
  return(counter / S)
}

# B(): function to calculate the fraction of basal species within a web
# input: adjacency matrix
# output: non-negative real number
basal <- function(adjacency_matrix) {
  S <- nrow(adjacency_matrix) # store the number of species
  counter <- 0 # initialize counter for number of basal species
  for (i in 1:S) { # loop through all species in web
    if (sum(adjacency_matrix[i,]) == 0) { # if species i eats no one, add it to counter
      counter <- counter + 1
    }
  }
  return(counter / S)
}