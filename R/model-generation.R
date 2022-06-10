#' Create populations based on given tree
#'
#' Creates a list of slendr populations according to a given tree. Starting from
#' the root, the tree is traversed recursively and a new population is created
#' for every left child of a given node.
#'
#' The simulation length parameter is used to scale the length of the edges of
#' the tree. The last population is created at 3/4 of the simulation length, so
#' that it will exist for a bit before the simulation ends. The root of the tree
#' is the ancestral population, which is created at time 1.
#'
#' @param tree An ape tree
#' @param population_size Integer number defining the populations size of all
#'   populations
#' @param simulation_length Integer number defining on how long the simulation
#'   of the populations will last, to scale the edge lengths of the tree
#'   accordingly
#' @return List of slendr populations
#' @examples
#' tree <- ape::rtree(4)
#' tree_populations(tree = tree, population_size = 1000, simulation_length = 50)
#' tree <- ape::rtree(6)
#' tree_populations(tree = tree, population_size = 200, simulation_length = 500)
tree_populations <- function(tree, population_size, simulation_length) {

  if (length(tree$tip.label) + tree$Nnode < 3) {
    stop("At least 2 populations must be created", call. = FALSE)
  }

  env <- new.env()
  env$populations <- vector(mode="list", length=tree$Nnode +
                              length(tree$tip.label))

  # scale defined by maximum edge length of the tree and given simulation length
  # only 3/4 of this scale used, so that all populations exist at least for 1/4
  # of the simulation
  scale <- floor(3/4*simulation_length/max(ape::node.depth.edgelength(tree)))

  # create ancestral population
  root <- length(tree$tip.label) + 1
  env$populations[[root]] <- population(paste0("pop", root), time = 1,
                                        N = population_size)

  # recurse through tree
  recursion(tree, root, env$populations[[root]], population_size, env, scale)

  # remove empty populations from list
  env$populations <- env$populations[-which(sapply(env$populations, is.null))]

  return(env$populations)
}

#' Create random populations
#'
#' Creates a random list of slendr populations of given length. A random tree is
#' created and the function \code{tree_populations} is called to create the
#' populations.
#'
#' @param n_populations The desired number of populations that should be created
#' @param population_size Integer number defining the populations size of all
#'   populations
#' @param simulation_length Integer number defining on how long the simulation
#'   of the populations will last, to scale the edge lengths of the tree
#'   accordingly
#' @return List of slendr populations
#' @examples
#' tree_populations(n_populations = 4, population_size = 1000,
#'   simulation_length = 50)
#' tree_populations(n_populations = 6, population_size = 200,
#'   simulations_length = 500)
random_populations <- function(n_populations, population_size,
                               simulation_length) {

  if (n_populations < 2) {
    stop("At least 2 populations must be created", call. = FALSE)
  }

  tree <- ape::rtree(n_populations)
  populations <- tree_populations(tree, population_size, simulation_length)
  return(populations)
}

#' Create model based on given tree
#'
#' Creates a slendr model based on a given tree and a given number of gene flow
#' events. The function \code{tree_populations} is used to get a list of
#' populations from given tree. Gene flow events between these populations are
#' created randomly.
#'
#' @param tree An ape tree
#' @param population_size Integer number defining the populations size of all
#'   populations
#' @param n_gene_flow Integer number defining the number of gene_flow events
#' @param rate_gene_flow Vector that specifies the min and max gene flow rate
#' @param simulation_length Integer number defining on how long the simulation
#'   of the populations will last, to scale the edge lengths of the tree
#'   accordingly
#' @return A slendr model
#' @examples
#' tree <- ape::rtree(4)
#' tree_model(tree = tree, population_size = 1000, n_gene_flow = 3,
#'   rate_gene_flow = 0.5, simulation_length = 100)
#' tree <- ape::rtree(6)
#' tree_model(tree = tree, population_size = 200, n_gene_flow = 4,
#'   rate_gene_flow = c(0.2, 0.9), simulation_length = 1000)
tree_model <- function(tree, population_size, n_gene_flow,
                       rate_gene_flow = c(0.01, 0.99), simulation_length) {

  populations <- tree_populations(tree, population_size, simulation_length)
  gf <- random_gene_flow(populations, n_gene_flow, rate_gene_flow,
                         simulation_length)
  model <- compile_model(populations = populations, gene_flow = gf,
                         generation_time = 1, sim_length = simulation_length)
  return(model)
}

#' Create random model
#'
#' Creates a slendr model based on a given number of populations and a given
#' number of gene flow events. A random ape tree is created according to the
#' specified number of populations. With this tree, the function
#' \code{tree_model} is called to create the model.
#'
#' @param n_populations The desired number of populations for the model
#' @param population_size Integer number defining the populations size of all
#'   populations
#' @param n_gene_flow Integer number defining the number of gene_flow events
#' @param rate_gene_flow Vector that specifies the min and max gene flow rate
#' @param simulation_length Integer number defining on how long the simulation
#'   of the populations will last, to scale the edge lengths of the tree
#'   accordingly
#' @return A slendr model
#' @examples
#' tree_model(n_populations = 4, population_size = 1000, n_gene_flow = 3,
#'   rate_gene_flow = 0.5, simulation_length = 100)
#' tree_model(n_populations = 6, population_size = 200, n_gene_flow = 4,
#'   rate_gene_flow = c(0.2, 0.9), simulation_length = 1000)
random_model <- function(n_populations, population_size, n_gene_flow,
                         rate_gene_flow = c(0.01, 0.99), simulation_length) {

  tree <- ape::rtree(n_populations)
  model <- tree_model(tree, population_size, n_gene_flow,
                      rate_gene_flow, simulation_length)
  return(model)
}

recursion <- function(tree, current_node, parent_population, population_size,
                      env, scale) {

  edge_list <- tree$edge
  children <- edge_list[edge_list[, 1] == current_node, 2]
  left <- children[1]
  right <- children[2]

  if (!is.na(left)){
    # adjust the edge length of node by the scale
    length <- ceiling(ape::node.depth.edgelength(tree)[left] * scale)
    # check that parent and child can not have the same length
    if(length == attr(parent_population, "history")[[1]]$time) {
      length <- length + 1
    }

    env$populations[[left]] <- population(paste0("pop", left), time = length,
                                          N = population_size,
                                          parent = parent_population)

    recursion(tree, left, env$populations[[left]], population_size, env, scale)
  }
  if (!is.na(right)){
    recursion(tree, right, parent_population, population_size, env, scale)
  }
}

random_gene_flow <- function(populations, n, rate, simulation_length) {
  gf <- vector(mode = "list", length = n)
  if (n != 0){
    for(i in 1:n) {

      number1 <- sample(1:length(populations), 1)
      number2 <- sample(1:length(populations), 1)
      # check that gene flow not between same population
      while(number1 == number2){
        number2 <- sample(1:length(populations), 1)
      }
      pop1 <- populations[[number1]]
      pop2 <- populations[[number2]]
      # find out which population was created later to use as possible start for
      # gene flow
      max_time <- max(attr(pop1, "history")[[1]]$time,
                    attr(pop2, "history")[[1]]$time)
      start_gf <- sample((max_time + 1):(simulation_length - 2), 1)
      end_gf <- sample((start_gf + 1):(simulation_length), 1)

      # given rate of gene flow can be either list of min/max or single value
      if (length(rate) == 2){
        if (rate[1] > rate[2]) {
          stop("No valid range for gene flow rate", call. = FALSE)
        }
        rate_gf <- runif(1, rate[1], rate[2])
      } else {
        rate_gf <- rate
      }
      gf[[i]] <- gene_flow(from = pop1, to = pop2, start = start_gf,
                           end = end_gf, rate_gf)
    }
  }
  return(gf)
}