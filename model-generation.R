library(phangorn)
library(ape)
library(phytools)


devtools::load_all(".")

#' Create a list of slendr populations according to a given tree. Starting from
#' the root, the tree is traversed recursively and a new population is created
#' for every left child of a given node.
#'
#' todo explain in more detail (z.b. 3/4)
#'
#' @param tree An ape tree
#' @param population_size Integer number defining the populations size of all
#'   populations
#' @param simulation_length Integer number defining on how long the simulation
#'   of the populations will last, to scale the edge lengths of the tree
#'   accordingly
#' @return List of slendr populations
#' @examples
#' tree <- rtree(4)
#' tree_populations(tree, 1000, 50)
#' tree <- rtree(6)
#' tree_populations(tree, 200, 500)
tree_populations <- function(tree, population_size, simulation_length) {
  # plot to check
  plot(tree, show.tip.label = F)
  nodelabels()
  tiplabels()

  env <- new.env()
  env$populations <- vector(mode="list", length=tree$Nnode +
                              length(tree$tip.label))
  ## TODO explain
  scale <- floor((3*simulation_length/4)/max(node.depth.edgelength(tree)))

  root <- length(tree$tip.label) + 1
  env$populations[[root]] <- population(paste0("pop", root), time = 1,
                                        N = population_size)

  recursion(tree, root, env$populations[[root]], population_size, env, scale)
  env$populations <- env$populations[-which(sapply(env$populations, is.null))]
  return(env$populations)
}

#' Create a random list of slendr populations of given length. A random tree is
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
#' tree_populations(4, 1000, 50)
#' tree_populations(6, 200, 500)
random_populations <- function(n_populations, population_size, simulation_length) {
  # todo: sanity check: if (npops < 1) stop(...), stop("blah blah", call. = FALSE)
  tree <- rtree(n_populations)
  populations <- tree_populations(tree, population_size, simulation_length)
  return(populations)
}

#' Create a slendr model based on a given tree and a given number of gene flow
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
#' tree <- rtree(4)
#' tree_model(tree, 1000, 3, 0.5, 100)
#' tree <- rtree(6)
#' tree_model(tree, 200, 4, c(0.2, 0.9), 1000)
tree_model <- function(tree, population_size, n_gene_flow,
                       rate_gene_flow = c(0.01, 0.99), simulation_length) {

  populations <- tree_populations(tree, population_size, simulation_length)
  gf <- random_gene_flow(populations, n_gene_flow, rate_gene_flow,
                         simulation_length)
  model <- compile_model(populations = populations, gene_flow = gf,
                         generation_time = 1, sim_length = simulation_length)
  return(model)
}

#' Create a slendr model based on a given number of populations and a given
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
#' tree_model(4, 1000, 3, 0.5, 100)
#' tree_model(6, 200, 4, c(0.2, 0.9), 1000)
random_model <- function(n_populations, population_size, n_gene_flow,
                         rate_gene_flow = c(0.01, 0.99), simulation_length) {

  tree <- rtree(n_populations)
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
    length <- ceiling(node.depth.edgelength(tree)[left] * scale)
    if(length == attr(parent_population, "history")[[1]]$time) {
      length <- length + 1
    }
    env$populations[[left]] <- population(paste0("pop", left), time = length,
                                          N = population_size,
                                          parent = parent_population)

    recursion(tree, left, env$populations[[left]], population_size, env, scale)
  }
  if (!is.na(right)){
    length <- floor(node.depth.edgelength(tree)[right] * scale) + 1
    #env$populations[[right]] <- population(paste0("pop", right), time = length,
    #N = population_size, parent = env$populations[[current_node]])
    recursion(tree, right, parent_population, population_size, env, scale)
  }
}

random_gene_flow <- function(populations, n, rate, simulation_length) {
  gf <- vector(mode = "list", length = n)
  if (n != 0){
    for(i in 1:n) {

      number1 <- sample(1:length(populations), 1)
      number2 <- sample(1:length(populations), 1)
      while(number1 == number2){
        number2 <- sample(1:length(populations), 1)
      }
      pop1 <- populations[[number1]]
      pop2 <- populations[[number2]]
      max_time <- max(attr(pop1, "history")[[1]]$time,
                    attr(pop2, "history")[[1]]$time)
      start_gf <- sample((max_time + 1):(simulation_length - 2), 1)
      end_gf <- sample((start_gf + 1):(simulation_length), 1)
      if (length(rate) == 2){
        rate_gf <- runif(1, rate[1], rate[2])
      } else {
        rate_gf <- rate
      }
      gf[[i]] <- gene_flow(from = pop1, to = pop2, start = start_gf, end = end_gf, rate_gf)
    }
  }
  return(gf)
}

# test calls
tree <- rtree(4)
pops <- tree_populations(tree, 1000, 50)
pops2 <- random_populations(3, 1000, 100)
model <- tree_model(tree, 1000, 3, c(0.2, 0.9), 100)
model2 <- random_model(3, 1000, 2, 0.5, 1000)
plot_model(model, proportions = TRUE)
plot_model(model2, proportions = TRUE)

