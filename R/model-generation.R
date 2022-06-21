# variable for the generation of the first split, defining the length of the
# ancestral population existing by itself
burn_in <- 2

#' Create populations based on given tree
#'
#' Creates a list of slendr populations according to a given tree. An ancestral
#' population is created, that exists for a specific burn-in time. Then the
#' populations specified by the tree are added. The root defines the first split
#' of a population from the ancestral population. It is created at the burn-in
#' time. The creation times for all other populations are sampled in the
#' beginning. They are randomly chosen in the interval of one generation after
#' the burn-in until three generations before the end of the simulation (to
#' allow gene flow). Starting from the root, the tree is traversed recursively
#' and a new population is created for every left child of a given node. The
#' right branch is the branch where the parent population continues. All
#' populations continue to exist from the time they were created.
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
tree_populations <- function(tree, population_size, simulation_length,
                             use_tree_lengths = FALSE) {

  if (length(tree$tip.label) + tree$Nnode < 3) {
    stop("At least 2 populations must be created", call. = FALSE)
  }
  if (use_tree_lengths) {
    # get a list of the depths of all nodes in the tree
    depths <- ape::node.depth.edgelength(tree)
    # calculate a scale to adjust the depths so that the one with the max depth
    # is created 3 generations before the simulation ends (to allow gene flow)
    scale <- floor((simulation_length-3)/max(depths))
    # save the creation times as the depths adjusted by the scale
    creation_times <- as.integer(depths * scale)
  } else {
    # sample creation times
    # first possible creation time is 1 generation after the burn-in time
    # last possible creation time is 3 generation before the simulation ends
    # (to leave time for one gene flow event)
    # the number of creation times needed is defined by the amount of nodes in the
    # tree - 1, because the creation time of the root is already specified by the
    # burn-in
    # replace is set to false as only one population can be created at given time
    creation_times <- sample ((burn_in + 1):(simulation_length - 3),
                              size = tree$Nnode + length(tree$tip.label) - 1,
                              replace = FALSE)
    # the creation time for the root is added
    creation_times <- append(creation_times, burn_in)

    # the idea now is that the node in the tree with the smallest depth has the
    # the creation times are sorted to know which is the earliest
    sorted_creation_times <- sort(creation_times)
    # get the order of the depths of the node
    order_depths <- order(ape::node.depth.edgelength(tree))
    # the order can be used to move the sorted creation times into the position
    # corresponding to the ordered node depths
    # the smallest creation time is now at the index of the population with the
    # smallest depth in the tree
    creation_times <- sorted_creation_times[order(order_depths)]
  }

  # create a new environment object (needed to adjust populations in recursion)
  env <- new.env()
  # create an empty list of populations within the new environment
  env$populations <- vector(mode="list", length=tree$Nnode +
                              length(tree$tip.label) + 1)

  # create the ancestral population in generation 1
  env$populations[[length(env$populations)]] <- population(paste0("pop", 0),
                                                           time = 1,
                                                           N = population_size)
  # create a new population for first split at root node
  # creation time is the burn-in variable, the parent is the ancestral population
  root <- length(tree$tip.label) + 1
  env$populations[[root]] <- population(paste0("pop", root), time = burn_in,
                                        N = population_size,
                                        parent = env$populations[[length(env$populations)]])

  # create list of parent populations for left and right child
  # left child should have root population as parent
  # right child should ancestral population as parent
  parent_populations <- list(env$populations[[root]],
                             env$populations[[length(env$populations)]])
  # recurse through tree
  recursion(tree = tree, current_node = root,
            parent_populations = parent_populations,
            population_size = population_size,
            env = env,
            creation_times = creation_times)

  # remove empty populations from list
  env$populations <- env$populations[-which(sapply(env$populations, is.null))]

  # return the list of populations
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

  # create a random tree with the number of leaves equal to the number of
  # populations
  tree <- ape::rtree(n_populations)
  # get the list of populations by calling tree_populations based on the tree
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
                       rate_gene_flow = c(0.01, 0.99), simulation_length,
                       use_tree_lengths = FALSE) {

  # get the list of populations by calling the tree_populations function
  populations <- tree_populations(tree, population_size, simulation_length,
                                  use_tree_lengths)
  # get the list of gene flow events by calling the random_gene_flow function
  gf <- random_gene_flow(populations, n_gene_flow, rate_gene_flow,
                         simulation_length)
  # compile the model based on populations, gene flow events and simulation length
  # the generation time is set to 1
  model <- compile_model(populations = populations, gene_flow = gf,
                         generation_time = 1, sim_length = simulation_length)
  # return the created model
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

  # create a random tree with the number of leaves equal to the number of
  # populations
  tree <- ape::rtree(n_populations)
  # create the model by calling the tree_model function with the created tree
  model <- tree_model(tree, population_size, n_gene_flow,
                      rate_gene_flow, simulation_length)
  # return the created model
  return(model)
}

recursion <- function(tree, current_node, parent_populations, population_size,
                      env, creation_times) {
  # find children of current node
  edge_list <- tree$edge
  children <- edge_list[edge_list[, 1] == current_node, 2]
  left <- children[1]
  right <- children[2]

  # if the left child is an internal node
  if (left > length(tree$tip.label)){
    # create a new population for the left child
    env$populations[[left]] <- population(paste0("pop", left),
                                          time = creation_times[[left]],
                                          N = population_size,
                                          parent = parent_populations[[1]])

    # update list of parent populations for left and right child
    # left child should have newly created population as parent
    # right child should have parent of newly created population as parent
    parent_populations <- list(env$populations[[left]], parent_populations[[1]])
    # continue to recurse through the tree
    recursion(tree = tree, current_node = left,
              parent_populations = parent_populations,
              population_size = population_size,
              env = env,
              creation_times = creation_times)
  }
  # if the right child is an internal node
  if (right > length(tree$tip.label)){
    # create a new population for the right child
    env$populations[[right]] <- population(paste0("pop", right),
                                          time = creation_times[[right]],
                                          N = population_size,
                                          parent = parent_populations[[2]])

    # update list of parent populations for left and right child
    # left child should have newly created population as parent
    # right child should have parent of newly created population as parent
    parent_populations <-  list(env$populations[[right]], parent_populations[[2]])
    # continue to recurse through the tree by only adjusting the current node
    recursion(tree = tree, current_node = right,
              parent_populations = parent_populations,
              population_size = population_size,
              env = env,
              creation_times = creation_times)
  }
}

random_gene_flow <- function(populations, n, rate, simulation_length) {

  # create empty gene flow list
  gf <- vector(mode = "list", length = n)
  if (n != 0){

    # find the times where the populations were created and store in list
    creation_times <- sapply(populations,
                             function(p) attr(p, "history")[[1]]$time)
    # for each gene flow event
    for(i in 1:n) {

      # sample start time of gene flow, can only be one generation after second
      # population is generated only one population exists
      start_gf <- sample((burn_in + 1):(simulation_length - 1), 1)
      # gene flow ends one generation after it started
      end_gf <- start_gf + 1

      # get a list of all populations that were created befor the gene flow
      # starts as only these can be considered
      possible <- list()
      j = 1
      for (p in populations) {
        # check if the creation time of the population is smaller than the time
        # where the gene flow starts
        if (p$time < start_gf) {
          # if so, add population to list of possible populations
          possible[[j]] <- p
          j = j + 1
        }
      }

      # sample two indices within the possible population list
      number1 <- sample(1:length(possible), 1)
      number2 <- sample(1:length(possible), 1)
      # check that gene flow would not be between the same population
      while(number1 == number2){
        number2 <- sample(1:length(possible), 1)
      }
      # get the populations corresponding to the sampled indices
      pop1 <- possible[[number1]]
      pop2 <- possible[[number2]]

      # given rate of gene flow can be either list of min/max or single value
      if (length(rate) == 2){
        if (rate[1] > rate[2]) {
          stop("No valid range for gene flow rate", call. = FALSE)
        }
        # sample rate within given range
        rate_gf <- runif(1, rate[1], rate[2])
      } else {
        # keep given rate
        rate_gf <- rate
      }

      # create new gene flow event and add to list
      gf[[i]] <- gene_flow(from = pop1, to = pop2, start = start_gf,
                           end = end_gf, rate = rate_gf)
    }
  }
  return(gf)
}