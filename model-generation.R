library(phangorn)
library(ape)
library(phytools)

devtools::load_all(".")


tree_populations <- function(tree, size, simulation_length) {
  # plot to check
  plot(tree, show.tip.label = F)
  nodelabels()
  tiplabels()

  env <- new.env()
  env$populations <- vector(mode="list", length=tree$Nnode + length(tree$tip.label))

  scale <- floor((3*simulation_length/4)/max(node.depth.edgelength(tree)))

  root <- length(tree$tip.label) + 1
  env$populations[[root]] <- population(paste0("pop", root), time = 1, N = size)

  recursion(tree, root, env$populations[[root]], size, env, scale)
  env$populations <- env$populations[-which(sapply(env$populations, is.null))]
  env$populations
}

random_populations <- function(n, size, simulation_length) {
  tree <- rtree(n)
  populations <- tree_populations(tree, size, simulation_length)
}

tree_model <- function(tree, pop_size, n_gene_flow, simulation_length, rate = c(0.01, 0.99)) {
  populations <- tree_populations(tree, pop_size, simulation_length)
  gf <- random_gene_flow(populations, n_gene_flow, simulation_length, rate)
  model <- compile_model(populations = populations, gene_flow = gf, generation_time = 1, sim_length = simulation_length)
  model
}

random_model <- function(n, pop_size, n_gene_flow, simulation_length, rate = c(0.01, 0.99)) {

  tree <- rtree(n)
  tree_model(tree, pop_size, n_gene_flow, simulation_length, rate)
}

recursion <- function(tree, current_node, parent_population, size, env, scale) {
  edge_list <- tree$edge
  children <- edge_list[edge_list[, 1] == current_node, 2]
  left <- children[1]
  right <- children[2]

  if (!(left %in% tree$tip.label) && !is.na(left)){
    length <- ceiling(node.depth.edgelength(tree)[left] * scale) + 1
    env$populations[[left]] <- population(paste0("pop", left), time = length, N = size, parent = parent_population)

    recursion(tree, left, env$populations[[left]], size, env, scale)
  }
  if (!(right %in% tree$tip.label) && !is.na(right)){
    length <- ceiling(node.depth.edgelength(tree)[right] * scale)
    print(right)
    #env$populations[[right]] <- population(paste0("pop", right), time = length, N = size, parent = env$populations[[current_node]])
    recursion(tree, right, parent_population, size, env, scale)
  }
}

random_gene_flow <- function(populations, n, simulation_length, rate) {
  gf <- vector(mode = "list", length = n)
  for(i in 1:n) {
    number1 <- sample(1:length(populations), 1)
    number2 <- sample(1:length(populations), 1)
    while(number1 == number2){
      number2 <- sample(1:length(populations), 1)
    }
    pop1 <- populations[[number1]]
    pop2 <- populations[[number2]]
    time1 <- attr(pop1, "history")[[1]]$time
    time2 <- attr(pop2, "history")[[1]]$time
    max_time <- max(time1, time2)
    start_gf <- sample((max_time + 1):(simulation_length - 2), 1)
    end_gf <- sample((start_gf + 1):(simulation_length), 1)
    if (length(rate) == 2){
      rate_gf <- runif(1, rate[1], rate[2])
    } else {
      rate_gf <- rate
    }
    gf[[i]] <- gene_flow(from = pop1, to = pop2, start = start_gf, end = end_gf, rate_gf)
  }
  return(gf)
}

# test calls
tree <- rtree(4)
pops <- tree_populations(tree, 1000, 50)
pops2 <- random_populations(3, 1000, 100)
model <- tree_model(tree, 1000, 3, 100, c(0.2, 0.9))
model2 <- random_model(3, 1000, 2, 1000, 0.5)
plot_model(model, proportions = TRUE)
plot_model(model2, proportions = TRUE)


