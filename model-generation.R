library(phangorn)
library(ape)
library(phytools)

devtools::load_all(".")


tree_populations <- function(tree, size, height) {
  # plot to check
  plot(tree, show.tip.label = F)
  nodelabels()
  tiplabels()

  env <- new.env()
  env$populations <- vector(mode="list", length=tree$Nnode + length(tree$tip.label))


  scale <- floor(height/sum(tree$edge.length))

  root <- length(tree$tip.label) + 1
  env$populations[[root]] <- population(paste0("pop", root), time = 1, N = size)

  recursion(tree, root, size, env, scale)

  env$populations
}

random_populations <- function(n, size, height) {
  tree <- rtree(n)
  populations <- tree_populations(tree, size, height)
}

tree_model <- function(tree, pop_size, n_gene_flow, height) {
  populations <- tree_populations(tree, pop_size, height)
  gf <- random_gene_flow(populations, n_gene_flow, height)
  model <- compile_model(populations = populations, gene_flow = gf, generation_time = 1, sim_length = height)
  model
}

random_model <- function(n, pop_size, n_gene_flow, height) {
  populations <- random_populations(n, pop_size, height)
  gf <- random_gene_flow(populations, n_gene_flow, height)
  model <- compile_model(populations = populations, gene_flow = gf, generation_time = 1, sim_length = height)
  model

  ### alternative
  # tree <- rtree(n)
  # tree_model(tree, pop_size, n_gene_flow, height)
}

recursion <- function(tree, current_node, size, env, scale) {
  edge_list <- tree$edge
  children <- edge_list[edge_list[, 1] == current_node, 2]
  left <- children[1]
  right <- children[2]

  if (!(left %in% tree$tip.label) && !is.na(left)){
    length <- ceiling(tree$edge.length[edge_list[, 1] == current_node & edge_list[, 2] == left] * scale) + attr(env$populations[[current_node]], "history")[[1]]$time
    env$populations[[left]] <- population(paste0("pop", left), time = length, N = size, parent = env$populations[[current_node]])

    recursion(tree, left, size, env, scale)
  }
  if (!(right %in% tree$tip.label)&& !is.na(right)){
    length <- ceiling(tree$edge.length[edge_list[, 1] == current_node & edge_list[, 2] == right] * scale) + attr(env$populations[[current_node]], "history")[[1]]$time
    env$populations[[right]] <- population(paste0("pop", right), time = length, N = size, parent = env$populations[[current_node]])

    recursion(tree, right, size, env, scale)
  }
}

random_gene_flow <- function(populations, n, height) {
  gf <- list()
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
    start_gf <- sample((max_time + 1):(height - 2), 1)
    end_gf <- sample((start_gf + 1):(height - 1), 1)
    gf <- append(gf, list(gene_flow(from = pop1, to = pop2, start = start_gf, end = end_gf, rate = 0.2)))
  }
  return(gf)
}

# test calls
tree <- rtree(3)
pops <- tree_populations(tree, 1000, 50)
pops2 <- random_populations(3, 1000, 100)
model <- tree_model(tree, 1000, 3, 100)
model2 <- random_model(3, 1000, 2, 100)
plot_model(model)
plot_model(model2)
