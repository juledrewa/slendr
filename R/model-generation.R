library(phangorn)
library(ape)
library(phytools)
library(slendr)


tree_model <- function(tree, size) {
  # plot to check
  plot(tree, show.tip.label = F)
  nodelabels()
  tiplabels()

  env <- new.env()
  env$populations <- vector(mode="list", length=tree$Nnode + length(tree$tip.label))

  root <- length(tree$tip.label) + 1
  env$populations[[root]] <- population(root, time = 1, N = size)

  recursion(tree, root, 1, size, env)
  env$populations
}

model <- function(n, size) {
  tree <- rtree(n)
  # plot to check
  plot(tree, show.tip.label = F)
  nodelabels()
  tiplabels()

  env <- new.env()
  env$populations <- vector(mode="list", length=tree$Nnode + length(tree$tip.label))

  root <- length(tree$tip.label) + 1
  env$populations[[root]] <- population(root, time = 1, N = size)

  recursion(tree, root, 1, size, env)
  env$populations
}

recursion <- function(tree, current_node, current_length, size, env) {
  edge_list <- tree$edge
  children <- edge_list[edge_list[, 1] == current_node, 2]
  left <- children[1]
  right <- children[2]
  if (!(left %in% tree$tip.label) && !is.na(left)){
    length <- tree$edge.length[edge_list[, 1] == current_node & edge_list[, 2] == left] * 100 + current_length
    env$populations[[left]] <- population(left, time = length, N = size, parent = env$populations[[current_node]])

    recursion(tree, left, length, size, env)
  }
  if (!(right %in% tree$tip.label)&& !is.na(right)){
    length <- tree$edge.length[edge_list[, 1] == current_node & edge_list[, 2] == right] * 100 + current_length
    env$populations[[right]] <- population(right, time = length, N = size, parent = env$populations[[current_node]])

    recursion(tree, right, length, size, env)
  }
}


# test calls
tree <- rtree(5)
tree_model(tree, 1000)
model(5, 1000)
