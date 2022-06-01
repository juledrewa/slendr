library(stringr)
library(ape)

devtools::load_all(".")
source("model-generation.R")

test_that("number of populations in tree_populations", {
  tree <- rtree(4)
  tree_pops <- tree_populations(tree, 1000, 50)
  expect_equal(length(tree_pops), 4)
})

test_that("number of populations in random_populations", {
  random_pops <- random_populations(5, 1000, 50)
  expect_equal(length(random_pops), 5)
})

test_that("number of populations in tree_model", {
  tree <- rtree(3)
  tree_model <- tree_model(tree, 1000, 3, c(0.2, 0.9), 100)
  expect_equal(length(tree_model$populations), 3)
  expect_equal(nrow(tree_model$splits), 3)
})

test_that("number of populations in random_model", {
  random_model <- random_model(6, 1000, 3, c(0.2, 0.9), 100)
  expect_equal(length(random_model$populations), 6)
  expect_equal(nrow(random_model$splits), 6)
})

test_that("number of gene flow events in tree_model", {
  tree <- rtree(3)
  tree_model <- tree_model(tree, 1000, 3, c(0.2, 0.9), 100)
  expect_equal(nrow(tree_model$geneflow), 3)
  tree_model <- tree_model(tree, 1000, 0, c(0.2, 0.9), 100)
  expect_equal(nrow(tree_model$geneflow), NULL)
  ######### is this correct???
})

test_that("number of gene flow events in random_model", {
  random_model <- random_model(6, 1000, 4, c(0.2, 0.9), 100)
  expect_equal(nrow(random_model$geneflow), 4)
})

test_that("simulation length in tree_model", {
  tree <- rtree(3)
  tree_model <- tree_model(tree, 1000, 3, c(0.2, 0.9), 100)
  expect_equal(tree_model$orig_length, 100)
})

test_that("simulation length in random_model", {
  random_model <- random_model(6, 1000, 4, c(0.2, 0.9), 1000)
  expect_equal(random_model$orig_length, 1000)
})