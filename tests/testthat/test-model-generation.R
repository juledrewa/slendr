test_that("number of populations in tree_populations", {
  tree <- ape::rtree(4)
  tree_pops <- tree_populations(tree = tree, population_size = 1000,
                                simulation_length = 50)
  expect_equal(length(tree_pops), 4)
  tree <- ape::rtree(2)
  tree_pops <- tree_populations(tree = tree, population_size = 750,
                                simulation_length = 100)
  expect_equal(length(tree_pops), 2)

  tree <- ape::rtree(1)
  expect_error(tree_populations(tree = tree, population_size = 50,
                                simulation_length = 1000),
               "At least 2 populations must be created")
})

test_that("number of populations in random_populations", {
  random_pops <- random_populations(n_populations = 5, population_size = 100,
                                    simulation_length = 500)
  expect_equal(length(random_pops), 5)
  random_pops <- random_populations(n_populations = 2, population_size = 10,
                                    simulation_length = 7500)
  expect_equal(length(random_pops), 2)

  expect_error(random_populations(n_populations = 0, population_size = 800,
                                  simulation_length = 1500),
               "At least 2 populations must be created")
  expect_error(random_populations(n_populations = 1, population_size = 70,
                                  simulation_length = 250),
               "At least 2 populations must be created")
})

test_that("number of populations in tree_model", {
  tree <- ape::rtree(3)
  tree_model <- tree_model(tree = tree, population_size = 900, n_gene_flow = 3,
                           rate_gene_flow = c(0.2, 0.8),
                           simulation_length = 600)
  expect_equal(length(tree_model$populations), 3)
  expect_equal(nrow(tree_model$splits), 3)
})

test_that("number of populations in random_model", {
  random_model <- random_model(n_populations = 6, population_size = 580,
                               n_gene_flow = 4,
                               rate_gene_flow = 0.6,
                               simulation_length = 400)
  expect_equal(length(random_model$populations), 6)
  expect_equal(nrow(random_model$splits), 6)
})

test_that("number of gene flow events in tree_model", {
  tree <- ape::rtree(3)
  tree_model <- tree_model(tree = tree, population_size = 8000, n_gene_flow = 3,
                           rate_gene_flow = c(0.3, 0.4),
                           simulation_length = 1000)
  expect_equal(nrow(tree_model$geneflow), 3)
  tree_model <- tree_model(tree = tree, population_size = 40, n_gene_flow = 0,
                           rate_gene_flow = 0.1,
                           simulation_length = 540)
  expect_equal(nrow(tree_model$geneflow), NULL)
})

test_that("number of gene flow events in random_model", {
  random_model <- random_model(n_populations = 6, population_size = 400,
                               n_gene_flow = 4,
                               rate_gene_flow = c(0.6, 0.8),
                               simulation_length = 600)
  expect_equal(nrow(random_model$geneflow), 4)

  random_model <- random_model(n_populations = 100, population_size = 3000,
                               n_gene_flow = 50,
                               rate_gene_flow = 0.9,
                               simulation_length = 2300)
  expect_equal(nrow(random_model$geneflow), 50)
})

test_that("rate gene flow in given range", {
  tree <- ape::rtree(6)
  tree_model <- tree_model(tree = tree, population_size = 60, n_gene_flow = 12,
                           rate_gene_flow = c(0.2, 0.9),
                           simulation_length = 1800)
  expect_true(max(tree_model$geneflow$rate) <= 0.9)
  expect_true(min(tree_model$geneflow$rate) >= 0.2)


  tree <- ape::rtree(7)
  tree_model <- tree_model(tree = tree, population_size = 520, n_gene_flow = 8,
                           rate_gene_flow = c(0.4, 0.7),
                           simulation_length = 6100)
  expect_true(max(tree_model$geneflow$rate) <= 0.7)
  expect_true(min(tree_model$geneflow$rate) >= 0.4)

  tree <- ape::rtree(9)
  expect_error(tree_model(tree = tree, population_size = 340, n_gene_flow = 6,
                          rate_gene_flow = c(0.7, 0.4),
                          simulation_length = 560),
               "No valid range for gene flow rate")
})

test_that("rate gene flow equal to given fixed value", {
  tree <- ape::rtree(5)
  tree_model <- tree_model(tree = tree, population_size = 700, n_gene_flow = 5,
                           rate_gene_flow = 0.5,
                           simulation_length = 420)
  expect_true(max(tree_model$geneflow$rate) == 0.5)
  expect_true(min(tree_model$geneflow$rate) == 0.5)

  tree <- ape::rtree(8)
  tree_model <- tree_model(tree = tree, population_size = 1000, n_gene_flow = 14,
                           rate_gene_flow = 0.2,
                           simulation_length = 4800)
  expect_true(max(tree_model$geneflow$rate) == 0.2)
  expect_true(min(tree_model$geneflow$rate) == 0.2)
})

test_that("simulation length in tree_model", {
  tree <- ape::rtree(3)
  tree_model <- tree_model(tree = tree, population_size = 290, n_gene_flow = 3,
                           rate_gene_flow = c(0.4, 0.7),
                           simulation_length = 100)
  expect_equal(tree_model$orig_length, 100)
})

test_that("simulation length in random_model", {
  random_model <- random_model(n_populations = 6, population_size = 25,
                               n_gene_flow = 4, rate_gene_flow = c(0.1, 0.8),
                               simulation_length = 1000)
  expect_equal(random_model$orig_length, 1000)
})

test_that("all populations created before simulation ends", {
  tree <- ape::rtree(8)
  tree_model <- tree_model(tree = tree, population_size = 800, n_gene_flow = 8,
                           rate_gene_flow = 0.5,
                           simulation_length = 450)
  expect_true(max(sapply(tree_model$populations,
                         function(p) attr(p, "history")[[1]]$time)) < 450)
})