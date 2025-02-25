test_that("all three variants of the Natural Earth data can be downloaded", {
  invisible(capture.output(small <- world(xrange = c(-15, 40), yrange = c(30, 60), crs = "EPSG:3035")))
  invisible(capture.output(small2 <- world(xrange = c(-15, 40), yrange = c(30, 60), crs = "EPSG:3035", scale = "small")))
  invisible(capture.output(medium <- world(xrange = c(-15, 40), yrange = c(30, 60), crs = "EPSG:3035", scale = "medium")))
  invisible(capture.output(large <- world(xrange = c(-15, 40), yrange = c(30, 60), crs = "EPSG:3035", scale = "large")))
  expect_equal(small, small2)
  expect_s3_class(small, "slendr_map")
  expect_s3_class(small2, "slendr_map")
  expect_s3_class(medium, "slendr_map")
  expect_s3_class(large, "slendr_map")
})

test_that("world dimensions are correctly specified", {
  expect_error(world(1:100, 1:100, landscape = "blank"))
  expect_s3_class(world(c(1, 100), c(1, 100), landscape = "blank"), "slendr_map")
})
