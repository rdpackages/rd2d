test_that("main user-facing functions are exported", {
  expect_true(is.function(rd2d))
  expect_true(is.function(rd2d.dist))
  expect_true(is.function(rdbw2d))
  expect_true(is.function(rdbw2d.dist))
})
