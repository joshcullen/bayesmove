test_that("MAP value is correctly extracted", {

  #simulate data
  ngibbs<- 1000
  y<- (-1000 * 501:1500)/(-500 + 501:1500) + rnorm(ngibbs, 0, 0.1)
  dat<- matrix(c(1, y), 1, 1001)
  dat<- data.frame(dat)
  names(dat)[1]<- "id"

  #run function
  MAP.est<- get_MAP(dat = dat, nburn = ngibbs/2)

  expect_equal(MAP.est, 1000)
  expect_type(MAP.est, "integer")
  expect_length(MAP.est, 1)
})
