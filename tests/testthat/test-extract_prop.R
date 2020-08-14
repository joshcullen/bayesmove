test_that("function can extract values for theta from LDA results", {

  #simulate data
  id<- rep(1:4, each = 250)
  SL<- sample(1:5, 1000, replace = T)
  TA<- sample(1:8, 1000, replace = T)
  tseg<- rep(rep(1:10, each = 25), 4)

  dat<- data.frame(id, tseg, SL, TA)

  #summarize by time segment
  obs<- summarize_tsegs(dat = dat, nbins = c(5,8))

  #cluster data with LDA
  res<- cluster_segments(dat = obs, gamma1 = 0.1, alpha = 0.1, ngibbs = 250,
                         nburn = 125, nmaxclust = 7, ndata.types = 2)

  #Extract proportions of behaviors per time segment
  theta.estim<- extract_prop(res = res, ngibbs = 250, nburn = 125, nmaxclust = 7)

  expect_is(theta.estim, "matrix")
  expect_equal(dim(theta.estim), c(40, 7))
  expect_equal(sum(theta.estim[1,]), 1)
  expect_equal(sum(theta.estim[18,]), 1)
  expect_equal(sum(theta.estim[27,]), 1)

})
