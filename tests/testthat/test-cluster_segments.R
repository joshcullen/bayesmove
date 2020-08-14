test_that("function clusters time segments", {

  #simulate data
  id<- rep(1:4, each = 250)
  SL<- sample(1:5, 1000, replace = T)
  TA<- sample(1:8, 1000, replace = T)
  tseg<- rep(rep(1:10, each = 25), 4)

  dat.seg<- data.frame(id, tseg, SL, TA)

  #summarize by time segment
  obs<- summarize_tsegs(dat = dat.seg, nbins = c(5,8))

  #cluster data with LDA
  res<- cluster_segments(dat = obs, gamma1 = 0.1, alpha = 0.1, ngibbs = 250,
                         nburn = 125, nmaxclust = 7, ndata.types = 2)


  expect_length(res, 4)
  expect_type(res$loglikel, "double")
  expect_is(res$theta, "matrix")
  expect_length(res$phi, 2)
  expect_type(res$phi, "list")
})
