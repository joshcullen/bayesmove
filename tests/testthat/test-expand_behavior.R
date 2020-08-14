test_that("behavior proportions augmented from segment to observation level", {

  #simulate data
  id<- rep(1:4, each = 250)
  date<- rep(seq(c(ISOdate(2020, 8, 14, tz = "UTC")), by = "hour", length.out = 250), 4)
  SL<- sample(1:5, 1000, replace = T)
  TA<- sample(1:8, 1000, replace = T)
  time1<- rep(1:250, 4)
  tseg<- rep(rep(1:10, each = 25), 4)

  dat<- data.frame(id, date, tseg, time1, SL, TA)

  # Select only id, tseg, SL, and TA columns
  dat2<- dat[,c("id","tseg","SL","TA")]

  #summarize by time segment
  obs<- summarize_tsegs(dat = dat2, nbins = c(5,8))

  #cluster data with LDA
  res<- cluster_segments(dat = obs, gamma1 = 0.1, alpha = 0.1, ngibbs = 250,
                         nburn = 125, nmaxclust = 7, ndata.types = 2)

  #Extract proportions of behaviors per time segment
  theta.estim<- extract_prop(res = res, ngibbs = 250, nburn = 125, nmaxclust = 7)

  #Create augmented matrix by replicating rows (tsegs) according to obs per tseg
  theta.estim.long<- expand_behavior(dat = dat, theta.estim = theta.estim, obs = obs,
                                 nbehav = 3, behav.names = c("Encamped","ARS","Transit"),
                                 behav.order = c(1,2,3))

  expect_equal(3*nrow(dat), nrow(theta.estim.long))
  expect_s3_class(theta.estim.long, "data.frame")
  expect_equal(length(unique(theta.estim.long$time1)), 250)
})
