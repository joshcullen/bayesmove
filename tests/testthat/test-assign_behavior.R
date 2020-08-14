test_that("behavior proportions are properly merged with original dataset", {

  #simulate data
  id<- rep(1:4, each = 250)
  date<- rep(seq(c(ISOdate(2020, 8, 14, tz = "UTC")), by = "hour", length.out = 250), 4)
  SL<- sample(1:5, 1000, replace = T)
  TA<- sample(1:8, 1000, replace = T)
  obs<- rep(1:250, 4)
  time1<- rep(1:250, 4)
  tseg<- rep(rep(1:10, each = 25), 4)

  dat.orig<- data.frame(id, date, obs, time1, SL, TA)
  dat.seg<- data.frame(id, date, tseg, obs, time1, SL, TA)
  dat.seg.list<- df_to_list(dat.seg, "id")

  # Select only id, tseg, SL, and TA columns
  dat2<- dat.seg[,c("id","tseg","SL","TA")]

  #summarize by time segment
  obs<- summarize_tsegs(dat = dat2, nbins = c(5,8))

  #cluster data with LDA
  res<- cluster_segments(dat = obs, gamma1 = 0.1, alpha = 0.1, ngibbs = 250,
                         nburn = 125, nmaxclust = 7, ndata.types = 2)

  #Extract proportions of behaviors per time segment
  theta.estim<- extract_prop(res = res, ngibbs = 250, nburn = 125, nmaxclust = 7)

  #Create augmented matrix by replicating rows (tsegs) according to obs per tseg
  theta.estim.long<- expand_behavior(dat = dat.seg, theta.estim = theta.estim, obs = obs,
                                     nbehav = 3, behav.names = c("Encamped","ARS","Transit"),
                                     behav.order = c(1,2,3))

  #Run function
  dat.out<- assign_behavior(dat.orig = dat.orig, dat.seg.list = dat.seg.list,
                            theta.estim.long = theta.estim.long,
                            behav.names = c("Encamped","ARS","Transit"))

  expect_named(dat.out[,ncol(dat.out), drop = F], "behav")
  expect_equal(sum(dat.out[1, 8:10]), 1)
  expect_equal(nrow(dat.out), nrow(dat.orig))
})
