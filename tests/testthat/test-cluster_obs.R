test_that("mixture model works", {

  #load data
  data(tracks.list)
  tracks.list<- lapply(tracks.list, function(x) x[1:250,])

  #convert from list to data frame
  tracks.list<- dplyr::bind_rows(tracks.list)

  #only retain id and discretized step length (SL) and turning angle (TA) columns
  tracks<- subset(tracks.list, select = c(SL, TA))


  set.seed(1)

  # Define model params
  alpha=0.1
  ngibbs=1000
  nburn=ngibbs/2
  nmaxclust=7

  dat.res<- cluster_obs(dat = tracks, alpha = alpha, ngibbs = ngibbs,
                             nmaxclust = nmaxclust, nburn = nburn)


  expect_length(dat.res, 6)
  expect_type(dat.res$loglikel, "double")
  expect_is(dat.res$theta, "matrix")
  expect_length(dat.res$phi, 2)
  expect_type(dat.res$phi, "list")
  expect_length(dat.res$z.MAP, nrow(tracks))
  expect_equal(rowSums(dat.res$z.posterior), rep(nburn, nrow(tracks)))
})
