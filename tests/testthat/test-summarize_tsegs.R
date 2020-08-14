test_that("observations are correctly counted by time segment", {

  #simulate data
  id<- rep(1:4, each = 250)
  SL<- sample(1:5, 1000, replace = T)
  TA<- sample(1:8, 1000, replace = T)
  tseg<- rep(rep(1:10, each = 25), 4)

  dat<- data.frame(id, tseg, SL, TA)

  #run function
  obs<- summarize_tsegs(dat = dat, nbins = c(5,8))


  expect_s3_class(obs, "data.frame")
  expect_length(obs$id, 40)
  expect_equal(sum(obs[1, which(grepl("y1", names(obs)))]), 25)
})
