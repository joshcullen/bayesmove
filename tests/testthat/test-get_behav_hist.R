test_that("function correctly wrangles phi matrix into data frame", {

  #simulate data
  SL<- matrix(runif(500*35), 500, 35)
  seq1<- seq(1, by = 7, length.out = 5)
  for (i in 0:6) {
    tmp<- t(apply(SL[,(seq1+i)], 1, function(x) x/sum(x)))
    SL[,(seq1+i)]<- tmp
  }

  TA<- matrix(runif(500*56), 500, 56)
  seq1<- seq(1, by = 7, length.out = 8)
  for (i in 0:6) {
    tmp<- t(apply(TA[,(seq1+i)], 1, function(x) x/sum(x)))
    TA[,(seq1+i)]<- tmp
  }

  phi<- list(SL, TA)

  res<- list(phi = phi)

  #to test for results from mixture model
  res2<- c(res, z = list(seq(1,1000)))


  #run function (for segmentation/LDA)
  behav.res<- get_behav_hist(dat = res, nburn = 250, ngibbs = 500, nmaxclust = 7,
                             var.names = c("Step Length","Turning Angle"))

  tmp_SL<- behav.res[behav.res$behav == 1 & behav.res$var == "Step Length", "prop"]
  tmp_TA<- behav.res[behav.res$behav == 1 & behav.res$var == "Turning Angle", "prop"]

  expect_equal(sum(dplyr::pull(tmp_SL, prop)), 1)
  expect_equal(sum(dplyr::pull(tmp_TA, prop)), 1)
  expect_equal(nrow(behav.res), 35+56)


  #run function (for mixture model)
  behav.res2<- get_behav_hist(dat = res2, nburn = 250, ngibbs = 500, nmaxclust = 7,
                              var.names = c("Step Length","Turning Angle"),
                              ord = 1:7, MAP.iter = 450)

  tmp_SL2<- behav.res2[behav.res2$behav == 1 & behav.res2$var == "Step Length", "prop"]
  tmp_TA2<- behav.res2[behav.res2$behav == 1 & behav.res2$var == "Turning Angle", "prop"]

  expect_equal(sum(dplyr::pull(tmp_SL2, prop)), 1)
  expect_equal(sum(dplyr::pull(tmp_TA2, prop)), 1)
  expect_equal(nrow(behav.res2), 35+56)
})
