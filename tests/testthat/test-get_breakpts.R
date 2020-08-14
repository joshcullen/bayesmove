test_that("correct breakpoints are extracted based on MAP estimate", {

  #simulate data
  id1 = id2 = id3 = list(4, c(2,4), 4, 4, c(4,8), c(4,8,17), c(4,8,17),
                         c(4,8,20), c(4,8,20,25), c(5,8,20,25))
  dat.list<- list(id1 = id1, id2 = id2, id3 = id3)

  MAP.est<- c(5, 8, 9)


  #run function
  brks<- get_breakpts(dat = dat.list, MAP.est = MAP.est)


  expect_s3_class(brks, "data.frame")
  expect_failure(expect_type(brks[,1], "double"))
  expect_type(brks[,2], "double")
  expect_type(brks[,3], "double")
  expect_type(brks[,4], "double")
  expect_equal(as.numeric(brks[1,-1]), c(4,8,NA,NA))
})
