test_that("time segments are assigned properly", {
  #simulate data

  tmp<- 1:105
  tmp<- tmp[-c(2,5,8,50,100)]
  time1<- rep(tmp, 10)
  id<- rep(1:10, each = 100)

  #simulate breakpoints
  brkpts<- sample(2:90, 70, replace = FALSE)
  brkpts<- matrix(brkpts, 10, 7, byrow = TRUE)
  brkpts<- t(apply(brkpts, 1, sort))
  brkpts<- data.frame(id = 1:10, brkpts)

  #create data frame
  dat<- data.frame(id, time1)

  #create list
  dat.list<- df_to_list(dat = dat, ind = "id")

  #run function
  dat1<- assign_tseg(dat = dat.list, brkpts = brkpts)

  #check number of unique tsegs per id
  n_uni<- dat1 %>%
    dplyr::group_by(id) %>%
    dplyr::summarise(n=dplyr::n_distinct(tseg)) %>%
    dplyr::pull(n)

  expect_equal(n_uni, rep(8, 10))
  expect_s3_class(dat1, "data.frame")
  expect_identical(nrow(dat1), length(id))
  expect_equal(sum(ifelse(is.na(dat1$tseg), 0, 1)), 1000)
  expect_identical(dat1$time1[1:100], 1:100)

})
