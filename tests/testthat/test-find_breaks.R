test_that("multiplication works", {
  #simuluate data
  var<- sample(1:3, size = 50, replace = TRUE)
  var<- rep(var, each = 20)
  id<- rep(1:10, each = 100)

  #create data frame
  dat<- data.frame(id, var)

  #true breakpoints
  true.breaks1<- which(diff(dat[dat$id == 1, 'var']) != 0) + 1
  true.breaks3<- which(diff(dat[dat$id == 3, 'var']) != 0) + 1
  true.breaks5<- which(diff(dat[dat$id == 5, 'var']) != 0) + 1
  true.breaks7<- which(diff(dat[dat$id == 7, 'var']) != 0) + 1
  true.breaks9<- which(diff(dat[dat$id == 9, 'var']) != 0) + 1

  #create list
  dat.list<- df_to_list(dat = dat, ind = "id")

  #run function using purrr::map()
  breaks<- purrr::map(dat.list, ~find_breaks(dat = ., ind = "var"))

  expect_is(breaks, "list")
  expect_is(dat$var, "integer")
  expect_equal(breaks[[1]], true.breaks1)
  expect_equal(breaks[[3]], true.breaks3)
  expect_equal(breaks[[5]], true.breaks5)
  expect_equal(breaks[[7]], true.breaks7)
  expect_equal(breaks[[9]], true.breaks9)
})
