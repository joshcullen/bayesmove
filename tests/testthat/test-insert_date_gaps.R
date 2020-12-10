test_that("NA gaps are properly inserted within a time series", {

  #simulate data
  date1<- seq.POSIXt(from = ISOdatetime(2020, 12, 01, 00, 00, 00),
                     to = ISOdatetime(2020, 12, 31, 23, 59, 59), by = "8 hours")

  #remove datetimes to create gaps
  set.seed(1)
  ind<- sort(sample(2:(length(date1) - 1), 19))
  date1.red<- date1[-ind]

  #create data frame
  dat<- data.frame(id = 1, tseg = 1, date = date1.red, prop = 0.5)
  dat<- rbind(dat, dat)
  dat$behavior<- factor(rep(1:2, each = length(date1.red)))


  #run function
  dat.NA<- insert_date_gaps(data = dat, tol = 9, units = "hours")

  n.gaps<- length(ind) - sum(diff(ind) == 1)


  expect_equal(nrow(dat.NA)/2 - length(date1.red), n.gaps)
  expect_true(is.na(dat.NA[ind[1]+1,"prop"]))



})



