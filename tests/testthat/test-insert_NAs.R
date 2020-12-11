test_that("NA gaps are properly inserted within a time series", {

  #simulate data
  date1<- seq.POSIXt(from = ISOdatetime(2020, 12, 01, 00, 00, 00),
                     to = ISOdatetime(2020, 12, 31, 23, 59, 59), by = "8 hours")

  #remove datetimes to create gaps
  set.seed(1)
  ind<- sort(sample(2:(length(date1) - 1), 19))
  date1.red<- date1[-ind]

  #create data frame
  dat<- data.frame(id = 1, date = date1.red,
                   dt = c(as.numeric(diff(date1.red)), NA),
                   covar1 = rnorm(length(date1.red)),
                   covar2 = rgamma(length(date1.red), 1, 1))

  #run function
  dat2<- insert_NAs(dat, int = 8, units = "hours")

  gaps<- which(is.na(dat2$covar1))


  expect_equal(length(gaps), length(ind))
  expect_identical(gaps, ind)
  expect_length(unique(na.omit(dat2$dt)), 1)

})
