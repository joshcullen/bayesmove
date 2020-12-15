test_that("in POSIXct format and only changed date and dt", {
  #simulate data
  set.seed(1)
  step<- rgamma(1000, c(1, 2.5, 10), c(1, 1, 1))
  angle<- runif(1000, -pi, pi)
  date<- seq(c(ISOdate(2020, 6, 17, tz = "UTC")), by = "hour", length.out = 1000)
  date<- date + lubridate::seconds(runif(length(date), -15, 15))  #introduce noise
  dt<- as.numeric(diff(date)) * 60  #convert time difference to seconds
  dt<- c(dt, NA)
  id<- rep(1:10, each = 100)

  #create data frame
  dat<- data.frame(id, date, dt, step, angle)

  #run function
  dat1<- round_track_time(dat = dat, id = "id", int = 3600, tol = 20, time.zone = "UTC",
                          units = "secs")

  expect_is(dat1, "data.frame")
  expect_is(dat1$date, "POSIXct")
  expect_equal(dat1$dt[1], 3600)
  expect_error(expect_identical(dat1$date, dat$date))
  expect_error(expect_identical(dat1$dt, dat$dt))
  expect_equal(dat1$step, dat$step)
})
