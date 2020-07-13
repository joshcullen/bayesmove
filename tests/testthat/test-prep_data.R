test_that("step lengths, turning angles, and time steps are calculated", {
  #simulate data
  lon<- c(1,1,3,4,4,5,7,9,10,13)
  lat<- c(2,1,2,2,3,5,8,1,1,2)
  date<- seq(c(ISOdate(2020, 6, 17, tz = "UTC")), by = "hour", length.out = 10)
  date<- date + lubridate::seconds(runif(length(date), -300, 300))  #introduce noise
  dat<- data.frame(id = 1, date, lon, lat)

  #run function
  dat1<- prep_data(dat = dat, coord.names = c("lon","lat"), id = "id")

  expect_is(dat1, "data.frame")
  expect_equal(dat1$x, dat$lon)
  expect_equal(dat1$y, dat$lat)
  expect_equal(dat1$date, dat$date)
  expect_lt(max(dat1$step, na.rm = TRUE), 8)
  expect_lte(max(abs(dat1$angle), na.rm = TRUE), pi)
  expect_condition(dat1$dt[10], NA)
  expect_gte(min(dat1$dt, na.rm = TRUE), 3000)
  expect_lte(max(dat1$dt, na.rm = TRUE), 4200)
})
