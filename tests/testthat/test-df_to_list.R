test_that("list is created properly from data frame", {
  #simulate data
  step<- rgamma(1000, c(1, 2.5, 10), c(1, 1, 1))
  angle<- runif(1000, -pi, pi)
  date<- seq(c(ISOdate(2020, 6, 17, tz = "UTC")), by = "hour", length.out = 1000)
  date<- date + lubridate::seconds(runif(length(date), -15, 15))  #introduce noise
  dt<- as.numeric(diff(date)) * 60  #convert time difference to seconds
  dt<- c(dt, NA)
  id<- rep(1:10, each = 100)

  #create data frame
  dat<- data.frame(id, date, dt, step, angle)
  dat<- round_track_time(dat = dat, id = "id", int = 3600, tol = 15, time.zone = "UTC",
                         units = "secs")

  #run function
  dat.list<- df_to_list(dat = dat, ind = "id")

  expect_is(dat.list, "list")
  expect_equal(length(dat.list), length(unique(id)))
  expect_equal(names(dat.list), unique(as.character(id)))
})
