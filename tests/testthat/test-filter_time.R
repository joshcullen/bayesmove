test_that("time is filtered and dt and date are updated", {
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
  dat<- round_track_time(dat = dat, id = "id", int = 3600, tol = 15, time.zone = "UTC")

  #create list
  dat.list<- df_to_list(dat = dat, ind = "id")

  #run function
  dat.list.filt<- filter_time(dat.list = dat.list, int = 3600)

  expect_is(dat.list.filt, "list")
  expect_is(dat.list.filt[[1]], "data.frame")
  expect_equal(unique(unlist(purrr::map(dat.list.filt, ~purrr::pluck(., "dt")))), 3600)
  expect_equal(max(dat.list.filt[[1]]$time1), nrow(dat.list.filt[[1]]))
  expect_is(dat.list.filt[[1]]$obs, "integer")

})
