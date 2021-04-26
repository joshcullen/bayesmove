test_that("segmentation model can effectively segment a dataset and return results", {

  #simulate data
  step<- rgamma(500, c(1, 2.5, 10), c(1, 1, 1))
  angle<- runif(500, -pi, pi)
  date<- seq(c(ISOdate(2020, 8, 13, tz = "UTC")), by = "hour", length.out = 500)
  date<- date + lubridate::seconds(runif(length(date), -15, 15))  #introduce noise
  dt<- as.numeric(diff(date)) * 60  #convert time difference to seconds
  dt<- c(dt, NA)
  id<- rep(paste0("id", 1), 500)


  #create data frame and round time
  dat<- data.frame(id, date, dt, step, angle)
  dat<- round_track_time(dat = dat, id = "id", int = 3600, tol = 15, time.zone = "UTC",
                         units = "secs")


  #define limits for each bin
  dist.lims<- quantile(step, c(0, 0.25, 0.5, 0.75, 0.95, 1))  #5 bins
  angle.lims<- c(-pi, -3*pi/4, -pi/2, -pi/4, 0, pi/4, pi/2, 3*pi/4, pi)  #8 bins

  #discretize step and angle
  dat1<- discrete_move_var(dat = dat, lims = list(dist.lims, angle.lims),
                           varIn = c("step", "angle"),
                           varOut = c("SL","TA"))


  #create list and filter by primary time step
  dat.list<- df_to_list(dat = dat1, ind = "id")
  dat.list.filt<- filter_time(dat.list = dat.list, int = 3600)


  #perform segmentation w/o pre-specifying breakpoints
  dat.list.filt1<- lapply(dat.list.filt,
                          function(x) subset(x, select = c(id, SL, TA)))
  progressr::with_progress({
    dat.res1<- segment_behavior(data = dat.list.filt1, ngibbs = 1000, nbins = c(5,8),
                                alpha = 1)
  })

  expect_length(dat.res1, 4)
  expect_type(dat.res1, "list")
  expect_type(dat.res1$brkpts, "list")
  expect_length(dat.res1$brkpts$id1, 1000)
  expect_s3_class(dat.res1$nbrks, "data.frame")
  expect_s3_class(dat.res1$LML, "data.frame")
  expect_s3_class(dat.res1$elapsed.time, "data.frame")
  expect_type(dat.res1$nbrks$Iter_1, "double")
  expect_type(dat.res1$LML$Iter_1, "double")
})
