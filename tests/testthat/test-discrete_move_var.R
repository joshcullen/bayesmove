test_that("discretization returns data frame of with vars of bin numbers", {
  #simulate data
  step<- rgamma(1000, c(1, 2.5, 10), c(1, 1, 1))
  angle<- runif(1000, -pi, pi)
  id<- rep(1:10, each = 100)

  #create data frame
  dat<- data.frame(id, step, angle)

  #define limits for each bin
  dist.lims<- c(quantile(step, c(0, 0.25, 0.5, 0.75, 0.95)), max(step))
  angle.lims<- c(-pi, -3*pi/4, -pi/2, -pi/4, 0, pi/4, pi/2, 3*pi/4, pi)

  #run function
  dat1<- discrete_move_var(dat = dat, lims = list(dist.lims, angle.lims),
                           varIn = c("step", "angle"),
                           varOut = c("SL","TA"))


  expect_equal(max(dat1$SL), 5)
  expect_equal(max(dat1$TA), 8)
  expect_is(dat1, "data.frame")

})
