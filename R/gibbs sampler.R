behav.gibbs.sampler=function(dat, ngibbs, nbins, alpha, breakpt) {

  start.time<- Sys.time()  #start timer
  set.seed(1)

  uni.id<- unique(dat$id)  #need this query if running multiple IDs in parallel
  dat<- subset(dat, select = -id)

  #useful stuff
  max.time<- nrow(dat)
  ndata.types<- length(nbins)

  #starting values
  if (is.null(breakpt)) breakpt<- floor(max.time/2)

  #to store results
  res.brks<- vector("list", ngibbs)
  res.LML<- matrix(NA, 1, (ngibbs+1))
  res.nbrks<- matrix(NA, 1, (ngibbs+1))
  store.param<- matrix(NA, ngibbs, 2)

  for (i in 1:ngibbs){
    vals<- samp.move(breakpt = breakpt, max.time = max.time, dat = dat,
                   alpha = alpha, nbins = nbins, ndata.types = ndata.types)

    breakpt<- vals[[1]]

    #store results
    res.brks[[i]]<- breakpt
    store.param[i,]<- c(length(breakpt), vals[[2]])

  }

  tmp<- store.param[,1]
  res.nbrks[1,]<- c(uni.id,tmp)
  colnames(res.nbrks)<- c('id', paste0("Iter_",1:ngibbs))

  tmp<- store.param[,2]
  res.LML[1,]<- c(uni.id,tmp)
  colnames(res.LML)<- c('id', paste0("Iter_",1:ngibbs))

  end.time<- Sys.time()
  elapsed.time<- difftime(end.time, start.time, units = "min")  #end timer

  list(breakpt=res.brks, nbrks=res.nbrks, LML=res.LML, elapsed.time=elapsed.time)
}
#----------------------------------------------------
segment.behavior=function(data, ngibbs, nbins, alpha, breakpt=NULL) {

  ## data must be list of data frames
  ## breakpt must be list of vectors

  tic()  #start timer
  mod<- future_map2(data, breakpt, ~behav.gibbs.sampler(dat = .x, ngibbs = ngibbs, nbins = nbins,
                                                         alpha = alpha, breakpt = .y),
                   .progress = TRUE)
  toc()  #provide elapsed time


  brkpts<- map(mod, 1)  #create list of all sets breakpoints by ID

  nbrks<- map_dfr(mod, 2) %>%
    t() %>%
    data.frame()  #create DF of number of breakpoints by ID
  names(nbrks)<- c('id', paste0("Iter_",1:ngibbs))

  LML<- map_dfr(mod, 3) %>%
    t() %>%
    data.frame()  #create DF of LML by ID
  names(LML)<- c('id', paste0("Iter_",1:ngibbs))

  elapsed.time<- map_dfr(mod, 4) %>%
    t() %>%
    data.frame()  #create DF of elapsed time
  names(elapsed.time)<- "time"


  list(brkpts = brkpts, nbrks = nbrks, LML = LML, elapsed.time = elapsed.time)
}
