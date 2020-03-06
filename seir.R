##' @param x a vector of nodes
##' @param size number of nodes to pick at random
sample2 <- function(x, size) {
  if(length(x)==1) {
    rep(x, size)
  } else{
    sample(x, size, replace=TRUE)
  }
}

seir.full <- function(size,
                      beta,
                      sigma,
                      gamma,
                      I0,
                      seed = NULL,
                      imax,
                      keep.intrinsic=FALSE){
  
  if (!is.null(seed)) set.seed(seed)
  
  V <- 1:size
  
  if (missing(I0)) {
    if (missing(I0)) stop("specify the initial conditions")
    
  }
  initial_infected <- 1:I0
  
  if (missing(imax)) imax <- size
  
  queue_v <- queue_t <- queue_infector <- rep(NA, I0)
  
  queue_v[1:I0] <- initial_infected 
  queue_t[1:I0] <- 0
  
  t_infected <- t_infectious <- t_recovered <- rep(NA, size)
  t_infected[initial_infected] <- 0
  
  t_gillespie <- NULL
  c_infected <- 0
  
  if (keep.intrinsic) {
    intrinsic_generation <- vector('list', length(V))
  } else {
    intrinsic_generation <- NULL
  }
  
  done <- rep(FALSE, size)
  infected_by <- rep(NA, size)
  
  stop <- FALSE
  
  while (!stop) {
    j.index <- which.min(queue_t)
    j <- queue_v[j.index]
    
    infected_by[j] <- queue_infector[j.index]
    t_infected[j] <- queue_t[j.index]
    
    t <- queue_t[j.index]; t_gillespie <- c(t_gillespie, t)
    
    latent <- rexp(1, sigma)
    t_infectious[j] <- t + latent
    
    c_infected <- c_infected +1
    
    n <- V[V != j]
    if (!keep.intrinsic) n <- n[!done[n]]
    
    rate <- length(n) * beta + gamma
    
    prob <- gamma/rate
    
    ncontact <- rnbinom(1, size=1, prob=prob)
    if (ncontact > 0) {
      queue_v <- c(queue_v, sample2(n, ncontact))
      queue_infector <- c(queue_infector, rep(j, ncontact))
    }
    
    time_between <- rexp(ncontact+1, rate=rate)
    c_time <- latent + cumsum(time_between)
    
    if (keep.intrinsic) intrinsic_generation[[j]] <- c_time[1:ncontact]
    
    if (ncontact > 0) {
      queue_t <- c(queue_t, t + c_time[1:ncontact])
    }
    
    t_recovered[j] <- t+tail(c_time,1)
    done[j] <- TRUE
    
    filter2 <- !done[queue_v]
    queue_v <- queue_v[filter2]
    queue_infector <- queue_infector[filter2]
    queue_t <- queue_t[filter2]
    
    stop <- (c_infected == length(V) || all(done[queue_v]) || c_infected == imax)
  }
  
  return(
    list(
      data=data.frame(
        time=t_gillespie[(I0):c_infected],
        infected=(I0):c_infected
      ),
      intrinsic_generation=intrinsic_generation,
      t_infected=t_infected,
      t_infectious=t_infectious,
      t_recovered=t_recovered,
      infected_by=infected_by
    )
  )
}
