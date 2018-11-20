##### UTILITY FUNCTIONS FOR TTE SIMULATION #####
## This is adapted from recurrent events with combination studies
## to TITE-PK model for dose escalation study.

## Calculates H(t) = \int_0^t h(t) dt (we always calculate from t=0)!

## do the same calculation on the basis of a pre-calculated lauc
## (two-column matrix with a row per time-point)
calc_lcumH_pk_full <- function(lbeta_coef, time, lauc) {
  res <- cbind(log(time) + lbeta_coef[1] - log(24*7*4),
               sweep(lauc, 2, lbeta_coef[2:3], "+"))
  return(apply(res, 1, log_sum_exp))
}



## weibull version
calc_lcumH_pk_full_weibull <- function(lalpha_coef, lbeta_coef, time, lauc) {
  ltime <- log(time) - log(24*7*4)
  lambda <- exp(lalpha_coef[1])
  res <- cbind(lambda * (ltime + lalpha_coef[2]),
               sweep(lauc, 2, lbeta_coef[1:2], "+"))
  return(apply(res, 1, log_sum_exp))
}

calc_lh_pk_full_weibull <- function(lalpha_coef, lbeta_coef, time, lpk) {
  ltime <- log(time) - log(24*7*4)
  lambda <- exp(lalpha_coef[1])
  res <- cbind(sum(lalpha_coef) + (lambda-1) * (ltime + lalpha_coef[2]),
               sweep(lpk, 2, lbeta_coef[1:2], "+"))
  return(apply(res, 1, log_sum_exp))
}

## assume we have a constant hazard at any time, then the integral is
## equal to b-a
calc_lcumH_const <- function(lbeta_coef, t0, time) {
  return(lbeta_coef + log(time - t0))
}

## calculates t for given u in 0...1 for F(t)=1-exp(-H(t)) such that u
## = F(t)
calc_t_Finv <- function(u, t0, tmax, calc_lcumH_fn) {
  ##optimize(function(time) (u - (1-exp( - exp(calc_lcumH_fn(t0=t0, time=time) ))) )^2, c(t0, tmax))$minimum
  l1mu <- log1p(-u) ## log(1-u)
  fn_direct <- function(time) {
    u - (1-exp( - exp(calc_lcumH_fn(t0=t0, time=time) )))
  }
  fn <- function(time) {
    l1mu + exp(calc_lcumH_fn(t0=t0, time=time))
  }
  ## check that there is a sign change
  test <- fn(c(t0, tmax))
  if(prod(test) > 0) {
    if(test[2] > 0)
      return(tmax)
    else
      return(t0)
  }
  uniroot(fn, c(t0, tmax), f.lower=test[1], f.upper=test[2])$root
}

## given as input a function discrete at points x with its derivative,
## this function returns an index vector which can be used to reduce
## the number of knots when using a linear approximation to the
## function
find_knots <- function(x, f, df, tol=1e-3) {
  K <- length(x)
  knots <- rep(TRUE, K)
  f_approx <- function(x0, f0, df0) { function(x) f0 + (x-x0) * df0 }
  
  k <- 1
  while(k < K) {
    ## setup Taylor approximation at current knot
    f_taylor <- f_approx(x[k], f[k], df[k])
    k <- k + 1
    while(k < K &  abs(f_taylor(x[k]) - f[k]) < tol) {
      knots[k] <- FALSE
      k <- k + 1
    }
  }
  knots
}

if(FALSE) {
  ## attempt to cache results, not used anymore
  simulate_pop <- function() {
    sim_cache <- list()
    fun <- simulate_patient_cached
    ## make sure the curried function searches in the right enclosing
    environment(fun) <- environment()
    fun
  }
  simulate_patient <- function(data, pk_fun1, pk_fun2, lalpha_coef, lbeta_coef, leta)
    simulate_patient_cached(data, pk_fun1, pk_fun2, lalpha_coef, lbeta_coef, leta, cache=FALSE)
}


## key simulation function. The function does a piece-wise linear
## approximation of the cumulative hazard H(t) prior to simulating.
simulate_patient <- function(data, pk_fun, lbeta_coef, leta, cache=FALSE) {
  ##cache <- FALSE
  time_cens <- data$time[data$dv==1 & data$mdv==0]
  #time_cens <- data$time[data$dv==0 & data$mdv==0]
  
  j <- data$id[1]
  jc <- as.character(j)
  cat("Patient", j, "\n")
  
  sim_hist <- NULL
  if(cache) {
    sim_hist <- sim_cache[[jc]]
  }
  
  if(!is.null(sim_hist)) {
    cat("Using cached lH for patient", j, "...\n")
    lauc_fn <- sim_hist$lauc_fn
  } else {
    cat("Calculating lH for patient", j, "...\n")
    
    ## fields needed for dosing information
    dosing_fields <- c("time", "lamt", "cmt", "addl", "tau", "evid", "mdv")
    
    dosing <- subset(data, evid==1, sel=dosing_fields)
    
    ## pre-calculate H(t) from 0 to censoring time in dt steps
    ## (note that the reached precision is higher as we linearly
    ## interpolate)
    ##dt <- 7*24
    dt <- 24/2
    time_grid <- seq(0, time_cens + dt, by=dt)
    lpk <- pk_subject(time_grid, pk_fun, dosing)
    
    ## find knots which approximate well the AUC of each drug
    f1 <- exp(lpk[,3,drop=FALSE])
    df1 <- c(diff(f1), Inf)/2
    knots1 <- find_knots(time_grid, f1, df1)

    lauc_fn <- approxfun(time_grid[knots1], lpk[knots1, 3 , drop=FALSE], rule=2)

    if(cache) {
      cat("Knots reduction:", mean(knots1), " / ")
      sim_cache[[jc]] <<- list(lauc_fn = lauc_fn)
    }
  }
  
  lbeta_coef_j <- lbeta_coef

  Hfun <- function(time) {
    time_month <- time / (24*7*4)
    exp(lbeta_coef_j + lauc_fn(time))
  }
  
  t0     <- 0
  H_ref  <- 0
  events <- c()
  ## get a random 0-1 number
  u <- runif(1)
  l1mu <- log1p(-u) ## log(1-u)
    
    ## log-space calculation
  F_shifted <- function(x) { l1mu + Hfun(x) - H_ref }
    ## linear-space calulation
    ##F_shifted <- function(x) { u - ( 1-exp(-(Hfun(x) - H_ref)) ) }
    ## version for gradient based minimization (not as robust as root-finding, but faster)
    ##F_shifted_sq <- function(x) { (u - ( 1-exp(-(Hfun(x) - H_ref)) ))^2 }
    
    ## invert F such that F(t_ev) == u
  left <- l1mu
  right <- F_shifted(time_cens)
    
  if(left * right < 0){
      ##t_ev <- optimize(F_shifted_sq, c(t0, time_cens))$minimum
      t_ev <- uniroot(F_shifted, c(t0, time_cens), f.lower=left, f.upper=right, tol=0.001 )$root

      
  } else {
      t_ev <- 2*time_cens
  }
    
  if(t_ev < time_cens) {
      ##cat("Registering event for patient", j, "at", t_ev, "\n")
      events <- c(events, t_ev)
      ## set t0 to current event & update reference lH
      H_ref <- Hfun(t_ev)
  } else {
      ##cat("Event is past censoring time", t_ev, ">", time_cens, "\n")
  }
    

  res <- data.frame(time = time_cens, dv = 0)
  if(length(events) > 0) {
    res <- rbind(data.frame(time=events, dv=1), data.frame(time=events, dv=0))
  }
  res
}


cloglog <- binomial(link="cloglog")$linkfun
inv_cloglog <- binomial(link="cloglog")$linkinv

convert <- function(Pref) {
    lHref   <- cloglog(Pref)
    beta    <- lHref 
    list(beta = beta)
}


### END UTILITY FUNCTION BLOCK
