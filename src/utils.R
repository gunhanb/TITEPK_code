relabel <- function(id) match(id, unique(id))

normPrior <- function(q, p=c(0.025, 1-0.025)) {
    z <- qnorm(p)
    s <- abs(diff(q) / diff(z))
    m <- q[1] + abs(z[1]) * s
    c(m=m, s=s)
}

inv_logit <- function(l) 1/(1 + exp(-l))

log_inv_logit <- function(l) -log1p(exp(-l))

logit <- function(p) log(p/(1-p))

cv2sd <- function(cv) sqrt(log(1 + cv^2))

breakPlots <- function(pl) {
    function(d) {
        pl %+% d
    }
}

reblock <- function(ID, n=9) ceiling(relabel(ID)/n)
swapData <- function(ggPl) function(newData) ggPl %+% newData


## numerically stable summation of logs on the natural scale
log_sum_exp <- function(x){ 
   xmax <- which.max(x) 
   log1p(sum(exp(x[-xmax]-x[xmax])))+x[xmax] 
} 

log_softmax <- function(gamma) gamma - log_sum_exp(gamma)

## compress the dosing records
pack_dosing <- function(dosing) {
    ind_dose <- c()
    addl <- c()
    tau <- c()
    D <- nrow(dosing)
    r <- rle(round(diff(c(dosing$time, Inf)), 2))
    l <- r$lengths
    d <- r$values
    ## ensure that we dose in the very same cmt and amt all the time!
    ## Use the label_dose_periods function to find periods where this
    ## stays constant.
    assert_that(all(dosing$cmt[1] == dosing$cmt))
    assert_that(all(dosing$amt[1] == dosing$amt))
    j <- 1
    i <- 1
    while(i <= length(l)) {
        ind_dose <- c(ind_dose, j)
        if(l[i] == 1) {
            tau <- c(tau, 0)
            addl <- c(addl, 0)
            j <- j + 1
            i <- i + 1
        } else if(l[i] > 1) {                
            tau <- c(tau, d[i])
            ## in case the next record is not a repeated one (l[i+1]
            ## == 1), then we include it with the current repeated one
            if(i < length(l) & l[i+1] == 1) {
                addl <- c(addl, l[i])
                j <- j + l[i] + 1
                i <- i + 2
            } else if(i == length(l)) {
                addl <- c(addl, l[i])
                j <- j + l[i] + 1
                i <- i + 1
            } else {
                ## in case the next record is also one which has
                ## repetitions, we do not include it into the current
                ## set
                addl <- c(addl, l[i]-1)
                j <- j + l[i]
                i <- i + 1
            }
        }
    }
    dosing$addl <- NULL
    dosing$tau <- NULL
    transform(dosing[ind_dose,], addl=addl, tau=tau)
}


unpack_dosing <- function(packed) {
    D <- nrow(packed)
    packed$DID <- 1:D
    unpacked <- packed[rep(1:D, times=packed$addl+1),]
    unpacked <- ddply(unpacked, .(DID), transform, time=time[1] + (0:addl[1]) * tau[1])
    transform(unpacked, DID=NULL, addl=0, tau=0)
}

label_dose_periods <- function(cmt, amt, ...) {
    dose_id <- do.call(paste, c(list(...), list(cmt, amt, sep="/")))
    rr <- rle(dose_id)
    L <- length(rr$lengths)
    rep(1:L, times=rr$lengths)
}

## function which is given a dosing history and returns at a given
## vector of times the system state
pk_subject <- function(time, pk_fun, dosing, ...) {
    ## time is expected to be a vector of times at which pk_fun values
    ## will be returned.  dosing must contain all the dosing
    ## information of a subject with columns time, lamt & cmt. pk_fun
    ## is called with an id column = 1, all needed dosing events, a
    ## respective evid column and an out vector which is marked 1 for
    ## the times requested

    dosing <- subset(dosing, evid==1 & mdv==1 & lamt > -30)
    
    names(dosing) <- paste("dose", names(dosing), sep="_")

    f <- c("dose_lamt", "dose_cmt", "dose_time", "dose_tau", "dose_addl")
    dosing <- dosing[,f]

    do.call(pk_fun, c(dosing, list(obs_time=time, init_time=0, init_lstate=rep(-35, 3), x_r=c(0), x_i=c(0L), ...)))
}
