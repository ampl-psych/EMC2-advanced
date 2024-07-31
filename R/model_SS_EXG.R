
#' Staircase algorithm
#'
#' @return A function to simulate staircase SSD data for
#' @export
staircase_fun <- function(dadm,p=.25,pars=NULL,
                      SSD0=.25,stairstep=.05,stairmin=0,stairmax=Inf)
  # random p of trials get NA, ready to be filled in by a staircase
{
  # if pars not supplied return an SSD column indicating all go trials
  if (is.null(pars)) return(rep(Inf,nrow(dadm)))

  # levels(dadm$lR) <- levels(dadm$lR)

  # if pars supplied potentially run staircase for each participant
  for(i in 1:length(unique(dadm$subjects))){
    dadm_i <- dadm[dadm$subjects==unique(dadm$subjects)[i],]
    pars_i <- pars[dadm$subjects==unique(dadm$subjects)[i],]
    nacc <- length(levels(dadm_i$lR))
    is1 <- dadm_i$lR==levels(dadm_i$lR)[1] # first accumulator
    ntrials <- sum(is1)
    if (!any(colnames(dadm_i)=="SSD")) { # pick staircase trials
      dadm_i$SSD <- rep(Inf, nrow(dadm_i))
      tmp <- matrix(dadm_i$SSD,nrow=nacc)
      tmp[,sample(1:ntrials,round(ntrials*p))] <- NA
      dadm_i$SSD <- as.vector(tmp)
    }
    out_i <- setNames(data.frame(matrix(NA,nrow=ntrials,ncol=3)),c("R","rt","SSD"))
    pick <- is.infinite(dadm_i$SSD)
    if (any(pick)) { # fill in no-stop trials
      out_i$SSD[pick[is1]] <- Inf
      pars_i[pick,"SSD"] <- Inf
      out_i[pick[is1],c("R","rt")] <- attributes(dadm_i)$model()$rfun(dadm_i$lR[pick],pars_i[pick,,drop=FALSE])
    }
    pick <- is.finite(dadm_i$SSD)
    if (any(pick)) { # fill in fixed SSD trials
      out_i$SSD[pick[is1]] <- dadm_i$SSD[pick][is1[pick]]
      pars_i[pick,"SSD"] <- dadm_i$SSD[pick]
      out_i[pick[is1],c("R","rt")] <- attributes(dadm_i)$model()$rfun(dadm_i$lR[pick],pars_i[pick,,drop=FALSE])
    }
    isna <- is.na(dadm_i$SSD)
    if (any(isna)) { # run staircase if any NAs to fill in
      nstair <- sum(isna)/nacc # number of staircase trials
      trials <- rep(0,ntrials*nacc) # used to pick out_i each staircase trial in pars_i
      trials[is.na(dadm_i$SSD)] <- rep(1:nstair,each=nacc)
      for (s in 1:nstair) { # run staircase
        current <-  trials == s # current staircase trial
        if (s==1)  dadm_i$SSD[current] <- out_i$SSD[current[is1]] <- SSD0 # initialize
        p_stairs <- pars_i[current,,drop=FALSE] # parameters for current staircase trial
        # simulate 1 trial, because is.na(SSD) in pars_i rfun returns dt, an nacc x 1 matrix
        dt <- attributes(dadm_i)$model()$rfun(dadm_i$lR[current],p_stairs)
        inhibit <- p_stairs[,"lI"]==1 # inhibition triggered
        # add SSD to stop and inhibition triggered
        dt[c(TRUE,inhibit),] <- dt[c(TRUE,inhibit),] + dadm_i$SSD[current][1]
        winner <- which.min(dt)
        iinhibit <- c(1,1+c(1:nacc)[inhibit]) # stop or inhibition triggered index
        if (s != nstair) { # set SSD for next staircase trial
          nexts <- trials == (s+1)
          if (any(iinhibit==winner))   # round as otherwise get spurious tiny differences
            dadm_i$SSD[nexts] <- round(dadm_i$SSD[current] + stairstep,3) else  # successful stop
              dadm_i$SSD[nexts] <- round(dadm_i$SSD[current] - stairstep,3)       # failed stop
            if ((dadm_i$SSD[nexts][1]<stairmin) | (dadm_i$SSD[nexts][1]>stairmax))
              dadm_i$SSD[nexts] <- dadm_i$SSD[current] # dont step
            out_i$SSD[c(trials == s+1)[is1]] <- dadm_i$SSD[nexts][1]
        }
        if (winner==1) { # stop wins
          if (any(inhibit)) { # inhibition triggered response
            if (all(is.infinite(dt[c(FALSE,inhibit),])))  # tf
              out_i[current[is1],c("R","rt")] <- c(NA,NA) else {
                pick <- which.min(dt[c(FALSE,inhibit),])
                out_i[current[is1],c("R","rt")] <-
                  c(c(1:nacc)[inhibit][pick],dt[c(FALSE,inhibit),][pick])
              }
          } # otherwise no response
        } else { # pick from all except stop
          out_i[current[is1],c("R","rt")] <- c(winner-1,dt[winner,])
        }
      }
    }
    out_i$R <- factor(out_i$R,levels=1:nacc,labels=levels(dadm_i$lR))
    if(i == 1){
      out <- out_i
    } else {
      out <- rbind(out, out_i)
    }
    if(length(unique(dadm$subjects)) == 1){
      return(out)
    }
  }
  return(out)
}


#### ExGaussian ----

dexGaussian <- function(rt,pars)
{
  isexp <- pars[,"sigma"] < 1e-4 # shifted exponential
  rt[isexp] <- dexp(rt[isexp]-pars[isexp,"mu"],1/pars[isexp,"tau"])
  isnorm <- !isexp & pars[,"tau"] < 0.05 * pars[,"sigma"] # normal
  rt[isnorm] <- dnorm(rt[isnorm], mean = pars[isnorm,"mu"], sd = pars[isnorm,"sigma"])
  isexg <- !(isexp | isnorm)
  if (any(isexg)) {
    s2 <- pars[isexg,"sigma"]^2
    z <- rt[isexg] - pars[isexg,"mu"] - (s2/pars[isexg,"tau"])
    rt[isexg] <- exp(
      log(pnorm(z/pars[isexg,"sigma"])) -
        log(pars[isexg,"tau"]) -
        (z + (s2/(2 *  pars[isexg,"tau"])))/pars[isexg,"tau"]
    )
  }
  rt
}

pexGaussian <- function(rt,pars)
  # cumulative density for single accumulator
{
  isexp <- pars[,"sigma"] < 1e-4 # shifted exponential
  rt[isexp] <- pexp(rt[isexp]-pars[isexp,"mu"],1/pars[isexp,"tau"])
  isnorm <- !isexp & pars[,"tau"] < 0.05 * pars[,"sigma"] # normal
  rt[isnorm] <- pnorm(rt[isnorm], mean = pars[isnorm,"mu"], sd = pars[isnorm,"sigma"])
  isexg <- !(isexp | isnorm)
  if (any(isexg)) {
    s2 <- pars[isexg,"sigma"]^2
    z <- rt[isexg] - pars[isexg,"mu"] - (s2/pars[isexg,"tau"])
    rt[isexg] <-
      pnorm((rt[isexg] - pars[isexg,"mu"])/pars[isexg,"sigma"]) -
      exp(log(pnorm(z/pars[isexg,"sigma"])) +
            ((pars[isexg,"mu"] + (s2/pars[isexg,"tau"]))^2 - (pars[isexg,"mu"]^2) -
               2 * rt[isexg] * (s2/pars[isexg,"tau"]))/(2 * s2))
  }
  rt
}

rexGaussian <- function(lR,pars,p_types=c("mu","sigma","tau"))
  # lR is an empty latent response factor lR with one level for each accumulator.
  # pars is a matrix of corresponding parameter values named as in p_types
  # pars must be sorted so accumulators and parameter for each trial are in
  # contiguous rows.
  #
  # test
  # pars=cbind(mu=c(.5,.6),sigma=c(.1,.1),tau=c(.2,.2)); lR=factor(c(1))
{
  if (!all(p_types %in% dimnames(pars)[[2]]))
    stop("pars must have columns ",paste(p_types,collapse = " "))
  dt <- matrix(
    rnorm(dim(pars)[1],mean=pars[,"mu"],sd=pars[,"sigma"]) +
      rexp(dim(pars)[1],rate=1/pars[,"tau"]),nrow=length(levels(lR)))
  R <- apply(dt,2,which.min)
  pick <- cbind(R,1:dim(dt)[2]) # Matrix to pick winner
  rt <- dt[pick]
  R <- factor(levels(lR)[R],levels=levels(lR))
  cbind.data.frame(R=R,rt=rt)
}

# Go cdf/pdf versions

dexGaussianG <- function(rt,pars)
{
  out <- numeric(length(rt))
  ok <- !is.na(rt)
  out[ok] <- dexGaussian(rt[ok],pars[ok,,drop=FALSE])
  out
}

pexGaussianG <- function(rt,pars)
{
  out <- numeric(length(rt))
  ok <- !is.na(rt)
  out[ok] <- pexGaussian(rt[ok],pars[ok,,drop=FALSE])
  out
}

# Stop cdf/pdf versions

dexGaussianS <- function(rt,pars)
{
  rt <- rt - pars[,"SSD"]
  dimnames(pars)[[2]][dimnames(pars)[[2]]=="muS"] <- "mu"
  dimnames(pars)[[2]][dimnames(pars)[[2]]=="sigmaS"] <- "sigma"
  dimnames(pars)[[2]][dimnames(pars)[[2]]=="tauS"] <- "tau"
  dexGaussian(rt,pars)
}


pexGaussianS <- function(rt,pars)
{
  rt <- rt - pars[,"SSD"]
  dimnames(pars)[[2]][dimnames(pars)[[2]]=="muS"] <- "mu"
  dimnames(pars)[[2]][dimnames(pars)[[2]]=="sigmaS"] <- "sigma"
  dimnames(pars)[[2]][dimnames(pars)[[2]]=="tauS"] <- "tau"
  pexGaussian(rt,pars)
}


# Stop signal random

rSSexGaussian <- function(lR,pars)
  # lR is an empty latent response factor lR with one level for each accumulator.
  # pars must contain an SSD column and an lI column. If SSD contains any
  # NAs then return the dt matrix (for use in staircase creation), else return
  # the usual Rrt data frame.
  # NB1: Go failures will only apply to accumulators where lI = TRUE
  #      and can still have a stop-triggered response on a go-failure trial.
{

  nacc <- length(levels(lR)) # Does not include stop runner
  ntrials <- dim(pars)[1]/nacc
  is1 <- lR==levels(lR)[1]
  acc <- 1:nacc

  # stop-triggered racers
  isST <- pars[,"lI"]==1
  accST <- acc[pars[1:nacc,"lI"]==1]

  # Default Inf finishing time so if not changed always looses
  dt <- matrix(Inf,nrow=nacc+1,ncol=ntrials)

  # Go failures
  isgf <- rep(pars[is1,"gf"] > runif(ntrials),each=nacc)
  # Expand to match go accumulators that don't fail
  isGO <- !isgf & !isST
  ngo <- sum(isGO)

  # Fill in go accumulators
  if (any(isGO)) dt[-1,][isGO] <-
    rnorm(ngo,mean=pars[isGO,"mu"],sd=pars[isGO,"sigma"]) +
    rexp(ngo,rate=1/pars[isGO,"tau"])

  # pick out stop trials and races with SSD that is not Inf (i.e., finite or
  # NA, the latter so a staircase can be filled in)
  isStrial <- !is.infinite(pars[is1,"SSD"])
  isSrace <- rep(isStrial,each=nacc)

  # pick out stop trials that are triggered
  isT <- pars[is1,"tf"][isStrial] < runif(sum(isStrial))

  # Logical to pick stop-triggered accumulators that are triggered
  isSTT <- logical(ntrials*nacc)

  # Pick out stop-triggered accumulators that are triggered
  isSTT[isSrace][rep(isT,each=nacc) & isST[isSrace]] <- TRUE
  nst <- sum(isSTT)

  # Fill in stop-triggered accumulators
  if (any(isSTT)) dt[-1,][isSTT] <-
    rnorm(nst,mean=pars[isSTT,"mu"],sd=pars[isSTT,"sigma"]) +
    rexp(nst,rate=1/pars[isSTT,"tau"])

  # pick out triggered stop racers
  isTS <- logical(ntrials)
  isTS[isStrial][isT] <- TRUE
  ns <- sum(isTS)

  # Fill in stop accumulators
  if (any(isTS)) dt[1,isTS] <-
    rnorm(ns,mean=pars[is1,"muS"][isTS],sd=pars[is1,"sigmaS"][isTS]) +
    rexp(ns,rate=1/pars[is1,"tauS"][isTS])

  # return dt to be used by a staircase algorithm
  if (any(is.na(pars[,"SSD"]))) return(dt)

  if (any(isTS)) dt[1,isTS] <- dt[1,isTS] + pars[is1,"SSD"][isTS]
  if (any(isSTT)) dt[-1,][isSTT] <- dt[-1,][isSTT] + pars[isSTT,"SSD"]

  # All SSD already filled in so return R and rt

  # R <- factor(rep(NA,ntrials),levels=levels(lR))
  R <- rt <- rep(NA,ntrials)

  # All accumulators Inf (usually when both go and tf)
  allinf <- apply(dt,2,\(x)all(is.infinite(x)))

  # get winner of stop and go where there is a race
  r <- c(1, 1 + acc)[apply(dt[,!allinf,drop=FALSE],2,which.min)]

  # stop wins
  stopwins <- r==1

  # First fill in cases where stop looses
  if (any(!stopwins)) {
    rgo <- r[!stopwins]-1
    R[!allinf][!stopwins] <- rgo
    pick <- cbind(rgo,c(1:sum(!stopwins))) # Matrix to pick winner
    rt[!allinf][!stopwins] <- dt[-1,!allinf,drop=FALSE][,!stopwins,drop=FALSE][pick]
  }

  # then if stop triggers extra accumulators find their winner
  if (any(isST) & any(stopwins)) {
    # stop triggered accumulators that are racing
    rst <- dt[-1,!allinf,drop=FALSE][accST,stopwins,drop=FALSE]
    # stop-triggered winners
    rtw <- apply(rst,2,which.min)
    # index for stop-triggered
    R[!allinf][stopwins] <- accST[rtw]
    pick <- cbind(rtw,1:ncol(rst))
    rt[!allinf][stopwins] <- rst[pick]
  }

  rt[is.na(R)] <- NA
  cbind.data.frame(R=factor(R,levels=1:nacc,labels=levels(lR)),rt=rt) #,SSD=SSD)
}



# Following functions moved to C++ model_SS_EXG.cpp

# pEXG <- function (q, mu = 5, sigma = 1, tau = 1, lower_tail = TRUE, log_p = FALSE)
#     # ex-Gaussian cumulative density
#     # Modified from gamlss.dist to make cdf in tau > 0.05 * sigma case robust,
#     # and robust to -Inf and Inf inputs, returns NA for bad sigma or tau and
#     # robust against small sigma cases.
#   {
#     if (sigma <= 0) return(rep(NA,length(q)))
#       if (tau <= 0) return(rep(NA,length(q)))
#
#       # if (sigma < 0.05*tau)
#       if (sigma < 1e-4)
#        return(pexp(q-mu,1/tau,log.p=log_p,lower.tail=lower_tail)) # shfited exponential
#
#     ly <- length(q)
#     sigma <- rep(sigma, length = ly)
#     mu <- rep(mu, length = ly)
#     tau <- rep(tau, length = ly)
#     index <- seq(along = q)
#     z <- q - mu - ((sigma^2)/tau)
#     cdf <- ifelse(is.finite(q),
#       ifelse(tau > 0.05 * sigma,
#          pnorm((q - mu)/sigma) - exp(log(pnorm(z/sigma)) + ((mu + (sigma^2/tau))^2 -
#            (mu^2) -  2 * q * ((sigma^2)/tau))/(2 * sigma^2)),
#          pnorm(q, mean = mu, sd = sigma)),
#         ifelse(q<0,0,1)
#     )
#     if (lower_tail == TRUE)
#       cdf <- cdf
#     else cdf <- 1 - cdf
#     if (log_p == FALSE)
#       cdf <- cdf
#     else cdf <- log(cdf)
#     cdf
#   }
#
# dEXG <- function (x, mu = 5, sigma = 1, tau = 1, log = FALSE)
#   # ex-Gaussian density
#   # gamlss.dist function, but returns NA for bad sigma or tau, and
#   # robust against small sigma cases.
# {
#     if (sigma <= 0) return(rep(NA,length(x)))
#     if (tau <= 0) return(rep(NA,length(x)))
#
#     # if (sigma < 0.05*tau)
#     if (sigma < 1e-4)
#       return(dexp(x-mu,1/tau,log=log)) # shfited exponential
#
#     ly <- length(x)
#     sigma <- rep(sigma, length = ly)
#     mu <- rep(mu, length = ly)
#     tau <- rep(tau, length = ly)
#     z <- x - mu - ((sigma^2)/tau)
#     logfy <- ifelse(tau > 0.05 * sigma,
#       -log(tau) - (z + (sigma^2/(2 *  tau)))/tau + log(pnorm(z/sigma)),
#       dnorm(x, mean = mu, sd = sigma, log = TRUE))
#     if (log == FALSE)
#       fy <- exp(logfy)
#     else fy <- logfy
#     fy
# }
#
# dEXGrace <- function(dt,mu,sigma,tau)
#   # Generates defective PDF for win by first runner, dt (decison time) is
#   # a matrix with length(mu) rows, one row for each runner, and one column
#   # for each decision time for which a defective density value will be
#   # returned.
# {
#   dt[1,] <- dEXG(dt[1,],mu[1],sigma[1],tau[1])
#   if (length(mu)>1) for (i in 2:length(mu))
#     dt[1,] <- dt[1,]*pEXG(dt[i,],mu[i],sigma[i],tau[i],lower_tail=FALSE)
#   dt[1,]
# }
#
#
# stopfn_exg <- function(t,mu,sigma,tau,SSD)
#   # Used by my.integrate, t = vector of times, SSD is a scalar stop-signal delay.
# {
#   dt <- matrix(rep(t+SSD,each=length(mu)),nrow=length(mu))
#   dt[1,] <- dt[1,]-SSD
#   dEXGrace(dt,mu,sigma,tau)
# }

pstopEXG <- function(parstop,n_acc,upper=Inf,
                     gpars=c("mu","sigma","tau"),spars=c("muS","sigmaS","tauS"))
{
  sindex <- seq(1,nrow(parstop),by=n_acc)
  ps <- parstop[sindex,spars,drop=FALSE]
  SSDs <- parstop[sindex,"SSD",drop=FALSE]
  ntrials <- length(SSDs)
  if (length(upper)==1) upper <- rep(upper,length.out=ntrials)
  pgo <- array(parstop[,gpars],dim=c(n_acc,ntrials,length(gpars)),
               dimnames=list(NULL,NULL,gpars))
  cells <- apply(cbind(SSDs,ps,upper,
    matrix(as.vector(aperm(pgo,c(2,1,3))),nrow=ntrials)),1,paste,collapse="")
  # cells <- character(ntrials)
  # for (i in 1:ntrials)
  #   cells[i] <- paste(SSDs[i],ps[i,],pgo[,i,],upper[i],collapse="")
  uniq <- !duplicated(cells)
  ups <- sapply(1:sum(uniq),function(i){
    my.integrate(f=stopfn_exg,lower=-Inf,SSD=SSDs[i],upper=upper[i],
                           mu=c(ps[i,"muS"],pgo[,i,"mu"]),
                           sigma=c(ps[i,"sigmaS"],pgo[,i,"sigma"]),
                           tau=c(ps[i,"tauS"],pgo[,i,"tau"]))
  })
  ups[as.numeric(factor(cells,levels=cells[uniq]))]
}

stopfn_exgST <- function(t,mu,sigma,tau,SSD,st=1)
  # Used by my.integrate, t = vector of times, SSD is a scalar stop-signal delay.
  # st is a vector of indices for the stop and stop-triggered accumulators
{
  dt <- matrix(rep(t+SSD,each=length(mu)),nrow=length(mu))
  dt[st,] <- dt[st,]-SSD
  dEXGrace(dt,mu,sigma,tau)
}

pstopEXGST <- function(parstop,n_acc,upper=Inf,st=1,
                     gpars=c("mu","sigma","tau"),spars=c("muS","sigmaS","tauS"))
{
  sindex <- seq(1,nrow(parstop),by=n_acc)
  ps <- parstop[sindex,spars,drop=FALSE]
  SSDs <- parstop[sindex,"SSD",drop=FALSE]
  ntrials <- length(SSDs)
  if (length(upper)==1) upper <- rep(upper,length.out=ntrials)
  pgo <- array(parstop[,gpars],dim=c(n_acc,ntrials,length(gpars)),
               dimnames=list(NULL,NULL,gpars))
  cells <- apply(cbind(SSDs,ps,upper,
    matrix(as.vector(aperm(pgo,c(2,1,3))),nrow=ntrials)),1,paste,collapse="")
  # cells <- character(ntrials)
  # for (i in 1:ntrials)
  #   cells[i] <- paste(SSDs[i],ps[i,],pgo[,i,],upper[i],collapse="")
  uniq <- !duplicated(cells)
  ups <- sapply(1:sum(uniq),function(i){
    my.integrate(f=stopfn_exgST,lower=-Inf,SSD=SSDs[i],upper=upper[i],
                           mu=c(ps[i,"muS"],pgo[,i,"mu"]),
                           sigma=c(ps[i,"sigmaS"],pgo[,i,"sigma"]),
                           tau=c(ps[i,"tauS"],pgo[,i,"tau"]),st=st)
  })
  ups[as.numeric(factor(cells,levels=cells[uniq]))]
}



#### Model list ----
#' Stop-signal exGaussian race
#'
#' @return A model list with all the necessary functions to sample
#' @export
SSexG <- function() {
  list(
    type="RACE",
    p_types=c(mu=log(.4),sigma=log(.05),tau=log(.1),
              muS=log(.3),sigmaS=log(.025),tauS=log(.05),tf=qnorm(0),gf=qnorm(0)),
    Ntransform=function(x,use=NULL) {
      # transform parameters back to real line
      isprobit <- dimnames(x)[[2]] %in% c("tf","gf")
      if (is.null(use)) {
        x[,!isprobit] <- exp(x[,!isprobit])
        x[,isprobit] <- pnorm(x[,isprobit])
      } else {
        ok <- dimnames(x)[[2]] %in% use
        x[,!isprobit & ok] <- exp(x[,!isprobit & ok])
        x[,isprobit & ok] <- pnorm(x[,isprobit & ok])
      }
      x
    },
    # p_vector transform
    transform = function(x) x,
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      if (any(names(dadm)=="SSD")) pars <- cbind(pars,SSD=dadm$SSD) else
        pars <- cbind(pars,SSD=rep(NA,dim(pars)[1]))
      pars <- cbind(pars,lI=as.numeric(dadm$lI))  # Only necessary for data generation.
      attr(pars,"ok") <- (pars[,"tau"] > 1e-3) & (pars[,"sigma"] > 1e-3) & (pars[,"mu"] > 1e-3) &
            (pars[,"tau"] < 1) & (pars[,"sigma"] < 1) &
            (pars[,"tauS"] > 1e-3) & (pars[,"sigmaS"] > 1e-3) & (pars[,"muS"] > 1e-3) &
            (pars[,"tauS"] < 1) & (pars[,"sigmaS"] < 1) &
            ((pars[,"tf"] > 1e-6) | pars[,"tf"] == 0) & ((pars[,"gf"] > 1e-6) | pars[,"gf"] == 0)
      pars
    },
    # Density function (PDF) for single go racer
    dfunG=function(rt,pars) dexGaussianG(rt,pars),
    # Probability function (CDF) for single go racer
    pfunG=function(rt,pars) pexGaussianG(rt,pars),
    # Density function (PDF) for single stop racer
    dfunS=function(rt,pars) dexGaussianS(rt,pars[,c("muS","sigmaS","tauS","SSD"),drop=FALSE]),
    # Probability function (CDF) for single stop racer
    pfunS=function(rt,pars) pexGaussianS(rt,pars[,c("muS","sigmaS","tauS","SSD"),drop=FALSE]),
    # Stop probability integral
    sfun=function(pars,n_acc,st=1,upper=Inf) pstopEXGST(pars,n_acc,upper=upper,st=st),
    # Random function for SS race
    rfun=function(lR=NULL,pars) {
      ok <- (pars[,"tau"] > 1e-3) & (pars[,"sigma"] > 1e-3) & (pars[,"mu"] > 1e-3) &
            (pars[,"tau"] < 1) & (pars[,"sigma"] < 1) &
            (pars[,"tauS"] > 1e-3) & (pars[,"sigmaS"] > 1e-3) & (pars[,"muS"] > 1e-3) &
            (pars[,"tauS"] < 1) & (pars[,"sigmaS"] < 1) &
            ((pars[,"tf"] > 1e-6) | pars[,"tf"] == 0) & ((pars[,"gf"] > 1e-6) | pars[,"gf"] == 0)

      if (is.null(lR)) ok else rSSexGaussian(lR,pars)
    },
    # Race likelihood combining pfun and dfun
    log_likelihood=function(p_vector,dadm,min_ll=log(1e-10))
      log_likelihood_race_ss(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
  )
}

#### Model list ----
#' Stop-signal exGaussian race with uni-valent stop-triggered responding
#'
#' @return A model list with all the necessary functions to sample
#' @export
SSexGuv <- function() {
  list(
    type="RACE",
    p_types=c(mu=log(.4),sigma=log(.05),tau=log(.1),
              muS=log(.3),sigmaS=log(.025),tauS=log(.05),tf=qnorm(0),gf=qnorm(0)),
    Ntransform=function(x,use=NULL) {
      # transform parameters back to real line
      isprobit <- dimnames(x)[[2]] %in% c("tf","gf")
      if (is.null(use)) {
        x[,!isprobit] <- exp(x[,!isprobit])
        x[,isprobit] <- pnorm(x[,isprobit])
      } else {
        ok <- dimnames(x)[[2]] %in% use
        x[,!isprobit & ok] <- exp(x[,!isprobit & ok])
        x[,isprobit & ok] <- pnorm(x[,isprobit & ok])
      }
      x
    },
    # p_vector transform
    transform = function(x) x,
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      if (any(names(dadm)=="SSD")) pars <- cbind(pars,SSD=dadm$SSD) else
        pars <- cbind(pars,SSD=rep(NA,dim(pars)[1]))
      attr(pars,"ok") <- (pars[,"tau"] > 1e-3) & (pars[,"sigma"] > 1e-3) & (pars[,"mu"] > 1e-3) &
            (pars[,"tau"] < 1) & (pars[,"sigma"] < 1) &
            (pars[,"tauS"] > 1e-3) & (pars[,"sigmaS"] > 1e-3) & (pars[,"muS"] > 1e-3) &
            (pars[,"tauS"] < 1) & (pars[,"sigmaS"] < 1) &
            ((pars[,"tf"] > 1e-6) | pars[,"tf"] == 0) & ((pars[,"gf"] > 1e-6) | pars[,"gf"] == 0)

      pars <- cbind(pars,lI=as.numeric(dadm$lI))  # Only necessary for data generation.
      pars
    },
    # Density function (PDF) for single go racer
    dfunG=function(rt,pars) dexGaussianG(rt,pars),
    # Probability function (CDF) for single go racer
    pfunG=function(rt,pars) pexGaussianG(rt,pars),
    # Density function (PDF) for single stop racer
    dfunS=function(rt,pars) dexGaussianS(rt,pars[,c("muS","sigmaS","tauS","SSD"),drop=FALSE]),
    # Probability function (CDF) for single stop racer
    pfunS=function(rt,pars) pexGaussianS(rt,pars[,c("muS","sigmaS","tauS","SSD"),drop=FALSE]),
    # Stop probability integral
    sfun=function(pars,n_acc,upper=Inf) pstopEXG(pars,n_acc,upper=upper),
    # Random function for SS race
    rfun=function(lR=NULL,pars) {
      ok <- (pars[,"tau"] > 1e-3) & (pars[,"sigma"] > 1e-3) & (pars[,"mu"] > 1e-3) &
            (pars[,"tau"] < 1) & (pars[,"sigma"] < 1) &
            (pars[,"tauS"] > 1e-3) & (pars[,"sigmaS"] > 1e-3) & (pars[,"muS"] > 1e-3) &
            (pars[,"tauS"] < 1) & (pars[,"sigmaS"] < 1) &
            ((pars[,"tf"] > 1e-6) | pars[,"tf"] == 0) & ((pars[,"gf"] > 1e-6) | pars[,"gf"] == 0)

      if (is.null(lR)) ok else rSSexGaussian(lR,pars)
    },
    # Race likelihood combining pfun and dfun
    log_likelihood=function(p_vector,dadm,min_ll=log(1e-10))
      log_likelihood_race_ss_uv(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
  )
}
