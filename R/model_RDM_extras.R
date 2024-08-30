#' The Racing Diffusion Model (RDM) with no C
#'
#' @return A list defining the cognitive model
#' @export
RDMnoC <- function(){
  list(
    type="RACE",
    p_types=c("v" = log(1),"B" = log(1),"A" = log(0),"t0" = log(0),"s" = log(1)),
    # Transform to natural scale
    Ntransform=function(x,use=NULL) {
      # transform parameters back to real line
      if (is.null(use)) {
        x <- exp(x)
      } else {
        x[,use] <- exp(x[,use])
      }
      x
    },
    # p_vector transform
    transform = function(x) x,
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      if (!is.null(attr(dadm,"adaptive"))) pars <- do_adaptive(pars,dadm)
       attr(pars,"ok") <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0) &
        ((pars[,"v"] > 1e-3) | pars[,"v"] == 0) &  (pars[,"B"] >= 0) &  (pars[,"s"] > 0)
      pars
    },
    # Random function for racing accumulators
    rfun=function(lR=NULL,pars) {
      ok <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0) &
        ((pars[,"v"] > 1e-3) | pars[,"v"] == 0) &  (pars[,"B"] >= 0) &  (pars[,"s"] > 0)
      if (is.null(lR)) ok else rRDM(lR,pars,ok=ok)
    },
    # Density function (PDF) for single accumulator
    dfun=function(rt,pars) dRDM(rt,pars),
    # Probability function (CDF) for single accumulator
    pfun=function(rt,pars) pRDM(rt,pars),
    # Race likelihood combining pfun and dfun
    log_likelihood=function(p_vector,dadm,min_ll=log(1e-10))
      log_likelihood_race(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
  )
}



#### Advantage ----

rARDM <- function(lR,pars,p_types=c("v","B","A","t0"),ok=rep(TRUE,dim(pars)[1]))
  # lR is an empty latent response factor lR with one level for each accumulator.
  # It is the length of the original data frame, whereas pars is rep'd
  # na=nr-1 times, where nr=length(levels(lR)).
  # pars is a matrix of parameter values named as in p_types already rep'd
  # pars must be sorted so na accumulators within responses fastest then
  # responses next, then each trial in contiguous rows. "s" parameter will be
  # used but can be omitted
{
  if (length(lR)==nrow(pars)) return(rRDM(lR,pars,p_types,ok)) else {
    nr <- length(levels(lR))
    na <- nr-1
    nt <- nrow(pars)/(na*nr)
    lR <- rep(lR,each=na)
    if (!all(p_types %in% dimnames(pars)[[2]]))
      stop("pars must have columns ",paste(p_types,collapse = " "))
    if (any(dimnames(pars)[[2]]=="s")) # rescale
      pars[,c("A","B","v")] <- pars[,c("A","B","v")]/pars[,"s"]
    pars[,"B"][pars[,"B"]<0] <- 0 # Protection for negatives
    pars[,"A"][pars[,"A"]<0] <- 0
    bad <- rep(NA, nt)
    out <- data.frame(R = bad, rt = bad)
    dt <- array(Inf,dim=c(na,nr,nt))
    t0 <- array(pars[,"t0"],dim=c(na,nr,nt))[1,,]
    pars <- pars[ok,]
    dt[ok] <- rWald(sum(ok),B=pars[,"B"],v=pars[,"v"],A=pars[,"A"])
    dt <- apply(dt,2:3,max)
    R <- apply(dt,2,which.min)
    pick <- cbind(R,1:dim(dt)[2]) # Matrix to pick winner
    # Any t0 difference with lR due to response production time (no effect on race)
    rt <- t0[pick] + dt[pick]
    ok <- matrix(ok,nrow=nr*na)[1,]
    out$R <- levels(lR)[R]
    out$R <- factor(out$R,levels=levels(lR))
    out$rt <- rt
  }
  out
}


#' Advantage Racing Diffusion Model, pbeta mapping
#'
#' @return A list defining the cognitive model
#' @export
ARDM <- function(){
  list(
    type="RACE",
    c_name = "ARDM",
    p_types=c("v" = log(1),"B" = log(1),"A" = log(0),"t0" = log(0),"s" = log(1),
              "AD"=log(1),"AS"=log(0),SD=0,SS=0,DD=0,DS=0),
    # Transform to natural scale
    Ntransform=function(x,use=NULL) {
      # transform parameters back to real line
      if (is.null(use)) {
        use <- !(colnames(x) %in% c("SD","SS","DD","DS"))
        x[,use] <- exp(x[,use])
      } else {
        use <- use[!(use %in% c("SD","SS","DD","DS"))] # never transform advantage shape
        x[,use] <- exp(x[,use])
      }
      x
    },
    # p_vector transform
    transform = function(x) x,
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      if (!is.null(attr(dadm,"adaptive"))) pars <- do_adaptive(pars,dadm)
      if (attr(dadm,"advantage")[[1]] == "parameter")
        pars <- advantage_parameters(pars,pars[,names(attr(dadm,"advantage"))[1]],
          length(levels(dadm$lR)),names(attr(dadm,"advantage"))[1]) else
        pars <- advantage_pars(pars,dadm[attr(dadm,"advantage")[[1]]],
          length(levels(dadm$lR)),names(attr(dadm,"advantage"))[1])
      attr(pars,"ok") <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0) &
        ((pars[,"v"] > 1e-3) | pars[,"v"] == 0) &  (pars[,"B"] >= 0) &  (pars[,"s"] > 0)
      pars
    },
    # Random function for racing accumulators
    rfun=function(lR=NULL,pars) {
      ok <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0) &
        ((pars[,"v"] > 1e-3) | pars[,"v"] == 0) &  (pars[,"B"] >= 0) &  (pars[,"s"] > 0)
      if (is.null(lR)) ok else rARDM(lR,pars,ok=ok)
    },
    # Density function (PDF) for single accumulator
    dfun=function(rt,pars) dRDM(rt,pars),
    # Probability function (CDF) for single accumulator
    pfun=function(rt,pars) pRDM(rt,pars),
    # Race likelihood combining pfun and dfun
    log_likelihood=function(p_vector,dadm,min_ll=log(1e-10))
      log_likelihood_race_advantage(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
  )
}

#' Advantage Racing Diffusion Model with no C, pbeta mapping
#'
#' @return A list defining the cognitive model
#' @export
ARDMnoC <- function(){
  list(
    type="RACE",
    p_types=c("v" = log(1),"B" = log(1),"A" = log(0),"t0" = log(0),"s" = log(1),
              "AD"=log(1),"AS"=log(0),SD=0,SS=0,DD=0,DS=0),
    # Transform to natural scale
    Ntransform=function(x,use=NULL) {
      # transform parameters back to real line
      if (is.null(use)) {
        use <- !(colnames(x) %in% c("SD","SS","DD","DS"))
        x[,use] <- exp(x[,use])
      } else {
        use <- use[!(use %in% c("SD","SS","DD","DS"))] # never transform advantage shape
        x[,use] <- exp(x[,use])
      }
      x
    },
    # p_vector transform
    transform = function(x) x,
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      if (!is.null(attr(dadm,"adaptive"))) pars <- do_adaptive(pars,dadm)
      if (attr(dadm,"advantage")[[1]] == "parameter")
        pars <- advantage_parameters(pars,pars[,names(attr(dadm,"advantage"))[1]],
          length(levels(dadm$lR)),names(attr(dadm,"advantage"))[1]) else
        pars <- advantage_pars(pars,dadm[attr(dadm,"advantage")[[1]]],
          length(levels(dadm$lR)),names(attr(dadm,"advantage"))[1])
      attr(pars,"ok") <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0) &
        ((pars[,"v"] > 1e-3) | pars[,"v"] == 0) &  (pars[,"B"] >= 0) &  (pars[,"s"] > 0)
      pars
    },
    # Random function for racing accumulators
    rfun=function(lR=NULL,pars) {
      ok <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0) &
        ((pars[,"v"] > 1e-3) | pars[,"v"] == 0) &  (pars[,"B"] >= 0) &  (pars[,"s"] > 0)
      if (is.null(lR)) ok else rARDM(lR,pars,ok=ok)
    },
    # Density function (PDF) for single accumulator
    dfun=function(rt,pars) dRDM(rt,pars),
    # Probability function (CDF) for single accumulator
    pfun=function(rt,pars) pRDM(rt,pars),
    # Race likelihood combining pfun and dfun
    log_likelihood=function(p_vector,dadm,min_ll=log(1e-10))
      log_likelihood_race_advantage(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
  )
}



#### Missing ----

#' MRDM: RDM parameterization with missing values
#'
#' @return A list defining the cognitive model
#' @export
MRDMnoC <- function(){
  list(
    type="RACE",
    p_types=c("v" = log(1),"B" = log(1),"A" = log(0),"t0" = log(0),"s" = log(1),
              pContaminant=qnorm(0)),
    p_types=c("v","B","A","t0","s","pContaminant"),
    # Transform to natural scale
    Ntransform=function(x,use=NULL) {
      # transform parameters back to real line
      doprobit <- dimnames(x)[[2]] == "pContaminant"
      x[,doprobit] <- pnorm(x[,doprobit])
      x[,!doprobit] <- exp(x[,!doprobit])
      x
    },
    # p_vector transform
    transform = function(x) x,
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      attr(pars,"ok") <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0)
      pars
    },
    # Random function for racing accumulators
    rfun=function(lR=NULL,pars) {
      ok <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0)
      if (is.null(lR)) ok else rRDM(lR,pars,ok=ok)
    },
    # Density function (PDF) for single accumulator
    dfun=function(rt,pars) dRDM(rt,pars),
    # Probability function (CDF) for single accumulator
    pfun=function(rt,pars) pRDM(rt,pars),
    # Race likelihood combining pfun and dfun
    log_likelihood=function(p_vector,dadm,min_ll=log(1e-10))
      log_likelihood_race_missing(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
  )
}

#' MRDM: RDM parameterization with missing values
#'
#' @return A list defining the cognitive model
#' @export
MRDM <- function(){
  list(
    type="RACE",
    c_name = "MRDM",
    p_types=c("v" = log(1),"B" = log(1),"A" = log(0),"t0" = log(0),"s" = log(1),
              pContaminant=qnorm(0)),
    p_types=c("v","B","A","t0","s","pContaminant"),
    # Transform to natural scale
    Ntransform=function(x,use=NULL) {
      # transform parameters back to real line
      doprobit <- dimnames(x)[[2]] == "pContaminant"
      x[,doprobit] <- pnorm(x[,doprobit])
      x[,!doprobit] <- exp(x[,!doprobit])
      x
    },
    # p_vector transform
    transform = function(x) x,
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      attr(pars,"ok") <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0)
      pars
    },
    # Random function for racing accumulators
    rfun=function(lR=NULL,pars) {
      ok <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0)
      if (is.null(lR)) ok else rRDM(lR,pars,ok=ok)
    },
    # Density function (PDF) for single accumulator
    dfun=function(rt,pars) dRDM(rt,pars),
    # Probability function (CDF) for single accumulator
    pfun=function(rt,pars) pRDM(rt,pars),
    # Race likelihood combining pfun and dfun
    log_likelihood=function(p_vector,dadm,min_ll=log(1e-10))
      log_likelihood_race_missing(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
  )
}



#### Timed Racing Diffusion Model (TRDM) ----


rTRDM <- function(lR,pars,p_types=c("v","B","t0","s","A",
                                    "vT","BT","t0T","sT","AT","pGuess"),
                  ok=rep(TRUE,dim(pars)[1]))
  # lR is an empty latent response factor lR with one level for each accumulator.
  # pars is a matrix of corresponding parameter values named as in p_types
  # pars must be sorted so accumulators and parameter for each trial are in
  # contiguous rows. "s" parameter will be used but can be omitted
  #
{
  if (!all(p_types %in% dimnames(pars)[[2]]))
    stop("pars must have columns ",paste(p_types,collapse = " "))

  if (any(dimnames(pars)[[2]]=="s_E")) # rescale
    pars[ok,c("A","B","v")] <- pars[ok,c("A","B","v")]/pars[ok,"s"]
  if (any(dimnames(pars)[[2]]=="sT")) # rescale
    pars[ok,c("AT","BT","vT")] <- pars[ok,c("AT","BT","vT")]/pars[ok,"sT"]

  pars[,"B"][pars[,"B"]<0] <- 0 # Protection for negatives
  pars[,"BT"][pars[,"BT"]<0] <- 0 # Protection for negatives

  # make containers for all trials
  nr <- length(levels(lR))
  bad <- rep(NA, length(lR)/nr)
  out <- data.frame(R = bad, rt = bad)
  dt <- matrix(Inf, nrow=nr,ncol=nrow(pars)/nr)
  tt <- rep(Inf,nrow(pars)/nr)
  t0 <- pars[,"t0"]
  pGuess <- matrix(pars[,"pGuess"],nrow=nr)
  nok <- sum(ok)
  pars <- pars[ok,]

  # first, sample responses from ok evidence accumulation process
  dt[ok] <- rWald(nok,B=pars[,"B"],v=pars[,"v"],A=pars[,"A"])
  R <- apply(dt,2,which.min)
  pick <- cbind(R,1:dim(dt)[2]) # Matrix to pick winner
  # Any t0 difference with lR due to response production time (no effect on race)
  rt <- matrix(t0,nrow=nr)[pick] + dt[pick]
  out$R <- levels(lR)[R]
  out$R <- factor(out$R,levels=levels(lR))
  out$rt <- rt

  # second, sample timing RTs then resample responses if timer terminated first
  is1 <- lR[ok] == levels(lR)[1]
  nok <- nok/nr
  tt[ok[is1]] <- pars[is1,"t0T"] +
                 rWald(nok, B=pars[is1,"BT"],v=pars[is1,"vT"],A=pars[is1,"AT"])
  timer_wins <- tt < out$rt

  # resample response only if there were >0 timer responses
  if(any(timer_wins)) {
    out$rt[timer_wins] <- tt[timer_wins]
    # then determine which response was given
    if (all(pGuess==-Inf)) {
      out$R[timer_wins] <- levels(lR)[sample(1:nr,sum(timer_wins),replace=TRUE,prob=rep(1/nr,nr))]
    } else {
      pGuess <- apply(pGuess[,timer_wins,drop=FALSE],2,function(x) {
        if (x[1]==-Inf) {
          x[-1] <- pnorm(x[-1])/(nr-1)
          x[1] <- 1-sum(x[-1])
        } else if (sum(x)!=1) stop("pGuess must sum to 1 over accumulators")
        x
      })
      out$R[timer_wins] <- levels(lR)[apply(pGuess,2,function(x){sample(1:nr,1,prob=x)})]
    }
    attr(out,"timer_wins") <- timer_wins
  }
  out
}

#' Timed Racing Diffusion Model with no C
#'
#' @return A list defining the cognitive model
#' @export
TRDMnoC <- function(){
  list(
    type="RACE",
    p_types=c("v" = log(1),"B" = log(1),A = log(0),"t0" = log(0),"s" = log(1),
      "vT" = log(1), "BT" = log(1), AT=log(0), "t0T" = log(.05), "sT" = log(1),
      pGuess = -Inf),
    # Transform to natural scale
    Ntransform=function(x,use=NULL) {
      # transform parameters back to real line
      if (is.null(use)) {
        x[,dimnames(x)[[2]]!="pGuess"] <- exp(x[,dimnames(x)[[2]]!="pGuess"])
      } else {
        use <- use[use!="pGuess"]
        x[,use] <- exp(x[,use])
      }
      x
    },
    # p_vector transform
    transform = function(x) x,
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      attr(pars,"ok") <- (pars[,"t0"] > .05) & (pars[,"s"] > 0) &
                        ((pars[,"v"] > 1e-3) | pars[,"v"] == 0) &
                        ((pars[,"B"] > 1e-3) | pars[,"B"] == 0) &
                        ((pars[,"A"] > 1e-6) | pars[,"A"] == 0) &
                         (pars[,"t0T"] > 0) & (pars[,"sT"] > 0) &
                        ((pars[,"vT"] > 1e-3) | pars[,"vT"] == 0) &
                        ((pars[,"BT"] > 1e-3) | pars[,"BT"] == 0) &
                        ((pars[,"AT"] > 1e-6) | pars[,"AT"] == 0)
      pars
    },
    # Random function for racing accumulators
    rfun=function(lR=NULL,pars) {
      ok <- (pars[,"t0"] > .05) & (pars[,"s"] > 0) &
           ((pars[,"v"] > 1e-3) | pars[,"v"] == 0) &
           ((pars[,"B"] > 1e-3) | pars[,"B"] == 0) &
           ((pars[,"A"] > 1e-6) | pars[,"A"] == 0) &
            (pars[,"t0T"] > 0) & (pars[,"sT"] > 0) &
           ((pars[,"vT"] > 1e-3) | pars[,"vT"] == 0) &
           ((pars[,"BT"] > 1e-3) | pars[,"BT"] == 0) &
           ((pars[,"AT"] > 1e-6) | pars[,"AT"] == 0)
      if (is.null(lR)) ok else rTRDM(lR,pars,ok=ok)
    },
    # Density function (PDF) for single evidence accumulator
    dfun=function(rt,pars) dRDM(rt,pars),
    # Probability function (CDF) for single evidence accumulator
    pfun=function(rt,pars) pRDM(rt,pars),
    # Density function (PDF) for single timing accumulator
    dfunT=function(rt,pars) {
       dRDM(rt,matrix(pars[,c("vT","BT","t0T","sT","AT")],ncol=5,
                      dimnames=list(NULL,c("v","B","t0","s","A"))))
    },
    # Probability function (CDF) for single timing accumulator
    pfunT=function(rt,pars) {
      pRDM(rt,matrix(pars[,c("vT","BT","t0T","sT","AT")],ncol=5,
                      dimnames=list(NULL,c("v","B","t0","s","A"))))
    },
    # Race likelihood combining pfun and dfun
    log_likelihood=function(p_vector,dadm,min_ll=log(1e-10))
      log_likelihood_race_trdm(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
  )
}

#' Timed Racing Diffusion Model
#'
#' @return A list defining the cognitive model
#' @export
TRDM <- function(){
  list(
    type="RACE",
    c_name = "TRDM",
    p_types=c("v" = log(1),"B" = log(1),A = log(0),"t0" = log(0),"s" = log(1),
      "vT" = log(1), "BT" = log(1), AT=log(0), "t0T" = log(.05), "sT" = log(1),
      pGuess = -Inf),
    # Transform to natural scale
    Ntransform=function(x,use=NULL) {
      # transform parameters back to real line
      if (is.null(use)) {
        x[,dimnames(x)[[2]]!="pGuess"] <- exp(x[,dimnames(x)[[2]]!="pGuess"])
      } else {
        use <- use[use!="pGuess"]
        x[,use] <- exp(x[,use])
      }
      x
    },
    # p_vector transform
    transform = function(x) x,
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      attr(pars,"ok") <- (pars[,"t0"] > .05) & (pars[,"s"] > 0) &
                        ((pars[,"v"] > 1e-3) | pars[,"v"] == 0) &
                        ((pars[,"B"] > 1e-3) | pars[,"B"] == 0) &
                        ((pars[,"A"] > 1e-6) | pars[,"A"] == 0) &
                         (pars[,"t0T"] > 0) & (pars[,"sT"] > 0) &
                        ((pars[,"vT"] > 1e-3) | pars[,"vT"] == 0) &
                        ((pars[,"BT"] > 1e-3) | pars[,"BT"] == 0) &
                        ((pars[,"AT"] > 1e-6) | pars[,"AT"] == 0)
      pars
    },
    # Random function for racing accumulators
    rfun=function(lR=NULL,pars) {
      ok <- (pars[,"t0"] > .05) & (pars[,"s"] > 0) &
           ((pars[,"v"] > 1e-3) | pars[,"v"] == 0) &
           ((pars[,"B"] > 1e-3) | pars[,"B"] == 0) &
           ((pars[,"A"] > 1e-6) | pars[,"A"] == 0) &
            (pars[,"t0T"] > 0) & (pars[,"sT"] > 0) &
           ((pars[,"vT"] > 1e-3) | pars[,"vT"] == 0) &
           ((pars[,"BT"] > 1e-3) | pars[,"BT"] == 0) &
           ((pars[,"AT"] > 1e-6) | pars[,"AT"] == 0)
      if (is.null(lR)) ok else rTRDM(lR,pars,ok=ok)
    },
    # Density function (PDF) for single evidence accumulator
    dfun=function(rt,pars) dRDM(rt,pars),
    # Probability function (CDF) for single evidence accumulator
    pfun=function(rt,pars) pRDM(rt,pars),
    # Density function (PDF) for single timing accumulator
    dfunT=function(rt,pars) {
       dRDM(rt,matrix(pars[,c("vT","BT","t0T","sT","AT")],ncol=5,
                      dimnames=list(NULL,c("v","B","t0","s","A"))))
    },
    # Probability function (CDF) for single timing accumulator
    pfunT=function(rt,pars) {
      pRDM(rt,matrix(pars[,c("vT","BT","t0T","sT","AT")],ncol=5,
                      dimnames=list(NULL,c("v","B","t0","s","A"))))
    },
    # Race likelihood combining pfun and dfun
    log_likelihood=function(p_vector,dadm,min_ll=log(1e-10))
      log_likelihood_race_trdm(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
  )
}

