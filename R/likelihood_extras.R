#### Advantage ----
log_likelihood_race_advantage <- function(p_vector,dadm,min_ll=log(1e-10))
  # Race model summed log likelihood
{

  pars <- get_pars_matrix(p_vector,dadm)

  nr <- length(levels(dadm$R))
  na <- nr-1
  nt <- nrow(dadm)/nr

  if (is.null(attr(pars,"ok")))
    ok <- !logical(dim(pars)[1]/(nr*na)) else
    ok <- apply(matrix(attr(pars,"ok"),nrow=(nr*na)),2,all)

  # This seems wasteful, can pre-compute, perhaps pass as dadm attr?
  winner <- rep(dadm$winner,each=na)
  rt <- rep(dadm$rt,each=na)

  # log pdf of winning response accumulators
  pdf <- matrix(attr(dadm,"model")()$dfun(rt=rt[winner],pars=pars[winner,]),
                  nrow=na,ncol=nt)
  # cdf (NOT log) of winning response accumulators
  if (na>1) cdfw <- matrix(attr(dadm,"model")()$pfun(rt=rt[winner],pars=pars[winner,]),
                            nrow=na,ncol=nt)
  # cdf (NOT log) of loosing response accumulators
  cdfl <- array(attr(dadm,"model")()$pfun(rt=rt[!winner],pars=pars[!winner,]),
               dim=c(na,na,nt))
  ia <- 1:na
  if (na==1) ll <- as.vector(log(pdf)) + log(1-as.vector(cdfl)) else {
    ll <- 0
    for (i in ia) {    # sum over accumulators for winning response
      ifirst <- ia[-i] # other accumulators, already done
                 # finish at t x already done
      ll <- ll + pdf[i,]*apply(cdfw[ifirst,,drop=FALSE],2,sum)
    }
    ll <- log(ll) + apply(log(1-apply(cdfl, 2:3, prod)),2,sum) # other responses survivors at t
  }
  ll[is.na(ll) | !ok] <- min_ll
  return(sum(pmax(min_ll,ll[attr(dadm,"expand_winner")])))
}


#### TRDM ----

log_likelihood_race_trdm <- function(p_vector,dadm,min_ll=log(1e-10))
  # TRDM model summed log likelihood
{
  pars <- get_pars_matrix(p_vector,dadm)

  if (is.null(attr(pars,"ok")))
    ok <- !logical(dim(pars)[1]) else ok <- attr(pars,"ok")

  nr <- length(levels(dadm$lR)) # number of *evidence* accumulators

  ll <- numeric(nrow(pars))
  pars <- pars[ok,]
  winner <- dadm$winner[ok]
  rt <- dadm$rt[ok]
  # unbiased guess if all -Inf, otherwise if first accumulator NA use lR=2,lR=3...
  # to set guess probability, else should already be a set of probabilities per
  # accumulator (in which case the values are not being sampled).
  if (all(pars[,"pGuess"]==-Inf)) pars[,"pGuess"] <- rep(1/nr,nrow(pars)) else {
    pars[,"pGuess"] <- as.vector(apply(matrix(pars[,"pGuess",drop=FALSE],nrow=nr),2,function(x) {
      x[-1] <- pnorm(x[-1])/(nr-1)
      x[1] <- 1-sum(x[-1])
      x
    }))
  }

  # pdf of timer (only take one of nr rows as only one timing accumulator)
  pdfTW <- pars[winner,"pGuess"] * attr(dadm,"model")()$dfunT(rt=rt[winner],pars=pars[winner,])
  # pdf of evidence winner
  pdfEW <- attr(dadm,"model")()$dfun(rt=rt[winner],pars=pars[winner,])
  # survivors of timing accumulator in case it looses
  surviveT <- 1-attr(dadm,"model")()$pfunT(rt=rt[winner],pars=pars[winner,])
  # survivors of all evidence accumulators (need winner survivor in case timer wins)
  surviveE <-  1-attr(dadm,"model")()$pfun(rt=rt,pars=pars)

  ll[ok] <- log(pdfEW*surviveT + pdfTW*surviveE[winner])
  if (nr>1) {
    if(nr==2) {
      ll[ok] <- ll[ok] + log(surviveE[!winner])
    } else {
      ll[ok] <- ll[ok] + log(apply(matrix(surviveE[!winner],nrow=nr-1),2,prod))
    }
  }
  ll[is.nan(ll) | is.na(ll) | !ok] <- min_ll
  return(sum(pmax(min_ll,ll[attr(dadm,"expand_winner")])))
}


#### Missing ----

my_integrate <- function(...,upper=Inf,big=10)
  # Avoids bug in integrate upper=Inf that uses only 1  subdivision
  # Use of  big=10 is arbitrary ...
{
  out <- try(integrate(...,upper=upper),silent=TRUE)
  if (!is(out,"try-error") && upper==Inf && out$subdivisions==1)
    out <- try(integrate(...,upper=big),silent=TRUE)
  out
}


# pars <- cbind(a=c(1,1),v=c(1,1),t0=c(.3,.3),z=c(.5,.5),d=c(0,0),sz=c(0,0),sv=c(0,0),st0=c(0,0),s=c(1,1))
# LC <- .4
# UC=2
#
# EMC2:::pDDM(c(LC,LC),"lower",pars)
# EMC2:::pDDM(UC,"lower",pars)
#
# EMC2:::pDDM(LC,"upper",pars)
# EMC2:::pDDM(UC,"upper",pars)


log_likelihood_ddm_missing <- function(p_vector,dadm,min_ll=log(1e-10))
  # DDM summed log likelihood, with protection against numerical issues
{

  pr_pt <- function(LT,UT,R,p)
  # p(untruncated response)/p(truncated response), > 1, multiplicative truncation correction
  {
    pr <- attr(dadm,"model")()$pfun(rep(Inf,length(R)),R,p)
    out <- rep(1,length(R))
    if (!any(is.na(pr))) {
      ok <- pr>0
      if (any(ok)) {
        pt <- attr(dadm,"model")()$pfun(UT[ok],R[ok],p[ok,,drop=FALSE]) -
            attr(dadm,"model")()$pfun(LT[ok],R[ok],p[ok,,drop=FALSE])
        if (!any(is.na(pt))) out[ok] <- pr[ok]/pt
      }
    }
    out[is.na(out) | is.nan(out) | !is.finite(out) | out < 1 ] <- 1
    out
  }

  pars <- get_pars_matrix(p_vector,dadm)
  like <- numeric(dim(dadm)[1])
  if (any(attr(pars,"ok"))) {
    rt <- dadm$rt[attr(pars,"ok")]
    R <- dadm$R[attr(pars,"ok")]
    p <- pars[attr(pars,"ok"),,drop=FALSE]

    # Calculate truncation?
    LT <- attr(dadm,"LT")
    UT <- attr(dadm,"UT")
    dotrunc <- (!is.null(LT) | !is.null(UT))
    if (is.null(LT)) LT <- 0
    if (is.null(UT)) UT <- Inf

    # Calculate censoring
    LC <- attr(dadm,"LC")
    UC <- attr(dadm,"UC")

    likeok <- rep(NA,sum(attr(pars,"ok")))

    # Response known
    # Fast
    nort <- rt==-Inf; nort[is.na(nort)] <- FALSE; nort <- nort & !is.na(R)
    if ( any(nort) ) {
      likeok[nort] <- pmax(0,attr(dadm,"model")()$pfun(rep(LC,sum(nort)),R[nort],p[nort,,drop=FALSE]))
    }
    # Slow
    nort <- rt==Inf; nort[is.na(nort)] <- FALSE; nort <- nort & !is.na(R)
    if ( any(nort) ) {
      likeok[nort] <- pmax(0,attr(dadm,"model")()$pfun(rep(Inf,sum(nort)),R[nort],p[nort,,drop=FALSE])-
                      attr(dadm,"model")()$pfun(rep(UC,sum(nort)),R[nort],p[nort,,drop=FALSE]))
    }
    # No direction
    nort <- is.na(rt) & !is.na(R)
    if ( any(nort) ) {
      likeok[nort] <- pmax(0,attr(dadm,"model")()$pfun(rep(LC,sum(nort)),R[nort],p[nort,,drop=FALSE]) +
                      (attr(dadm,"model")()$pfun(rep(Inf,sum(nort)),R[nort],p[nort,,drop=FALSE])-
                         attr(dadm,"model")()$pfun(rep(UC,sum(nort)),R[nort],p[nort,,drop=FALSE])))
    }

    # Response unknown.
    # Fast
    nort <- rt==-Inf; nort[is.na(nort)] <- FALSE; nort <- nort & is.na(R)
    if ( any(nort) ) {
      likeok[nort] <- pmax(0,attr(dadm,"model")()$pfun(rep(LC,sum(nort)),"lower",p[nort,,drop=FALSE]) +
                      attr(dadm,"model")()$pfun(rep(LC,sum(nort)),"upper",p[nort,,drop=FALSE]))
    }
    # Slow
    nort <- rt==Inf; nort[is.na(nort)] <- FALSE; nort <- nort & is.na(R)
    if ( any(nort) ) {
      likeok[nort] <- pmax(0,(attr(dadm,"model")()$pfun(rep(Inf,sum(nort)),"lower",p[nort,,drop=FALSE])-
                       attr(dadm,"model")()$pfun(rep(UC,sum(nort)),"lower",p[nort,,drop=FALSE])) +
                      (attr(dadm,"model")()$pfun(rep(Inf,sum(nort)),"upper",p[nort,,drop=FALSE])-
                       attr(dadm,"model")()$pfun(rep(UC,sum(nort)),"upper",p[nort,,drop=FALSE])))
    }
    # no direction
    nort <- is.na(rt) & is.na(R)
    likeok[nort] <- 0
    nort <- nort & (p[,"pContaminant"] == 0) # Otherwise not identifiable
    if ( any(nort) ) {
      likeok[nort] <- pmax(0,
        attr(dadm,"model")()$pfun(rep(LC,sum(nort)),"lower",p[nort,,drop=FALSE]) +
        attr(dadm,"model")()$pfun(rep(LC,sum(nort)),"upper",p[nort,,drop=FALSE]) +
        (attr(dadm,"model")()$pfun(rep(Inf,sum(nort)),"lower",p[nort,,drop=FALSE])-
         attr(dadm,"model")()$pfun(rep(UC,sum(nort)),"lower",p[nort,,drop=FALSE])) +
        (attr(dadm,"model")()$pfun(rep(Inf,sum(nort)),"upper",p[nort,,drop=FALSE])-
         attr(dadm,"model")()$pfun(rep(UC,sum(nort)),"upper",p[nort,,drop=FALSE])))
    }

    # Truncation where not censored or censored and response known
    ok <- is.na(likeok) & !is.na(R)
    mult <- rep(1,length(likeok))
    if ( dotrunc & any(ok) )
      mult[ok] <- pr_pt(rep(LT,sum(ok)),rep(UT,sum(ok)),R[ok],p[ok,,drop=FALSE])

    # Usual non-missing update x truncation ratio
    ok <- is.na(likeok)
    likeok[ok] <- mult[ok]*attr(dadm,"model")()$dfun(rt[ok],R[ok],p[ok,,drop=FALSE])

    # Non-process (contaminant) miss.
    ispContaminant <- p[,"pContaminant"]>0
    if ( any(ispContaminant) ) {
      pc <- p[,"pContaminant"]
      isMiss <- is.na(R)
      likeok[isMiss] <- pc[isMiss] + (1-pc[isMiss])*likeok[isMiss]
      likeok[!isMiss] <- (1-pc[!isMiss])*likeok[!isMiss]
    }

    like[attr(pars,"ok")] <- likeok
  }

  like[attr(pars,"ok")][is.na(like[attr(pars,"ok")])] <- 0
  sum(pmax(min_ll,log(like[attr(dadm,"expand")])))
}


f <- function(t,p,dfun,pfun) {
    # Called by integrate to get race density for vector of times t given
    # matrix of parameters where first row is the winner.

    out <- dfun(t,
                matrix(rep(p[1,],each=length(t)),nrow=length(t),dimnames=list(NULL,dimnames(p)[[2]])))
    if (dim(p)[1]>1) for (i in 2:dim(p)[1])
      out <- out*(1-pfun(t,
                         matrix(rep(p[i,],each=length(t)),nrow=length(t),dimnames=list(NULL,dimnames(p)[[2]]))))
    out
}

pr_pt <- function(LT,UT,ps,dadm)
    # log(p(untruncated response)/p(truncated response)), >= log(1), multiplicative truncation correction
  {
    pr <- my_integrate(f,lower=0,upper=Inf,p=ps,
                       dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun)
    if (inherits(pr, "try-error") || suppressWarnings(is.nan(pr$value))) return(NA)
    if (pr$value==0) return(0)
    pt <- my_integrate(f,lower=LT,upper=UT,p=ps,
                       dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun)
    if (inherits(pt, "try-error") || suppressWarnings(is.nan(pt$value)) || pt$value==0) return(NA)
    out <- pmax(0,pmin(pr$value,1))/pmax(0,pmin(pt$value,1))
    if (is.infinite(out)) return(NA)
    out
}

pLU <- function(LT,LC,UC,UT,ps,dadm)
    # Probability from LT-LC + UC-UT
  {
    pL <- my_integrate(f,lower=LT,upper=LC,p=ps,
                       dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun)
    if (inherits(pL,"try-error") || suppressWarnings(is.nan(pL$value))) return(NA)
    pU <- my_integrate(f,lower=UC,upper=UT,p=ps,
                       dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun)
    if (inherits(pU,"try-error") || suppressWarnings(is.nan(pU$value))) return(NA)
    pmax(0,pmin(pL$value,1))+pmax(0,pmin(pU$value,1))
}


log_likelihood_race_missing <- function(p_vector,dadm,min_ll=log(1e-10))
  # Race model summed log likelihood for models allowing missing values
{

  pars <- get_pars_matrix(p_vector,dadm)

  if (any(names(dadm)=="RACE")) # Some accumulators not present
    pars[as.numeric(dadm$lR)>as.numeric(as.character(dadm$RACE)),] <- NA

  if (is.null(attr(pars,"ok")))
    ok <- !logical(dim(pars)[1]) else ok <- attr(pars,"ok")

  lds <- numeric(dim(dadm)[1]) # log pdf (winner) or survivor (losers)
  lds[dadm$winner] <- log(attr(dadm,"model")()$dfun(rt=dadm$rt[dadm$winner],
                                                    pars=pars[dadm$winner,]))
  n_acc <- length(levels(dadm$R))
  if (n_acc>1) lds[!dadm$winner] <-
    log(1-attr(dadm,"model")()$pfun(rt=dadm$rt[!dadm$winner],pars=pars[!dadm$winner,]))
  lds[is.na(lds) | !ok] <- -Inf

  # Calculate truncation?
  LT <- attr(dadm,"LT")
  UT <- attr(dadm,"UT")
  dotrunc <- (!is.null(LT) | !is.null(UT))
  if (is.null(LT)) LT <- 0
  if (is.null(UT)) UT <- Inf

  # Calculate censoring
  LC <- attr(dadm,"LC")
  UC <- attr(dadm,"UC")

  # Response known
  # Fast
  nort <- dadm$rt==-Inf; nort[is.na(nort)] <- FALSE; nort <- nort & !is.na(dadm$R)
  if ( any(nort) ) {
    mpars <- array(pars[nort,,drop=FALSE],dim=c(n_acc,sum(nort)/n_acc,ncol(pars)),
                   dimnames = list(NULL,NULL,colnames(pars)))
    winner <- matrix(dadm$winner[nort],nrow=n_acc)
    tofix <- dadm$winner & nort
    for (i in 1:dim(mpars)[2]) {
      if (dim(mpars)[[1]]==1) pi <- t(as.matrix(mpars[,i,])) else
        pi <- mpars[,i,][order(!winner[,i]),]
      tmp <- my_integrate(f,lower=LT,upper=LC,p=pi,
                          dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun)
      if ( !inherits(tmp, "try-error") && suppressWarnings(!is.nan(log(tmp$value))) )
        lds[tofix][i] <- log(pmax(0,pmin(tmp$value,1)))
    }
  }
  # Slow
  nort <- dadm$rt==Inf; nort[is.na(nort)] <- FALSE; nort <- nort & !is.na(dadm$R)
  if ( any(nort) ) {
    mpars <- array(pars[nort,,drop=FALSE],dim=c(n_acc,sum(nort)/n_acc,ncol(pars)),
                   dimnames = list(NULL,NULL,colnames(pars)))
    winner <- matrix(dadm$winner[nort],nrow=n_acc)
    tofix <- dadm$winner & nort
    for (i in 1:dim(mpars)[2]) {
      if (dim(mpars)[[1]]==1) pi <- t(as.matrix(mpars[,i,])) else
        pi <- mpars[,i,][order(!winner[,i]),]
      tmp <- my_integrate(f,lower=UC,upper=UT,p=pi,
                          dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun)
      if (!inherits(tmp, "try-error") && suppressWarnings(!is.nan(log(tmp$value))))
        lds[tofix][i] <- log(pmax(0,pmin(tmp$value,1)))
    }
  }
  # No direction
  nort <- is.na(dadm$rt) & !is.na(dadm$R)
  if ( any(nort) ) {
    mpars <- array(pars[nort,,drop=FALSE],dim=c(n_acc,sum(nort)/n_acc,ncol(pars)),
                   dimnames = list(NULL,NULL,colnames(pars)))
    winner <- matrix(dadm$winner[nort],nrow=n_acc)
    tofix <- dadm$winner & nort
    for (i in 1:dim(mpars)[2]) {
      if (dim(mpars)[[1]]==1) pi <- t(as.matrix(mpars[,i,])) else
        pi <- mpars[,i,][order(!winner[,i]),]
      tmp <- my_integrate(f,lower=LT,upper=LC,p=pi,
                          dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun)
      if ( !inherits(tmp,"try-error") && suppressWarnings(!is.nan(tmp$value)) ) {
        p <- tmp$value
        if (dim(mpars)[[1]]==1) pi <- mpars[,i,,drop=FALSE] else
          pi <- mpars[,i,][order(!winner[,i]),]
        tmp <- my_integrate(f,lower=UC,upper=UT,p=pi,
                            dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun)
        if ( !inherits(tmp,"try-error") && suppressWarnings(!is.nan(tmp$value)) )
          p <- p + tmp$value else p <- 0
      } else p <- 0
      lds[tofix][i] <- log(pmax(0,pmin(p,1)))
    }
  }

  # Response unknown.
  # Fast
  nort <- dadm$rt==-Inf; nort[is.na(nort)] <- FALSE; nort <- nort & is.na(dadm$R)
  if ( any(nort) ) {
    mpars <- array(pars[nort,,drop=FALSE],dim=c(n_acc,sum(nort)/n_acc,ncol(pars)),
                   dimnames = list(NULL,NULL,colnames(pars)))
    winner <- matrix(dadm$winner[nort],nrow=n_acc)
    tofixfast <- dadm$winner & nort
    for (i in 1:dim(mpars)[2]) {
      if (dim(mpars)[[1]]==1) pi <- t(as.matrix(mpars[,i,])) else
        pi <- mpars[,i,]
      pc <- my_integrate(f,lower=LT,upper=LC,p=pi,
                         dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun)
      if (inherits(pc, "try-error") || suppressWarnings(is.nan(pc$value)))
        p <- NA else p <- pmax(0,pmin(pc$value,1))
      if (!is.na(p)) {
        if (p != 0 && !(LT==0 & UT==Inf))  cf <- pr_pt(LT,UT,mpars[,i,],dadm) else cf <- 1
        if (!is.na(cf)) p <- p*cf
      }
      if (!is.na(p) & n_acc>1) for (j in 2:n_acc) {
        pc <- my_integrate(f,lower=LT,upper=LC,p=mpars[,i,][c(j,c(1:n_acc)[-j]),],
                           dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun)
        if (inherits(pc, "try-error") || suppressWarnings(is.nan(pc$value))) {
          p <- NA; break
        }
        if (pc$value != 0 & !(LT==0 & UT==Inf))
          cf <- pr_pt(LT,UT,mpars[,i,][c(j,c(1:n_acc)[-j]),],dadm) else cf <- 1
          if (!is.na(cf)) p <- p + pc$value*cf
      }
      lp <- log(p)
      if (!is.nan(lp) & !is.na(lp)) lds[tofixfast][i] <- lp else lds[tofixfast][i] <- -Inf
    }
  } else tofixfast <- NA
  # Slow
  nort <- dadm$rt==Inf; nort[is.na(nort)] <- FALSE; nort <- nort & is.na(dadm$R)
  if ( any(nort) ) {
    mpars <- array(pars[nort,,drop=FALSE],dim=c(n_acc,sum(nort)/n_acc,ncol(pars)),
                   dimnames = list(NULL,NULL,colnames(pars)))
    winner <- matrix(dadm$winner[nort],nrow=n_acc)
    tofixslow <- dadm$winner & nort
    for (i in 1:dim(mpars)[2]) {
      if (dim(mpars)[[1]]==1) pi <- t(as.matrix(mpars[,i,])) else
        pi <- mpars[,i,]
      pc <- my_integrate(f,lower=UC,upper=UT,p=pi,
                         dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun)
      if (inherits(pc, "try-error") || suppressWarnings(is.nan(pc$value)))
        p <- NA else p <- pmax(0,pmin(pc$value,1))
      if (!is.na(p)) {
        if (p != 0 && !(LT==0 & UT==Inf))  cf <- pr_pt(LT,UT,mpars[,i,],dadm) else cf <- 1
        if (!is.na(cf)) p <- p*cf
      }
      if (!is.na(p) & n_acc>1) for (j in 2:n_acc) {
        pc <- my_integrate(f,lower=UC,upper=UT,p=mpars[,i,][c(j,c(1:n_acc)[-j]),],
                           dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun)
        if (inherits(pc, "try-error") || suppressWarnings(is.nan(pc$value))) {
          p <- NA; break
        }
        if (pc$value != 0 & !(LT==0 & UT==Inf))
          cf <- pr_pt(LT,UT,mpars[,i,][c(j,c(1:n_acc)[-j]),],dadm) else cf <- 1
          if (!is.na(cf)) p <- p + pc$value*cf
      }
      lp <- log(p)
      if (!is.nan(lp) & !is.na(lp)) lds[tofixslow][i] <- lp else lds[tofixslow][i] <- -Inf
    }
  } else tofixslow <- NA
  # no direction
  nort <- is.na(dadm$rt) & is.na(dadm$R)
  nort <- nort & (pars[,"pContaminant"] == 0) # Otherwise not identifiable
  if ( any(nort) ) {
    mpars <- array(pars[nort,,drop=FALSE],dim=c(n_acc,sum(nort)/n_acc,ncol(pars)),
                   dimnames = list(NULL,NULL,colnames(pars)))
    winner <- matrix(dadm$winner[nort],nrow=n_acc)
    tofix <- dadm$winner & nort
    for (i in 1:dim(mpars)[2]) {
      if (dim(mpars)[[1]]==1) pi <- t(as.matrix(mpars[,i,])) else
        pi <- mpars[,i,]
      pc <- pLU(LT,LC,UC,UT,pi,dadm)
      if (is.na(pc)) p <- NA else {
        if (pc!=0 & !(LT==0 & UT==Inf)) cf <- pr_pt(LT,UT,pi,dadm) else cf <- 1
        if (!is.na(cf)) p <- pc*cf else p <- NA
        if (!is.na(p) & n_acc>1) for (j in 2:n_acc) {
          pc <- pLU(LT,LC,UC,UT,mpars[,i,][c(j,c(1:n_acc)[-j]),],dadm)
          if (is.na(pc)) {
            p <- NA; break
          }
          if (pc!=0 & !(LT==0 & UT==Inf))
            cf <- pr_pt(LT,UT,mpars[,i,][c(j,c(1:n_acc)[-j]),],dadm) else cf <- 1
            if (is.na(cf)) {
              p <- NA; break
            } else p <- p +  pc*cf
        }
      }
      lp <- log(pmax(0,pmin(p,1)))
      if (!is.nan(lp) & !is.na(lp)) lds[tofix][i] <- lp
    }
  }


  # Truncation where not censored or censored and response known
  ok <- is.finite(lds[attr(dadm,"unique_nort") & dadm$winner])
  alreadyfixed <- is.na(dadm$R[attr(dadm,"unique_nort") & dadm$winner])
  ok <- ok & !alreadyfixed
  if ( dotrunc & any(ok) ) {
    tpars <- pars[attr(dadm,"unique_nort"),,drop=FALSE]
    tpars <- array(tpars[,,drop=FALSE],dim=c(n_acc,nrow(tpars)/n_acc,ncol(tpars)),
                   dimnames = list(NULL,NULL,colnames(tpars)))[,,,drop=FALSE]
    winner <- matrix(dadm$winner[attr(dadm,"unique_nort")],nrow=n_acc)[,,drop=FALSE]
    cf <- rep(NA,length(ok))
    for (i in 1:length(ok)) if (ok[i]) {
      if (dim(tpars)[[1]]==1) pi <- t(as.matrix(tpars[,i,])) else
        pi <- tpars[,i,][order(!winner[,i]),]
      cf[i] <- pr_pt(LT,UT,pi,dadm)
    }
    cf <- rep(log(cf),each=n_acc)[attr(dadm,"expand_nort")]
    fix <- dadm$winner & !is.na(cf) & !is.nan(cf) & is.finite(cf)
    if (any(fix)) lds[fix] <- lds[fix] + cf[fix]
    badfix <- dadm$winner & (is.na(cf) | is.nan(cf) | is.infinite(cf))
    if (!all(is.na(tofixfast))) badfix <- badfix & !tofixfast
    if (!all(is.na(tofixslow))) badfix <- badfix & !tofixslow
    if (any(badfix)) lds[badfix] <- -Inf
  }

  # Non-process (contaminant) miss.
  ispContaminant <- pars[,"pContaminant"]>0
  if ( any(ispContaminant) ) {
    p <- exp(lds[dadm$winner])
    pc <- pars[dadm$winner,"pContaminant"]
    isMiss <- is.na(dadm$R[dadm$winner])
    p[isMiss] <- pc[isMiss] + (1-pc[isMiss])*p[isMiss]
    p[!isMiss] <- (1-pc[!isMiss])*p[!isMiss]
    lds[dadm$winner] <- log(p)
  }

  if (n_acc>1) {
    ll <- lds[dadm$winner]
    if (n_acc==2) {
      ll <- ll + lds[!dadm$winner]
    } else {
      ll <- ll + apply(matrix(lds[!dadm$winner],nrow=n_acc-1),2,sum)
    }
  } else ll <- lds
  ll[is.na(ll) | is.nan(ll)] <- -Inf

  return(sum(pmax(min_ll,ll[attr(dadm,"expand_winner")])))
}


log_likelihood_race_missing_LBAU <- function(p_vector,dadm,min_ll=log(1e-10))
  # Race model summed log likelihood for an LBA allowing negative rates and
  # missing values
{

  f <- function(t,p,dfun,pfun) {
    # Called by integrate to get race density for vector of times t given
    # matrix of parameters where first row is the winner.

    out <- dfun(t,
                matrix(rep(p[1,],each=length(t)),nrow=length(t),dimnames=list(NULL,dimnames(p)[[2]])))
    if (dim(p)[1]>1) for (i in 2:dim(p)[1])
      out <- out*(1-pfun(t,
                         matrix(rep(p[i,],each=length(t)),nrow=length(t),dimnames=list(NULL,dimnames(p)[[2]]))))
    out
  }

  pr_pt <- function(LT,UT,ps,dadm)
    # p(untrucated response)/p(truncated response), > 1, multiplicative truncation correction
  {
    pr <- my_integrate(f,lower=0,upper=Inf,p=ps,
                       dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun)
    if (inherits(pr, "try-error") || suppressWarnings(is.nan(pr$value))) return(NA)
    if (pr$value==0) return(0)
    pt <- my_integrate(f,lower=LT,upper=UT,p=ps,
                       dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun)
    if (inherits(pt, "try-error") || suppressWarnings(is.nan(pt$value)) || pt$value==0) return(NA)
    out <- pmax(0,pmin(pr$value,1))/pmax(0,pmin(pt$value,1))
    if (is.infinite(out)) return(NA)
    out
  }

  pLU <- function(LT,LC,UC,UT,ps,dadm)
    # Probability from LT-LC + UC-UT
  {
    pL <- my_integrate(f,lower=LT,upper=LC,p=ps,
                       dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun)
    if (inherits(pL,"try-error") || suppressWarnings(is.nan(pL$value))) return(NA)
    pU <- my_integrate(f,lower=UC,upper=UT,p=ps,
                       dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun)
    if (inherits(pU,"try-error") || suppressWarnings(is.nan(pU$value))) return(NA)
    pmax(0,pmin(pL$value,1))+pmax(0,pmin(pU$value,1))
  }


  pars <- get_pars_matrix(p_vector,dadm)

  if (any(names(dadm)=="RACE")) # Some accumulators not present
    pars[as.numeric(dadm$lR)>as.numeric(as.character(dadm$RACE)),] <- NA

  if (is.null(attr(pars,"ok")))
    ok <- !logical(dim(pars)[1]) else ok <- attr(pars,"ok")

  lds <- numeric(dim(dadm)[1]) # log pdf (winner) or survivor (losers)
  lds[dadm$winner] <- log(attr(dadm,"model")()$dfun(rt=dadm$rt[dadm$winner],
                                                    pars=pars[dadm$winner,]))
  n_acc <- length(levels(dadm$R))
  if (n_acc>1) lds[!dadm$winner] <-
    log(1-attr(dadm,"model")()$pfun(rt=dadm$rt[!dadm$winner],pars=pars[!dadm$winner,]))
  lds[is.na(lds) | !ok] <- -Inf

  # Calculate truncation?
  LT <- attr(dadm,"LT")
  UT <- attr(dadm,"UT")
  dotrunc <- (!is.null(LT) | !is.null(UT))
  if (is.null(LT)) LT <- 0
  if (is.null(UT)) UT <- Inf

  # Calculate censoring
  LC <- attr(dadm,"LC")
  UC <- attr(dadm,"UC")

  # Response known
  # Fast
  nort <- dadm$rt==-Inf; nort[is.na(nort)] <- FALSE; nort <- nort & !is.na(dadm$R)
  if ( any(nort) ) {
    mpars <- array(pars[nort,,drop=FALSE],dim=c(n_acc,sum(nort)/n_acc,ncol(pars)),
                   dimnames = list(NULL,NULL,colnames(pars)))
    winner <- matrix(dadm$winner[nort],nrow=n_acc)
    tofix <- dadm$winner & nort
    for (i in 1:dim(mpars)[2]) {
      if (dim(mpars)[[1]]==1) pi <- t(as.matrix(mpars[,i,])) else
        pi <- mpars[,i,][order(!winner[,i]),]
      tmp <- my_integrate(f,lower=LT,upper=LC,p=pi,
                          dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun)
      if ( !inherits(tmp, "try-error") && suppressWarnings(!is.nan(log(tmp$value))) )
        lds[tofix][i] <- log(pmax(0,pmin(tmp$value,1)))
    }
  }
  # Slow
  nort <- dadm$rt==Inf; nort[is.na(nort)] <- FALSE; nort <- nort & !is.na(dadm$R)
  if ( any(nort) ) {
    mpars <- array(pars[nort,,drop=FALSE],dim=c(n_acc,sum(nort)/n_acc,ncol(pars)),
                   dimnames = list(NULL,NULL,colnames(pars)))
    winner <- matrix(dadm$winner[nort],nrow=n_acc)
    tofix <- dadm$winner & nort
    for (i in 1:dim(mpars)[2]) {
      if (dim(mpars)[[1]]==1) pi <- t(as.matrix(mpars[,i,])) else
        pi <- mpars[,i,][order(!winner[,i]),]
      tmp <- my_integrate(f,lower=UC,upper=UT,p=pi,
                          dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun)
      # if (!inherits(tmp, "try-error") && suppressWarnings(!is.nan(log(tmp$value))))
      #   lds[tofix][i] <- log(pmax(0,pmin(tmp$value,1)))
      if (!inherits(tmp, "try-error") && suppressWarnings(!is.nan(log(tmp$value))))
        lds[tofix][i] <- pmax(0,pmin(tmp$value,1))
      lds[tofix][i] <- log(lds[tofix][i] + (1-lds[tofix][i])*prod(pnorm(0,mpars[,i,"v"],mpars[,i,"sv"])))
    }
  }
  # No direction
  nort <- is.na(dadm$rt) & !is.na(dadm$R)
  if ( any(nort) ) {
    mpars <- array(pars[nort,,drop=FALSE],dim=c(n_acc,sum(nort)/n_acc,ncol(pars)),
                   dimnames = list(NULL,NULL,colnames(pars)))
    winner <- matrix(dadm$winner[nort],nrow=n_acc)
    tofix <- dadm$winner & nort
    for (i in 1:dim(mpars)[2]) {
      if (dim(mpars)[[1]]==1) pi <- t(as.matrix(mpars[,i,])) else
        pi <- mpars[,i,][order(!winner[,i]),]
      tmp <- my_integrate(f,lower=LT,upper=LC,p=pi,
                          dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun)
      if ( !inherits(tmp,"try-error") && suppressWarnings(!is.nan(tmp$value)) ) {
        p <- tmp$value
        if (dim(mpars)[[1]]==1) pi <- mpars[,i,,drop=FALSE] else
          pi <- mpars[,i,][order(!winner[,i]),]
        tmp <- my_integrate(f,lower=UC,upper=UT,p=pi,
                            dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun)
        if ( !inherits(tmp,"try-error") && suppressWarnings(!is.nan(tmp$value)) )
          p <- p + tmp$value else p <- 0
      } else p <- 0
      lds[tofix][i] <- log(pmax(0,pmin(p,1)))
    }
  }

  # Response unknown.
  # Fast
  nort <- dadm$rt==-Inf; nort[is.na(nort)] <- FALSE; nort <- nort & is.na(dadm$R)
  if ( any(nort) ) {
    mpars <- array(pars[nort,,drop=FALSE],dim=c(n_acc,sum(nort)/n_acc,ncol(pars)),
                   dimnames = list(NULL,NULL,colnames(pars)))
    winner <- matrix(dadm$winner[nort],nrow=n_acc)
    tofixfast <- dadm$winner & nort
    for (i in 1:dim(mpars)[2]) {
      if (dim(mpars)[[1]]==1) pi <- t(as.matrix(mpars[,i,])) else
        pi <- mpars[,i,]
      pc <- my_integrate(f,lower=LT,upper=LC,p=pi,
                         dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun)
      if (inherits(pc, "try-error") || suppressWarnings(is.nan(pc$value)))
        p <- NA else p <- pmax(0,pmin(pc$value,1))
      if (!is.na(p)) {
        if (p != 0 && !(LT==0 & UT==Inf))  cf <- pr_pt(LT,UT,mpars[,i,],dadm) else cf <- 1
        if (!is.na(cf)) p <- p*cf
      }
      if (!is.na(p) & n_acc>1) for (j in 2:n_acc) {
        pc <- my_integrate(f,lower=LT,upper=LC,p=mpars[,i,][c(j,c(1:n_acc)[-j]),],
                           dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun)
        if (inherits(pc, "try-error") || suppressWarnings(is.nan(pc$value))) {
          p <- NA; break
        }
        if (pc$value != 0 & !(LT==0 & UT==Inf))
          cf <- pr_pt(LT,UT,mpars[,i,][c(j,c(1:n_acc)[-j]),],dadm) else cf <- 1
          if (!is.na(cf)) p <- p + pc$value*cf
      }
      lp <- log(p)
      if (!is.nan(lp) & !is.na(lp)) lds[tofixfast][i] <- lp else lds[tofixfast][i] <- -Inf
    }
  } else tofixfast <- NA
  # Slow
  nort <- dadm$rt==Inf; nort[is.na(nort)] <- FALSE; nort <- nort & is.na(dadm$R)
  if ( any(nort) ) {
    mpars <- array(pars[nort,,drop=FALSE],dim=c(n_acc,sum(nort)/n_acc,ncol(pars)),
                   dimnames = list(NULL,NULL,colnames(pars)))
    winner <- matrix(dadm$winner[nort],nrow=n_acc)
    tofixslow <- dadm$winner & nort
    for (i in 1:dim(mpars)[2]) {
      if (dim(mpars)[[1]]==1) pi <- t(as.matrix(mpars[,i,])) else
        pi <- mpars[,i,]
      pc <- my_integrate(f,lower=UC,upper=UT,p=pi,
                         dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun)
      if (inherits(pc, "try-error") || suppressWarnings(is.nan(pc$value)))
        p <- NA else p <- pmax(0,pmin(pc$value,1))
      if (!is.na(p)) {
        if (p != 0 && !(LT==0 & UT==Inf))  cf <- pr_pt(LT,UT,mpars[,i,],dadm) else cf <- 1
        if (!is.na(cf)) p <- p*cf
      }
      if (!is.na(p) & n_acc>1) for (j in 2:n_acc) {
        pc <- my_integrate(f,lower=UC,upper=UT,p=mpars[,i,][c(j,c(1:n_acc)[-j]),],
                           dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun)
        if (inherits(pc, "try-error") || suppressWarnings(is.nan(pc$value))) {
          p <- NA; break
        }
        if (pc$value != 0 & !(LT==0 & UT==Inf))
          cf <- pr_pt(LT,UT,mpars[,i,][c(j,c(1:n_acc)[-j]),],dadm) else cf <- 1
          if (!is.na(cf)) p <- p + pc$value*cf
      }
      lp <- log(p + (1-p)*prod(pnorm(0,mpars[,i,"v"],mpars[,i,"sv"])))
      if (!is.nan(lp) & !is.na(lp)) lds[tofixslow][i] <- lp else lds[tofixslow][i] <- -Inf
    }
  } else tofixslow <- NA
  # no direction
  nort <- is.na(dadm$rt) & is.na(dadm$R)
  nort <- nort & (pars[,"pContaminant"] == 0) # Otherwise not identifiable
  if ( any(nort) ) {
    mpars <- array(pars[nort,,drop=FALSE],dim=c(n_acc,sum(nort)/n_acc,ncol(pars)),
                   dimnames = list(NULL,NULL,colnames(pars)))
    winner <- matrix(dadm$winner[nort],nrow=n_acc)
    tofix <- dadm$winner & nort
    for (i in 1:dim(mpars)[2]) {
      if (dim(mpars)[[1]]==1) pi <- t(as.matrix(mpars[,i,])) else
        pi <- mpars[,i,]
      pc <- pLU(LT,LC,UC,UT,pi,dadm)
      if (is.na(pc)) p <- NA else {
        if (pc!=0 & !(LT==0 & UT==Inf)) cf <- pr_pt(LT,UT,pi,dadm) else cf <- 1
        if (!is.na(cf)) p <- pc*cf else p <- NA
        if (!is.na(p) & n_acc>1) for (j in 2:n_acc) {
          pc <- pLU(LT,LC,UC,UT,mpars[,i,][c(j,c(1:n_acc)[-j]),],dadm)
          if (is.na(pc)) {
            p <- NA; break
          }
          if (pc!=0 & !(LT==0 & UT==Inf))
            cf <- pr_pt(LT,UT,mpars[,i,][c(j,c(1:n_acc)[-j]),],dadm) else cf <- 1
            if (is.na(cf)) {
              p <- NA; break
            } else p <- p +  pc*cf
        }
      }
      lp <- log(pmax(0,pmin(p,1)))
      if (!is.nan(lp) & !is.na(lp)) lds[tofix][i] <- lp
    }
  }


  # Truncation where not censored or censored and response known
  ok <- is.finite(lds[attr(dadm,"unique_nort") & dadm$winner])
  alreadyfixed <- is.na(dadm$R[attr(dadm,"unique_nort") & dadm$winner])
  ok <- ok & !alreadyfixed
  if ( dotrunc & any(ok) ) {
    tpars <- pars[attr(dadm,"unique_nort"),,drop=FALSE]
    tpars <- array(tpars[,,drop=FALSE],dim=c(n_acc,nrow(tpars)/n_acc,ncol(tpars)),
                   dimnames = list(NULL,NULL,colnames(tpars)))[,,,drop=FALSE]
    winner <- matrix(dadm$winner[attr(dadm,"unique_nort")],nrow=n_acc)[,,drop=FALSE]
    cf <- rep(NA,length(ok))
    for (i in 1:length(ok)) if (ok[i]) {
      if (dim(tpars)[[1]]==1) pi <- t(as.matrix(tpars[,i,])) else
        pi <- tpars[,i,][order(!winner[,i]),]
      cf[i] <- pr_pt(LT,UT,pi,dadm)
    }
    cf <- rep(log(cf),each=n_acc)[attr(dadm,"expand_nort")]
    fix <- dadm$winner & !is.na(cf) & !is.nan(cf) & is.finite(cf)
    if (any(fix)) lds[fix] <- lds[fix] + cf[fix]
    badfix <- dadm$winner & (is.na(cf) | is.nan(cf) | is.infinite(cf))
    if (!all(is.na(tofixfast))) badfix <- badfix & !tofixfast
    if (!all(is.na(tofixslow))) badfix <- badfix & !tofixslow
    if (any(badfix)) lds[badfix] <- -Inf
  }

  # Non-process (contaminant) miss.
  ispContaminant <- pars[,"pContaminant"]>0
  if ( any(ispContaminant) ) {
    p <- exp(lds[dadm$winner])
    pc <- pars[dadm$winner,"pContaminant"]
    isMiss <- is.na(dadm$R[dadm$winner])
    p[isMiss] <- pc[isMiss] + (1-pc[isMiss])*p[isMiss]
    p[!isMiss] <- (1-pc[!isMiss])*p[!isMiss]
    lds[dadm$winner] <- log(p)
  }

  if (n_acc>1) {
    ll <- lds[dadm$winner]
    if (n_acc==2) {
      ll <- ll + lds[!dadm$winner]
    } else {
      ll <- ll + apply(matrix(lds[!dadm$winner],nrow=n_acc-1),2,sum)
    }
  } else ll <- lds
  ll[is.na(ll) | is.nan(ll)] <- -Inf

  return(sum(pmax(min_ll,ll[attr(dadm,"expand_winner")])))
}

# log_likelihood_mt <- function(p_vector,dadm,min_ll=log(1e-10))
#   # Multiple threshold (BE or TC) summed log likelihood
#   # attr(dadm,"dL")
# {
#
#   pars <- get_pars_matrix(p_vector,dadm)
#
#   Dnams <- dimnames(pars)[[2]][substr(dimnames(pars)[[2]],1,2)=="DT"]
#   tmats <- array(cbind(rep(0,dim(pars)[1]),pars[,c(Dnams,"b")]),dim=c(2,nrow(pars)/2,2+length(Dnams)),
#                 dimnames=list(NULL,NULL,c("DT0",Dnams,"b")))
#   pnams <- dimnames(pars)[[2]][!(dimnames(pars)[[2]] %in% c(Dnams,"b"))]
#   pmats <- array(pars[,pnams],dim=c(2,nrow(pars)/2,length(pnams)),
#                 dimnames=list(NULL,NULL,pnams))
#   nt <- dim(pmats)[2]
#   like <- numeric(nt)
#   for (i in 1:nt) {
#     i2 <- i*2
#     # Get look up table for current rating (dadm$R[i2])
#     pick <- attr(dadm,"dL")[as.numeric(dadm$R[i2]),]
#     tmp <- try(n1PDF_MTR_1(rt=dadm$rt[i*2], pars = pmats[c(pick[1],pick[3]),i,],
#       dl = tmats[pick[3],i,pick[4]], du = tmats[pick[5],i,pick[6]], b = tmats[pick[1],i,pick[2]]),silent=TRUE)
#     if (inherits(tmp,"try-error")) like[i] <- 0 else like[i] <- tmp
#   }
#   ll <- log(ll)
#   ll[is.na(ll) | is.nan(ll) | ll == Inf] <- 0
#   return(sum(pmax(min_ll,like[attr(dadm,"expand_winner")])))
# }

#### Multiple threshold ----

log_likelihood_mt <- function(p_vector,dadm,min_ll=log(1e-10),n_cores=10)
  # Multiple threshold (BE or TC) summed log likelihood
  # attr(dadm,"dL")
{

  mt <- function(i,dadm,tmats) {
    i2 <- i*2
    # Get look up table for current rating (dadm$R[i2])
    pick <- attr(dadm,"dL")[as.numeric(dadm$R[i2]),]
    tmp <- try(n1PDF_MTR_1(rt=dadm$rt[i*2], pars = pmats[c(pick[1],pick[3]),i,],
                           dl = tmats[pick[3],i,pick[4]], du = tmats[pick[5],i,pick[6]],
                           b = tmats[pick[1],i,pick[2]]),silent=TRUE)
    if (inherits(tmp,"try-error")) 0 else tmp
  }

  pars <- get_pars_matrix(p_vector,dadm)

  Dnams <- dimnames(pars)[[2]][substr(dimnames(pars)[[2]],1,2)=="DT"]
  tmats <- array(cbind(rep(0,dim(pars)[1]),pars[,c(Dnams,"b")]),dim=c(2,nrow(pars)/2,2+length(Dnams)),
                 dimnames=list(NULL,NULL,c("DT0",Dnams,"b")))
  pnams <- dimnames(pars)[[2]][!(dimnames(pars)[[2]] %in% c(Dnams,"b"))]
  pmats <- array(pars[,pnams],dim=c(2,nrow(pars)/2,length(pnams)),
                 dimnames=list(NULL,NULL,pnams))
  nt <- dim(pmats)[2]
  ll <- log(unlist(mclapply(1:nt,mt,dadm=dadm,tmats=tmats,mc.cores=n_cores)))
  ll[is.na(ll) | is.nan(ll) | ll == Inf] <- -Inf
  return(sum(pmax(min_ll,ll[attr(dadm,"expand_winner")])))
}

#### Stop signal ----
  # ptrials <- rep(trials,each=n_acc)        # trial number for each row
  #     # Used with trials (length = n_trials)
  # istGOwin <- # Go wins (either on GO or STOP trials)
  #   apply(matrix(dadm$winner[ispGOacc],nrow=n_accG),2,any)
  # tGOwin <- trials[istGOwin]
  #   # Used pars and dadm  (length = nrow(dadm))
  # ispGOwin <- ptrials %in% tGOwin # GO wins rows
  #
  #   # Used with trials (length = n_trials)
  # istStop <- ispStop[isp1]                # stop trials
  # tStop <- c(trials)[istStop]             # stop trials
  # istStopGOwin <- (tGOwin %in% tStop)     # Stop trials where GO wins
  # istStopGOwin <- istStopGOwin


  # # Go failure index, no response and either go trial or stop trial with one or
  # # more stop-triggered runners
  # istgf <- is.na(dadm$R[isp1]) & (!istStop | (istStop & (n_accS>0)))
  #
  # # Stop trials where stop wins and not gf
  # tStopStopWin <- tStop[!(tStop %in% c(tGOwin,trials[istgf]))]
  #
  #
  # # 5 cases, first two in this branch
  # # 1) go wins, no stop runner
  # # 2) go wins with stop runner (same as case 1 + stop survivor)




log_likelihood_race_ss <- function(p_vector,dadm,min_ll=log(1e-10))
{

  pars <- get_pars_matrix(p_vector,dadm)

  # Set up indices:
  # "is" = logical, "isp" pars/dadm index, "t" trials index, "ist" logical on trial index
  # "n_" number of integer
  n_acc <- length(levels(dadm$lR))                   # total number of accumulators

  if (is.null(attr(pars,"ok")))
    ok <- !logical(dim(pars)[1]) else ok <- attr(pars,"ok")
  if (!any(ok)) return(min_ll*length(attr(dadm, "expand_winner")))

  # # spurious go winners on no-response trials
  # dadm$winner[is.na(dadm$R)] <- FALSE
if (any(is.infinite(dadm$rt))) stop("BUGGER!")

  # Counts
  n_accG <- sum(as.numeric(dadm[1:n_acc,"lI"])==2)   # go accumulators
  n_accST <- sum(as.numeric(dadm[1:n_acc,"lI"])==1)  # stop-triggered accumulators
  GOR <- levels(dadm$lR)[as.numeric(dadm[1:n_acc,"lI"])==2]
  if (n_accST>0) STR <- levels(dadm$lR)[as.numeric(dadm[1:n_acc,"lI"])==1]

  # Likelihood for all trials and for ok trials
  allLL <- rep(-Inf,nrow(dadm)/n_acc)
  allok <- ok[c(dadm$lR==levels(dadm$lR)[1])] # used to put back into allLL

  # remove bad trials
  pars <- pars[ok,,drop=FALSE]
  dadm <- dadm[ok,,drop=FALSE]
  isp1 <- dadm$lR==levels(dadm$lR)[1]      # 1st accumulator rows
  ispGOacc <- dadm$lI==levels(dadm$lI)[2] # Go accumulator rows
  ispStop <- is.finite(dadm$SSD) # stop-trial rows
  gf <- pars[isp1,"gf"]
  tf <- pars[isp1,"tf"]

  # No response
  ispNR <- is.na(dadm$R)
  if (any(ispNR)) { # Test as some models always respond
    ispgoNR <- ispNR & !ispStop
    tgoNR <- c(1:sum(isp1))[ispgoNR[isp1]]
    if (any(ispgoNR)) allLL[allok][tgoNR] <- log(gf[tgoNR]) # Go trials
    ispstopNR <- ispNR & ispStop
    tstopNR <- (rep(1:(nrow(dadm)/n_acc),each=n_acc))[ispstopNR & isp1]
    if (any(ispstopNR)) { # Stop trials
      if (n_accST>0) pStop <- 0 else
        pStop <- pmax(0,attr(dadm,"model")()$sfun(pars[ispStop & ispGOacc & ispNR,,drop=FALSE],n_acc=n_accG))
      allLL[allok][tstopNR] <- log(gf[tstopNR] + (1-gf[tstopNR])*(1-tf[tstopNR])*pStop)
    }
   }

  # remove no response trials
  allr <- !ispNR[isp1] # used to put back into allLL
  pars <- pars[!ispNR,,drop=FALSE]
  dadm <- dadm[!ispNR,,drop=FALSE]
  isp1 <- dadm$lR==levels(dadm$lR)[1]      # 1st accumulator rows
  ispGOacc <- dadm$lI==levels(dadm$lI)[2] # Go accumulator rows
  ispStop <- is.finite(dadm$SSD) # stop-trial rows with a response
  gf <- pars[isp1,"gf"]
  tf <- pars[isp1,"tf"]

  n_trials <- nrow(dadm)/n_acc # number of trial
  trials <- 1:n_trials         # trial number
  ptrials <- rep(trials,each=n_acc)
  accST <- c(1:n_acc)[pars[1:n_acc,"lI"]==1]

  like <- numeric(n_trials)
  lds <- numeric(nrow(dadm)) # log density and survivor, used for both go and stop trials

  # Go trials with response
  if (any(!ispStop)) {
    ispGOwin <-  !ispStop & dadm$winner # Winner go accumulator rows
    tGO <- trials[!ispStop[isp1]]
    # Winner density
    lds[ispGOwin] <- log(attr(dadm,"model")()$dfunG(
      rt=dadm[ispGOwin,"rt"],pars=pars[ispGOwin,,drop=FALSE]))
    like[tGO] <- lds[ispGOwin]
    if (n_accG >1) {  # Looser survivor go accumulator(s)
      ispGOloss <- !ispStop & !dadm$winner & ispGOacc # Looser go accumulator rows
      lds[ispGOloss] <- log(1-attr(dadm,"model")()$pfunG(
         rt=dadm$rt[ispGOloss],pars=pars[ispGOloss,,drop=FALSE]))
      like[tGO] <- like[tGO] + apply(matrix(lds[ispGOloss],nrow=n_accG-1),2,sum)
    }
    like[tGO] <- (1-gf[tGO])*exp(like[tGO])
  }

  # Stop trials with a response
  if (any(ispStop)) {
    # Stop looses
    ispSwin <-       ispStop & dadm$winner              # Winner go of ST accumulator rows
    ispSlossGOacc <- ispStop & !dadm$winner & ispGOacc  # Loosing go accumulator rows
    ispSlossSTacc <- ispStop & !dadm$winner & !ispGOacc # Loosing ST accumulator rows

    # pStop at observed rt if ST present (calculate before correcting rt)
    if (n_accST>0) {
      tST <- trials[ispStop[isp1] & as.numeric(dadm$lI)[dadm$winner] == 1]
      ispST <- ptrials %in% tST
      if (any(ispST)) {
        pStop <- numeric(n_trials)
        upper <- dadm$rt[dadm$winner][tST]
        pStop[tST] <- pmax(0,attr(dadm,"model")()$sfun(pars[ispST,,drop=FALSE],
          upper=upper,n_acc=n_acc,st=c(1,1+accST)))
      }
    }

    # For following race calculations correct rt with SSD for ST accumulators
    if (any(ispSlossSTacc)) dadm[ispSlossSTacc,"rt"] <-
        dadm[ispSlossSTacc,"rt"]-dadm[ispSlossSTacc,"SSD"]

    # Fill in lds and sums over survivors for race
    lds[ispSwin] <- log(attr(dadm,"model")()$dfunG(
      rt=dadm[ispSwin,"rt"],pars=pars[ispSwin,,drop=FALSE]))
    if (n_acc >1) {  # Survivor for looser go and/or ST accumulator(s)
      lds[ispSlossGOacc | ispSlossSTacc] <- log(1-attr(dadm,"model")()$pfunG(
          rt=dadm[ispSlossGOacc | ispSlossSTacc,"rt"],
          pars=pars[ispSlossGOacc | ispSlossSTacc,,drop=FALSE]))
      # Sum survivor over loosing ST and GO accumulators
      SSTGO <- tapply(lds[ispSlossGOacc | ispSlossSTacc],
        cbind.data.frame(trials=ptrials[ispSlossGOacc | ispSlossSTacc],
                         lI=dadm[ispSlossGOacc | ispSlossSTacc,"lI"]),sum)
      SSTGO[is.na(SSTGO)] <- 0 # cases where no ST or GO survivor
    } else SSTGO <- matrix(0,ncol=2)

    # Stop accumulator survivor
    sStop <- log(1-attr(dadm,"model")()$pfunS(
      rt=dadm[ispSwin,"rt"],pars=pars[ispSwin,,drop=FALSE]))

    # Get like
    tS <- trials[ispStop[isp1]]
    # Sub-select from tS
    istSgo <- dadm$R[isp1][tS] %in% GOR
    # Stop looses or not present, can produce ST but only if wins absolute race
    like[tS] <- # no failures, works for all responses
      (1-gf[tS])*(1-tf[tS])*exp(lds[ispSwin]+SSTGO[,1]+SSTGO[,2]+sStop)
    if (any(istSgo)) like[tS][istSgo] <- like[tS][istSgo] + # tf (no stop runner), no gf, only produces GO responses
          (1-gf[tS][istSgo])*(tf[tS][istSgo])*exp(lds[ispSwin][istSgo]+SSTGO[,2][istSgo])
    # If both tf and gf then no response, handled previously

    ### Stop wins at some time before rt, and so must be ST response, add to ST race winners
    if (n_accST>0) {
      istSst <- dadm$R[isp1][tS] %in% STR
      like[tS][istSst] <- like[tS][istSst] + # gf (no go runner) no tf, only produces ST responses
          (gf[tS][istSst])*(1-tf[tS][istSst])*exp(lds[ispSwin][istSst]+SSTGO[,1][istSst]+sStop[istSst])
      SST <- numeric(n_trials)
      # Winner and looser ST density already computed.
      SST[tST] <- lds[dadm$winner][tST]
      if (n_accST>1) SST[tST] <- SST[tST] + apply(matrix(
        lds[(ptrials %in% tST) & !dadm$winner & !ispGOacc],nrow=n_accST-1),2,sum)
      like[tS][istSst] <- like[tS][istSst] + pStop[tS][istSst]*exp(SST[tST])
    }
  }
  allLL[allok][allr] <- log(like)
  sum(pmax(min_ll,allLL[attr(dadm,"expand_winner")]))
}


log_likelihood_race_ss_uv <- function(p_vector,dadm,min_ll=log(1e-10))
{

  # univalent responding, assumes go and stop-triggered in the same order and
  # same length.

  pars <- get_pars_matrix(p_vector,dadm)

  if (is.null(attr(pars,"ok")))
    ok <- !logical(dim(pars)[1]) else ok <- attr(pars,"ok")
  if (!any(ok)) return(min_ll*nrow(dadm)/n_acc)

  # Likelihood for all trials and for ok trials
  alllike <- numeric(nrow(dadm)/length(levels(dadm$lR)))
  isp1 <- dadm$lR==levels(dadm$lR)[1]
  isnr <- is.na(dadm$R[isp1])
  alllike[isnr] <- pars[isnr & isp1,"gf"] # only way to get NA
  ok <- ok & !isnr

  # remove bad trials
  goodt <- ok[isp1]
  like <- alllike[goodt]
  pars <- pars[ok,,drop=FALSE]
  dadm <- dadm[ok,,drop=FALSE]
  isp1 <- dadm$lR==levels(dadm$lR)[1]

  # failures
  gf <- pars[isp1,"gf"]
  tf <- pars[isp1,"tf"]

  # number of accumulators
  n_acc <- length(levels(dadm$lR))/2 # must be half go half stop triggered
  n_trials <- nrow(dadm)/(2*n_acc) # number of trials
  ispGOacc <- dadm$lI==levels(dadm$lI)[2] # Go accumulator rows

  # Compute on indices to minimize object copying

  # Stop indices
    # Used pars and dadm  (length = nrow(dadm))
  ispStop <- is.finite(dadm$SSD) # stop-trial rows
    # Used with trials (length = n_trials)
  istStop <- c(1:n_trials)[ispStop[isp1]] # stop trials

  # Stop-triggered equivalent of go winner (assumes 1st go = 1st stop triggered etc.)
  isSTwinner <- logical(nrow(dadm))
  isSTwinner[!ispGOacc] <- dadm$winner[ispGOacc]

  # Go
  lds <- numeric(nrow(dadm)/2) # log pdf (winner) or survivor (losers)
  lds[dadm$winner] <- log(attr(dadm,"model")()$dfunG(
      rt=dadm[dadm$winner,"rt"],pars=pars[dadm$winner & ispGOacc,]))
  if (n_acc>1) {
    lds[!dadm$winner & ispGOacc] <- log(1-attr(dadm,"model")()$pfunG(
      rt=dadm$rt[dadm$winner],pars=pars[!dadm$winner & ispGOacc,,drop=FALSE]))
    if (n_acc==2) {
      lds[dadm$winner] <- lds[dadm$winner] + lds[!dadm$winner & ispGOacc]
    } else {
      lds[dadm$winner] <- lds[dadm$winner] + apply(matrix(lds[!dadm$winner & ispGOacc],nrow=n_acc-1),2,sum)
    }
  }
  like <- lds[dadm$winner]  # exponentiated to make like below
  if ( any(istStop)) { # stop trials
    # Add p stop loss to go win on stop trial
    like[istStop] <- like[istStop] + log(1-attr(dadm,"model")()$pfunS(
      rt=dadm$rt[isp1][istStop],pars=pars[ispStop & isp1,,drop=FALSE]))

    # Stop wins
    pStop <- pmax(0,attr(dadm,"model")()$sfun(pars[ispStop & ispGOacc,,drop=FALSE],n_acc))

    # Get like for stop-triggered accumulators
    lds <- numeric(sum(ispStop & !ispGOacc)) # log pdf (winner) or survivor (losers)
    lds[isSTwinner[ispStop & !ispGOacc]] <- log(attr(dadm,"model")()$dfunG(
      rt=dadm$rt[isp1][istStop]-dadm$SSD[isp1][istStop],
      pars=pars[isSTwinner & ispStop & !ispGOacc,]))
    lds[!isSTwinner[ispStop & !ispGOacc]] <- log(1-attr(dadm,"model")()$pfunG(
      rt=dadm$rt[isp1][istStop]--dadm$SSD[isp1][istStop],
      pars=pars[!isSTwinner & ispStop & !ispGOacc,,drop=FALSE]))
    if (n_acc==2) { # must have at least 2 accumulators
      lds[isSTwinner[ispStop & !ispGOacc]] <-
        lds[isSTwinner[ispStop & !ispGOacc]] + lds[!isSTwinner[ispStop & !ispGOacc]]
    } else {
      lds[isSTwinner[ispStop & !ispGOacc]] <- lds[isSTwinner[ispStop & !ispGOacc]] +
        apply(matrix(lds[!isSTwinner[ispStop & !ispGOacc]],nrow=n_acc-1),2,sum)
    }

    # Put three cases together
    like <- exp(like)
    like[istStop] <- tf[istStop]*like[istStop] +        # Trigger failure
      (1-tf[istStop])*((1-pStop)*like[istStop] +        # No trigger failure  but stop looses
       pStop*exp(lds[isSTwinner[ispStop & !ispGOacc]])) # No trigger failure + stop-triggered
  } else like <- exp(like)
  like[like<0 | is.na(like) | is.nan(like)] <- 0

  # Expand results
  alllike[goodt] <- (1-gf)*like
  sum(pmax(min_ll,log(alllike)[attr(dadm,"expand_winner")]))
}

