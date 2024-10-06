#' The Diffusion Decision Model
#'
#' @return A model list with all the necessary functions to sample
#' @export
DDMnoC <- function(){
  list(
    type="DDM",
    p_types=c("v" = 1,"a" = log(1),"sv" = log(0),"t0" = log(0),"st0" = log(0),
              "s" = log(1),"Z" = qnorm(0.5),"SZ" = qnorm(0),"DP" = qnorm(0.5)),
    ptransform=list(transform=c(v = "identity",a = "exp",sv = "exp",t0 = "exp",
      st0 = "exp",s = "exp",Z = "pnorm",SZ = "pnorm",DP = "pnorm"),lower=c(t0=.05)),
    bound=list(minmax=cbind(v=c(-20,20),a=c(0,10),Z=c(.001,.999),t0=c(0,Inf),
                sv=c(.001,10),s=c(0,Inf),SZ=c(.001,.999),st0=c(0,.5),DP=c(0,1)),
                exceptions=c(sv=0,SZ=0,st0=0)),
    # The "TZD" parameterization defined relative to the "rtdists" package is:
    # natural scale
    #   v = rtdists rate v (positive favors upper)
    # log scale
    #   t0 > 0: lower bound of non-decision time
    #   st0 > 0: rtdists width of non-decision time distribution
    #   a > 0: rtdists upper threshold, a
    #   sv > 0: rtdists v standard deviation sv
    #   s > 0: rtdists moment-to-moment standard deviation, s
    # probit scale
    #   0 < Z < 1: rtdists start point z = Z*a
    #   0 < SZ < 1: rtdists start-point variability, sz = 2*SZ*min(c(a*Z,a*(1-Z))
    #   0 < DP < 1: rtdists d = t0(upper)-t0(lower) = (2*DP-1)*t0  #
    #
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      pars <- cbind(pars,z=pars[,"a"]*pars[,"Z"],
                    sz = 2*pars[,"SZ"]*pars[,"a"]*apply(cbind(pars[,"Z"],1-pars[,"Z"]),1,min))
      pars <- cbind(pars, d = pars[,"t0"]*(2*pars[,"DP"]-1))
      if (!is.null(attr(dadm,"adaptive"))) pars <- do_adaptive(pars,dadm)
      pars
    },
    # p_vector transform
    transform = function(p) p,
    # Random function
    rfun=function(lR=NULL,pars,ok) rDDM(lR,pars,ok),
    # Density function (PDF)
    dfun=function(rt,R,pars) dDDM(rt,R,pars),
    # Probability function (CDF)
    pfun=function(rt,R,pars) pDDM(rt,R,pars),
    log_likelihood=function(p_vector,dadm,min_ll=log(1e-10)){
      log_likelihood_ddm(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
    }
  )
}



#### Missing ----

#' MDDM: Missing value DDM with contamination
#'
#' @return A model list with all the necessary functions to sample
#' @export
MDDM <- function(){
  list(
    type="DDM",
    p_types=c("v" = 1,"a" = log(1),"sv" = log(0),"t0" = log(0),"st0" = log(0),
              "s" = log(1),"Z" = qnorm(0.5),"SZ" = qnorm(0),"DP" = qnorm(0.5),
              pContaminant = qnorm(0)),
    ptransform=list(transform=c(v = "identity",a = "exp",sv = "exp",t0 = "exp",
      st0 = "exp",s = "exp",Z = "pnorm",SZ = "pnorm",DP = "pnorm"),lower=c(t0=.05),
      pContaminant="pnorm"),
    bound=list(minmax=cbind(v=c(-20,20),a=c(0,10),Z=c(.001,.999),t0=c(0,Inf),
                sv=c(.001,10),s=c(0,Inf),SZ=c(.001,.999),st0=c(0,.5),DP=c(0,1),
                pContaminant=c(0,1)),
                exceptions=c(sv=0,SZ=0,st0=0,pContaminant=0)),
    # The "TZD" parameterization defined relative to the "rtdists" package is:
    # natural scale
    #   v = rtdists rate v (positive favors upper)
    # log scale
    #   t0 > 0: lower bound of non-decision time
    #   st0 > 0: rtdists width of non-decision time distribution
    #   a > 0: rtdists upper threshold, a
    #   sv > 0: rtdists v standard deviation sv
    #   s > 0: rtdists moment-to-moment standard deviation, s
    # probit scale
    #   0 < Z < 1: rtdists start point z = Z*a
    #   0 < SZ < 1: rtdists start-point variability, sz = 2*SZ*min(c(a*Z,a*(1-Z))
    #   0 < DP < 1: rtdists d = t0(upper)-t0(lower) = (2*DP-1)*t0  #
    #
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      pars <- cbind(pars,z=pars[,"a"]*pars[,"Z"],
                    sz = 2*pars[,"SZ"]*pars[,"a"]*apply(cbind(pars[,"Z"],1-pars[,"Z"]),1,min))
      pars <- cbind(pars, d = pars[,"t0"]*(2*pars[,"DP"]-1))
      pars
    },
    # p_vector transform, sets s as a scaling parameter
    transform = function(p) p,
    # Random function
    rfun=function(lR=NULL,pars,ok) rDDM(lR,pars,ok),
    # Density function (PDF)
    dfun=function(rt,R,pars) dDDM(rt,R,pars),
    # Probability function (CDF)
    pfun=function(rt,R,pars) pDDM(rt,R,pars),
    log_likelihood=function(p_vector,dadm,min_ll=log(1e-10)){
      log_likelihood_ddm_missing(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
    }
  )
}


#### Advantage ----

#' Advantage DDM, allowing for sv increasing with stimulus sum
#'
#' @return A list defining the cognitive model
#' @export
ADDMnoC <- function(){
  list(
    type="DDM",
    p_types=c("v" = 1,"a" = log(1),"sv" = log(0),"t0" = log(0),"st0" = log(0),
              "s" = log(1),"Z" = qnorm(0.5),"SZ" = qnorm(0),"DP" = qnorm(0.5),
              "AD"=log(1),ssv=log(0),SD=0,SS=0,DD=0,DS=0),
    ptransform=list(transform=c(v = "identity",a = "exp",sv = "exp",t0 = "exp",
      st0 = "exp",s = "exp",Z = "pnorm",SZ = "pnorm",DP = "pnorm",
      AD="log",ssv="log"),lower=c(t0=.05)),
    bound=list(minmax=cbind(v=c(-20,20),a=c(0,10),Z=c(.001,.999),t0=c(0,Inf),
                sv=c(.001,10),s=c(0,Inf),SZ=c(.001,.999),st0=c(0,.5),DP=c(0,1)),
                exceptions=c(sv=0,SZ=0,st0=0,AD=c(0,Inf),ssv=c(0,Inf))),
    # The "TZD" parameterization defined relative to the "rtdists" package is:
    # natural scale
    #   v = rtdists rate v (positive favors upper)
    # log scale
    #   t0 > 0: lower bound of non-decision time
    #   st0 > 0: rtdists width of non-decision time distribution
    #   a > 0: rtdists upper threshold, a
    #   sv > 0: rtdists v standard deviation sv
    #   s > 0: rtdists moment-to-moment standard deviation, s
    # probit scale
    #   0 < Z < 1: rtdists start point z = Z*a
    #   0 < SZ < 1: rtdists start-point variability, sz = 2*SZ*min(c(a*Z,a*(1-Z))
    #   0 < DP < 1: rtdists d = t0(upper)-t0(lower) = (2*DP-1)*t0  #
    #
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      pars <- cbind(pars,z=pars[,"a"]*pars[,"Z"], d = pars[,"t0"]*(2*pars[,"DP"]-1),
                    sz = 2*pars[,"SZ"]*pars[,"a"]*apply(cbind(pars[,"Z"],1-pars[,"Z"]),1,min))
      if (!is.null(attr(dadm,"adaptive"))) pars <- do_adaptive(pars,dadm)
      pars <- advantage_pars_ddm(pars,dadm[attr(dadm,"advantage")[[1]]])
      pars
    },
    # p_vector transform, sets s as a scaling parameter
    transform = function(p) p,
    # Random function
    rfun=function(lR=NULL,pars,ok) rDDM(lR,pars,ok),
    # Density function (PDF)
    dfun=function(rt,R,pars) dDDM(rt,R,pars),
    # Probability function (CDF)
    pfun=function(rt,R,pars) pDDM(rt,R,pars),
    log_likelihood=function(p_vector,dadm,min_ll=log(1e-10)){
      log_likelihood_ddm(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
    }
  )
}

#' #' Advantage DDM, allowing for sv increasing with stimulus sum and 3 parameter shape mapping
#' #'
#' #' @return A list defining the cognitive model
#' #' @export
#' AiDDMnoC <- function(){
#'   list(
#'     type="DDM",
#'     p_types=c("v" = 1,"a" = log(1),"sv" = log(0),"t0" = log(0),"st0" = log(0),
#'               "s" = log(1),"Z" = qnorm(0.5),"SZ" = qnorm(0),"DP" = qnorm(0.5),
#'               "AD"=log(1),SD=0,S1=0,ssv=log(0),SW=qnorm(.5)),
#'     # The "TZD" parameterization defined relative to the "rtdists" package is:
#'     # natural scale
#'     #   v = rtdists rate v (positive favors upper)
#'     # log scale
#'     #   t0 > 0: lower bound of non-decision time
#'     #   st0 > 0: rtdists width of non-decision time distribution
#'     #   a > 0: rtdists upper threshold, a
#'     #   sv > 0: rtdists v standard deviation sv
#'     #   s > 0: rtdists moment-to-moment standard deviation, s
#'     # probit scale
#'     #   0 < Z < 1: rtdists start point z = Z*a
#'     #   0 < SZ < 1: rtdists start-point variability, sz = 2*SZ*min(c(a*Z,a*(1-Z))
#'     #   0 < DP < 1: rtdists d = t0(upper)-t0(lower) = (2*DP-1)*t0  #
#'     #
#'     # Trial dependent parameter transform
#'     Ttransform = function(pars,dadm) {
#'       pars <- cbind(pars,z=pars[,"a"]*pars[,"Z"], d = pars[,"t0"]*(2*pars[,"DP"]-1),
#'                     sz = 2*pars[,"SZ"]*pars[,"a"]*apply(cbind(pars[,"Z"],1-pars[,"Z"]),1,min))
#'       if (!is.null(attr(dadm,"adaptive"))) pars <- do_adaptive(pars,dadm)
#'       pars <- advantage_pars_ddm(pars,dadm[attr(dadm,"advantage")[[1]]],inv=TRUE)
#'       attr(pars,"ok") <-
#'         !( abs(pars[,"v"])> 20 | pars[,"a"]> 10 | pars[,"a"] < 0 | pars[,"sv"] > 10 |
#'              pars[,"SZ"]> .999 |pars[,"Z"]> .999 | pars[,"Z"] < .001 | pars[,"ssv"] < 0 |
#'              pars[,"t0"] < .05 | pars[,"st0"] > .5 | pars[,"st0"] < 0 | pars[,"DP"] < 0 | pars[,"DP"] > 1)
#'       if (pars[1,"sv"] !=0) attr(pars,"ok") <- attr(pars,"ok") & pars[,"sv"] > .001
#'       if (pars[1,"SZ"] !=0) attr(pars,"ok") <- attr(pars,"ok") & pars[,"SZ"] > .001
#'       pars
#'     },
#'     # p_vector transform, sets s as a scaling parameter
#'     transform = function(p) p,
#'     # Random function
#'     rfun=function(lR=NULL,pars) {
#'      ok <- !( abs(pars[,"v"])> 20 | pars[,"a"]> 10 | pars[,"a"] < 0 | pars[,"sv"] > 10 |
#'              pars[,"SZ"]> .999 |pars[,"Z"]> .999 | pars[,"Z"] < .001 |
#'              pars[,"t0"] < .05 | pars[,"st0"] > .5 | pars[,"st0"] < 0 | pars[,"DP"] < 0 | pars[,"DP"] > 1)
#'       if (pars[1,"sv"] !=0) attr(pars,"ok") <- attr(pars,"ok") & pars[,"sv"] > .001
#'       if (pars[1,"SZ"] !=0) attr(pars,"ok") <- attr(pars,"ok") & pars[,"SZ"] > .001
#'       if (is.null(lR)) ok else rDDM(lR,pars,precision=2.5,ok)
#'     },
#'     # Density function (PDF)
#'     dfun=function(rt,R,pars) dDDM(rt,R,pars,precision=2.5),
#'     # Probability function (CDF)
#'     pfun=function(rt,R,pars) pDDM(rt,R,pars,precision=2.5),
#'     log_likelihood=function(p_vector,dadm,min_ll=log(1e-10)){
#'       log_likelihood_ddm(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
#'     }
#'   )
#' }
