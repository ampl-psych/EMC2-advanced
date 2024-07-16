advantage_pars <- function(pars,SV,nr,pname="v")  {
  na <- nr-1                     # accumulators within responses
  nt <- nrow(pars)/nr            # trials
  SV <- pbeta(SV[[1]],exp((pars[,"SS"]+pars[,"SD"])/2),
              exp((pars[,"SS"]-pars[,"SD"])/2))
  q <- matrix(SV,nrow=nr)
  # sum and difference arrays, dimensions so as.vector fits with rpars
  s <- d <- array(dim=c(na,nr,ncol(q)))
  nrs <- 1:nr
  for (i in 1:na) {
    for (j in 1:nr) {
      notj <- nrs[-j]
      # Full looping algorithm in C would be over 1:nt
      d[i,j,] <- apply(q[c(notj[i],j),],2,diff)
      s[i,j,] <- apply(q[c(notj[i],j),],2,sum)
    }
  }
  # Beta difference mapping
  if (!all(c("DS","DD") %in% dimnames(pars)[[2]])) d <- (1+as.vector(d))/2 else
    d <- pbeta((1+as.vector(d))/2,exp((pars[,"DS"]+pars[,"DD"])/2),
               exp((pars[,"DS"]-pars[,"DD"])/2))
  # expand pars
  rpars <- apply(pars[,dimnames(pars)[[2]]!=pname],2,rep,each=na)
  v <- rep(pars[,pname],each=na) +
    rep(pars[,"AS"],each=na)*as.vector(s) +
    rep(pars[,"AD"],each=na)*(2*d-1)
  rpars <- cbind(v,rpars)
  dimnames(rpars)[[2]][1] <- pname
  rpars
}

advantage_pars_ddm <- function(pars,SV)  {
  # Beta stimulus value mapping
  SV1 <- pbeta(SV[[1]],exp((pars[,"SS"]+pars[,"SD"])/2),
               exp((pars[,"SS"]-pars[,"SD"])/2))
  SV2 <- pbeta(SV[[2]],exp((pars[,"SS"]+pars[,"SD"])/2),
               exp((pars[,"SS"]-pars[,"SD"])/2))
  # Beta difference mapping
  if (!all(c("DS","DD") %in% dimnames(pars)[[2]])) d <- (1+SV2-SV1)/2 else
    d <- pbeta((1+SV2-SV1)/2,exp((pars[,"DS"]+pars[,"DD"])/2),
               exp((pars[,"DS"]-pars[,"DD"])/2))
  pars[,"v"] <- pars[,"v"] + pars[,"AD"]*(2*d-1)
  pars[,"sv"] <- pars[,"sv"] + pars[,"ssv"]*(SV2^2+SV1^2)^(1/2)
  pars
}
