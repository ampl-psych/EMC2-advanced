rm(list=ls())
devtools::load_all()
#library(EMC2)

#set.seed(123) # miss.sigmaS = 0.4000000
#set.seed(234) # miss.sigmaS = 0.4000000, miss.tauS = 0.400000
set.seed(678) # miss.sigmaS = 0.4000000, miss.tauS = 0.400000, miss.muS = 0.40000
# set.seed(971) # good


#### Setup ----
# We use the usual contrast for average and difference. Note that as we are
# working with the exGaussian we expect d to be negative because match will be
# less than mismatch
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))

# staircase is a replacement for the model random function which runs a
# staircase algorithm on a dadm given a corresponding set of parameters (pars),
# returning a data frame with R, rt and SSD for each trial
#
# If SSD is not supplied as a covariate it randomly chooses a proportion p of
# trials to be stop trials.
#
# If SSD is supplied as a covariate it should contain Inf (for trials with no
# stop signal), and for the other (stop signal) trials either positive reals
# (where SSD is fixed) or NA, in which case values are filled into the dadm by
# the staircase.
#
# Because of the way EMC2 treats design functions this function also needs to
# return a vector when pars not supplied (as is the case when it is run by
# design_model when initially making a dadm).
#
# If make_data is given a design with a staircase function it will run it with
# pars and so make an SSD column, if not present already, or if present it will
# replace any NA values with staircase generated SSDs along with the
# corresponding R and rt values.

staircase <- function(dadm,p=.25,pars=NULL,
                      SSD0=.25,stairstep=.05,stairmin=0,stairmax=Inf)
  # random p of trials get NA, ready to be filled in by a staircase
{
  # if pars not supplied return an SSD column indicating all go trials
  if (is.null(pars)) return(rep(Inf,nrow(dadm)))

  # levels(dadm$lR) <- levels(dadm$lR)

  # if pars supplied potentially run staircase
  nacc <- length(levels(dadm$lR))
  is1 <- dadm$lR==levels(dadm$lR)[1] # first accumulator
  ntrials <- sum(is1)
  if (!any(colnames(dadm)=="SSD")) { # pick staircase trials
    dadm$SSD <- rep(Inf, nrow(dadm))
    tmp <- matrix(dadm$SSD,nrow=nacc)
    tmp[,sample(1:ntrials,round(ntrials*p))] <- NA
    dadm$SSD <- as.vector(tmp)
  }
  out <- setNames(data.frame(matrix(NA,nrow=ntrials,ncol=3)),c("R","rt","SSD"))
  pick <- is.infinite(dadm$SSD)
  if (any(pick)) { # fill in no-stop trials
    out$SSD[pick[is1]] <- Inf
    pars[pick,"SSD"] <- Inf
    out[pick[is1],c("R","rt")] <- attributes(dadm)$model()$rfun(dadm$lR[pick],pars[pick,,drop=FALSE])
  }
  pick <- is.finite(dadm$SSD)
  if (any(pick)) { # fill in fixed SSD trials
    out$SSD[pick[is1]] <- dadm$SSD[pick][is1[pick]]
    pars[pick,"SSD"] <- dadm$SSD[pick]
    out[pick[is1],c("R","rt")] <- attributes(dadm)$model()$rfun(dadm$lR[pick],pars[pick,,drop=FALSE])
  }
  isna <- is.na(dadm$SSD)
  if (any(isna)) { # run staircase if any NAs to fill in
    nstair <- sum(isna)/nacc # number of staircase trials
    trials <- rep(0,ntrials*nacc) # used to pick out each staircase trial in pars
    trials[is.na(dadm$SSD)] <- rep(1:nstair,each=nacc)
    for (i in 1:nstair) { # run staircase
      current <-  trials == i # current staircase trial
      if (i==1)  dadm$SSD[current] <- out$SSD[current[is1]] <- SSD0 # initialize
      p <- pars[current,,drop=FALSE] # parameters for current staircase trial
      # simulate 1 trial, because is.na(SSD) in pars rfun returns dt, an nacc x 1 matrix
      dt <- attributes(dadm)$model()$rfun(dadm$lR[current],p)
      inhibit <- p[,"lI"]==1 # inhibition triggered
      # add SSD to stop and inhibition triggered
      dt[c(TRUE,inhibit),] <- dt[c(TRUE,inhibit),] + dadm$SSD[current][1]
      winner <- which.min(dt)
      iinhibit <- c(1,1+c(1:nacc)[inhibit]) # stop or inhibition triggered index
      if (i != nstair) { # set SSD for next staircase trial
        nexts <- trials == (i+1)
        if (any(iinhibit==winner))   # round as otherwise get spurious tiny differences
          dadm$SSD[nexts] <- round(dadm$SSD[current] + stairstep,3) else  # successful stop
            dadm$SSD[nexts] <- round(dadm$SSD[current] - stairstep,3)       # failed stop
          if ((dadm$SSD[nexts][1]<stairmin) | (dadm$SSD[nexts][1]>stairmax))
            dadm$SSD[nexts] <- dadm$SSD[current] # dont step
          out$SSD[c(trials == i+1)[is1]] <- dadm$SSD[nexts][1]
      }
      if (winner==1) { # stop wins
        if (any(inhibit)) { # inhibition triggered response
          if (all(is.infinite(dt[c(FALSE,inhibit),])))  # tf
            out[current[is1],c("R","rt")] <- c(NA,NA) else {
              pick <- which.min(dt[c(FALSE,inhibit),])
              out[current[is1],c("R","rt")] <-
                c(c(1:nacc)[inhibit][pick],dt[c(FALSE,inhibit),][pick])
            }
        } # otherwise no response
      } else { # pick from all except stop
        out[current[is1],c("R","rt")] <- c(winner-1,dt[winner,])
      }
    }
  }
  out$R <- factor(out$R,levels=1:nacc,labels=levels(dadm$lR))
  out
}

#### Standard stop-signal paradigm (no stop-triggered) ----

# Now no need for a TA factor, usual binary matchfun
GA <- c("g1","g2")
matchfun <- function(d) as.numeric(d$GA)==as.numeric(d$lR)

# Latent inhibit function has to have levels as otherwise only ever TRUE
lI <- function(d) factor(d$lR %in% c("g1","g2"),levels=c(F,T),labels=c("st","go"))

designSSexG0 <- make_design(model=SSexG,
                            factors=list(subjects=1,GA=GA),
                            functions = list(staircase=staircase,lI=lI),
                            Rlevels=c("g1","g2"),matchfun=matchfun,
                            contrasts=list(lM=ADmat),
                            covariates=c("SSD"),
                            formula=list(mu~lM,sigma~1,tau~1,muS~1,sigmaS~1,tauS~1,tf~1,gf~1)
)

# Only Go mu
gf=0; tf=0
p_vector <- sampled_p_vector(designSSexG0,doMap = FALSE)
p_vector[1:9] <- c(log(.5),log(.85),
                   log(.05),log(.1),
                   log(.3),log(.05),log(.1),
                   qnorm(tf),qnorm(gf))

# does not work
#mapped_par(p_vector, design = designSSexG0)

dat <- make_data(p_vector,designSSexG0,n_trials = 10000)
#plot_defective_density(dat,layout=c(2,2))

# GO TRIAL DATA
gdat <- dat[!is.finite(dat$SSD),]
# Probability of go response
mean(gdat$R %in% c("g1","g2"))
# Go accuracy
tapply(as.numeric(gdat$GA)==as.numeric(gdat$R),gdat[,c("GA")],mean, na.rm = T )

# STOP TRIAL DATA
sdat <- dat[is.finite(dat$SSD),]
# SSD
#plot(sdat$SSD,type="l")
# Probability of stopping
mean(is.na(sdat$R))
# signal-respond accuracy, higher than go accuracy as quicker wins more often
tapply(as.numeric(sdat$GA)==as.numeric(sdat$R),sdat[,c("GA")],mean,na.rm=TRUE)


# Profiles good.
# undebug(EMC2:::log_likelihood_race_ss)
# EMC2:::log_likelihood_race_ss(p_vector,dadm)
#EMC2:::log_likelihood_race_ss(p_vector,dadm)

par(mfrow=c(2,5))
for (i in names(p_vector)) if ((!(i %in% c("gf","tf"))) | (i=="gf" & gf>0) | (i=="tf" & tf>0))
  print(profile_plot(i,p_vector,p_vector[i]-.4,p_vector[i]+.4,design = designSSexG0,cores=6, data = dat))


