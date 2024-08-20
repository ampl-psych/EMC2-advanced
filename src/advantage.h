#ifndef advantage_h
#define advantage_h

#include "utility_functions.h"
#include <Rcpp.h>
using namespace Rcpp;

// advantage_pars <- function(pars,SV,nr,pname="v")  {
//   na <- nr-1                     # accumulators within responses
//   nt <- nrow(pars)/nr            # trials
// # Beta stimulus value mapping
//   SV <- pbeta(SV[[1]],exp((pars[,"SS"]+pars[,"SD"])/2),
//                      exp((pars[,"SS"]-pars[,"SD"])/2))
//     q <- matrix(SV,nrow=nr)
// # sum and difference arrays, dimensions so as.vector fits with rpars
//     s <- d <- array(dim=c(na,nr,ncol(q)))
//     nrs <- 1:nr
//   for (i in 1:na) {
//       for (j in 1:nr) {
//         notj <- nrs[-j]
// # Full looping algorithm in C would be over 1:nt
//         d[i,j,] <- apply(q[c(notj[i],j),],2,diff)
//         s[i,j,] <- apply(q[c(notj[i],j),],2,sum)
//       }
//     }
// # expand pars
//   rpars <- apply(pars[,dimnames(pars)[[2]]!=pname],2,rep,each=na)
// # add in advantage parameter (usually v)
//     v <- rep(pars[,pname],each=na) + rep(pars[,"AS"],each=na)*as.vector(s) +
//       rep(pars[,"AD"],each=na)*as.vector(d)
//     rpars <- cbind(v,rpars)
//     dimnames(rpars)[[2]][1] <- pname
//   rpars
// }


// [[Rcpp::export]]
NumericMatrix advantagepars(NumericMatrix pars, NumericVector SV, int nr, String pname){
  int na = nr -1;
  int nt = pars.nrow()/nr;
  // Contains is just a custom ''which'' function, which returns a logical vector
  // with only true at the match.
  LogicalVector is_SS = contains(colnames(pars), "SS");
  LogicalVector is_SD = contains(colnames(pars), "SD");
  LogicalVector is_AS = contains(colnames(pars), "AS");
  LogicalVector is_AD = contains(colnames(pars), "AD");
  LogicalVector is_DD = contains(colnames(pars), "DD");
  LogicalVector is_DS = contains(colnames(pars), "DS");
  // I have yet to find a better way to go through a matrix and get matching columns
  // Works better with dataframes that are treated as lists.
  // But this way is definitely fastests (since it's only one explicit loop,
  // instead of four one for each search). Just more code.
  NumericVector SS(pars.nrow());
  NumericVector SD(pars.nrow());
  NumericVector AS(pars.nrow());
  NumericVector AD(pars.nrow());
  NumericVector DD(pars.nrow());
  NumericVector DS(pars.nrow());
  for(int i = 0; i < pars.ncol(); i ++){
    if(is_SS[i] == TRUE){
      SS = pars(_, i);
    }
    if(is_SD[i] == TRUE){
      SD = pars(_, i);
    }
    if(is_AS[i] == TRUE){
      AS = pars(_, i);
    }
    if(is_AD[i] == TRUE){
      AD = pars(_, i);
    }
    if(is_DD[i] == TRUE){
      DD = pars(_, i);
    }
    if(is_DS[i] == TRUE){
      DS = pars(_, i);
    }
  }
  // R:: has a ton of distributional functions
  // Unlike base R they do require single values instead of vectors.
  for(int j = 0; j < SV.length(); j ++){
    SV[j] = R::pbeta(SV[j], exp((SS[j]+SD[j])/2), exp((SS[j]-SD[j])/2), true, false);
  }
  // Set up our empty matrices
  NumericMatrix q(nr, nt, SV.begin());
  // These two are now matrices instead of arrays, with the first dimension
  // equal to the first*second dimension of the array.
  NumericMatrix s_mat(na*nr, nt);
  NumericMatrix d_mat(na*nr, nt);
  IntegerVector notm;
  IntegerVector nrs = seq(0, nr-1);
  for (int k = 0; k < nt; k++) {
    for(int l = 0; l < na; l++){
      for(int m = 0; m < nr; m++){
        notm = nrs[nrs != m];
        d_mat(m * na + l, k) =  q(m, k) - q(notm[l], k);
        s_mat(m * na + l, k) = q(notm[l], k) + q(m, k);
      }
    }
  }

  NumericMatrix pars_out(pars.nrow()*na, pars.ncol());
  colnames(pars_out) = colnames(pars);

  // Find out for which parameter I need to do the advantage transformation
  LogicalVector is_ADpar = contains(colnames(pars), pname);
  NumericVector pars_curr(pars.nrow());
  NumericVector d = as<NumericVector>(d_mat);
  NumericVector DS2 = rep_each(DS, na);
  NumericVector DD2 = rep_each(DD, na);

  for(int j = 0; j < DS2.length(); j ++){
    if(sum(is_DS) > 0){
      d[j] = R::pbeta((d[j] + 1)/2, exp((DS2[j]+DD2[j])/2), exp((DS2[j]-DD2[j])/2), true, false);
    }else{
      d[j] = (d[j] + 1)/2;
    }
  }
  NumericVector s = as<NumericVector>(s_mat);
  for(int n = 0; n < pars.ncol(); n++){
    pars_curr = pars(_, n);
    pars_out(_,n) = rep_each(pars_curr, na);
    if(is_ADpar[n] == TRUE){ // If it's the advantage parameter also do:
      // NumericVector tmp = pars_out(_,n);
      pars_out(_,n) = pars_out(_,n) + rep_each(AS, na) * s + rep_each(AD, na) * (2*d -1);
      // NumericVector tmp2 = tmp + rep_each(AS, na) * s + rep_each(AD, na) * (2*d -1);
      // Rcout << tmp2;
    }
  }
  return(pars_out);
}


#endif

// advantage_pars_ddm <- function(pars,SV)  {
// # Beta stimulus value mapping
//   SV1 <- pbeta(SV[[1]],exp((pars[,"SS"]+pars[,"SD"])/2),
//                exp((pars[,"SS"]-pars[,"SD"])/2))
//   SV2 <- pbeta(SV[[2]],exp((pars[,"SS"]+pars[,"SD"])/2),
//                exp((pars[,"SS"]-pars[,"SD"])/2))
// # Beta difference mapping
//   if (!all(c("DS","DD") %in% dimnames(pars)[[2]])) d <- (1+SV2-SV1)/2 else
//     d <- pbeta((1+SV2-SV1)/2,exp((pars[,"DS"]+pars[,"DD"])/2),
//                exp((pars[,"DS"]-pars[,"DD"])/2))
//     pars[,"v"] <- pars[,"v"] + pars[,"AD"]*(2*d-1)
//     pars[,"sv"] <- pars[,"sv"] + pars[,"ssv"]*(SV2^2+SV1^2)^(1/2)
//     pars
// }
//

// NumericMatrix advantagepars(NumericMatrix pars, NumericVector SV1, NumericVector SV2){
//   // Contains is just a custom ''which'' function, which returns a logical vector
//   // with only true at the match.
//   LogicalVector is_SS = contains(colnames(pars), "SS");
//   LogicalVector is_SD = contains(colnames(pars), "SD");
//   LogicalVector is_AS = contains(colnames(pars), "AS");
//   LogicalVector is_AD = contains(colnames(pars), "AD");
//   LogicalVector is_DD = contains(colnames(pars), "DD");
//   LogicalVector is_DS = contains(colnames(pars), "DS");
//   // I have yet to find a better way to go through a matrix and get matching columns
//   // Works better with dataframes that are treated as lists.
//   // But this way is definitely fastests (since it's only one explicit loop,
//   // instead of four one for each search). Just more code.
//   NumericVector SS(pars.nrow());
//   NumericVector SD(pars.nrow());
//   NumericVector AS(pars.nrow());
//   NumericVector AD(pars.nrow());
//   NumericVector DD(pars.nrow());
//   NumericVector DS(pars.nrow());
//   for(int i = 0; i < pars.ncol(); i ++){
//     if(is_SS[i] == TRUE){
//       SS = pars(_, i);
//     }
//     if(is_SD[i] == TRUE){
//       SD = pars(_, i);
//     }
//     if(is_AS[i] == TRUE){
//       AS = pars(_, i);
//     }
//     if(is_AD[i] == TRUE){
//       AD = pars(_, i);
//     }
//     if(is_DD[i] == TRUE){
//       DD = pars(_, i);
//     }
//     if(is_DS[i] == TRUE){
//       DS = pars(_, i);
//     }
//   }
//   // R:: has a ton of distributional functions
//   // Unlike base R they do require single values instead of vectors.
//   for(int j = 0; j < SV1.length(); j ++){
//     SV1[j] = R::pbeta(SV1[j], exp((SS[j]+SD[j])/2), exp(SS[j]-SD[j]/2), true, false);
//     SV2[j] = R::pbeta(SV2[j], exp((SS[j]+SD[j])/2), exp(SS[j]-SD[j]/2), true, false);
//   }
//
//   // Set up our empty matrices
//   NumericMatrix q(nr, nt, SV.begin());
//   // These two are now matrices instead of arrays, with the first dimension
//   // equal to the first*second dimension of the array.
//   NumericMatrix s(na*nr, nt);
//   NumericMatrix d(na*nr, nt);
//   IntegerVector notm;
//   IntegerVector nrs = seq(0, nr-1);
//   for (int k = 0; k < nt; k++) {
//     for(int l = 0; l < na; l++){
//       for(int m = 0; m < nr; m++){
//         notm = nrs[nrs != m];
//         d(m * na + l, k) =  q(m, k) - q(notm[l], k);
//         s(m * na + l, k) = q(notm[l], k) + q(m, k);
//       }
//     }
//   }
//   NumericMatrix pars_out(pars.nrow()*na, pars.ncol());
//   colnames(pars_out) = colnames(pars);
//
//   // Find out for which parameter I need to do the advantage transformation
//   LogicalVector is_ADpar = contains(colnames(pars), pname);
//   NumericVector pars_curr(pars.nrow());
//
//   // Convert our matrices to NumericVector (same as as.vector() in R)
//   d.attr("dim") = Dimension(d.nrow() * d.ncol(), 1);
//   d = (1+d)/2;
//   if(sum(is_DS) > 0){
//     for(int j = 0; j < DS.length(); j ++){
//       d[j] = R::pbeta(d[j], exp((DS[j]+DD[j])/2), exp(DS[j]-DD[j]/2), true, false);
//     }
//   }
//
//   s.attr("dim") = Dimension(s.nrow() * s.ncol(), 1);
//   for(int n = 0; n < pars.ncol(); n++){
//     pars_curr = pars(_, n);
//     pars_out(_,n) = rep_each(pars_curr, na);
//     if(is_ADpar[n] == TRUE){ // If it's the advantage parameter also do:
//       pars_out(_,n) = pars_out(_,n) + rep_each(AS, na) * s + rep_each(AD, na) * d;
//     }
//   }
//   return(pars_out);
// }


