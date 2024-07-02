#ifndef lnr_h
#define lnr_h

#include <Rcpp.h>
#include "utility_functions.h"
#include "dynamic.h"
#include "advantage.h"

using namespace Rcpp;


NumericVector transform_lnr(NumericVector x){
  return(x);
}

NumericMatrix Ntransform_lnr(NumericMatrix x, CharacterVector use, DataFrame data, List adaptive) {
  NumericMatrix out(clone(x));
  LogicalVector col_idx = contains(colnames(x), "m");
  LogicalVector use_idx = contains_multiple(colnames(x), use);

  for(int i = 0; i < x.ncol(); i ++){
    if(use_idx[i] == TRUE){
      if(col_idx[i] == FALSE){
        out (_, i) = exp(out(_, i));
      };
    }
  };
  if(adaptive.length() > 0){
    out = map_adaptive(adaptive, out, data);
  }
  List advantage = data.attr("advantage");
  if(advantage.length() > 0){
    String par_name = advantage.names();
    String S_name = advantage[0];
    NumericVector SV = data[S_name];
    CharacterVector R = data["R"];
    out = advantagepars(out, SV,unique(R).length(),par_name);
  }
  return(out);
}

NumericVector plnr_c(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll){
  // 0 = m, 1 = s, 2 = t0
  int n = sum(idx);
  NumericVector out(n);
  int k = 0;
  for(int i = 0; i < rts.length(); i++){
    if(idx[i] == TRUE){
      if(!NumericVector::is_na(pars(i,0)) & (rts[i] - pars(i,2) > 0) & (pars(i, 1) > 0) &  (pars(i, 2) > 0.05)){
        out[k] = R::plnorm(rts[i] - pars(i,2), pars(i, 0), pars(i, 1), TRUE, FALSE);
      } else{
        out[k] = min_ll;
      }
      k++;
    }
  }

  return(out);
}

NumericVector dlnr_c(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll){
  int n = sum(idx);
  NumericVector out(n);
  int k = 0;
  for(int i = 0; i < rts.length(); i++){
    if(idx[i] == TRUE){
      if(!NumericVector::is_na(pars(i,0)) & (rts[i] - pars(i,2) > 0) & (pars(i, 1) > 0) &  (pars(i, 2) > 0.05)){
        out[k] = R::dlnorm(rts[i] - pars(i,2), pars(i, 0), pars(i, 1), FALSE);
      } else{
        out[k] = min_ll;
      }
      k++;
    }

  }

  return(out);
}


#endif
