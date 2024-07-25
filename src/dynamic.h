#ifndef dynamic_h
#define dynamic_h

#include "utility_functions.h"
#include <Rcpp.h>
using namespace Rcpp;

// run_delta <- function(q0,p,target,return_extras=FALSE) {
//   if (any(dimnames(target)[[2]]=="winner")) {
//     winner <- as.numeric(target[,1]) # R always first
//     target <- target[,-1,drop=FALSE]
//   } else winner <- rep(1,nrow(target))
//     if (!return_extras) out <- matrix(ncol=ncol(target),nrow=nrow(target)) else
//       out <- array(dim=c(ncol(target),nrow(target),2),
//                    dimnames=list(dimnames(target)[[2]],NULL,c("q","pe")))
//       for (i in 1:ncol(target)) {
//         if (!return_extras)
//           out[,i] <- run_delta_i(q0,p,target[,i],winner,return_extras) else
//             out[i,,] <- run_delta_i(q0,p,target[,i],winner,return_extras)
//       }
//       out
// }

// run_delta_i <- function(q0,alpha,target,winner,return_extras=FALSE) {
//
//   q <- pe <- numeric(length(target))
//   q[1] <- q0
//   alpha <- pnorm(alpha)
//   for (i in 2:length(q)) {
//     if (is.na(target[i-1]) | winner[i-1] != 1) {
//       pe[i-1] <- NA
//       q[i] <- q[i-1]
//     } else {
//       pe[i-1] <- target[i-1]-q[i-1]
//       q[i] <- q[i-1] + alpha[i-1]*pe[i-1]
//     }
//   }
//   if (return_extras) {
//     pe <- target - q
//     cbind(q=q,pe=pe)
//   } else q
// }

NumericVector run_delta_i_dyn(double q0, double alpha, NumericVector target, NumericVector winner){
  NumericVector q(target.length());
  q[0] = q0;
  alpha = R::pnorm(alpha, 0, 1, TRUE, FALSE);
  double pe;
  for(int t = 1; t < q.length(); t ++){
    if(NumericVector::is_na(target[t-1]) || winner[t-1] != 1){
      q[t] = q[t-1];
    } else{
      pe = target[t-1] - q[t-1];
      q[t] = q[t-1] + alpha*pe;
    }
  }
  return(q);
}

NumericVector run_delta_dyn(double q0, double p, NumericMatrix target, bool isRL){
  NumericVector winner;
  NumericVector s_col;
  if(isRL){
    s_col = target(_, 1);
  }
  CharacterVector tmp = colnames(target);
  if(is_true(any(contains(colnames(target), "winner")))){
    winner = target(_, 0);
    if(isRL){
      target = target(_, Range(2, target.ncol() - 1));
    } else{
      target = target(_, Range(1, target.ncol() - 1));
    }
  } else{
    winner.fill(1);
  }
  NumericMatrix out(target.nrow(), target.ncol());
  for(int r = 0; r < target.ncol(); r ++){
    out(_, r) = run_delta_i_dyn(q0, p, target(_, r), winner);
  }
  NumericVector out_v(target.nrow());
  if(isRL){
    for(int v = 0; v < target.nrow(); v ++){
      out_v[v] = out(v, s_col[v]-1);
    }
  } else{
    out_v = as<NumericVector>(out);
  }
  return(out_v);
}

// Here base is the base parameter, and dp are the additional parameters. They are all doubles.
// Dyntype and maptype inform us what transformation to be done.
// Data and cnames, tell us what columns in the data are informing our transformations
NumericVector dynfuns_c(double base, NumericVector dp, String dyntype, String maptype, DataFrame data, CharacterVector cnames,
                        bool do_lR) {
  bool isRL = is_true(any(contains(cnames, "winner")));
  NumericVector out(data.nrow());
  NumericVector lR = data["lR"];
  NumericMatrix covariates(data.nrow(), cnames.length());
  LogicalVector cnames_idx = contains_multiple(data.names(), cnames);
  int k = 0;
  for(int j = 0; j < data.ncol(); j ++){
    if(cnames_idx[j] == TRUE){
      covariates(_, k) = as<NumericVector>(data[j]);
      k++;
    }
  }
  if(do_lR){
    covariates = submat_rcpp(covariates, lR == 1);
  } else{
    covariates = submat_rcpp(covariates, lR < 9999); // I'm doing something wrong here but this fixes it????
  }
  colnames(covariates) = cnames;
  NumericVector res(covariates.nrow());
  //
  // dyntypes
  if(dyntype == "ld"){
    res = -1* as<NumericVector>(covariates);
  }
  if(dyntype == "li"){
    res = as<NumericVector>(covariates);
  }
  if(dyntype == "ed"){
    res = exp(-exp(dp[1])*as<NumericVector>(covariates));
  }
  if(dyntype == "ei"){
    res = 1-exp(-exp(dp[1])*as<NumericVector>(covariates));
  }
  if(dyntype == "pd"){
    res = pow((1+as<NumericVector>(covariates)), -exp(dp[1]));
  }
  if(dyntype == "pi"){
    res = 1-pow((1+as<NumericVector>(covariates)), -exp(dp[1]));
  }
  if(dyntype == "p2"){
    res = dp[0]*as<NumericVector>(covariates) + dp[1]* pow(as<NumericVector>(covariates), 2);
  }
  if(dyntype == "p3"){
    res = dp[0]*as<NumericVector>(covariates) + dp[1]* pow(as<NumericVector>(covariates), 2) + dp[2]* pow(as<NumericVector>(covariates), 3);
  }
  if(dyntype == "p4"){
    res = dp[0]*as<NumericVector>(covariates) + dp[1]* pow(as<NumericVector>(covariates), 2) + dp[2]* pow(as<NumericVector>(covariates), 3) + dp[3]* pow(as<NumericVector>(covariates), 4);
  }
  if(dyntype == "d1"){
    res = run_delta_dyn(dp[1], dp[2], covariates, isRL);
  }
  // Map types
  if(maptype == "lin"){
    res = base + dp[0]*res;
  }
  if(maptype == "plin"){
    res = exp(base) + exp(dp[0])*res;
  }
  if(maptype == "add"){
    res = base + res;
  }
  if(maptype == "ucent"){
    res = base + dp[0]*(res - .5);
  }
  if(do_lR){
    for(int i = 0; i < max(lR); i ++){
      out[lR == (i + 1)] = res;
    }
  } else{
    out = res;
  }
  return(out);
}

NumericMatrix map_dyn(List dynamic, DataFrame data, NumericVector p, CharacterVector curr_names,
                      LogicalVector isin) {
  NumericVector dp;
  List cur_dynamic;
  CharacterVector dyn_names = dynamic.names();
  NumericVector p_curr = p[curr_names];
  NumericMatrix out(data.nrow(), curr_names.length());
  NumericVector input(data.nrow());
  for(int q = 0; q < isin.length(); q ++){
    if(isin[q] == TRUE){
      String curr_name = curr_names[q];
      cur_dynamic = dynamic[curr_name];
      String dyntype = cur_dynamic["dyntype"];
      String maptype = cur_dynamic["maptype"];
      CharacterVector dpnames = cur_dynamic["dpnames"];
      CharacterVector cnames = cur_dynamic["covnames"];
      bool do_lR = cur_dynamic["lR1"];
      dp = p[dpnames];
      out(_, q) = dynfuns_c(p_curr[q],dp, dyntype, maptype, data, cnames, do_lR);
    } else{
      input.fill(p_curr[q]);
      out(_, q) = input;
    }
  }
  return(out);
}


NumericVector vector_pow(NumericVector x1, NumericVector x2){
  NumericVector out(x1.length());
  for(unsigned int i = 0; i < out.length(); i ++){
    out[i] = pow(x1[i], x2[i]);
  }
  return(out);
}


// Here base is the base parameter, and dp are the additional parameters. They are all vectors.
// Because here they have already been mapped back to the design (so one parameter per observation * number of accumulators)
// Dyntype and maptype inform us what transformation to be done.
// Data and cnames, tell us what columns in the data are informing our transformations
NumericVector adaptfuns_c(NumericVector base, NumericMatrix dp, String dyntype, String maptype, DataFrame data, CharacterVector cnames, bool do_lR) {
  NumericVector out(data.nrow());
  NumericVector lR = data["lR"];
  // NumericVector covariates = data[cnames];
  //
  //
  NumericMatrix covariates(data.nrow(), cnames.length());
  LogicalVector cnames_idx = contains_multiple(data.names(), cnames);
  int k = 0;
  for(int j = 0; j < data.ncol(); j ++){
    if(cnames_idx[j] == TRUE){
      covariates(_, k) = as<NumericVector>(data[j]);
      k++;
    }
  }
  if(do_lR){
    covariates = submat_rcpp(covariates, lR == 1);
    dp = submat_rcpp(dp, lR == 1);
  } else{
    covariates = submat_rcpp(covariates, lR < 9999); // I'm doing something wrong here but this fixes it????
  }
  NumericVector res(covariates.nrow());

  // dyntypes
  if(dyntype == "ld"){
    res = -1*as<NumericVector>(covariates);
  }
  if(dyntype == "li"){
    res = as<NumericVector>(covariates);
  }
  if(dyntype == "ed"){
    res = exp(-exp(dp(_, 1))*as<NumericVector>(covariates));
  }
  if(dyntype == "ei"){
    res = 1-exp(-exp(dp(_, 1))*as<NumericVector>(covariates));
  }
  if(dyntype == "pd"){
    res = vector_pow((1+as<NumericVector>(covariates)), -exp(dp(_, 1)));
  }
  if(dyntype == "pi"){
    res = 1-vector_pow((1+as<NumericVector>(covariates)), -exp(dp(_, 1)));
  }
  if(dyntype == "p2"){
    res = dp(_, 0)*as<NumericVector>(covariates) + dp(_, 1)*pow(as<NumericVector>(covariates), 2);
  }
  if(dyntype == "p3"){
    res = dp(_, 0)*as<NumericVector>(covariates) + dp(_, 1)*pow(as<NumericVector>(covariates), 2) + dp(_, 2)*pow(as<NumericVector>(covariates), 3);
  }
  if(dyntype == "p4"){
    res = dp(_, 0)*as<NumericVector>(covariates) + dp(_, 1)*pow(as<NumericVector>(covariates), 2) + dp(_, 2)*pow(as<NumericVector>(covariates), 3) + dp(_, 3)*pow(as<NumericVector>(covariates), 4);
  }
  // Map types
  if(maptype == "lin"){
    res = base + dp(_, 0)*res;
  }
  if(maptype == "plin"){
    res = exp(base) + exp(dp(_, 0))*res;
  }
  if(maptype == "add"){
    res = base + res;
  }
  if(maptype == "ucent"){
    res = base + dp(_, 0)*(res - .5);
  }
  if(do_lR){
    for(int i = 0; i < max(lR); i ++){
      out[lR == (i + 1)] = res;
    }
  } else{
    out = res;
  }

  return(out);
}



NumericMatrix map_adaptive(List adaptive, NumericMatrix p, DataFrame data) {
  List cur_dynamic;
  CharacterVector adapt_names = adaptive.names();
  CharacterVector curr_names = colnames(p);
  LogicalVector isin = contains_multiple(curr_names, adapt_names);
  for(int q = 0; q < isin.length(); q ++){
    if(isin[q] == TRUE){
      String curr_name = curr_names[q];
      cur_dynamic = adaptive[curr_name];
      String dyntype = cur_dynamic["dyntype"];
      String maptype = cur_dynamic["maptype"];
      CharacterVector cnames = cur_dynamic["covnames"];
      CharacterVector aptypes = cur_dynamic["aptypes"];
      LogicalVector is_adapt = contains_multiple(curr_names, aptypes);
      NumericMatrix dp(p.nrow(), sum(is_adapt));
      int z = 0;
      for(int k = 0; k < isin.length();k++){
        if(is_adapt[k] == TRUE){
          dp(_, z) = p(_,k);
          z++;
        }
      }
      bool do_lR = cur_dynamic["lR1"];
      p(_, q) = adaptfuns_c(p(_, q),dp, dyntype, maptype, data, cnames, do_lR);
    }
  }
  return(p);
}


#endif
