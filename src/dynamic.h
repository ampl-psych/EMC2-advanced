#ifndef dynamic_h
#define dynamic_h

#include "utility_functions.h"
#include <Rcpp.h>
using namespace Rcpp;

NumericVector dynfuns_c(double base, NumericVector dp, String dyntype, String maptype, DataFrame data, String cnames) {
  NumericVector out(data.nrow());
  NumericVector lR = data["lR"];
  NumericVector tmp = data[cnames];
  bool do_lR = TRUE;
  if(do_lR){
    tmp = tmp[lR == 1];
  }
  NumericVector res(tmp.length());

  // dyntypes
  if(dyntype == "ld"){
    res = -1*tmp;
  }
  if(dyntype == "li"){
    res = tmp;
  }
  if(dyntype == "ed"){
    res = exp(-exp(dp[1])*tmp);
  }
  if(dyntype == "ei"){
    res = 1-exp(-exp(dp[1])*tmp);
  }
  if(dyntype == "pd"){
    res = pow((1+tmp), -exp(dp[1]));
  }
  if(dyntype == "pi"){
    res = 1-pow((1+tmp), -exp(dp[1]));
  }
  if(dyntype == "p2"){
    res = dp[0]*tmp + dp[1]* pow(tmp, 2);
  }
  if(dyntype == "p3"){
    res = dp[0]*tmp + dp[1]* pow(tmp, 2) + dp[2]* pow(tmp, 3);
  }
  if(dyntype == "p4"){
    res = dp[0]*tmp + dp[1]* pow(tmp, 2) + dp[2]* pow(tmp, 3) + dp[3]* pow(tmp, 4);
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
      String cnames = cur_dynamic["covnames"];
      dp = p[dpnames];
      out(_, q) = dynfuns_c(p_curr[q],dp, dyntype, maptype, data, cnames);
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

NumericVector adaptfuns_c(NumericVector base, NumericMatrix dp, String dyntype, String maptype, DataFrame data, String cnames, bool do_lR) {
  NumericVector out(data.nrow());
  NumericVector lR = data["lR"];
  NumericVector tmp = data[cnames];
  // do_lR = TRUE;
  NumericVector base_test = base;
  if(do_lR){
    tmp = tmp[lR == 1];
    dp = submat_rcpp(dp, lR == 1);
    base = base[lR == 1];
  } else{
    tmp = tmp[lR < 99999]; // I'm doing something wrong here but this fixes it????
  }
  NumericVector res(tmp.length());

  // dyntypes
  if(dyntype == "ld"){
    res = -1*tmp;
  }
  if(dyntype == "li"){
    res = tmp;
  }
  if(dyntype == "ed"){
    res = exp(-exp(dp(_, 1))*tmp);
  }
  if(dyntype == "ei"){
    res = 1-exp(-exp(dp(_, 1))*tmp);
  }
  if(dyntype == "pd"){
    res = vector_pow((1+tmp), -exp(dp(_, 1)));
  }
  if(dyntype == "pi"){
    res = 1-vector_pow((1+tmp), -exp(dp(_, 1)));
  }
  if(dyntype == "p2"){
    res = dp(_, 0)*tmp + dp(_, 1)*pow(tmp, 2);
  }
  if(dyntype == "p3"){
    res = dp(_, 0)*tmp + dp(_, 1)*pow(tmp, 2) + dp(_, 2)*pow(tmp, 3);
  }
  if(dyntype == "p4"){
    res = dp(_, 0)*tmp + dp(_, 1)*pow(tmp, 2) + dp(_, 2)*pow(tmp, 3) + dp(_, 3)*pow(tmp, 4);
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
      String cnames = cur_dynamic["covnames"];
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


























