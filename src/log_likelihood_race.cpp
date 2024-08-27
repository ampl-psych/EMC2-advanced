#include <Rcpp.h>
#include "utility_functions.h"
#include "dynamic.h"
#include "model_lnr.h"
#include "model_LBA.h"
#include "model_RDM.h"
#include "model_DDM.h"
#include "advantage.h"
#include "numeric_integration.h"
using namespace Rcpp;


NumericVector c_expand(NumericVector x1, NumericVector expand){
  const int n_out = expand.length();
  NumericVector out(n_out);
  int curr_idx;
  for(int i = 0; i < n_out; i++){
    curr_idx = expand[i] - 1; //expand created in 1-based R
    out[i] = x1[curr_idx];
  }
  return(out);
}

LogicalVector c_bool_expand(LogicalVector x1, NumericVector expand){
  const int n_out = expand.length();
  LogicalVector out(n_out);
  int curr_idx;
  for(int i = 0; i < n_out; i++){
    curr_idx = expand[i] - 1; //expand created in 1-based R
    out[i] = x1[curr_idx];
  }
  return(out);
}


CharacterVector c_add_charvectors(CharacterVector x, CharacterVector y) {
  CharacterVector z(x.size() + y.size());
  std::copy(x.begin(), x.end(), z.begin());
  std::copy(y.begin(), y.end(), z.begin() + x.size());
  return(z);
}

NumericVector c_add_vectors(NumericVector x1, NumericVector x2){
  if(is_na(x2)[0] ){
    return(x1);
  }
  NumericVector output(x1.size() + x2.size());
  std::copy(x1.begin(), x1.end(), output.begin());
  std::copy(x2.begin(), x2.end(), output.begin() + x1.size());
  CharacterVector all_names(x1.size() + x2.size());
  CharacterVector x1_names = x1.names();
  CharacterVector x2_names = x2.names();
  std::copy(x1_names.begin(), x1_names.end(), all_names.begin());
  std::copy(x2_names.begin(), x2_names.end(), all_names.begin() + x1.size());
  output.names() = all_names;
  return output;
}


// LL generic functions
// [[Rcpp::export]]
NumericMatrix c_map_p(NumericVector p_vector, CharacterVector p_types, List designs, int n_trials, List dynamic,
                      DataFrame data){
  NumericMatrix pars(n_trials, p_types.length());
  NumericVector p_mult_design;
  for(int i = 0; i < p_types.length(); i++){

    NumericMatrix curr_design = designs[i];

    CharacterVector curr_names = colnames(curr_design);
    if(dynamic.length() > 0){
      LogicalVector isin = contains_multiple(curr_names, dynamic.names());
      if(sum(isin) > 0){
        NumericMatrix p_mat = map_dyn(dynamic, data, p_vector, curr_names, isin);
        for(int k = 0; k < curr_design.ncol(); k ++){
          pars(_, i) = pars(_, i) + p_mat(_, k) * curr_design(_, k);
        };
      } else{
        for(int j = 0; j < curr_design.ncol(); j ++){
          String curr_name(curr_names[j]);
          pars(_, i) = pars(_, i) + p_vector[curr_name] * curr_design(_, j);
        };
      }
    } else{
      for(int j = 0; j < curr_design.ncol(); j ++){
        String curr_name(curr_names[j]);
        p_mult_design =  p_vector[curr_name] * curr_design(_, j);
        p_mult_design[is_nan(p_mult_design)] = 0;
        pars(_, i) = pars(_, i) + p_mult_design;
      };

    }
  };
  colnames(pars) = p_types;
  return(pars);
}

NumericMatrix get_pars_matrix(NumericVector p_vector, NumericVector constants,
                       NumericVector (*transform)(NumericVector),
                       NumericMatrix (*Ntransform)(NumericMatrix, CharacterVector, DataFrame, List),
                       CharacterVector p_types, List designs, int n_trials,
                       List dynamic, DataFrame data){
  NumericVector p_vector_updtd(clone(p_vector));
  p_vector_updtd = c_add_vectors(p_vector_updtd, constants);
  p_vector_updtd = transform(p_vector_updtd);
  List adaptive = data.attr("adaptive");
  if(adaptive.length() > 0){
    p_types = c_add_charvectors(p_types, adaptive.attr("aptypes"));
  }
  NumericMatrix pars = c_map_p(p_vector_updtd, p_types, designs, n_trials, dynamic, data);
  CharacterVector use = data.attr("transform_names");
  pars = Ntransform(pars, use, data, adaptive);
  return(pars);
}

double c_log_likelihood_DDM(NumericMatrix pars, DataFrame data,
                          const int n_trials, NumericVector expand,
                          double min_ll, List group_idx){
  const int n_out = expand.length();
  NumericVector rts = data["rt"];
  NumericVector R = data["R"];
  NumericVector lls(n_trials);
  NumericVector lls_exp(n_out);
  lls = log(d_DDM_c(rts, R, group_idx, pars));
  lls_exp = c_expand(lls, expand); // decompress
  lls_exp[is_na(lls_exp)] = min_ll;
  lls_exp[is_infinite(lls_exp)] = min_ll;
  lls_exp[lls_exp < min_ll] = min_ll;

  return(sum(lls_exp));
}

double c_log_likelihood_race(NumericMatrix pars, DataFrame data,
                             NumericVector (*dfun)(NumericVector, NumericMatrix, LogicalVector, double),
                             NumericVector (*pfun)(NumericVector, NumericMatrix, LogicalVector, double),
                             const int n_trials, LogicalVector winner, NumericVector expand,
                             double min_ll){
  const int n_out = expand.length();
  NumericVector lds(n_trials);
  NumericVector rts = data["rt"];
  CharacterVector R = data["R"];
  NumericVector lds_exp(n_out);
  const int n_acc = unique(R).length();
  if(sum(contains(data.names(), "RACE")) == 1){
    NumericVector lR = data["lR"];
    NumericVector NACC = data["RACE"];
    CharacterVector vals_NACC = NACC.attr("levels");
    for(int x = 0; x < pars.nrow(); x++){
      if(lR[x] > atoi(vals_NACC[NACC[x]-1])){
        pars(x,0) = NA_REAL;
      }
    }
  }
  NumericVector win = log(dfun(rts, pars, winner, exp(min_ll))); //first for compressed
  lds[winner] = win;
  if(n_acc > 1){
    NumericVector loss = log(1- pfun(rts, pars, !winner, exp(min_ll))); //cdfs
    loss[is_na(loss)] = min_ll;
    loss[loss == log(1 - exp(min_ll))] = min_ll;
    lds[!winner] = loss;
  }
  lds[is_na(lds)] = min_ll;
  lds_exp = c_expand(lds, expand); // decompress
  if(n_acc > 1){
    LogicalVector winner_exp = c_bool_expand(winner, expand);
    NumericVector ll_out = lds_exp[winner_exp];
    if(n_acc == 2){
      NumericVector lds_los = lds_exp[!winner_exp];
      ll_out = ll_out + lds_los;
    } else{
      NumericVector lds_los = lds_exp[!winner_exp];
      for(int z = 0; z < ll_out.length(); z++){
        ll_out[z] = ll_out[z] + sum(lds_los[seq( z * (n_acc -1), (z+1) * (n_acc -1) -1)]);
      }
    }
    ll_out[is_na(ll_out)] = min_ll;
    ll_out[is_infinite(ll_out)] = min_ll;
    ll_out[ll_out < min_ll] = min_ll;
    return(sum(ll_out));
  } else{
    lds_exp[is_na(lds_exp)] = min_ll;
    lds_exp[is_infinite(lds_exp)] = min_ll;
    lds_exp[lds_exp < min_ll] = min_ll;
    return(sum(lds_exp));
  }
}

NumericVector convert_weights(NumericVector weights, int nr){
  NumericVector out(nr);
  for(int k = 0; k < weights.length(); k ++){
    if(k != 0){
      out[k] = R::pnorm(weights[k], 0, 1, TRUE, FALSE)/(nr-1);
    }
  }
  out[0] = 1-sum(out);
  return(out);
}

double c_log_likelihood_race_trdm(NumericMatrix pars, DataFrame data,
                             const int n_trials, LogicalVector winner, NumericVector expand,
                             double min_ll){
  expand = data.attr("expand_winner"); //overwrite expand
  CharacterVector R = data["R"];
  NumericVector rt = data["rt"];
  const int nr = unique(R).length();
  NumericVector pGuess_data = pars(_, pars.ncol() - 1);
  NumericVector pGuess(sum(winner));
  if(all(is_infinite(pGuess_data))){
    pGuess.fill(nr);
    // so annoying I can't fill it with 1/x...
    pGuess = 1/pGuess;
  } else{
    NumericVector weights(nr);
    int q = 0;
    for(int i = 0; i < (pGuess_data.length()/nr); i ++){
      for(int j = 0; j < nr; j ++){
        weights[j] = pGuess_data[i*nr + j];
      }
      weights = convert_weights(weights, nr);
      for(int j = 0; j < nr; j ++){
        if(winner[i*nr + j] == TRUE){
          pGuess[q] = weights[j]; // only p guess for winner is used
          // and pGuess is only length winner
          q ++;
        }
      }
    }
  }
  NumericVector ll(pars.nrow());
  LogicalVector winner_all(pars.nrow(), true);
  NumericMatrix pars_T = pars(_, Range(5,9));
  NumericMatrix pars_E = pars(_, Range(0,4));
  colnames(pars_T) = colnames(pars_E);
  NumericVector pdfTW = pGuess * drdm_c(rt, pars_T, winner, exp(min_ll));
  NumericVector pdfEW = drdm_c(rt, pars_E, winner, exp(min_ll));
  NumericVector surviveT = 1 - prdm_c(rt, pars_T, winner, exp(min_ll));
  NumericVector surviveE = 1 - prdm_c(rt, pars_E, winner_all, exp(min_ll));
  ll = log(pdfEW*surviveT + pdfTW * as<NumericVector>(surviveE[winner]));
  if(nr > 1){
    if(nr == 2){
      ll = ll + log(as<NumericVector>(surviveE[!winner]));
    } else{
      NumericVector surviveE_nw = surviveE[!winner];
      NumericMatrix surviveE_nw_m(nr-1, surviveE_nw.length()/(nr-1), surviveE_nw.begin());
      NumericVector prodsE (surviveE_nw_m.ncol());
      prodsE.fill(1);
      for(int i = 0; i < surviveE_nw_m.nrow(); i ++){
        prodsE = prodsE * surviveE_nw_m(i,_);
      }
      ll = ll + log(prodsE);
    }
  }
  ll = c_expand(ll, expand);
  ll[is_na(ll)] = min_ll;
  ll[is_infinite(ll)] = min_ll;
  ll[ll < min_ll] = min_ll;
  return(sum(ll));
}

double c_log_likelihood_race_advantage(NumericMatrix pars, DataFrame data,
                                       NumericVector (*dfun)(NumericVector, NumericMatrix, LogicalVector, double),
                                       NumericVector (*pfun)(NumericVector, NumericMatrix, LogicalVector, double),
                                       const int n_trials, LogicalVector winner, NumericVector expand,
                                       double min_ll){
  expand = data.attr("expand_winner"); //overwrite expand
  CharacterVector R = data["R"];
  NumericVector rt = data["rt"];
  const int nr = unique(R).length();
  const int na = nr-1;
  const int nt = data.nrow()/nr;
  winner = rep_each(winner, na);
  rt = rep_each(rt, na);

  // This is a way to fill a matrix with a numericVector
  // Similar to as.matrix or matrix() in R.
  NumericVector pdf = dfun(rt, pars, winner, exp(min_ll));
  NumericMatrix pdf_m(na, nt, pdf.begin()); //first for compressed
  NumericMatrix cdfw(na, nt); // Need to initialize this even if I don't need it or it won't compile
  if(na > 1){
    NumericVector cdfw_vec = pfun(rt, pars, winner, exp(min_ll));
    cdfw = NumericMatrix(na, nt, cdfw_vec.begin());
  }
  // Initialize our matrix (this is an array in R)
  // first dimension is equal to first*second dimension of the R array.
  NumericMatrix cdfl(na*na, nt, pfun(rt, pars, !winner, exp(min_ll)).begin());
  // Rcout << cdfl;
  NumericVector ll(nt);
  if(na == 1){ // Rcpp is a lot less fussy about what is a matrix and what is a vector
    // Even though I initialized it as a matrix and log expects a vector
    // A matrix with a dimension of 1 on one side is just treated as a vector.
    ll = log(pdf) + log(1 - cdfl);
  } else{
    // This takes care of the apply loop: ll <- ll + pdf[i,]*apply(cdfw[ifirst,,drop=FALSE],2,sum)
    // Here I explicitly loop since apply doesn't (really) exist in Rcpp
    NumericVector cdfw_sum(nt);
    for(int i = 0; i < na; i++){
      cdfw_sum.fill(0); // clear our sums to be zero
      for(int j = 0; j < na; j ++){
        if(j != i){
          cdfw_sum = cdfw_sum + cdfw(j,_); // Add current row to existing sum
        }
      }
      ll = ll + pdf_m(i,_) * cdfw_sum; // at this point cdfw_sum == apply(cdfw[ifirst,,drop=FALSE],2,sum)
    }
    ll = log(ll);
    // This is the more tricky part
    // Initialize our empty product matrix
    // cdfl is now an array.
    // This code in Rcpp is actually a triple loop:
    // Apply(cdfl, 2:3, prod)
    // For the second and third dimension take the product across the first.
    NumericMatrix prod_out(na, nt);
    for(int k = 0; k < na; k++){ // Loop across the first
      for(int l = 0; l < nt; l ++){ // and second dimension
        prod_out(k, l) = 1;
        for(int m = 0; m < na*na; m ++){
          // this makes sure we only take elements of the column that were interested in.
          // This does the product loop (just done through prod in R)
          if(m >= k*na && m < (k+1)*na){
            prod_out(k, l) = prod_out(k, l) * cdfl(m, l);
          }
        }

      }
    }
    // apply(log(1-apply(cdfl, 2:3, prod)),2,sum)
    // Now we need to do the apply((log(1-x)), 2, sum) part
    // Where x is prod_out that we already calculated
    NumericVector log_sums_prod_out(nt);
    for(int n = 0; n < na; n ++){
      log_sums_prod_out = log_sums_prod_out + log(1-prod_out(n,_));
    }
    ll = ll + log_sums_prod_out;
  }
  ll = c_expand(ll, expand);
  ll[is_na(ll)] = min_ll;
  ll[is_infinite(ll)] = min_ll;
  ll[ll < min_ll] = min_ll;
  return(sum(ll));
}

double c_log_likelihood_race_missing(NumericMatrix pars, DataFrame data,
                                     NumericVector (*dfun)(NumericVector, NumericMatrix, LogicalVector, double),
                                     NumericVector (*pfun)(NumericVector, NumericMatrix, LogicalVector, double),
                                     const int n_trials, LogicalVector winner, NumericVector expand,
                                     double min_ll){
  const int n_out = expand.length();
  NumericVector lds(n_trials);
  NumericVector rts = data["rt"];
  CharacterVector R = data["R"];
  NumericVector pCont = pars(_ , pars.ncol() - 1);
  NumericVector lds_exp(n_out);
  int n_acc = unique(R).length();
  if (any(is_na(R))) {
    n_acc -= 1;
  }
  if(sum(contains(data.names(), "NACC")) == 1){
    NumericVector lR = data["lR"];
    NumericVector NACC = data["NACC"];
    for(int x = 0; x < pars.nrow(); x++){
      if(lR[x] > NACC[x]){
        pars(x,0) = NA_REAL;
      }
    }
  }
  NumericVector win = log(dfun(rts, pars, winner, exp(min_ll))); //first for compressed
  lds[winner] = win;
  if(n_acc > 1){
    NumericVector loss = log(1- pfun(rts, pars, !winner, exp(min_ll))); //cdfs
    loss[is_na(loss)] = min_ll;
    loss[loss == log(1 - exp(min_ll))] = min_ll;
    lds[!winner] = loss;
  }
  lds[is_na(lds) | (winner & is_infinite(rts))] = min_ll;
  lds[(!winner) & (is_infinite(rts) | is_na(rts))] = 0;


  // Calculate truncation?
  double LT = 0;
  double UT = R_PosInf;

  Nullable<double> LT_ = data.attr("LT");
  Nullable<double> UT_ = data.attr("UT");
  bool dotrunc = LT_.isNotNull() | UT_.isNotNull();

  if (LT_.isNotNull()) {
    LT = as<double>(LT_);
  }
  if (UT_.isNotNull()) {
    UT = as<double>(UT_);
  }

  // Calculate censoring
  double LC;
  double UC;
  Nullable<double> LC_ = data.attr("LC");
  Nullable<double> UC_ = data.attr("UC");

  if (LC_.isNotNull()) {
    LC = as<double>(LC_);
  }
  if (UC_.isNotNull()) {
    UC = as<double>(UC_);
  }

  // Response known
  // Fast
  LogicalVector neginf = rts == R_NegInf; // also sets NA to FALSE
  LogicalVector nortfast = neginf & !is_na(R);

  if (is_true(any(nortfast))) {
    NumericMatrix mparsfast(sum(nortfast),pars.ncol());
    for (int i = 0, j = 0; i < nortfast.length(); i++) {
      if (nortfast[i]) {
        mparsfast(j,_) = pars(i,_);
        j++;
      }
    }

    LogicalVector winnerfastvec = winner[nortfast];
    LogicalMatrix winnerfast(n_acc, winnerfastvec.length() / n_acc, winnerfastvec.begin());

    LogicalVector tofixfast = (winner & nortfast);
    NumericVector ldstofixfast(sum(tofixfast));

    for (int i = 0; i < sum(tofixfast); i++) {
      NumericMatrix pifast(n_acc, pars.ncol());
      if (n_acc == 1) {
        pifast(0,_) = mparsfast(i,_);
      } else {
        for (int j = 0; j < n_acc; j++) {
          pifast(j,_) = mparsfast(i * n_acc + j,_);
        }
      }
      NumericVector tmp = f_integrate(pifast, winnerfast(_,i), dfun, pfun, exp(min_ll), LT, LC);
      ldstofixfast[i] = std::log(std::max(0.0, std::min(tmp[0], 1.0)));
    }
    lds[tofixfast] = ldstofixfast;
  }

  // Slow
  LogicalVector posinf = rts == R_PosInf; // also sets NA to FALSE
  LogicalVector nortslow = posinf & !is_na(R);

  if (is_true(any(nortslow))) {
    NumericMatrix mparsslow(sum(nortslow),pars.ncol());
    for (int i = 0, j = 0; i < nortslow.length(); i++) {
      if (nortslow[i]) {
        mparsslow(j,_) = pars(i,_);
        j++;
      }
    }

    LogicalVector winnerslowvec = winner[nortslow];
    LogicalMatrix winnerslow(n_acc, winnerslowvec.length() / n_acc, winnerslowvec.begin());

    LogicalVector tofixslow = (winner & nortslow);
    NumericVector ldstofixslow(sum(tofixslow));

    for (int i = 0; i < sum(tofixslow); i++) {
      NumericMatrix pislow(n_acc, pars.ncol());
      if (n_acc == 1) {
        pislow(0,_) = mparsslow(i,_);
      } else {
        for (int j = 0; j < n_acc; j++) {
          pislow(j,_) = mparsslow(i * n_acc + j,_);
        }
      }
      NumericVector tmp = f_integrate(pislow, winnerslow(_,i), dfun, pfun, exp(min_ll), UC, UT);
      ldstofixslow[i] = std::log(std::max(0.0, std::min(tmp[0], 1.0)));
    }
    lds[tofixslow] = ldstofixslow;
  }

  // No direction
  LogicalVector nortno = is_na(rts) & !is_na(R);

  if (is_true(any(nortno))) {
    NumericMatrix mparsno(sum(nortno), pars.ncol());

    for (int i = 0, j = 0; i < nortno.length(); i++) {
      if (nortno[i]) {
        mparsno(j,_) = pars(i,_);
        j++;
      }
    }

    LogicalVector winnernovec = winner[nortno];
    LogicalMatrix winnerno(n_acc, winnernovec.length() / n_acc, winnernovec.begin());

    LogicalVector tofixno = (winner & nortno);
    NumericVector ldstofixno(sum(tofixno));

    for (int i = 0; i < sum(tofixno); i++) {
      NumericMatrix pino(n_acc, pars.ncol());
      if (n_acc  == 1) {
        pino(0,_) = mparsno(i,_);
      } else {
        for (int j = 0; j < n_acc; j++) {
          pino(j,_) = mparsno(i * n_acc + j,_);
        }
      }
      NumericVector tmpslow = f_integrate(pino, winnerno(_,i), dfun, pfun, exp(min_ll), UC, UT);
      NumericVector tmpfast = f_integrate(pino, winnerno(_,i), dfun, pfun, exp(min_ll), LT, LC);
      double tmp = tmpslow[0] + tmpfast[0];

      ldstofixno[i] = std::log(std::max(0.0, std::min(tmp, 1.0)));
    }
    lds[tofixno] = ldstofixno;
  }

  // Response unknown
  // Fast
  LogicalVector nortfastu = (rts == R_NegInf) & is_na(R);
  LogicalVector tofixfast(winner.length());

  if (is_true(any(nortfastu))) {
    NumericMatrix mpars(sum(nortfastu), pars.ncol());

    for (int i = 0, j = 0; i < nortfastu.length(); i++) {
      if (nortfastu[i]) {
        mpars(j,_) = pars(i,_);
        j++;
      }
    }

    LogicalVector winnerfastuvec = winner[nortfastu];
    LogicalMatrix winnerfastu(n_acc, winnerfastuvec.length() / n_acc, winnerfastuvec.begin());

    tofixfast = (winner & nortfastu);
    NumericVector ldstofixfast(sum(tofixfast));

    for (int i = 0; i < sum(tofixfast); i++) {
      NumericMatrix pi(n_acc, pars.nrow());

      for (int j = 0; j < n_acc; j++) {
        pi(j,_) = mpars(i * n_acc + j,_);
      }

      LogicalVector idx(n_acc, 0);
      idx[0] = 1;

      NumericVector pc = f_integrate(pi, idx, dfun, pfun, exp(min_ll), LT, LC);
      double p;
      if (pc[2] != 0 || traits::is_nan<REALSXP>(pc[0])) {
        p = NA_REAL;
      } else{
        p = std::max(0.0 ,std::min(pc[0],1.0));
      }

      double cf;
      if (p != 0 && !(LT==0 && UT==R_PosInf)) {
        cf = pr_pt(pi, idx, dfun, pfun, exp(min_ll), LT, UT);
      } else {
        cf = 1;
      }

      if (!traits::is_na<REALSXP>(cf)) {
        p *= cf;
      }

      if (!traits::is_na<REALSXP>(p) && n_acc > 1) {
        for (int j = 1; j < n_acc; j++) {
          idx.fill(0);
          idx[j] = 1;
          pc = f_integrate(pi, idx, dfun, pfun, exp(min_ll), LT, LC);
          if (pc[2] != 0 || traits::is_nan<REALSXP>(pc[0])) {
            p = NA_REAL;
            break;
          }
          if (pc[0] != 0.0 && !(LT == 0 && UT == R_PosInf)) {
            cf = pr_pt(pi, idx, dfun, pfun, exp(min_ll), LT, UT);
          } else{
            cf = 1;
          }
          if (!traits::is_na<REALSXP>(cf)) {
            p += pc[0] * cf;
          }
        }
      }
      double lp = std::log(p);
      if (!traits::is_na<REALSXP>(lp)) {
        ldstofixfast[i] = lp;
      } else{
        ldstofixfast[i] = R_NegInf;
      }
    }
    lds[tofixfast] = ldstofixfast;
  }

  // Slow
  LogicalVector nortslowu = (rts == R_PosInf) & is_na(R);
  LogicalVector tofixslow(winner.length());

  if (is_true(any(nortslowu))) {
    NumericMatrix mpars(sum(nortslowu), pars.ncol());

    for (int i = 0, j = 0; i < nortslowu.length(); i++) {
      if (nortslowu[i]) {
        mpars(j,_) = pars(i,_);
        j++;
      }
    }

    LogicalVector winnerslowuvec = winner[nortslowu];
    LogicalMatrix winnerslowu(n_acc, winnerslowuvec.length() / n_acc, winnerslowuvec.begin());

    tofixslow = (winner & nortslowu);
    NumericVector ldstofixslow(sum(tofixslow));

    for (int i = 0; i < sum(tofixslow); i++) {
      NumericMatrix pi(n_acc, pars.nrow());

      for (int j = 0; j < n_acc; j++) {
        pi(j,_) = mpars(i * n_acc + j,_);
      }

      LogicalVector idx(n_acc);
      idx[0] = 1;

      NumericVector pc = f_integrate(pi, idx, dfun, pfun, exp(min_ll), UC, UT);
      double p;
      if (pc[2] != 0 || traits::is_nan<REALSXP>(pc[0])) {
        p = NA_REAL;
      } else{
        p = std::max(0.0,std::min(pc[0],1.0));
      }

      double cf;
      if (p != 0 && !(LT==0 && UT==R_PosInf)) {
        cf = pr_pt(pi, idx, dfun, pfun, exp(min_ll), LT, UT);
      } else {
        cf = 1;
      }

      if (!traits::is_na<REALSXP>(cf)) {
        p *= cf;
      }

      if (!traits::is_na<REALSXP>(p) && n_acc > 1) {
        for (int j = 1; j < n_acc; j++) {
          idx.fill(0);
          idx[j] = 1;
          pc = f_integrate(pi, idx, dfun, pfun, exp(min_ll), UC, UT);
          if (pc[2] != 0 || traits::is_nan<REALSXP>(pc[0])) {
            p = NA_REAL;
            break;
          }
          if (pc[0] != 0.0 && !(LT == 0 && UT == R_PosInf)) {
            cf = pr_pt(pi, idx, dfun, pfun, exp(min_ll), LT, UT);
          } else{
            cf = 1;
          }
          if (!traits::is_na<REALSXP>(cf)) {
            p += pc[0] * cf;
          }
        }
      }
      double lp = std::log(p);
      if (!traits::is_na<REALSXP>(lp)) {
        ldstofixslow[i] = lp;
      } else{
        ldstofixslow[i] = R_NegInf;
      }
    }
    lds[tofixslow] = ldstofixslow;
  }

  // No direction
  LogicalVector nortnou = is_na(rts) & is_na(R);
  nortnou = nortnou & (pCont == 0); // Otherwise non-identifiable

  if (is_true(any(nortnou))) {
    NumericMatrix mpars(sum(nortnou), pars.ncol());

    for (int i = 0, j = 0; i < nortnou.length(); i++) {
      if (nortnou[i]) {
        mpars(j,_) = pars(i,_);
        j++;
      }
    }

    LogicalVector winnernouvec = winner[nortnou];
    LogicalMatrix winnernou(n_acc, winnernouvec.length() / n_acc, winnernouvec.begin());

    LogicalVector tofix = (winner & nortnou);
    NumericVector ldstofix(sum(tofix));

    for (int i = 0; i < sum(tofix); i++) {
      NumericMatrix pi(n_acc, pars.ncol());

      for (int j = 0; j < n_acc; j++) {
        pi(j,_) = mpars(i * n_acc + j,_);
      }

      LogicalVector idx(n_acc);
      idx[0] = 1;

      double pc = pLU(pi, idx, dfun, pfun, exp(min_ll), LT, LC, UC, UT);
      double p;
      double cf;
      if (traits::is_na<REALSXP>(pc)) {
        p = NA_REAL;
      } else{
        if (pc != 0.0 && !(LT == 0 && UT == R_PosInf)) {
          cf = pr_pt(pi, idx, dfun, pfun, exp(min_ll), LT, UT);
        } else{
          cf = 1;
        }

        if (!traits::is_na<REALSXP>(cf)) {
          p = pc*cf;
        } else {
          p = NA_REAL;
        }

        if (!traits::is_na<REALSXP>(p) && n_acc > 1) {
          for (int j = 1; j < n_acc; j++) {
            idx.fill(0);
            idx[j] = 1;

            pc = pLU(pi, idx, dfun, pfun, exp(min_ll), LT, LC, UC, UT);
            if (traits::is_na<REALSXP>(pc)) {
              p = NA_REAL;
              break;
            }
            if (pc != 0 && !(LT == 0 && UT == R_PosInf)) {
              cf = pr_pt(pi, idx, dfun, pfun, exp(min_ll), LT, UT);
            } else {
              cf = 1;
            }

            if (traits::is_na<REALSXP>(cf)) {
              p = NA_REAL;
              break;
            } else{
              p += pc * cf;
            }
          }
        }
      }
      double lp = std::log(p);
      if (!traits::is_na<REALSXP>(lp)) {
        ldstofix[i] = lp;
      } else{
        ldstofix[i] = R_NegInf;
      }
    }
    lds[tofix] = ldstofix;
  }

  // Truncation where not censored or censored and response known
  LogicalVector unique_nort = data.attr("unique_nort");
  NumericVector uniquewinlike = lds[unique_nort & winner];
  LogicalVector ok = is_finite(uniquewinlike);
  CharacterVector uniquewinresp = R[unique_nort & winner];
  LogicalVector alreadyfixed = is_na(uniquewinresp);
  ok = ok & !alreadyfixed;

  if (dotrunc & is_true(any(ok))) {
    NumericVector expand_nort = data.attr("expand_nort");
    NumericMatrix tpars(sum(unique_nort),pars.ncol());
    for (int i=0, j=0; i < unique_nort.length(); i++) {
      if (unique_nort[i]) {
        tpars(j,_) = pars(i,_);
        j++;
      }
    }
    LogicalVector winnertruncvec = winner[unique_nort];
    LogicalMatrix winnertrunc(n_acc, winnertruncvec.length() / n_acc, winnertruncvec.begin());
    NumericVector cf = rep(NA_REAL, ok.length());

    for (int i = 0; i < ok.length(); i++) {
      if (ok[i]) {
        NumericMatrix pi(n_acc, pars.ncol());

        for (int j = 0; j < n_acc; j++) {
          pi(j,_) = tpars(i * n_acc + j , _ );
        }

        cf[i] = pr_pt(pi, winnertrunc(_,i), dfun, pfun, exp(min_ll), LT, UT);
      }
    }
    NumericVector cf_log = rep_each(log(cf), n_acc);
    NumericVector cf_exp = c_expand(cf_log, expand_nort);
    LogicalVector fix = winner & !is_na(cf_exp) & is_finite(cf_exp);
    if (is_true(any(fix))) {
      lds[fix] = lds[fix] + cf_exp[fix];
    }
    LogicalVector badfix = winner & (is_na(cf_exp) | is_infinite(cf_exp));
    if (all(!is_na(tofixfast))) {
      badfix = badfix & !tofixfast;
    }
    if (all(!is_na(tofixslow))) {
      badfix = badfix & !tofixslow;
    }
    if (is_true(any(badfix))) {
      lds[badfix] = R_NegInf;
    }
  }

  // Non-process (contaminant) miss.
  LogicalVector isPCont = (pCont > 0);
  if (is_true(any(isPCont))) {
    NumericVector p = exp(lds[winner]);
    NumericVector pc = pCont[winner];
    CharacterVector Rwin = R[winner];
    LogicalVector isMiss = is_na(Rwin);
    for (int i = 0; i < p.length(); i++) {
      if (isMiss[i]) {
        p[i] = pc[i] + (1  - pc[i]) * p[i];
      } else {
        p[i] = (1 - pc[i]) * p[i];
      }
    }
    NumericVector ldswinner = log(p);
    lds[winner] = ldswinner;
  }

  // original code
  lds_exp = c_expand(lds, expand); // decompress

  if(n_acc > 1){
    LogicalVector winner_exp = c_bool_expand(winner, expand);
    NumericVector ll_out = lds_exp[winner_exp];
    if(n_acc == 2){
      NumericVector lds_los = lds_exp[!winner_exp];
      ll_out = ll_out + lds_los;
    } else{
      NumericVector lds_los = lds_exp[!winner_exp];
      for(int z = 0; z < ll_out.length(); z++){
        ll_out[z] = ll_out[z] + sum(lds_los[seq( z * (n_acc -1), (z+1) * (n_acc -1) -1)]);
      }
    }
    ll_out[is_na(ll_out)] = min_ll;
    ll_out[is_infinite(ll_out)] = min_ll;
    ll_out[ll_out < min_ll] = min_ll;
    return(sum(ll_out));
  } else{
    lds_exp[is_na(lds_exp)] = min_ll;
    lds_exp[is_infinite(lds_exp)] = min_ll;
    lds_exp[lds_exp < min_ll] = min_ll;
    return(sum(lds_exp));
  }
}

// [[Rcpp::export]]
NumericVector calc_ll(NumericMatrix p_matrix, DataFrame data, NumericVector constants,
                      List designs, String type, CharacterVector p_types,
                      double min_ll, List group_idx){
  const int n_particles = p_matrix.nrow();
  const int n_trials = data.nrow();
  NumericVector lls(n_particles);
  NumericVector p_vector(p_matrix.ncol());
  CharacterVector p_names = colnames(p_matrix);
  NumericMatrix pars(n_trials, p_types.length());
  p_vector.names() = p_names;
  NumericVector expand = data.attr("expand");
  List dynamic = data.attr("dynamic");
  if(type == "DDM"){
    for(int i = 0; i < n_particles; i++){
      p_vector = p_matrix(i, _);
      pars = get_pars_matrix(p_vector, constants, transform_DDM, Ntransform_DDM, p_types, designs, n_trials, dynamic, data);
      lls[i] = c_log_likelihood_DDM(pars, data, n_trials, expand, min_ll, group_idx);
    }
  } else if(type == "TRDM"){
    LogicalVector winner = data["winner"];
    for(int i = 0; i < n_particles; i++){
      p_vector = p_matrix(i, _);
      pars = get_pars_matrix(p_vector, constants, transform_trdm, Ntransform_trdm, p_types, designs, n_trials, dynamic, data);
      lls[i] = c_log_likelihood_race_trdm(pars, data, n_trials, winner, expand, min_ll);
    }
  } else{
    LogicalVector winner = data["winner"];
    // Love me some good old ugly but fast c++ pointers
    NumericVector (*dfun)(NumericVector, NumericMatrix, LogicalVector, double);
    NumericVector (*pfun)(NumericVector, NumericMatrix, LogicalVector, double);
    NumericVector (*transform)(NumericVector);
    NumericMatrix (*Ntransform)(NumericMatrix, CharacterVector, DataFrame, List);
    // NumericMatrix (*Ttransform)(NumericMatrix);
    if(type == "LBA" || type == "ALBA" || type == "MLBA"){
      dfun = dlba_c;
      pfun = plba_c;
      transform = transform_lba;
      Ntransform = Ntransform_lba;
    } else if(type == "RDM" || type == "ARDM" || type == "MRDM"){
      dfun = drdm_c;
      pfun = prdm_c;
      transform = transform_rdm;
      Ntransform = Ntransform_rdm;
    } else{
      dfun = dlnr_c;
      pfun = plnr_c;
      transform = transform_lnr;
      Ntransform = Ntransform_lnr;
    }
    if(type == "RDM" || type == "LBA" || type == "LNR"){
      for(int i = 0; i < n_particles; i++){
        p_vector = p_matrix(i, _);
        pars = get_pars_matrix(p_vector, constants, transform, Ntransform, p_types, designs, n_trials, dynamic, data);
        lls[i] = c_log_likelihood_race(pars, data, dfun, pfun, n_trials, winner, expand, min_ll);
      }
    } else if(type == "MRDM" || type == "MLBA" || type == "MLNR"){
      for(int i = 0; i < n_particles; i++){
        p_vector = p_matrix(i, _);
        pars = get_pars_matrix(p_vector, constants, transform, Ntransform, p_types, designs, n_trials, dynamic, data);
        lls[i] = c_log_likelihood_race_missing(pars, data, dfun, pfun, n_trials, winner, expand, min_ll);
      }
    } else{ // apparently we're a Advantage model and not a standard model
      for(int i = 0; i < n_particles; i++){
        p_vector = p_matrix(i, _);
        pars = get_pars_matrix(p_vector, constants, transform, Ntransform, p_types, designs, n_trials, dynamic, data);
        lls[i] = c_log_likelihood_race_advantage(pars, data, dfun, pfun, n_trials, winner, expand, min_ll);
      }
    }

  }
  return(lls);
}

