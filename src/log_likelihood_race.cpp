#include <Rcpp.h>
#include "utility_functions.h"
#include "dynamic.h"
#include "model_lnr.h"
#include "model_LBA.h"
#include "model_RDM.h"
#include "model_DDM.h"
#include "advantage.h"
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

NumericMatrix get_pars(NumericVector p_vector, NumericVector constants,
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


// log_likelihood_race_advantage <- function(p_vector,dadm,min_ll=log(1e-10))
// # Race model summed log likelihood
// {
//
//   pars <- get_pars(p_vector,dadm)
//
//   if (is.null(attr(pars,"ok")))
//     ok <- !logical(dim(pars)[1]) else ok <- attr(pars,"ok")
//
//       nr <- length(levels(dadm$R))
//       na <- nr-1
//     nt <- nrow(dadm)/nr
//
// # This seems wasteful, can pre-compute, perhaps pass as dadm attr?
//     winner <- rep(dadm$winner,each=na)
//       rt <- rep(dadm$rt,each=na)
//
// # log pdf of winning response accumulators
//       pdf <- matrix(attr(dadm,"model")()$dfun(rt=rt[winner],pars=pars[winner,]),
//                     nrow=na,ncol=nt)
// # cdf (NOT log) of winning response accumulators
//         if (na>1) cdfw <- matrix(attr(dadm,"model")()$pfun(rt=rt[winner],pars=pars[winner,]),
//             nrow=na,ncol=nt)
// # cdf (NOT log) of loosing response accumulators
//           cdfl <- array(attr(dadm,"model")()$pfun(rt=rt[!winner],pars=pars[!winner,]),
//                         dim=c(na,na,nt))
//           ia <- 1:na
//         if (na==1) ll <- as.vector(log(pdf)) + log(1-as.vector(cdfl)) else {
//           ll <- 0
//           for (i in ia) {    # sum over accumulators for winning response
//           ifirst <- ia[-i] # other accumulators, already done
// # finish at t x already done
//           ll <- ll + pdf[i,]*apply(cdfw[ifirst,,drop=FALSE],2,sum)
//           }
//           ll <- log(ll) + apply(log(1-apply(cdfl, 2:3, prod)),2,sum) # other responses survivors at t
//         }
//         ll[is.na(ll) | !ok] <- min_ll
//         return(sum(pmax(min_ll,ll[attr(dadm,"expand_winner")])))
// }


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
  // Rcout << ll;
  // Rcout << "\n \n";
  ll = c_expand(ll, expand);
  ll[is_na(ll)] = min_ll;
  ll[is_infinite(ll)] = min_ll;
  ll[ll < min_ll] = min_ll;
  return(sum(ll));
}

// NumericVector calc_ll(NumericMatrix p_matrix, DataFrame data, NumericVector constants,
//                       List designs, String type, CharacterVector p_types,
//                       double min_ll, List group_idx){
//   const int n_particles = p_matrix.nrow();
//   const int n_trials = data.nrow();
//   NumericVector lls(n_particles);
//   NumericVector p_vector(p_matrix.ncol());

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
      pars = get_pars(p_vector, constants, transform_DDM, Ntransform_DDM, p_types, designs, n_trials, dynamic, data);
      lls[i] = c_log_likelihood_DDM(pars, data, n_trials, expand, min_ll, group_idx);
    }
  } else if(type == "TRDM"){
    LogicalVector winner = data["winner"];
    for(int i = 0; i < n_particles; i++){
      p_vector = p_matrix(i, _);
      pars = get_pars(p_vector, constants, transform_trdm, Ntransform_trdm, p_types, designs, n_trials, dynamic, data);
      lls[i] = c_log_likelihood_race_trdm(pars, data,n_trials, winner, expand, min_ll);
    }
  } else{
    LogicalVector winner = data["winner"];
    // Love me some good old ugly but fast c++ pointers
    NumericVector (*dfun)(NumericVector, NumericMatrix, LogicalVector, double);
    NumericVector (*pfun)(NumericVector, NumericMatrix, LogicalVector, double);
    NumericVector (*transform)(NumericVector);
    NumericMatrix (*Ntransform)(NumericMatrix, CharacterVector, DataFrame, List);
    // NumericMatrix (*Ttransform)(NumericMatrix);
    if(type == "LBA" || type == "ALBA"){
      dfun = dlba_c;
      pfun = plba_c;
      transform = transform_lba;
      Ntransform = Ntransform_lba;
    } else if(type == "RDM" || type == "ARDM"){
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
        pars = get_pars(p_vector, constants, transform, Ntransform, p_types, designs, n_trials, dynamic, data);
        lls[i] = c_log_likelihood_race(pars, data, dfun, pfun, n_trials, winner, expand, min_ll);
      }
    } else{ // apparently we're a Advantage model and not a standard model
      for(int i = 0; i < n_particles; i++){
        p_vector = p_matrix(i, _);
        pars = get_pars(p_vector, constants, transform, Ntransform, p_types, designs, n_trials, dynamic, data);
        lls[i] = c_log_likelihood_race_advantage(pars, data, dfun, pfun, n_trials, winner, expand, min_ll);
      }
    }

  }
  return(lls);
}

