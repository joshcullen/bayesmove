// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <ctime>
#include <fstream>
using namespace Rcpp;

/***************************************************************************************************************************/
/*********************************                      UTILS          *****************************************************/
/***************************************************************************************************************************/

//' Internal function that generates nmat matrix to help with multinomial draws
//'
//' @param z An integer vector.
//' @param dat An integer vector.
//' @param ncateg An integer.
//' @param nbehav An integer.
//' @param nobs An integer.
//'
//'
// [[Rcpp::export]]
IntegerMatrix SummarizeDat(IntegerVector z, IntegerVector dat, int ncateg,
                           int nbehav, int nobs){
  IntegerMatrix nmat(nbehav,ncateg);
  LogicalVector na1=!is_na(dat);

  for (int i = 0; i < nobs; i++){
    if (na1[i]){
      nmat(z[i],dat[i])=nmat(z[i],dat[i])+1;
    }
  }
  return nmat;
}
