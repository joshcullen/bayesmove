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

// This function helps with multinomial draws
int cat1(double value, NumericVector prob) {
  int res=prob.length()-1;
  double probcum = 0;
  
  for (int i = 0; i < prob.length(); i++) {
    probcum = probcum + prob(i);
    if (value < probcum) {
      res = i;
      break;
    }
  }
  return res;
}

//' This function samples z's from a categorical distribution
// [[Rcpp::export]]
IntegerVector rmultinom1(NumericMatrix prob, NumericVector randu) {
  
  IntegerVector z(prob.nrow());
  
  for(int i=0; i<prob.nrow();i++){
    z[i]=cat1(randu[i],prob(i,_));
  }
  return z;
}

//' This function samples z's from a multinomial distribution
// [[Rcpp::export]]
IntegerVector rmultinom2(NumericVector prob, int n, NumericVector randu, int nmaxclust) {
  IntegerVector ZAgg(nmaxclust);
  int tmp;
  
  for(int i=0; i<n;i++){
    tmp=cat1(randu[i],prob);
    ZAgg[tmp]=ZAgg[tmp]+1;
  }
  return ZAgg;
}

//' This function samples z1 aggregate
// [[Rcpp::export]]
List SampleZAgg(int ntsegm,int b1,IntegerMatrix y1, int nmaxclust,
                NumericMatrix lphi1, NumericMatrix ltheta,NumericVector zeroes){
  //convert array into arma::cube
  NumericVector vecArray=clone(zeroes);
  arma::cube Z1Agg(vecArray.begin(),ntsegm, b1, nmaxclust);
  
  IntegerVector tmp(nmaxclust);
  NumericVector lprob(nmaxclust);
  NumericVector prob(nmaxclust);
  for (int i=0; i<ntsegm; i++){
    for (int j=0; j<b1; j++){
      if (y1(i,j)>0){
        //calculate probability
        for (int k=0; k<nmaxclust; k++){
          lprob[k]=lphi1(k,j)+ltheta(i,k);
        }
        lprob=lprob-max(lprob);
        prob=exp(lprob);
        prob=prob/sum(prob);
        
        //sample from multinomial and store results
        tmp=rmultinom2(prob,y1(i,j),runif(y1(i,j)),nmaxclust);
        for (int k=0; k<nmaxclust; k++){
          Z1Agg(i,j,k)=tmp[k];
        }
      }
    }
  }
  List L = List::create(Named("Z1Agg") =Z1Agg);
  
  return(L);
}

//' This function calculates the inverted cumsum
// [[Rcpp::export]]
IntegerMatrix CumSumInv(int ntsegm, int nmaxclust, IntegerMatrix z){
  IntegerMatrix res(ntsegm,z.ncol());
  IntegerVector soma(ntsegm);
  
  for (int j=z.ncol()-1; j> -1; j--){
    if (j==z.ncol()-1){
      soma=z(_,j);
    }            
    if (j!=z.ncol()-1){
      soma=soma+z(_,j);
    }
    res(_,j)=soma;
  }
  return(res);
}