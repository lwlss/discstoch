#include <Rcpp.h> 
using namespace Rcpp;

/*** R
aft = function(x,y,t){
  af = approx(x,y,t,method="constant")
  return(af$y)
}
*/

// [[Rcpp::export]]
double simDt_cpp(int K=1000, double r=1.0, int N0=1, int NSwitch=100, double t0=0, double t1=1) {
  if (NSwitch>N0){
    int eventN0=NSwitch-N0;
    NumericVector unifs=runif(eventN0);
    IntegerVector nn = seq(N0,NSwitch);
    NumericVector clist = as<NumericVector>(nn);
    NumericVector dts=-log(unifs)/(r*clist[seq(1,(eventN0))]*(1-clist[seq(1,(eventN0))]/K));
    NumericVector ats1=cumsum(dts);
    NumericVector atsadd(eventN0,t0);
    NumericVector ats=ats1+atsadd;
    ats.push_front(t0);
    double tmax=max(ats);
    if (tmax>=t1){
      NumericVector aft1;
      Environment G = Rcpp::Environment::global_env();
      Function aft = G["aft"];
      aft1=aft(ats,clist,t1);
      return aft1[0];
    }else{
      return K*NSwitch*exp(r*(t1-tmax))/(K+NSwitch*(exp(r*(t1-tmax))-1));
    }
  }else{
    return K*N0*exp(r*(t1-t0))/(K+N0*(exp(r*(t1-t0))-1));
  }
}