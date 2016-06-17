#include <Rcpp.h> 
using namespace Rcpp;

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
      NumericVector which=ats[ats < t1];
      int ind=which.length();
      double ats_t1=(clist[ind-1]-clist[ind])/(ats[ind-1]-ats[ind])*(t1-ats[ind])+clist[ind];
      return ats_t1;
    }else{
      return K*NSwitch*exp(r*(t1-tmax))/(K+NSwitch*(exp(r*(t1-tmax))-1));
    }
  }else{
    return K*N0*exp(r*(t1-t0))/(K+N0*(exp(r*(t1-t0))-1));
  }
}

// [[Rcpp::export]]
List pf_cpp(Function simx0, int n, int t0, 
            NumericVector times, NumericVector deltas, 
            Function dataLik, Function stepFun, 
            int i, NumericVector area){
  NumericMatrix w(n,1);
  NumericMatrix xmat=simx0(n,t0);
  for (int j=0; j<n; j=j+1){
    NumericVector xm=stepFun(xmat(j,_), times(i-1), deltas(i-1));
    xmat(j,_)=xm;
    NumericVector wj = dataLik(xmat(j,_),times(i),area(i-1));
    w(j,_)=wj;
  }
  return List::create(_["xmat"]=xmat, _["w"]=w);
}

// [[Rcpp::export]]
NumericMatrix mcmc_cpp(int p, double tune, int iters, int thin, Function mLLik, 
                       NumericVector th, NumericVector pmin, NumericVector pmax){
  double ll=-1e99;
  NumericMatrix thmat(iters,p);
  for (int i=0; i<iters; i=i+1){
    NumericVector thi=th;
    Rcout << i << " ";
    for(int j=0; j<thin; j=j+1){
      NumericVector thprob=pmin/2;
      while(sum((thprob<pmin)|(thprob>pmax))>0){
        thprob=thi*exp(rnorm(p,0,tune));
      }
      NumericVector llprob=mLLik(thprob);
      NumericVector pr=log(runif(1));
      if (pr[0]<(llprob[0]-ll)){
        thi=thprob;
        ll=llprob[0];
      }
    }
    thmat(i,_)=thi;
  }
  return thmat;
}
