#include <iostream>
#include "BS.h"
#include <cmath>


BS::BS(){
  S_ = 0;
  r_ = 0;
  sigma_ = 0;
}

BS::BS(double S, double r, double sigma){
  S_ = S;
  r_ = r;
  sigma_ = sigma;
}

BS::~BS(){
}

double &BS::S(){
  return S_;
}

double & BS::r(){
  return r_;
}

double & BS::sigma(){
  return sigma_;
}

double BS::S() const{
  return S_;
}

double BS::r() const{
  return r_;
}

double BS::sigma() const{
  return sigma_;
}

void BS::asset(std::vector<double> *path, double T, int N,MTRand  *myGen){
  
  double St = S_;
  double Zs,U;
  double sigmaW = sigma_*sqrt(T/double(N));
  double varW = T/double(N);
  double drift = (r_ - pow(sigma_,2)/2)*varW;
  
  (*path)[0] = S_;
  
  for (int i = 1; i <= N; i++){
    U = (*myGen)();
    Zs = normsinv(U);
    
    St = St*exp(drift + sigmaW*Zs);
    (*path)[i] = St;
  }
}


