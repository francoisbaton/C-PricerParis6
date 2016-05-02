#include <iostream>
#include "Heston.h"

Heston::Heston(){
  theta_  = 0;
  k_ = 0;
  sigma_  = 0;
  V0_ = 0;
  X0_ = 0;
  rho_ = 0;
}

Heston::Heston(double theta,double k,double sigma,double V0,double X0,double rho){
  theta_  = theta;
  k_ =  k ;
  sigma_  = sigma;
  V0_ = V0;
  X0_ = X0;
  rho_ = rho;  
}

Heston::~Heston(){
  
}

double Heston::theta() const{
  return theta_;
}

double Heston::k() const{
  return k_;
}

double Heston::sigma() const{
  return sigma_;
}
double Heston::V0() const{
  return V0_;
}

double Heston::X0() const{
  return X0_;
}

double Heston::rho() const{
  return rho_;
}

double & Heston::theta(){
  return theta_;
}

double & Heston::k(){
  return k_;
}
double & Heston::sigma(){
  return sigma_;
}
double & Heston::V0(){
  return V0_;
}
double & Heston::X0(){
  return X0_;
}
double & Heston::rho(){
  return rho_;
}

void Heston::assetEuler (std::vector<double> *pathX,std::vector<double> *pathV, double T, int N,MTRand  *Gen){
  
  double sigmaBis = sqrt(T/double(N));
  double varBis = T/double(N);
  double corr = sqrt(1-pow(rho_,2));
  double Zx,Zv,U;
  double logXt, Vt;
  
  (*pathX)[0] = X0_;
  (*pathV)[0] = V0_;
  
  logXt = log(X0_);
  Vt = (*pathV)[0];
  
  for (int i = 1; i <= N; i++){
    
    if (Vt > 0){
      U = (*Gen)();
      Zv = normsinv(U);
      U = (*Gen)();
      Zx = normsinv(U);
      Zx = rho_*Zv + corr*Zx;
      
      logXt = logXt - Vt*varBis/double(2) + sqrt(Vt)*sigmaBis*Zx;
      Vt = Vt + k_*(theta_ - Vt)*varBis + sigma_*sqrt(Vt)*sigmaBis*Zv;
    }
    else{
      Vt = Vt + k_*theta_*varBis;
    }
    
    (*pathX)[i] = exp(logXt);
    (*pathV)[i] = Vt;
  }
}


void Heston::assetJK (std::vector<double> *pathX,std::vector<double> *pathV, double T, int N,MTRand  *Gen){
  double sigmaW = sqrt(T/double(N));
  double varW = T/double(N);
  double varV = pow(sigma_,2);
  double corr = sqrt(1-pow(rho_,2));
  double Zx,Zv,U;
  double logXt, Vt, VtNext;
  
  
  (*pathX)[0] = X0_;
  (*pathV)[0] = V0_;
  
  logXt = log(X0_);
  Vt = (*pathV)[0];
  
  for (int i = 1; i <= N; i++){
    
    if (Vt > 0){
      
      U = (*Gen)();
      Zv = normsinv(U);
      U = (*Gen)();
      Zx = normsinv(U);

      Zx = rho_*Zv + corr*Zx;
      
      VtNext = (Vt + k_*theta_*varW + sigma_*sqrt(Vt)*sigmaW*Zv + varV*varW*(pow(Zv,2)-1)/double(4))/(double(1) + k_*varW);
      
      if (VtNext > 0){
        logXt = logXt - varW*(VtNext + Vt)/double(4) + rho_*sqrt(Vt)*sigmaW*Zv + (sqrt(VtNext) + sqrt(Vt))*sigmaW*(Zx - rho_*Zv)/double(2) + 
              sigma_*rho_*varW*(pow(Zv,2)-1)/double(4);
      }
      else{
        logXt = logXt - varW*(Vt)/double(4) + rho_*sqrt(Vt)*sigmaW*Zv + (sqrt(Vt))*sigmaW*(Zx - rho_*Zv)/double(2) + 
        sigma_*rho_*varW*(pow(Zv,2)-double(1))/double(4);        
      }
    }
    else{
      U = (*Gen)();
      Zv = normsinv(U);
      
      VtNext = Vt + k_*(theta_)*varW;
      if (VtNext > 0){
        U = (*Gen)();
        Zx = normsinv(U);
        Zx = rho_*Zv + corr*Zx;
        
        logXt = logXt - varW*(VtNext)/double(4) + (sqrt(VtNext))*sigmaW*(Zx - rho_*Zv)/double(2) + 
        sigma_*rho_*varW*(pow(Zv,2)-double(1))/double(4);        
      }
      else{
        logXt = logXt + sigma_*rho_*varW*(pow(Zv,2)-double(1))/double(4);
      }
    }
    
    (*pathX)[i] = exp(logXt);
    (*pathV)[i] = VtNext;
    Vt = VtNext;
  }
}


double Heston::inversePsi(double u,double p, double beta){
  if (u <= p){
    return 0;
  }
  else{
    return log((double(1)-p)/(double(1)-u))/beta;
  }
}


void Heston::assetQE (std::vector<double> *pathX,std::vector<double> *pathV, double T, int N,MTRand  *Gen){
  double varW = T/double(N);
  double varV = pow(sigma_,2);
  double Zx,Zv, U;
  double logXt, Vt, VtNext;
  double m, s2, psi, psiC;
  double a,b2, beta, p;
  
  double g1, g2;
  double K0,K1,K2,K3,K4;
  
  g1 = double(1)/double(2);
  g2 = double(1)/double(2);
  
  K0 = -rho_*k_*theta_*varW/sigma_;
  K1 = g1*varW*(k_*rho_/sigma_ - double(1)/double(2)) - rho_/sigma_;
  K2 = g2*varW*(k_*rho_/sigma_ - double(1)/double(2)) + rho_/sigma_;
  K3 = g1*varW*(double(1)-pow(rho_,2));
  K4 = g2*varW*(double(1)-pow(rho_,2));
  
  
  psiC = 1.5;
  
  double ekD = exp(-k_*varW);
  double ekD2 = pow(double(1)-ekD,2);
  
  
  (*pathX)[0] = X0_;
  (*pathV)[0] = V0_;
  
  logXt = log(X0_);
  Vt = (*pathV)[0];
  
  for (int i = 1; i <= N; i++){
    
    //Constructin de V(t+delta)
    U = (*Gen)();
    
    
    m = theta_ + (Vt - theta_)*ekD;
    s2 = Vt*varV*ekD*(double(1)-ekD)/k_ + theta_*varV*ekD2/(double(2)*k_);
    psi = s2/pow(m,2);
    
    if (psi <= psiC){
      Zv = normsinv(U);
      b2  = double(2)/psi - double(1) + sqrt(double(2)/psi)*sqrt(double(2)/psi - double(1));
      a = m / (double(1) + b2);
      
      VtNext = a*pow(sqrt(b2) + Zv,2);
    }
    else{
      p = (psi - double(1))/(psi + double(1));
      beta = (double(1)-p)/m;
      VtNext = inversePsi(U, p, beta);
    }
    
    //Construction de log(X(t + Delta)
    U = (*Gen)();
    Zx = normsinv(U);
    
    logXt = logXt + K0 + K1*Vt + K2*VtNext + sqrt(K3*Vt + K4*VtNext)*Zx;
    
    (*pathX)[i] = exp(logXt);
    (*pathV)[i] = VtNext;
    Vt = VtNext;
    
  }
  
}
