#include <iostream>
#include "MonteCarlo.h"





MonteCarlo::MonteCarlo() : MC(){
}

MonteCarlo::MonteCarlo(int M) : MC(M){
}

MonteCarlo::~MonteCarlo(){
  
}


/*void MonteCarlo::Sn(double &mean, double &varIc,CIR *modeleCir,void (CIR::*asset)(std::vector<double> *,double , int,MTRand  *),double T, int N){
  
  MTRand  myGen(2);
  mean = 0;
  varIc = 0; //largeur de l'intervalle de confiance
  double tmp;
  
  for (int i = 0; i < M_; i++){
    
   tmp = modeleCir->supX(asset, T, N, &myGen);
   mean = mean +tmp;
   varIc = varIc + pow(tmp,2);
  }
  
  mean = mean / double(M_);
  varIc = varIc/double(M_) - pow(mean,2);
  varIc = 2*1.96*sqrt(varIc)/sqrt(double(M_));
  
}



void MonteCarlo::weakError(double &mean,double &varIc,CIR *modeleCir,void (CIR::*asset)(std::vector<double> *,double , int,MTRand  *),double T, int N){
 
  MTRand  myGen(2);
  mean = 0;
  varIc = 0; //largeur de l'intervalle de confiance
  double tmp;
  
  for (int i = 0; i < M_; i++){
    
    tmp = modeleCir->estimateur(asset, T, N, &myGen);
    //tmp = modeleCir->estimateurRomberg(asset, T, N, &myGen);
    mean = mean +tmp;
    varIc = varIc + pow(tmp,2);
  }
  
  mean = mean / double(M_);
  varIc = varIc/double(M_) - pow(mean,2);
  varIc = 2*1.96*sqrt(varIc)/sqrt(double(M_));  
}

*/

void MonteCarlo::computePriceOption(double &price, double &varIc,int N,Option *opt, Heston* modele,void (Heston::*asset)(std::vector<double> *,std::vector<double> *, double, int ,MTRand *)){
  
  MTRand  myGen(2);
  price = 0;
  varIc = 0; //largeur de l'intervalle de confiance
  std::vector<double> pathX(N+1);
  std::vector<double> pathV(N+1);
  double T = opt->T();

  double tmp;
  
  for (int i = 0; i < M_; i++){
    
    (modele->*asset)(&pathX, &pathV, T, N, &myGen);
    
    tmp = opt->payoff(&pathX);
    price = price + tmp;
    varIc = varIc + pow(tmp,2);
  }
  
  price = price / double(M_);
  varIc = varIc/double(M_) - pow(price,2);
  varIc = 2*1.96*sqrt(varIc)/sqrt(double(M_));
  
}



void MonteCarlo::computePriceOptionBS(double &price, double &varIc,int N,Option *opt, BS* modele){
  price = 0;
  varIc = 0; //largeur de l'intervalle de confiance
  std::vector<double> pathS(N+1);
  double T = opt->T();
  double mean = 0;
  double std = 0;

  MTRand  myGen(2);
  double tmp;
   
  for (int i = 0; i < M_; i++){
    
    modele->asset(&pathS,T,N,&myGen);
    
    tmp = opt->payoff(&pathS);
    mean = mean + tmp;
    std = std + pow(tmp,2);
  }
  
  price = mean / double(M_);
  varIc = std/double(M_) - pow(price,2);
  varIc = 2*1.96*sqrt(varIc)/sqrt(double(M_));
  
}
  

