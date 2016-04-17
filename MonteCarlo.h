#ifndef MonteCarlo_MonteCarlo_h
#define MonteCarlo_MonteCarlo_h
#include "mtrand.h"
#include <cmath>
#include "Option.h"
#include "Heston.h"
#include "BS.h"
#include "MC.h"

class MonteCarlo : public MC
{
  
public:
  MonteCarlo();
  
  MonteCarlo(int M);
  
  ~MonteCarlo();
  
  int M() const;
  int & M();
  
  
  // Calcul E[Sup|X_2n - X_n|]
  //virtual void Sn(double &mean,double &varIc,CIR *modeleCir,void (CIR::*asset)(std::vector<double> *,double , int,MTRand  *),double T, int N);
  
  
  // Calcul E[f(X_n)]
  //virtual void weakError(double &mean,double &varIc,CIR *modeleCir,void (CIR::*asset)(std::vector<double> *,double , int,MTRand  *),double T, int N);

  //Calcul le prix d'une option selon une méthode de discretisation passée en paramètre
  virtual void computePriceOption(double &price, double &varIc,int N, Option *opt, Heston* modele,void (Heston::*asset)(std::vector<double> *pathX,std::vector<double> *pathV, double T, int N,MTRand  *myGen));
  
  //Calcul le prix Black-Scholes
   virtual void computePriceOptionBS(double &price, double &varIc,int N,Option *opt, BS* modele);

};


#endif
