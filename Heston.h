#ifndef MonteCarlo_Heston_h
#define MonteCarlo_Heston_h
#include <vector>
#include "mtrand.h"
#include <cmath>


class Heston{
private:
  double theta_;
  double k_;
  double sigma_;
  double V0_;
  double X0_;
  double rho_;
  
  double inversePsi(double u, double p, double beta);
  
public:
  Heston();
  
  Heston(double theta,double k,double sigma,double V0,double X0,double rho);
  
  ~Heston();
  
  double theta() const;
  double k() const;
  double sigma() const;
  double V0() const;
  double X0() const;
  double rho() const;
  
  double & theta();
  double & k();
  double & sigma();
  double & V0();
  double & X0();
  double & rho();
  

  void assetEuler (std::vector<double> *pathX,std::vector<double> *pathV, double T, int N,MTRand  *gen);
  void assetJK (std::vector<double> *pathX,std::vector<double> *pathV, double T, int N,MTRand  *gen);
  void assetQE (std::vector<double> *pathX,std::vector<double> *pathV, double T, int N,MTRand  *gen);
                
};

#endif
