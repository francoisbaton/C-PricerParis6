#ifndef MonteCarlo_BS_h
#define MonteCarlo_BS_h
#include <vector>
#include "mtrand.h"

class BS{
private:
  double S_;
  double r_;
  double sigma_;
  
public:
  
  BS();
  
  BS(double S, double r, double sigma);
  
  ~BS();
  
  
  double & S();
  double & r();
  double & sigma();
  
  double S() const;
  double r() const;
  double sigma() const;
  
  void asset (std::vector<double> *path, double T, int N,MTRand  *myGen);
  
};

#endif
