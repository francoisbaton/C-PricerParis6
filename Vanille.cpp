#include <iostream>
#include "Vanille.h"
#include <climits>
#include "mtrand.h"
#include <math.h>

//Constructeur par defaut

Vanille::Vanille():Option(){
}

Vanille::Vanille(double strike, double T){
  T_ = T;
  K_ = strike;
  
}


Vanille::~Vanille(){
}

//Méthode de calcul de payoff d'une option panier
double Vanille::payoff (const std::vector<double> *path) {
  
  int N = int((*path).size());
	double s_T = (*path)[N-1];

	if ((s_T-K_) >=0) {
		return s_T-K_;
	}
	else {
		return 0;
	}
}



std::complex<double> Vanille:: C(Heston* modele,double r,double t, double phi,double b,double u){
  std::complex<double> valC;
  double k = modele->k();
  double theta = modele->theta();
  double a = k*theta;
  
  double sigma = modele->sigma();
  double rho = modele->rho();
  
  //Construction de d
  std::complex<double> d;
  std::complex<double> tmp;
  std::complex<double> tmp2;
  
  tmp.real() = -b;
  tmp.imag() = rho*sigma*phi;
  
  tmp = tmp*tmp;
  
  d = tmp;
  
  tmp.real() = -double(1)*pow(phi,2);
  tmp.imag() = double(2)*u*phi;
  tmp = pow(sigma,2)*tmp;
  d = d - tmp;
  d = sqrt(d);
  
  //Construction de g
  std::complex<double> g;
  tmp = b + d;
  tmp.imag() = tmp.imag() - rho*sigma*phi;
  g  = tmp;
  
  tmp = b-d;
  tmp.imag() = tmp.imag() - rho*sigma*phi;
  g = g/tmp;
  
  
  //Construction de (1-exp(d*t))/(1-g)
  tmp = double(1) - exp(d*t)*g;
  tmp2 = double(1) - g;
  tmp = tmp/tmp2;
  
  tmp2.imag() = r*phi*t - a/(pow(sigma,2))*rho*sigma*phi*t;
  tmp2.real() = 0;
  
  
  valC = a/(pow(sigma,2))*((b+d)*t - double(2)*log(tmp));
  valC = valC + tmp2;
  

  return valC;
  
}


std::complex<double> Vanille :: D(Heston* modele,double r,double t, double phi,double b,double u){
  std::complex<double> valD;
  
  double sigma = modele->sigma();
  double rho = modele->rho();
  
  //Construction de d
  std::complex<double> d;
  std::complex<double> tmp;
  std::complex<double> tmp2;
  
  tmp.real() = -b;
  tmp.imag() = rho*sigma*phi;
  
  tmp = tmp*tmp;
  
  d = tmp;
  
  tmp.real() = -double(1)*pow(phi,2);
  tmp.imag() = double(2)*u*phi;
  tmp = pow(sigma,2)*tmp;
  d = d - tmp;
  d = sqrt(d);
  
  //Construction de g
  std::complex<double> g;
  tmp = b + d;
  tmp.imag() = tmp.imag() - rho*sigma*phi;
  g  = tmp;
  
  tmp = b-d;
  tmp.imag() = tmp.imag() - rho*sigma*phi;
  g = g/tmp;
  
  
  //Construction de (1-exp(dt))/(1-g*exp(dt))
  tmp = double(1) - exp(d*t);
  tmp2 = double(1) - g*exp(d*t);
  tmp = tmp/tmp2;
  
  tmp2 = b + d;
  tmp2.imag() = tmp2.imag() - rho*sigma*phi;
  tmp2 = tmp2/(pow(sigma,2));
  
  valD = tmp*tmp2;
  
  
  return valD;
  
}

std::complex<double> Vanille :: F(Heston* modele,double r,double x,double v,double t, double phi,double b,double u){
  std::complex<double> valF;
  
  valF = exp(C(modele,r,T_ - t, phi,b,u) + D(modele,r,T_ - t, phi,b,u)*v);
  std::complex<double> tmp;
  
  tmp.real() = 0.00;
  tmp.imag() = phi*x;
  tmp = exp(tmp);
  
  valF = valF * tmp;
  
  return valF;
  
}

double Vanille :: G(Heston* modele,double r,double x,double v,double t, double phi,double b,double u){
  double valG;
  
  std::complex<double> tmp;
  std::complex<double> tmp2;
  tmp.real() = 0.00;
  tmp.imag() = -phi*log(K_);
  tmp = exp(tmp);

  tmp = tmp*F(modele,r,x,v,t, phi,b,u);
  
  tmp2.real() = 0.00;
  tmp2.imag() = phi;
  
  tmp = tmp/tmp2;
  valG = tmp.real();
  
  return valG;
}

double Vanille::integrande(Heston* modele,double r,double t, double phi){
  double valIntegrande;
  
  double S0 = modele->X0();
  double V0 = modele->V0();
  
  double u1 = 0.5;
  double u2 = -u1;
  double sigma = modele->sigma();
  double rho = modele->rho();
  
  double k = modele->k();
  double b1 = k - rho*sigma;
  double b2 = k ;
  
  valIntegrande = S0*G(modele,r,log(S0),V0,t, phi,b1,u1) - K_*exp(-r*(T_-t))*G(modele,r,log(S0),V0,t, phi,b2,u2);
  
  return valIntegrande;
}


double Vanille::computePrice (Heston* modele,double r, double t){
 double price;
  double integrale = 0;
 double tmp;
  double phi;
 
 double S0 = modele->X0();
 
 //Approcimation de l'integrale
 
 double a = 0.00;
 double b = 1.00;
 
 int n = 10000;
 double h = (b-a)/double(n);
 a = 0.0007;
  
 
 //Calcul de l'integrale par la méthode des trapezes
 tmp = (integrande(modele, r, t, a) + integrande(modele, r, t, double(1)/a)/(a*a))/double(2);
  if (!std::isnan(tmp)){
    integrale = integrale + tmp;
  }
 tmp = (integrande(modele, r, t, b) + integrande(modele, r, t, double(1)/b)/(b*b))/double(2);
 if (!std::isnan(tmp)){
   integrale = integrale + tmp;
 }
 for (int i = 1; i < n; i++){
   phi = a + double(i)*h;
   tmp = integrande(modele, r, t, phi)+ integrande(modele, r, t, double(1)/phi)/(phi*phi);
   if (!std::isnan(tmp)){
     integrale = integrale + tmp;
   }
 }
 integrale = h*integrale/M_PI;
 
 price = (S0 - K_*exp(-r*(T_-t)))/double(2) + integrale;
 
 return price;
 }

