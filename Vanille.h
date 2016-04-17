#ifndef MonteCarlo_Vanille_h
#define MonteCarlo_Vanille_h
#include "Option.h"
#include <complex>


/**
 * \struct Vanille
 * \brief Description de l'option Vanille
 */

class Vanille :public Option
{
private :
  std::complex<double> D(Heston* modele,double r,double t, double phi,double b,double u);
  std::complex<double> C(Heston* modele,double r,double t, double phi,double b,double u);
  std::complex<double> F(Heston* modele,double r,double x,double v,double t, double phi,double b,double u);
  double G(Heston* modele,double r,double x,double v,double t, double phi,double b,double u);
  double integrande(Heston* modele,double r,double t, double phi);
  
  void evalParam(double x, Heston* modele,std::complex<double> &h1, std::complex<double> &h2);
public:
	/**
   * Constructeur par defaut
   */ 
  Vanille();
  
	/**
   * Constructeur avec paramètre
   */   
  Vanille(double strike, double T);
  
  
  /**
   * Destructeur
   */
  ~Vanille();
  
  /**
   * Calcul la valeur du payoff sur la trajectoire passée en paramètre
   *
   * @param path (input) est une matrice de taille 1 x (N+1) contenant une
   * trajectoire du modèle
   * @return phi(trajectoire)
   */
  virtual double payoff (const std::vector<double> *path);
  
  /**
   * Calcul le prix de l'option selon une formule semi fermée
   */
  virtual double computePrice (Heston* modele,double r, double t);
  
  
  
  
  
};


#endif
