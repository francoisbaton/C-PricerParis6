//
//  Option.h
//  MonteCarlo
//
//  Created by Yoann on 22/02/14.
//  Copyright (c) 2014 __MyCompanyName__. All rights reserved.
//

#ifndef MonteCarlo_Option_h
#define MonteCarlo_Option_h

#include <vector>
#include "Heston.h"

class Option {
  
protected:
  double T_; /*! maturité */
  double K_; /*! Strike */
  
  
public:
  
  /**
   * Constructeur par defaut
   */ 
  Option();
  
  /**
   * Destructeur
   */
  virtual ~Option() = 0;
  
  
  /**
   * Accesseur à la valeur de la maturité
   * @return T
   */
  double T() const;
  
  /**
   * Mutateur de la valeur de la maturité
   * @return T
   */
  double & T();
  
  /**
   * Accesseur au strike
   * @return K
   */
  double K() const;
  
  /**
   * Mutateur du strike
   * @return K
   */
  double & K();
  
  
  /**
   * Calcul la valeur du payoff sur la trajectoire passée en paramètre
   *
   * @param path (input) est une matrice de taille 1 x (N+1) contenant une
   * trajectoire du modèle
   * @return prix
   */
  virtual double payoff (const std::vector<double> *path) = 0;

  /**
   * Calcul le prix de l'option selon une formule semi fermée
   */
  virtual double computePrice(Heston *modele,double r, double t) = 0;
  
};


#endif
