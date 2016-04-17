//
//  Option.cpp
//  MonteCarlo
//
//  Created by Yoann on 22/02/14.
//  Copyright (c) 2014 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "Option.h"


//Constructeur par défaut
Option::Option()
{
  T_ = 0;
  K_ = 0;
}


//Destructeur
Option::~Option() {}

//Accesseur à la valeur de la maturité
double Option::T() const{
	return T_;
}

//Mutateur de la valeur de la maturité
double & Option::T() {
	return T_;
}

//Accesseur au strike
double Option::K() const{
	return K_;
}

double & Option::K() {
	return K_;
}