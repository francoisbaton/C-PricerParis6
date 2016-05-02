#include <iostream>
#include "Vanille.h"
#include "mtrand.h"
#include "MonteCarlo.h"
#include <fstream>
#include <cmath>
#include "Heston.h"
#include <climits>
#include <ctime>
#include "BS.h"
#include "MC.h"


//Test de la classe Heston
int main (int argc, const char * argv[]){
  

  
  double sigma,k,rho,T,X0,V0,theta;

  
  sigma = 0.1;
  k = 2;
  rho = 0.5;
  T = 5;
  theta = 0.01;
  V0 = 0.01;
  X0 = 100;
  
  //creation d'un modele de Heston
  Heston modele(theta,k,sigma,V0,X0,rho);
  
  
  double delta = double(1)/double(2);
  int N = int(T/delta);
  
  MTRand myGen(37);
  
  
  std::vector<double> pathX(N+1);
  std::vector<double> pathV(N+1);
  
  
  //Test de la classe Option
  double K = 100;
  Vanille* opt = new Vanille(K,T);
  double r = 0;
  double t = 0;
  
  
  double priceFFHeston;
  priceFFHeston = opt->computePrice(&modele,r,t);
  std::cout << "prix Formule Fermée : " << priceFFHeston << "\n \n";
  
  
  //Test de la classe MonteCarlo
  int M = 1000000;
  MC* MC = new MonteCarlo(M);
  
  double priceHeston,varIc;
  time_t tbegin,tend;


  tbegin=time(NULL); 
  MC->computePriceOption(priceHeston, varIc,N, opt, &modele, &Heston::assetEuler);
  tend=time(NULL); 
  std::cout << "Prix MC Euler:  " << priceHeston << "\n";
  std::cout << "IC = [ " << priceHeston - varIc/2 << " , " << priceHeston + varIc/2 << " ]\n";
  std::cout << "Largeur " << varIc << "\n";
  std::cout << "Diff = " << priceFFHeston - priceHeston << "\n";
  std::cout << "Time = " << difftime(tend,tbegin) << "\n \n";
  
  tbegin=time(NULL);
  MC->computePriceOption(priceHeston, varIc,N, opt, &modele, &Heston::assetJK);
  tend=time(NULL);
  std::cout << "Prix MC JK:  " << priceHeston << "\n";
  std::cout << "IC = [ " << priceHeston - varIc/2 << " , " << priceHeston + varIc/2 << " ]\n";
  std::cout << "Largeur " << varIc << "\n";
  std::cout << "Diff = " << priceFFHeston - priceHeston << "\n";
  std::cout << "Time = " << difftime(tend,tbegin) << "\n \n";
  
  tbegin=time(NULL);
  MC->computePriceOption(priceHeston, varIc,N, opt, &modele, &Heston::assetQE);
  tend=time(NULL);
  std::cout << "Prix MC QE:  " << priceHeston << "\n";
  std::cout << "IC = [ " << priceHeston - varIc/2 << " , " << priceHeston + varIc/2 << " ]\n";
  std::cout << "Largeur " << varIc << "\n";
  std::cout << "Diff = " << priceFFHeston - priceHeston << "\n";
  std::cout << "Time = " << difftime(tend,tbegin) << "\n \n";
  
  
  //Difference entre prixBS et prix Heston
  double priceBS = 0;
  std::ofstream fichierDiff("comparisonBSHeston.txt"); 
  double S0 = 10;
  fichierDiff << "x=c(";
  
  for (int i = 0; i < 40; i++){
    fichierDiff << S0 << ",";
    S0 = S0 + 5;
  }
  fichierDiff << S0 << ")";
  
  fichierDiff << "\n" << "\n";
  
  fichierDiff << "y=c(";
  
  S0 = 10;
  
  for (int i = 0; i < 40; i++){
    modele.X0() = S0;
    
    modele.rho() = 0.5;
    priceHeston = opt->computePrice(&modele,r,t);
    modele.rho() = 0;
    priceBS = opt->computePrice(&modele,r,t);
    
    fichierDiff << priceHeston - priceBS << ",";
    S0 = S0 + 5;
  }
  modele.rho() = 0.5;
  priceHeston = opt->computePrice(&modele,r,t);
  modele.rho() = 0;
  priceBS = opt->computePrice(&modele,r,t);

  fichierDiff << priceHeston - priceBS << ")";
  
  fichierDiff << "\n" << "\n";
  fichierDiff << "plot(x,y,type = \"l\")";
  
  fichierDiff.close();
  modele.rho() = 0.5;
  
  modele.X0() = 100;
  
  std::cout << "Generation du fichier DataError.txt \n";
  
  int len = 10;
  double vectDelta[10] = {1, 1/double(2), 1/double(4), 1/double(6), 1/double(8), 1/double(10), 1/double(12),1/double(16),1/double(24), 1/double(32)};
  
  priceFFHeston = opt->computePrice(&modele,r,t);
  std::cout << "prix FF : " << priceFFHeston << "\n \n";

  std::ofstream fichierBis("dataError.txt"); 
  
  fichierBis << "x=c(";
  for (int i = 0; i < len-1; i++){
    delta = vectDelta[i];
    fichierBis << delta << ",";
  }
  delta = vectDelta[len-1];
  fichierBis << delta << ")";
  
  fichierBis << "\n" << "\n";
  
  std::cout << "Euler\n";
  fichierBis << "yEuler = c(";
  for (int i = 0; i < len-1; i++){
    delta = vectDelta[i];
    std::cout << "Delta = " << delta << "\n";
    MC->computePriceOption(priceHeston, varIc,int(T/delta), opt, &modele, &Heston::assetEuler);

    
    fichierBis << priceFFHeston - priceHeston << ",";
  }
  delta = vectDelta[len-1];
  std::cout << "Delta = " << delta << "\n \n";
  MC->computePriceOption(priceHeston, varIc,int(T/delta), opt, &modele, &Heston::assetEuler);
  
  fichierBis << priceFFHeston - priceHeston << ")";
  fichierBis << "\n" << "\n";
  
  std::cout << "JK \n";
  fichierBis << "yJK = c(";
  for (int i = 0 ; i < len -1; i++){
    delta = vectDelta[i];
    std::cout << "Delta = " << delta << "\n";
    MC->computePriceOption(priceHeston, varIc,int(T/delta), opt, &modele, &Heston::assetJK);
    
    
    fichierBis << priceFFHeston - priceHeston << ",";
    
  }
  delta = vectDelta[len-1];
  std::cout << "Delta = " << delta << "\n \n";
  MC->computePriceOption(priceHeston, varIc,int(T/delta), opt, &modele, &Heston::assetJK);
  
  fichierBis << priceFFHeston - priceHeston << ")";
  fichierBis << "\n" << "\n";
  
  std::cout << "QE\n";
  fichierBis << "yQE = c(";
  
  for (int i = 0; i < len-1; i++){
    delta = vectDelta[i];
    std::cout << "Delta = " << delta << "\n";
    MC->computePriceOption(priceHeston, varIc,int(T/delta), opt, &modele, &Heston::assetQE);
    
    
    fichierBis << priceFFHeston - priceHeston << ",";
    
  }
  delta = vectDelta[len-1];
  std::cout << "Delta = " << delta << "\n \n";
  MC->computePriceOption(priceHeston, varIc,int(T/delta), opt, &modele, &Heston::assetQE);
  
  fichierBis << priceFFHeston - priceHeston << ")";
  fichierBis << "\n" << "\n";
    
  for (int i = 0; i<len-1; i++){
    delta = vectDelta[i];
    std::cout << "Delta = " << delta << "\n";
    MC->computePriceOption(priceHeston, varIc,int(T/delta), opt, &modele, &Heston::assetQE);
    
    
    fichierBis << priceFFHeston - priceHeston << ",";
  }
  delta = vectDelta[len-1];
  std::cout << "Delta = " << delta << "\n \n";
  MC->computePriceOption(priceHeston, varIc,int(T/delta), opt, &modele, &Heston::assetQE);
  
  fichierBis << priceFFHeston - priceHeston << ")";
  fichierBis << "\n" << "\n";

  fichierBis << "y = cbind(yEuler,yJK,yQE) \n";
  fichierBis << "y = abs(y) \n";
  fichierBis << "matplot(x,y,type = c(\"l\",\"l\",\"l\",\"l\"), col = c(\"blue\",\"red\",\"green\",\"magenta\"),lwd=c(1,1,1,1), lty=c(1,1,1,1),xlab = \"Delta\",ylab = \"|e|\") \n";
  fichierBis << "legend(\"topright\",c(\"Euler\",\"JK\",\"QE\"),col = c(\"blue\",\"red\",\"green\",\"magenta\"),lwd=c(1,1,1,1), lty=c(1,1,1,1)) \n";
  
  fichierBis.close();

  return 0;
}
