#include <iostream>
#include "MC.h"





MC::MC(){
  M_ = 0;
}

MC::MC(int M){
  M_ = M;
}

MC::~MC(){
  
}

int MC::M() const{
  return M_;
}

int & MC::M(){
  return M_;
}

