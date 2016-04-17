#include <iostream>


// mtrand.cpp, see include file mtrand.h for information

#include "mtrand.h"
// non-inline function definitions and static member definitions cannot
// reside in header file because of the risk of multiple declarations

// initialization of static private members
//unsigned long MTRand_int32::state[n] = {0x0UL};
//int MTRand_int32::p = 0;
//bool MTRand_int32::init = false;

void MTRand_int32::gen_state() { // generate new state vector
  for (int i = 0; i < (n - m); ++i)
    state[i] = state[i + m] ^ twiddle(state[i], state[i + 1]);
  for (int i = n - m; i < (n - 1); ++i)
    state[i] = state[i + m - n] ^ twiddle(state[i], state[i + 1]);
  state[n - 1] = state[m - 1] ^ twiddle(state[n - 1], state[0]);
  p = 0; // reset position
}

void MTRand_int32::seed(unsigned long s) {  // init by 32 bit seed
  state[0] = s & 0xFFFFFFFFUL; // for > 32 bit machines
  for (int i = 1; i < n; ++i) {
    state[i] = 1812433253UL * (state[i - 1] ^ (state[i - 1] >> 30)) + i;
    // see Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier
    // in the previous versions, MSBs of the seed affect only MSBs of the array state
    // 2002/01/09 modified by Makoto Matsumoto
    state[i] &= 0xFFFFFFFFUL; // for > 32 bit machines
  }
  p = n; // force gen_state() to be called for next random number
}

void MTRand_int32::seed(const unsigned long* array, int size) { // init by array
  seed(19650218UL);
  int i = 1, j = 0;
  for (int k = ((n > size) ? n : size); k; --k) {
    state[i] = (state[i] ^ ((state[i - 1] ^ (state[i - 1] >> 30)) * 1664525UL))
    + array[j] + j; // non linear
    state[i] &= 0xFFFFFFFFUL; // for > 32 bit machines
    ++j; j %= size;
    if ((++i) == n) { state[0] = state[n - 1]; i = 1; }
  }
  for (int k = n - 1; k; --k) {
    state[i] = (state[i] ^ ((state[i - 1] ^ (state[i - 1] >> 30)) * 1566083941UL)) - i;
    state[i] &= 0xFFFFFFFFUL; // for > 32 bit machines
    if ((++i) == n) { state[0] = state[n - 1]; i = 1; }
  }
  state[0] = 0x80000000UL; // MSB is 1; assuring non-zero initial array
  p = n; // force gen_state() to be called for next random number
}

double normal(double mean, double sigma,MTRand  *myGen){
  double u = (*myGen)();
  double v = (*myGen)();
  
  double rac = sqrt(-2*log(u))*cos(2*M_PI*v);
  return mean + sqrt(sigma)*rac;
}

/*
 * The inverse standard normal distribution.
 *
 *   Author:      Peter John Acklam <pjacklam@online.no>
 *   URL:         http://home.online.no/~pjacklam
 *
 * This function is based on the MATLAB code from the address above,
 * translated to C, and adapted for our purposes.
 */
#define  A1  (-3.969683028665376e+01)
#define  A2   2.209460984245205e+02
#define  A3  (-2.759285104469687e+02)
#define  A4   1.383577518672690e+02
#define  A5  (-3.066479806614716e+01)
#define  A6   2.506628277459239e+00

#define  B1  (-5.447609879822406e+01)
#define  B2   1.615858368580409e+02
#define  B3  (-1.556989798598866e+02)
#define  B4   6.680131188771972e+01
#define  B5  (-1.328068155288572e+01)

#define  C1  (-7.784894002430293e-03)
#define  C2  (-3.223964580411365e-01)
#define  C3  (-2.400758277161838e+00)
#define  C4  (-2.549732539343734e+00)
#define  C5   4.374664141464968e+00
#define  C6   2.938163982698783e+00

#define  D1   7.784695709041462e-03
#define  D2   3.224671290700398e-01
#define  D3   2.445134137142996e+00
#define  D4   3.754408661907416e+00

#define P_LOW   0.02425
/* P_high = 1 - p_low*/
#define P_HIGH  0.97575

long double normsinv(long double p)
{
  long double x;
  long double q, r, u, e;
  if ((0 < p )  && (p < P_LOW)){
    q = sqrt(-2*log(p));
    x = (((((C1*q+C2)*q+C3)*q+C4)*q+C5)*q+C6) / ((((D1*q+D2)*q+D3)*q+D4)*q+1);
  }
  else{
    if ((P_LOW <= p) && (p <= P_HIGH)){
      q = p - 0.5;
      r = q*q;
      x = (((((A1*r+A2)*r+A3)*r+A4)*r+A5)*r+A6)*q /(((((B1*r+B2)*r+B3)*r+B4)*r+B5)*r+1);
    }
    else{
      if ((P_HIGH < p)&&(p < 1)){
        q = sqrt(-2*log(1-p));
        x = -(((((C1*q+C2)*q+C3)*q+C4)*q+C5)*q+C6) / ((((D1*q+D2)*q+D3)*q+D4)*q+1);
      }
    }
  }
  
   if(( 0 < p)&&(p < 1)){
   e = 0.5 * erfc(-x/sqrt(2)) - p;
   u = e * sqrt(2*M_PI) * exp(x*x/2);
   x = x - u/(1 + x*u/2);
   }
  
  return x;
}


/*double Fnormal01(double x){
  double a = -10;
  
  int n = 100;
  double h = (x-a)/double(n);
  
  double tmp = 0;
  double integrale = (exp(-pow(a,2)/2) + exp(-pow(x,2)/2))/2;
  
  for (int i = 1; i < n; i++){
    tmp = exp(-pow(a + double(i)*h,2)/double(2));
    integrale = integrale + tmp;
  } 
  integrale = h*integrale/sqrt(M_2_PI);
  
  return integrale;
}*/

double Fnormal01(double z)
{
  if (z > 6.0)
    return 1.0; // éviter les valeurs illicites
  if (z < -6.0)
    return 0.0;
  double b1 = 0.31938153;
  double b2 = -0.356563782;
  double b3 = 1.781477937;
  double b4 = -1.821255978;
  double b5 = 1.330274429;
  double p = 0.2316419;
  double c2 = 0.3989423;
  double a=fabs(z);
  double t = 1.0/(1.0+a*p);
  double b = c2*exp((-z)*(z/2.0));
  double n = ((((b5*t+b4)*t+b3)*t+b2)*t+b1)*t;
  n = 1.0-b*n;
  if ( z < 0.0 )
    n = 1.0 - n;
  return n;
};


