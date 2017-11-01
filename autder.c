#include <stdio.h>
#include <stdlib.h>
#include <math.h>
  /*
    Automatic differentiation of a function f(t) = (1/c(t))^alpha
    c_coeff is padded to length of order
  */

void autder(double* f, double* c_coeff, int N_coeff, double alpha, int order)
{
  // Allocate memory for temporary coefficients
  double* a_ = (double*)malloc(order*sizeof(double));
  
  
  a_[0] = 1/c[0];
  
  // calculate temporary coefficients
  for (int n = 1; n <= order; n++)
  {
    for (int j = 1; j <= n; j++)
    {
      
      a_[n] -= c[j]*a_[n-j];
    }
    a_[n] /= c[0];
  }
  
  f[0] = 1;
  
  // calculate tayler coefficients of order, "order" :D
  for (int n = 1; n <= order; n++)
  {
    for (int j = 0; j < n; j++)
    {
      f[n] += (n*alpha - j*(alpha + 1))*a_[n-j]*f[j];
    }
    f[n] /= (n*a_[0]);
  }
  
  free(a_);
  
  return;
}