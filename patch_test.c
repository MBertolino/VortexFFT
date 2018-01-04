#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "orig_functions.h"
#include "fft_functions.h"
#include "misc.h"
#include <unistd.h>

#define TWOPI 6.2831853071795864769


int main(int argc, char **argv) {
  
  // Patch parameters
  double alpha, theta;
  alpha = 0.7; // Interpolation between 2D Euler and Quasi-geostrophic
  theta = -1.0;
  
  // Runge-kutta parameters
  long double tol_rk45_time, tol_rk45_space, h;
  tol_rk45_time = 1.e-8;
  tol_rk45_space = 1.e-10;
  h = 1.e-4;
  
  // Interpolation and vector field pointers
  double *d, *kappa, *mu, *beta, *gamma, *t, *n, *norm;
  double *x, *k1, *k2, *k3, *k4, *k5, *k6;
 
  for (int i = 2; i < 1024*4096; i *=2)
  {
    printf("i = %d\n", i);
    // Generate patch
    x = (double*)malloc(2*i*sizeof(double));
    for (int j = 0; j < i; j++)
    {
      x[2*j] = cos(TWOPI*j/(double)i); //cos(TWOPI*j/(double)M) + 0.45*sin(TWOPI*5*j/(double)M); // 
      x[2*j+1] = sin(TWOPI*j/(double)i); //sin(TWOPI*j/(double)M) + 0.3*cos(TWOPI*3*j/(double)M);
    }
    // Allocate
    allocate(&d, &kappa, &mu, &beta, &gamma, &t, &n, &norm, &k1, &k2, &k3, &k4, &k5, &k6, i);
    
    // Interpolate
    interpolate(x, 0, i, t, n, d, kappa, mu, beta, gamma);
    
    // Compare algorithms
    //compare_algo(x, mu, beta, gamma, t, n, i, i, h, tol_rk45_space, alpha, theta);
    compare_algo_time(x, mu, beta, gamma, t, n, i, i, h, tol_rk45_space, alpha, theta);
   
    // Free memory
    free_step(d, kappa, mu, beta, gamma, t, n, norm, k1, k2, k3, k4, k5, k6);
    free(x);
  }

  return 0;
}
