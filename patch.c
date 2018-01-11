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
  int M, M2, N;
  double alpha, theta, area1, area2;
  M = atoi(argv[1]); // Number of points in each circle
  M2 = atoi(argv[1]);
  N = M + M2;
  alpha = 0.7; // Interpolation between 2D Euler and Quasi-geostrophic
  theta = -1.0;
  
  // Runge-kutta parameters
  long double tol_rk45_time, tol_rk45_space, h;
  tol_rk45_time = 1.e-8;
  tol_rk45_space = 1.e-10;
  h = 1.e-4;
  
  // Time parameters
  int T;
  double dt, time;
  T = 1000;
  dt = 1.e-3;
  time = 0.0;
  
  // Interpolation and vector field pointers
  double *d, *kappa, *mu, *beta, *gamma, *t, *n, *norm;
  double *k1, *k2, *k3, *k4, *k5, *k6;

  // Generate patch
  double* x = (double*)malloc(2*N*sizeof(double));
  for (int j = 0; j < M; j++) {
    x[2*j] = cos(TWOPI*j/(double)M) + 1.1; //cos(TWOPI*j/(double)M) + 0.45*sin(TWOPI*5*j/(double)M); // 
    x[2*j+1] = sin(TWOPI*j/(double)M); //sin(TWOPI*j/(double)M) + 0.3*cos(TWOPI*3*j/(double)M);
    
    x[2*(j+M)] = cos(TWOPI*j/(double)M) - 1.1;
    x[2*(j+M)+1] = sin(TWOPI*j/(double)M);
  }

  // Evolve the patch
  for (int k = 0; k <= T; k++)
  {
    
    if (k % 1 == 0)
    {
      // Print to file
      print_to_file(x, M, N, k);
      printf(" \n \n--------------------------\n \n");
      printf("k = %d\n", k);
      printf("N = %d\n", N);
      printf("dt = %e\n", dt);
      printf("time = %1.15lf\n", time);
    }

    // Allocate
    allocate(&d, &kappa, &mu, &beta, &gamma, &t, &n, &norm, &k1, &k2, &k3, &k4, &k5, &k6, N);
    
    // Interpolate
    interpolate(x, 0, M, t, n, d, kappa, mu, beta, gamma);
    //interpolate(x, M, N, t, n, d, kappa, mu, beta, gamma);
     
    // Evolve patches
    dt = runge_kutta45(x, k1, k2, k3, k4, k5, k6,\
                       tol_rk45_time, dt, M, N,\
                  mu, beta, gamma, t, n, alpha, tol_rk45_space, h, theta, norm);
    time += dt;
    
    // Compute area
    //area1 = compute_area(x, 0, M, t, n, mu, beta, gamma);
    //area2 = compute_area(x, M, N, t, n, mu, beta, gamma);
    //area1 = area_fft(x, 0, M);
    //area2 = area_fft(x, M, N);
    //printf("area1 = %e, area2 = %e \n", area1, area2); 
    
    // Reallocate the points
    //interpolate(x, 0, M, t, n, d, kappa, mu, beta, gamma);
    //interpolate(x, M, N, t, n, d, kappa, mu, beta, gamma);
    //points_reloc(&x, t, n, &N, kappa, mu, beta, gamma, &M, &M2, 2);
    
    // Free memory
    free_step(d, kappa, mu, beta, gamma, t, n, norm, k1, k2, k3, k4, k5, k6);
  }
  free(x);

  return 0;
}

