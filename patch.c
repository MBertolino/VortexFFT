#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "functions.h"

#define TWOPI 6.2831853071795864769

int main() {
 // Number of pointssssssssssssssssssssssssssss

  
  // Number of points
  int N = 20; // Points
  int P = 10;  // Interpolation points
  int n_dim = 2;
  int T = 100;
  double eps = 0.1;
  double h = 0.1;
  double alpha = 0.5; // Interpolation between 2D Euler and Quasi-geostrophic
  /*
    
    More to look at is declaring variables in the RK-schemes.
    
    Next: general clean up of the code, adding some more comments explaining and clearifying steps.
    Evolve the system and check correctness of solutions.
    
  */
  
  double* p = (double*)malloc(P*sizeof(double));
  double* eta = (double*)malloc(P*sizeof(double));
  double** t = (double**)malloc(N*sizeof(double*)); 
  double** n = (double**)malloc(N*sizeof(double*));
  double* mu = (double*)malloc(N*sizeof(double*));
  double* beta = (double*)malloc(N*sizeof(double*));
  double* gamma = (double*)malloc(N*sizeof(double*));
  for (int i = 0; i < N; i++) {
    t[i] = (double*)malloc(n_dim*sizeof(double));
    n[i] = (double*)malloc(n_dim*sizeof(double));
  }
  double* d = (double*)malloc(N*sizeof(double));
  double* kappa = (double*)malloc(N*sizeof(double));
  double kappa_den[2];

  
  // Allocate coordinates
  double** x = (double**)malloc(N*P*sizeof(double*));
  double** dxdt = (double**)malloc(N*P*sizeof(double*));
  
  for (int i = 0; i < N*P; i++) {
    x[i] = (double*)malloc(n_dim*sizeof(double));
    dxdt[i] = (double*)malloc(n_dim*sizeof(double));
  }

  // Generate circle. (j*P) to avoid the interpolated nodes)
  for (int j = 0; j < N; j++) {
    x[j*P][0] = cos(TWOPI*j/(double)N);
    x[j*P][1] = sin(TWOPI*j/(double)N);
  }
  
  // Step in time
  for (int dt = 0; dt < T; dt++) {
    // Interpolate
    interpolate(x, N, P, n_dim, t, n, p, eta, d, kappa, kappa_den, mu, beta, gamma);
    printf("dt = %d\n", dt);
     
    // Calculate derivatives
    for (int j = 0; j < N*P; j++)
      compute_derivative(dxdt[j], x, mu, beta, gamma, t, n, N*P, alpha, h, eps, j);
    
    printf("dxdt = %lf\n", dxdt[1][0]);
    // Time integrate with RK4
    
    // Redistribute the nodes
    // points_reloc();
  }
  
  /*
  // Print to file
  char str[80] = "../circle.csv";
  FILE* f = fopen(str, "wb");
  for (int i = 0; i < N*P; i++) {
    fprintf(f, "%lf,%lf\n", x[i], y[i]);
  }
  fclose(f);
  */
  
  // Free memory
  for (int i = 0; i < N*P; i++)
    free(x[i]);
  free(x);
  for (int i = 0; i < N; i++) {
    free(t[i]);
    free(n[i]);
  }
  free(t);
  free(n);
  free(d);
  free(p);
  free(eta);
  free(kappa);
  
  return 0;
}
