#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "functions.h"

#define TWOPI 6.2831853071795864769

int main() {
 // Number of pointssssssssssssssssssssssssssss

  
  // Number of points
  int N = 50; // Points
  int P = 100;  // Interpolation points
  int n_dim = 2;
  int T = 100;
  
  // Interpolation between 2D Euler and Quasi-geostrophic
  double alpha = 0.5;
  
  // Allocate coordinates
  double** x = (double**)malloc(N*P*sizeof(double*));
  double** dxdt = (double**)malloc(N*P*sizeof(double*));
  for (int i = 0; i < N*P; i++) {
    x[i] = (double*)malloc(n_dim*sizeof(double));
    dxdt[i] = (double*)malloc(n_dim*sizeof(double));
  }
  double derivative;
  
  // Generate circle. (j*P) to avoid the interpolated nodes)
  for (int j = 0; j < N; j++) {
    x[j*P][0] = cos(TWOPI*j/(double)N);
    x[j*P][1] = sin(TWOPI*j/(double)N);
  }
  
  // Step in time
  for (int t = 0; t < T; t++) {
    
    // Send in mu, beta, gamma, eta, t and n etc to interpolate
    
    // Interpolate
    interpolate(x, N, P, n_dim);
    
    // Calculate derivatives
    for (int j = 0; j < N*P; j++) {
      // Now we multiply both integrals by (t_x[i] + mu[i]*n_x[i]) which is wrong!! See formula (28)
      derivative = compute_derivative(x, mu, beta, gamma, t, n, N*P, alpha);
    }
    
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
  
  return 0;
}
