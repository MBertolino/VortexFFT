#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "functions.h"

#define TWOPI 6.2831853071795864769

int main() {
  
  // Number of points
  int N = 50; // Points
  int P = 100;  // Interpolation points
  int n_dim = 2;
  
  // Interpolation between 2D Euler and Quasi-geostrophic
  alpha = 0.5;
  
  // Coordinates
  double** x = (double**)malloc(n_dim*sizeof(double*));
  double** dxdt = (double**)malloc(n_dim*sizeof(double*));
  for (int i = 0; i < n_dim; i++) {
    x[i] = (double*)malloc(P*N*sizeof(double));
    dxdt[i] = (double*)malloc(P*N*sizeof(double));
  }
  double derivative;
  
  // Step in time
  for (int t = 0; t < T; t++) {
    
    // Send in mu, beta, gamma, eta, t and n etc to interpolate
    
    // Interpolate
    interpolate(x, N, P);
    
    // Calculate derivatives
    for (int j = 0; j < N*P; j++) {
      // Now we multiply both integrals by (t_x[i] + mu[i]*n_x[i]) which is wrong!! See formula (28)
      derivative = compute_derivative(x, p, mu, beta, gamma, t, n, N*P, alpha);
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
  free(x);
  free(y);

  return 0;
}