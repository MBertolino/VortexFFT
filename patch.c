#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "functions.h"

#define TWOPI 6.2831853071795864769

int main() {
  int N = 50; // Number of points
  
  // Number of interpolation points
  int P = 100;
  
  // Coordinates
  double* x = (double*)malloc(P*N*sizeof(double));
  double* y = (double*)malloc(P*N*sizeof(double));
  
  /*
  // Step in time
  for (int t = 0; t < T; t++) {
    
    // Interpolate
    interpolate(x, y, N, P);
    
    // Evolve the contour integrals
    for (int i = 0; i < N*P; i++) {
      if ((x[j] == x[i]) && (y[j] == y[i])) {
        evolve_integral(p, mu, beta, gamma);
        
      } else if ((x[j] == x[i+1]) && (y[j] == y[i+1])) {
        evolve_integral(1-p, mu + 2*beta + 3*gamma, -beta -3*gamma, gamma);
      
      } else if (sqrt((x[j] - x[i])*(x[j] - x[i]) + (x[j] - x[i])*(x[j] - x[i])) > 100) {
        evolve_integral_g();
      
      } else if (sqrt((x[j] - x[i])*(x[j] - x[i]) + (x[j] - x[i])*(x[j] - x[i])) < 0.01) {
        evolve_integral_g
      }
    }
  }
  */
  
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