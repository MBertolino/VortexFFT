#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TWOPI 6.2831853071795864769

int main() {
  int N = 4; // Number of points
  
  // Cubic interpolation params
  int P = 10; // Number of interpolation points
  double* p = (double*)malloc(P*sizeof(double));
  double* eta = (double*)malloc((N-1)*sizeof(double));
  for (int i = 0; i < P; i++)
    p[i] = (double)i/P;
  double* t_x = (double*)malloc((N-1)*sizeof(double));
  double* t_y = (double*)malloc((N-1)*sizeof(double));
  double* n_x = (double*)malloc((N-1)*sizeof(double));
  double* n_y = (double*)malloc((N-1)*sizeof(double));
  double* mu = (double*)malloc((N-1)*sizeof(double));
  
  // Coordinates
  double* x = (double*)malloc(P*N*sizeof(double));
  double* y = (double*)malloc(P*N*sizeof(double));
  
  // Generate circle (j*P to avoid the interpolated nodes)
  for (int j = 0; j < N; j++) {
    x[j*P] = cos(TWOPI*j/N);
    y[j*P] = sin(TWOPI*j/N);
  } 
  
  // Construct the cubic interpolation
 
    
  for (int j = 0; j < N-1; j++) {
    t_x[j] = x[j+1] - x[j];
    t_y[j] = y[j+1] - y[j];
    n_x[j] = -t_y[j];
    n_y[j] = t_x[j];
    d = sqrt((x[(j+1)*p] - x[j*p])*(x[(j+1)*p] - x[j*p]) + (y[(j+1)*p] - y[j*p])*(y[(j+1)*p] - y[j*p]))
    kappa_x[j] = 2*(t_x[j-1]*t_y[j] - t_y[j-1]*t_x[j])/abs(d[j-1]*d[j-1]*t
    
    
      for (int i = 0; i < P; i++) {
        eta[i] = mu[j]*p[i] + beta[j]*p[i]*p[i] + gamma[j]*p[i]*p[i]*p[i];
        x[j*P + i] = x[j*P] + p[i]*t_x[j] + eta_x[i]*n[j]
      }
    
  }
  
  
  // Free memory
  free(x);
  free(y);

  return 0;
}