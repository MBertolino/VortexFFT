#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "functions.h"

#define TWOPI 6.2831853071795864769

int main() {
  
  // Number of points
  int N = 4500*2; // Points
  int P = 2;  // Interpolation points
  int n_dim = 2;
  int T = 100;
  double eps = 0.1;
  double h = 0.1;
  double alpha = 0.5; // Interpolation between 2D Euler and Quasi-geostrophic
  double theta = -1.0;

  // Allocate coordinates
  double** x = (double**)malloc(N*P*sizeof(double*));
  double** dxdt = (double**)malloc(N*sizeof(double*));
  
  double* p = (double*)malloc(P*sizeof(double));
  double* eta = (double*)malloc(P*sizeof(double));
  double* d = (double*)malloc(N*sizeof(double));
  double* kappa = (double*)malloc(N*sizeof(double));
  double kappa_den[2];  
  double* mu = (double*)malloc(N*sizeof(double*));
  double* beta = (double*)malloc(N*sizeof(double*));
  double* gamma = (double*)malloc(N*sizeof(double*));
  double** t = (double**)malloc(N*sizeof(double*)); 
  double** n = (double**)malloc(N*sizeof(double*));
  for (int i = 0; i < N; i++)
  {
    t[i] = (double*)malloc(n_dim*sizeof(double));
    n[i] = (double*)malloc(n_dim*sizeof(double));
  }
  double* mu_loc = (double*)malloc(N*P*sizeof(double));
  double* beta_loc = (double*)malloc(N*P*sizeof(double));
  double* gamma_loc = (double*)malloc(N*P*sizeof(double));
  double** t_loc = (double**)malloc(N*P*sizeof(double*));
  double** n_loc = (double**)malloc(N*P*sizeof(double*));
  for (int i = 0; i < N*P; i++)
  {
    x[i] = (double*)malloc(n_dim*sizeof(double));
    t_loc[i] = (double*)malloc(n_dim*sizeof(double));
    n_loc[i] = (double*)malloc(n_dim*sizeof(double));
  }
  
  for (int j = 0; j < N; j++)
    dxdt[j] = (double*)malloc(n_dim*sizeof(double));
  
  // Generate circle. (j*P) to avoid the interpolated nodes)
  for (int j = 0; j < N; j++) {
    x[j*P][0] = cos(TWOPI*j/(double)N);
    x[j*P][1] = sin(TWOPI*j/(double)N);
  }
  
    
  // Step in time
  T = 1;
  for (int dt = 0; dt < T; dt++) {
    // Interpolate
    interpolate(x, N, P, n_dim, t, n, p, eta, d, kappa, kappa_den, mu, beta, gamma);
    printf("dt = %d\n", dt);
    local_coeffs(N*P, x, t_loc, n_loc, mu_loc, beta_loc, gamma_loc);
    
    // Calculate derivatives
    for (int j = 0; j < N; j++)
    {
      compute_derivative(dxdt[j], x, mu_loc, beta_loc, gamma_loc, t_loc, n_loc, N, P, alpha, h, eps, j);
      dxdt[j][0] = dxdt[j][0]*theta/(TWOPI);
      dxdt[j][1] = dxdt[j][1]*theta/(TWOPI);
      //printf("dxdt[%d][%d] = %lf\n", j, 0, dxdt[j][0]);
    }
    printf("\n");
     //for (int j = 0; j < N; j++)
      //printf("dxdt[%d][%d] = %lf\n", j, 1, dxdt[j][1]);
    printf("dxdt[0][0] = %lf\n", dxdt[0][0]);
    printf("dxdt[0][1] = %lf\n", dxdt[0][1]);
    // Time integrate with RK4
    
    
    // Redistribute the nodes
    // points_reloc();
  }
  
  /*
  // Print to file
  char str[80] = "../circle.csv";
  FILE* f = fopen(str, "wb");
  for (int i = 0; i < N*P; i++) {
    fprintf(f, "%lf,%lf\n", x[i][0], x[i][1]);
  }
  fclose(f);
  */
  
  // Free memory
  for (int i = 0; i < N*P; i++)
  { 
    free(x[i]);
    free(t_loc[i]);
    free(n_loc[i]);
  }
  free(x);
  free(t_loc);
  free(n_loc);
  for (int i = 0; i < N; i++) {
    free(dxdt[i]);
    free(t[i]);
    free(n[i]);
  }
  free(dxdt);
  free(t);
  free(n);
  free(d);
  free(p);
  free(eta);
  free(kappa);
  free(mu);
  free(beta);
  free(gamma);
  
  free(mu_loc);
  free(beta_loc);
  free(gamma_loc);

  
  return 0;
}
