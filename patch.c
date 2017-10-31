#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TWOPI 6.2831853071795864769

int main() {
  int N = 10; // Number of points
  
  // Cubic interpolation params
  int P = 100; // Number of interpolation points
  
  double* p = (double*)malloc(P*sizeof(double));
  double* eta = (double*)malloc(P*sizeof(double));
  for (int i = 0; i < P; i++)
    p[i] = (double)i/P;
  double* t_x = (double*)malloc(N*sizeof(double));
  double* t_y = (double*)malloc(N*sizeof(double));
  double* n_x = (double*)malloc(N*sizeof(double));
  double* n_y = (double*)malloc(N*sizeof(double));
  double* d = (double*)malloc(N*sizeof(double));
  double* kappa = (double*)malloc(N*sizeof(double));
  double kappa_den_x, kappa_den_y;
  double mu, beta, gamma;
  
  // Coordinates
  double* x = (double*)malloc(P*N*sizeof(double));
  double* y = (double*)malloc(P*N*sizeof(double));
  
  // Generate circle (j*P to avoid the interpolated nodes)
  for (int j = 0; j < N; j++) {
    x[j*P] = cos(TWOPI*j/(double)N);
    y[j*P] = sin(TWOPI*j/(double)N);
  }
  
  // Calculate t an n
  for (int j = 0; j < N-1; j++) {
    t_x[j] = x[(j+1)*P] - x[j*P];
    t_y[j] = y[(j+1)*P] - y[j*P];
    n_x[j] = -t_y[j];
    n_y[j] = t_x[j];
    d[j] = sqrt(t_x[j]*t_x[j] + t_y[j]*t_y[j]);
  }  
  // Special case j = N-1
  t_x[N-1] = x[0] - x[(N-1)*P];
  t_y[N-1] = y[0] - y[(N-1)*P];
  n_x[N-1] = -t_y[N-1];
  n_y[N-1] = t_x[N-1];
  d[N-1] = sqrt(t_x[N-1]*t_x[N-1] + t_y[N-1]*t_y[N-1]);
  
  // kappa local curvature
  kappa_den_x = (d[N-1]*d[N-1]*t_x[0] + d[0]*d[0]*t_x[N-1]);
  kappa_den_y = (d[N-1]*d[N-1]*t_y[0] + d[0]*d[0]*t_y[N-1]);
  kappa[0] = 2*(t_x[N-1]*t_y[0] - t_y[N-1]*t_x[0])\
    /sqrt(kappa_den_x*kappa_den_x + kappa_den_y*kappa_den_y);
  
  for (int j = 1; j < N; j++) {
    // kappa local curvature
    kappa_den_x = (d[j-1]*d[j-1]*t_x[j] + d[j]*d[j]*t_x[j-1]);
    kappa_den_y = (d[j-1]*d[j-1]*t_y[j] + d[j]*d[j]*t_y[j-1]);
    kappa[j] = 2*(t_x[j-1]*t_y[j] - t_y[j-1]*t_x[j])\
      /sqrt(kappa_den_x*kappa_den_x + kappa_den_y*kappa_den_y);
  }
  
  // Construct the cubic interpolation
  for (int j = 0; j < N; j++) {
    
    // Cubic interpolation coefficients
    mu = -(double)1/3*d[j]*kappa[j] - (double)1/6*d[j]*kappa[j+1];
    beta = 0.5*d[j]*kappa[j];
    gamma = (double)1/6*d[j]*(kappa[j+1] - kappa[j]);
    
    for (int i = 0; i < P; i++) {
      eta[i] = mu*p[i] + beta*p[i]*p[i] + gamma*p[i]*p[i]*p[i];
      x[j*P + i] = x[j*P] + p[i]*t_x[j] + eta[i]*n_x[j];
      y[j*P + i] = y[j*P] + p[i]*t_y[j] + eta[i]*n_y[j];
    }
  }
  
  // Print to file
  char str[80] = "../circle.csv";
  FILE* f = fopen(str, "wb");
  for (int i = 0; i < N*P; i++) {
    fprintf(f, "%lf,%lf\n", x[i], y[i]);
  }
  fclose(f);
  
  // Free memory
  free(t_x);
  free(t_y);
  free(n_x);
  free(n_y);
  free(d);
  free(p);
  free(eta);
  free(kappa);
  
  free(x);
  free(y);

  return 0;
}