#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "contiguous_functions.h"
#include <unistd.h>

#define TWOPI 6.2831853071795864769

int main() {
  
  // Number of points
  int M = 64; // Number of points in each circle
  int N = 2*M;
  int n_dim = 2;
  int size = N*n_dim;
  int T = 2000;
  double tol_rk45_time = 1.e-8;
  long double tol_rk45_space = 1.e-8;
  long double h = 1.e-3;
  double alpha = 0.5; // Interpolation between 2D Euler and Quasi-geostrophic
  double theta = -1.0;
  double dt = 1.e-3;//1.*h;
  double F, tpi, time;
  int* pN;
  double** px;
  int zeros;
  int N_old;
  time = 0;
    
  // Allocate coordinates
  double* x = (double*)malloc(size*sizeof(double));
  double* dxdt = (double*)malloc(size*sizeof(double));
  double* dxdt_k1 = (double*)malloc(size*sizeof(double));
  double* dxdt_k2 = (double*)malloc(size*sizeof(double));
  double* dxdt_k3 = (double*)malloc(size*sizeof(double));
  double* dxdt_k4 = (double*)malloc(size*sizeof(double));
  double* dxdt_k5 = (double*)malloc(size*sizeof(double));
  double* dxdt_k6 = (double*)malloc(size*sizeof(double));
  double* dxdt_RK4 = (double*)malloc(size*sizeof(double));
  double* dxdt_RK5 = (double*)malloc(size*sizeof(double));
  double kappa_den[2];  

  px = &x;
  pN = &N;
  zeros = size*sizeof(double);
  memset(dxdt, 0, zeros);
  memset(dxdt_k1, 0, zeros);
  memset(dxdt_k2, 0, zeros);
  memset(dxdt_k3, 0, zeros);
  memset(dxdt_k4, 0, zeros);
  memset(dxdt_k5, 0, zeros);
  memset(dxdt_k6, 0, zeros);
  memset(dxdt_RK4, 0, zeros);
  memset(dxdt_RK5, 0, zeros);
   
  // Generate circle
  for (int j = 0; j < M; j++) {
    x[2*j] = cos(TWOPI*j/(double)M) - 1.01;
    x[2*j+1] = sin(TWOPI*j/(double)M);
    
    x[2*j+2*M] = cos(TWOPI*j/(double)M) + 1.01;
    x[2*j+1+2*M] = sin(TWOPI*j/(double)M);
   // printf("x[%d] = %e, x[%d] = %e, \n", 2*j, x[2*j], 2*j+2*M, x[2*j+2*M] );
  }
  double area1, area2;
 
  // Print to file  
  char str[80] = "../circle_";
  char str2[80] = "";
  sprintf(str2, "%d", 1);
  strcat(str, str2);
  strcat(str, "a.csv");
  FILE* f = fopen(str, "wb");
  for (int i = 0; i < N; i++) {
    fprintf(f, "%lf,%lf\n", x[2*i], x[2*i+1]);
  }
  fclose(f);
  
  // Evolve
  for (int k = 0; k <= T; k++) {
    size = N*n_dim;
    double* d = (double*)malloc(N*sizeof(double));
    double* kappa = (double*)malloc(N*sizeof(double));
    double* mu = (double*)malloc(N*sizeof(double));
    double* beta = (double*)malloc(N*sizeof(double));
    double* gamma = (double*)malloc(N*sizeof(double));
    double* t = (double*)malloc(size*sizeof(double)); 
    double* n = (double*)malloc(size*sizeof(double));

  
  
    // Interpolate
    interpolate(x, 0, M, n_dim, t, n, d, kappa, kappa_den, mu, beta, gamma);
    interpolate(x, M, N, n_dim, t, n, d, kappa, kappa_den, mu, beta, gamma);
    
    // Compute area
    area1 = compute_area(x, 0, M, t, n, mu, beta, gamma);
    area2 = compute_area(x, M, N, t, n, mu, beta, gamma);
    printf("area1 = %lf\n", area1);
    printf("area2 = %lf\n\n", area2);
      
    // Evolve patches
    dt = runge_kutta45(x, dxdt_k1, dxdt_k2, dxdt_k3, dxdt_k4, dxdt_k5,\
                  dxdt_k6, dxdt_RK4, dxdt_RK5, tol_rk45_time, dt, M, N,\
                  mu, beta, gamma, t, n, alpha, tol_rk45_space, h);
    time += dt;
    printf("time = %1.15lf\n", time);
    

    
    printf("--------------------------\n");
    
    //Print to file
    if (k%100 == 0) {
      // Print to file
      char str[80] = "../circle_";
      char str2[80] = "";
      sprintf(str2, "%d", k);
      strcat(str, str2);
      strcat(str, ".csv");
      FILE* f = fopen(str, "wb");
      for (int i = 0; i < N; i++) {
        fprintf(f, "%lf,%lf\n", x[2*i], x[2*i+1]);
      }
      fclose(f);
    }
    
    // Redistribute the nodes

    N_old  = N;
    px = &x;
    //points_reloc(px, t, n, pN, kappa, mu, gamma, beta);
    //printf("N = %d\n", N);
    //printf(" \n");
   	
    free(t);
    free(n);
    free(d);
    free(kappa);
    free(mu);
    free(beta);
    free(gamma); 
  }

  
  // Free memory

  free(x);


  free(dxdt);
  free(dxdt_k1);
  free(dxdt_k2);
  free(dxdt_k3);
  free(dxdt_k4);
  free(dxdt_k5);
  free(dxdt_k6);
  free(dxdt_RK4);
  free(dxdt_RK5);

  
  return 0;
}
