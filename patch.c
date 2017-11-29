#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "functions.h"
#include <unistd.h>

#define TWOPI 6.2831853071795864769

int main() {
  
  // Number of points
  int M = 256; // Number of points in each circle
  int N = M;
  int n_dim = 2;
  int size = N*n_dim;
  int T = 2000;
  double tol_rk45_time = 1.e-8;
  long double tol_rk45_space = 1.e-8;
  long double h = 1.e-3;
  double alpha = 0.2; // Interpolation between 2D Euler and Quasi-geostrophic
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
  double* dxdt_fft = (double*)malloc(size*sizeof(double));
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
  
 // memset(dxdt, 0, zeros);
  //memset(dxdt_fft, 0, zeros);
  /*memset(dxdt_k1, 0, zeros);
  memset(dxdt_k2, 0, zeros);
  memset(dxdt_k3, 0, zeros);
  memset(dxdt_k4, 0, zeros);
  memset(dxdt_k5, 0, zeros);
  memset(dxdt_k6, 0, zeros);
  memset(dxdt_RK4, 0, zeros);
  memset(dxdt_RK5, 0, zeros);
  */
   
  // Generate circle
  for (int j = 0; j < M; j++) {
    x[2*j] = 5*cos(TWOPI*j/(double)M)+0.2*cos(6*TWOPI*j/(double)M);// - 1.1;
    x[2*j+1] = sin(TWOPI*j/(double)M)-0.2*sin(4*TWOPI*j/(double)M);
    
    //x[2*j+2*M] = cos(TWOPI*j/(double)M) + 1.1;
    //x[2*j+1+2*M] = sin(TWOPI*j/(double)M);
    // printf("x[%d] = %e, x[%d] = %e, \n", 2*j, x[2*j], 2*j+2*M, x[2*j+2*M] );
  }
  double area1, area2;
 
  // Print to file  
  char str[80] = "../circle_";
  char str2[80] = "";
  sprintf(str2, "%d", 1);
  strcat(str, str2);
  strcat(str, "a.txt");
  FILE* f = fopen(str, "wb");
  for (int i = 0; i < N; i++) {
    fprintf(f, "%lf,%lf\n", x[2*i], x[2*i+1]);
  }
  fclose(f);
  
  // Evolve
  for (int k = 0; k <= T; k++) {
    double* d = (double*)malloc(N*sizeof(double));
    double* kappa = (double*)malloc(N*sizeof(double));
    double* mu = (double*)malloc(N*sizeof(double));
    double* beta = (double*)malloc(N*sizeof(double));
    double* gamma = (double*)malloc(N*sizeof(double));
    double* t = (double*)malloc(size*sizeof(double)); 
    double* n = (double*)malloc(size*sizeof(double));
    //double* dxdt = (double*)malloc(size*sizeof(double));
    //double* dxdt_fft = (double*)malloc(N*sizeof(double));
    double* dxdt_k1 = (double*)malloc(size*sizeof(double));
    double* dxdt_k2 = (double*)malloc(size*sizeof(double));
    double* dxdt_k3 = (double*)malloc(size*sizeof(double));
    double* dxdt_k4 = (double*)malloc(size*sizeof(double));
    double* dxdt_k5 = (double*)malloc(size*sizeof(double));
    double* dxdt_k6 = (double*)malloc(size*sizeof(double));
    double* dxdt_RK4 = (double*)malloc(size*sizeof(double));
    double* dxdt_RK5 = (double*)malloc(size*sizeof(double));
    
    zeros = size*sizeof(double);
    memset(dxdt_k1, 0, zeros);
    memset(dxdt_k2, 0, zeros);
    memset(dxdt_k3, 0, zeros);
    memset(dxdt_k4, 0, zeros);
    memset(dxdt_k5, 0, zeros);
    memset(dxdt_k6, 0, zeros);
    memset(dxdt_RK4, 0, zeros);
    memset(dxdt_RK5, 0, zeros);
    
    printf("test1\n");
    // Interpolate
    interpolate(x, 0, M, n_dim, t, n, d, kappa, kappa_den, mu, beta, gamma);
    //interpolate(x, M, N, n_dim, t, n, d, kappa, kappa_den, mu, beta, gamma);
    
    /*
    // Compare FFT and Mancho
    for (int j = 0; j < N; j+=2)
    {
      compute_fft(&dxdt_fft[j], &dxdt_fft[j+1], x, N, alpha, j);
      compute_derivative(&dxdt[j], &dxdt[j+1], x, mu, beta, gamma, t, n, M, N, alpha, h, tol_rk45_space, j);
    }
    
    for (int j = 0; j < 8; j++)
    {
      printf("dxdt_fft[%d] = %lf\n", j, dxdt_fft[j]);
      printf("dxdt_ama[%d] = %lf\n\n", j, dxdt[j]);
    }
    printf("Done\n");
    sleep(5);
    */
      
    // Evolve patches
    dt = runge_kutta45(x, dxdt_k1, dxdt_k2, dxdt_k3, dxdt_k4, dxdt_k5,\
                  dxdt_k6, dxdt_RK4, dxdt_RK5, tol_rk45_time, dt, M, N,\
                  mu, beta, gamma, t, n, alpha, tol_rk45_space, h, &time);
 
    printf("time = %1.15lf\n", time);    
    printf("--------------------------\n");
    
    //Print to file
    if (k%1 == 0) {
      // Print to file
      char str[80] = "../circle_";
      char str2[80] = "";
      sprintf(str2, "%d", k);
      strcat(str, str2);
      strcat(str, ".txt");
      FILE* f = fopen(str, "wb");
      for (int i = 0; i < N; i++) {
        fprintf(f, "%lf,%lf\n", x[2*i], x[2*i+1]);
      }
      fclose(f);
    }
    
    // Compute area
    area1 = compute_area(x, 0, M, t, n, mu, beta, gamma);
    // area2 = compute_area(x, M, N, t, n, mu, beta, gamma);
    printf("area1 = %lf\n", area1);
    // printf("area2 = %lf\n\n", area2);

    // Interpolate
    interpolate(x, 0, M, n_dim, t, n, d, kappa, kappa_den, mu, beta, gamma);
    N_old  = N;
    px = &x;
    
    // Redistribute the nodes
    points_reloc(px, t, n, pN, kappa, mu, gamma, beta);
    printf("N = %d, N_old = %d \n", N, N_old);
    printf(" \n");
    size = N*n_dim;
    M = N;
   	
    // Free memory
    free(t);
    free(n);
    free(d);
    free(kappa);
    free(mu);
    free(beta);
    free(gamma); 

    //free(dxdt);
    //free(dxdt_fft);
    free(dxdt_k1);
    free(dxdt_k2);
    free(dxdt_k3);
    free(dxdt_k4);
    free(dxdt_k5);
    free(dxdt_k6);
    free(dxdt_RK4);
    free(dxdt_RK5);
  }
  free(x);

  return 0;
}
