#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "functions.h"
#include <unistd.h>

#define TWOPI 6.2831853071795864769

int main() {
  
  // Number of points
  int M = 4; // Number of points in each circle
  int N = 2*M;
  int n_dim = 2;
  int T = 11;
  long double eps = 1.e-8;
  long double h = 1.e-8;
  double alpha = 0.5; // Interpolation between 2D Euler and Quasi-geostrophic
  double theta = -1.0;
  double dt = 1.e-3;//1.*h;
  double F, tpi;
    
  // Allocate coordinates
  double** x = (double**)malloc(N*sizeof(double*));
  double** x_temp = (double**)malloc(N*sizeof(double*));
  double** dxdt = (double**)malloc(N*sizeof(double*));
  double** dxdt_k1 = (double**)malloc(N*sizeof(double*));
  double** dxdt_k2 = (double**)malloc(N*sizeof(double*));
  double** dxdt_k3 = (double**)malloc(N*sizeof(double*));
  double** dxdt_k4 = (double**)malloc(N*sizeof(double*));
  double** dxdt_k5 = (double**)malloc(N*sizeof(double*));
  double** dxdt_k6 = (double**)malloc(N*sizeof(double*));
  double** dxdt_RK4 = (double**)malloc(N*sizeof(double*));
  double** dxdt_RK5 = (double**)malloc(N*sizeof(double*));
  
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
  for (int i = 0; i < N; i++)
  {
    x[i] = (double*)malloc(n_dim*sizeof(double));
    x_temp[i] = (double*)malloc(n_dim*sizeof(double));
  }
  
  for (int j = 0; j < N; j++)
  {
    dxdt[j] = (double*)malloc(n_dim*sizeof(double));
    dxdt[j][0] = 0.;
    dxdt[j][1] = 0.;
    dxdt_k1[j] = (double*)malloc(n_dim*sizeof(double));
    dxdt_k2[j] = (double*)malloc(n_dim*sizeof(double));
    dxdt_k3[j] = (double*)malloc(n_dim*sizeof(double));
    dxdt_k4[j] = (double*)malloc(n_dim*sizeof(double));
    dxdt_k5[j] = (double*)malloc(n_dim*sizeof(double));
    dxdt_k6[j] = (double*)malloc(n_dim*sizeof(double));
    dxdt_RK4[j] = (double*)malloc(n_dim*sizeof(double));
    dxdt_RK5[j] = (double*)malloc(n_dim*sizeof(double));
  }
  
  for (int j = 0; j < N; j++)
    {
      dxdt_k1[j][0] = 0.;
      dxdt_k1[j][1] = 0.;
      dxdt_k2[j][0] = 0.;
      dxdt_k2[j][1] = 0.;
      dxdt_k3[j][0] = 0.;
      dxdt_k3[j][1] = 0.;
      dxdt_k4[j][0] = 0.;
      dxdt_k4[j][1] = 0.;
      dxdt_k5[j][0] = 0.;
      dxdt_k5[j][1] = 0.;
      dxdt_k6[j][0] = 0.;
      dxdt_k6[j][1] = 0.;
      dxdt_RK4[j][0] = 0.;
      dxdt_RK4[j][1] = 0.;
      dxdt_RK5[j][0] = 0.;
      dxdt_RK5[j][1] = 0.;
    }
    
  
  // Generate circle
  for (int j = 0; j < M; j++) {
    x[j][0] = cos(TWOPI*j/(double)M) - 1.1;
    x[j][1] = sin(TWOPI*j/(double)M);
    
    x[j+M][0] = cos(TWOPI*j/(double)M) + 1.1;
    x[j+M][1] = sin(TWOPI*j/(double)M);
  }
  double area1;
  double area2;
  area1 = compute_area(x, 0, M);
  area2 = compute_area(x, M, N);
  printf("area1 = %lf\n", area1);
  printf("area2 = %lf\n", area2);

  // Print to file  
  char str[80] = "../circle_";
  char str2[80] = "";
  sprintf(str2, "%d", 1);
  strcat(str, str2);
  strcat(str, "a.csv");
  FILE* f = fopen(str, "wb");
  for (int i = 0; i < N; i++) {
    fprintf(f, "%lf,%lf\n", x[i][0], x[i][1]);
  }
  fclose(f);
  
  // Step in time
  tpi = theta/(TWOPI);	
	F = dt*tpi;
  for (int k = 0; k < T; k++) {
    if (k % 10 == 0)
      printf("k = %d\n", k);
    
    // Interpolate
    interpolate(x, 0, M, n_dim, t, n, d, kappa, kappa_den, mu, beta, gamma);
    interpolate(x, M, N, n_dim, t, n, d, kappa, kappa_den, mu, beta, gamma);
    for (int j = 0; j < N; j++)
    {
      x_temp[j][0] = 0.;
      x_temp[j][1] = 0.;
    }
    
    /*
    // Runge-Kutta
    // Step 1 in RK
    for (int j = 0; j < N; j++)
    {
      compute_derivative(dxdt_k1[j], x, mu, beta, gamma, t, n, M, N, alpha, h, eps, j);
      dxdt_k1[j][0] = F*dxdt_k1[j][0];
      dxdt_k1[j][1] = F*dxdt_k1[j][1];
      //printf("dxdt dot x = %e\n", dxdt_k1[j][0]*x[j][0] + dxdt_k1[j][1]*x[j][1]);
    }
    for (int j = 0; j < N; j++)
    {
      x_temp[j][0] = x[j][0] + 0.5*dxdt_k1[j][0];
      x_temp[j][1] = x[j][1] + 0.5*dxdt_k1[j][1];
    }
    
    // Step 2 in RK
    for (int j = 0; j < N; j++)
    {
      compute_derivative(dxdt_k2[j], x_temp, mu, beta, gamma, t, n, M, N, alpha, h, eps, j);
      dxdt_k2[j][0] = F*dxdt_k2[j][0];
      dxdt_k2[j][1] = F*dxdt_k2[j][1];
    }
    for (int j = 0; j < N; j++)
    {
      x_temp[j][0] = x[j][0] + 0.5*dxdt_k2[j][0];
      x_temp[j][1] = x[j][1] + 0.5*dxdt_k2[j][1];
    }
    
    // Step 3 in RK
    for (int j = 0; j < N; j++)
    {
      compute_derivative(dxdt_k3[j], x_temp, mu, beta, gamma, t, n, M, N, alpha, h, eps, j);
      dxdt_k3[j][0] = F*dxdt_k3[j][0];
      dxdt_k3[j][1] = F*dxdt_k3[j][1];
    }
    for (int j = 0; j < N; j++)
    {
      x_temp[j][0] = x[j][0] + dxdt_k3[j][0];
      x_temp[j][1] = x[j][1] + dxdt_k3[j][1];
    }
    
    // Step 4 in RK
    for (int j = 0; j < N; j++)
    {
      compute_derivative(dxdt_k4[j], x_temp, mu, beta, gamma, t, n, M, N, alpha, h, eps, j);
      dxdt_k4[j][0] = F*dxdt_k4[j][0];
      dxdt_k4[j][1] = F*dxdt_k4[j][1];
    }
    for (int j = 0; j < N; j++)
    {
      x[j][0] = x[j][0] + (dxdt_k1[j][0] + 2*dxdt_k2[j][0] + 2*dxdt_k3[j][0] + dxdt_k4[j][0])/6;
      x[j][1] = x[j][1] + (dxdt_k1[j][1] + 2*dxdt_k2[j][1] + 2*dxdt_k3[j][1] + dxdt_k4[j][1])/6;
    }
    */
      
    // Compute area
    area1 = compute_area(x, 0, M);
    area2 = compute_area(x, M, N);
    printf("area1 = %lf\n", area1);
    printf("area2 = %lf\n", area2);
    
    //Print to file
    if (k%1 == 0) {
      // Print to file
      char str[80] = "../circle_";
      char str2[80] = "";
      sprintf(str2, "%d", k);
      strcat(str, str2);
      strcat(str, ".csv");
      FILE* f = fopen(str, "wb");
      for (int i = 0; i < N; i++) {
        fprintf(f, "%lf,%lf\n", x[i][0], x[i][1]);
      }
      fclose(f);
    }
    
    // Redistribute the nodes
    // N = points_reloc(x, t, n, N, kappa, mu, gamma, beta);
  }

  
  // Free memory
  for (int i = 0; i < N; i++)
  { 
    free(x[i]);
    free(x_temp[i]);
  }
  free(x);
  free(x_temp);
  for (int i = 0; i < N; i++) {
    free(dxdt[i]);
    free(dxdt_k1[i]);
    free(dxdt_k2[i]);
    free(dxdt_k3[i]);
    free(dxdt_k4[i]);
    free(dxdt_k5[i]);
    free(dxdt_k6[i]);
    free(dxdt_RK4[i]);
    free(dxdt_RK5[i]);
    free(t[i]);
    free(n[i]);
  }
  free(dxdt);
  free(dxdt_k1);
  free(dxdt_k2);
  free(dxdt_k3);
  free(dxdt_k4);
  free(dxdt_k5);
  free(dxdt_k6);
  free(dxdt_RK4);
  free(dxdt_RK5);
  free(t);
  free(n);
  free(d);
  free(kappa);
  free(mu);
  free(beta);
  free(gamma);
  
  return 0;
}
