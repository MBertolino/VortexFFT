#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "functions.h"

#define TWOPI 6.2831853071795864769

int main() {
  
  // Number of points
  int N = 450; // Points
  int n_dim = 2;
  int T = 2;
  long double eps = 0.00001;
  long double h = 0.001;
  double alpha = 0.5; // Interpolation between 2D Euler and Quasi-geostrophic
  double theta = -1.0;
  double dt = 0.1;

  // Allocate coordinates
  double** x = (double**)malloc(N*sizeof(double*));
  double** x_temp = (double**)malloc(N*sizeof(double*));
  double** dxdt_k1 = (double**)malloc(N*sizeof(double*));
  double** dxdt_k2 = (double**)malloc(N*sizeof(double*));
  double** dxdt_k3 = (double**)malloc(N*sizeof(double*));
  double** dxdt_k4 = (double**)malloc(N*sizeof(double*));
  
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
    dxdt_k1[j] = (double*)malloc(n_dim*sizeof(double));
    dxdt_k2[j] = (double*)malloc(n_dim*sizeof(double));
    dxdt_k3[j] = (double*)malloc(n_dim*sizeof(double));
    dxdt_k4[j] = (double*)malloc(n_dim*sizeof(double));
  }
  
  // Generate circle. (j*P) to avoid the interpolated nodes)
  for (int j = 0; j < N; j++) {
    x[j][0] = cos(TWOPI*j/(double)N);
    x[j][1] = sin(TWOPI*j/(double)N);
  }
  
  // Interpolates
  interpolate(x, N, n_dim, t, n, d, kappa, kappa_den, mu, beta, gamma);
  
  printf("mu[0] = %lf\n", mu[0]);
  
  // Print to file  
  char str[80] = "../circle_";
  char str2[80] = "";
  sprintf(str2, "%d", 1);
  strcat(str, str2);
  strcat(str, ".csv");
  FILE* f = fopen(str, "wb");
  for (int i = 0; i < N; i++) {
    fprintf(f, "%lf,%lf\n", x[i][0], x[i][1]);
  }
  fclose(f);
  
  // Step in time
  T = 1;
  for (int k = 0; k < T; k++) {
    printf("k = %d\n", k);
    
    for (int j = 0; j < N; j++)
    {
      //printf("j = %d\n", j);
      compute_derivative(dxdt_k1[j], x, mu, beta, gamma, t, n, N, alpha, h, eps, j);
      dxdt_k1[j][0] = dxdt_k1[j][0]*theta/(TWOPI);
      dxdt_k1[j][1] = dxdt_k1[j][1]*theta/(TWOPI);
      //printf("\n");
    }
    printf("dxdt[0][0] = %lf\n", dxdt_k1[0][0]);
    printf("dxdt[0][1] = %lf\n", dxdt_k1[0][1]);
    
    /*
    // Runge-Kutta
    // Step 1 in RK
    for (int j = 0; j < N; j++)
    {
      compute_derivative(dxdt_k1[j], x, mu, beta, gamma, t, n, N, P, alpha, h, eps, j);
      dxdt_k1[j][0] = dt*dxdt_k1[j][0]*theta/(TWOPI);
      dxdt_k1[j][1] = dt*dxdt_k1[j][1]*theta/(TWOPI);
    }
    for (int j = 0; j < N; j++)
    {
      x_temp[j][0] = x[j][0] + 0.5*dxdt_k1[j][0];
      x_temp[j][1] = x[j][1] + 0.5*dxdt_k1[j][1];
    }
    
    // Step 2 in RK
    for (int j = 0; j < N; j++)
    {      
      compute_derivative(dxdt_k2[j], x_temp, mu, beta, gamma, t, n, N, P, alpha, h, eps, j);
      printf("test\n");
      dxdt_k2[j][0] = dt*dxdt_k2[j][0]*theta/(TWOPI);
      dxdt_k2[j][1] = dt*dxdt_k2[j][1]*theta/(TWOPI);
    }
    for (int j = 0; j < N; j++)
    {
      x_temp[j][0] = x[j][0] + 0.5*dxdt_k2[j][0];
      x_temp[j][1] = x[j][1] + 0.5*dxdt_k2[j][1];
    }
    
    // Step 3 in RK
    for (int j = 0; j < N; j++)
    {
      compute_derivative(dxdt_k3[j], x_temp, mu, beta, gamma, t, n, N, P, alpha, h, eps, j);
      dxdt_k3[j][0] = dt*dxdt_k3[j][0]*theta/(TWOPI);
      dxdt_k3[j][1] = dt*dxdt_k3[j][1]*theta/(TWOPI);
    }
    for (int j = 0; j < N; j++)
    {
      x_temp[j][0] = x[j][0] + dxdt_k2[j][0];
      x_temp[j][1] = x[j][1] + dxdt_k2[j][1];
    }
    
    // Step 4 in RK
    for (int j = 0; j < N; j++)
    {
      compute_derivative(dxdt_k4[j], x_temp, mu, beta, gamma, t, n, N, P, alpha, h, eps, j);
      dxdt_k4[j][0] = dxdt_k4[j][0]*theta/(TWOPI);
      dxdt_k4[j][1] = dxdt_k4[j][1]*theta/(TWOPI);
    }
    
    for (int j = 0; j < N; j++)
    {
      x[j][0] = x[j][0] + (dxdt_k1[j][0] + 2*dxdt_k2[j][0] + 2*dxdt_k3[j][0] + dxdt_k4[j][0])/6;
      x[j][1] = x[j][1] + (dxdt_k1[j][1] + 2*dxdt_k2[j][1] + 2*dxdt_k3[j][1] + dxdt_k4[j][1])/6;
    }
    /*
    
    /*
    // Print to file  
    char str[80] = "../circle_";
    char str2[80] = "";
    sprintf(str2, "%d", k);
    strcat(str, str2);
    strcat(str, ".csv");
    FILE* f = fopen(str, "wb");
    for (int i = 0; i < N*P; i++) {
      fprintf(f, "%lf,%lf\n", x[i][0], x[i][1]);
    }
    fclose(f);
    */
    
    //intf("\n");
     //for (int j = 0; j < N; j++)
      //printf("dxdt_k1[%d][%d] = %lf\n", j, 1, dxdt_k1[j][1]);
    //intf("dxdt_k1[0][0] = %lf\n", dxdt_k1[0][0]);
    //intf("dxdt_k1[0][1] = %lf\n", dxdt_k1[0][1]);
    // Time integrate with RK4
    
    
    // Redistribute the nodes
    // points_reloc();
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
    free(dxdt_k1[i]);
    free(dxdt_k2[i]);
    free(dxdt_k3[i]);
    free(dxdt_k4[i]);
    free(t[i]);
    free(n[i]);
  }
  free(dxdt_k1);
  free(dxdt_k2);
  free(dxdt_k3);
  free(dxdt_k4);
  free(t);
  free(n);
  free(d);
  free(kappa);
  free(mu);
  free(beta);
  free(gamma);
  
  return 0;
}
