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
  int T = 50000;
  double tol_rk45_time = 1.e-8;
  long double tol_rk45_space = 1.e-8;
  long double h = 1.e-3;
  double alpha = 0.5; // Interpolation between 2D Euler and Quasi-geostrophic
  double theta = -1.0;
  double dt = 1.e-3;//1.*h;
  double F, tpi, time;
  int* pN;
  double*** px;
  
  int N_old;
  time = 0;
    
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
  double kappa_den[2];  
  for (int i = 0; i < N; i++)
  {
    x[i] = (double*)malloc(n_dim*sizeof(double));
    x_temp[i] = (double*)malloc(n_dim*sizeof(double));
  }
  px = &x;
  pN = &N;
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
    x[j][0] = cos(TWOPI*j/(double)M);// - 1.01;
    x[j][1] = sin(TWOPI*j/(double)M);
    
    //x[j+M][0] = cos(TWOPI*j/(double)M) + 1.01;
    //x[j+M][1] = sin(TWOPI*j/(double)M);
  }
  double area1;
  //double area2;
  area1 = compute_area(x, 0, M);
  //area2 = compute_area(x, M, N);
  //printf("area1 = %lf\n", area1);
  //printf("area2 = %lf\n", area2);

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
    double* d = (double*)malloc(N*sizeof(double));
    double* kappa = (double*)malloc(N*sizeof(double));
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

  
  
  printf("k = %d\n", k);
    M = N;//M = (int)N/2;
    // Interpolate
    interpolate(x, 0, M, n_dim, t, n, d, kappa, kappa_den, mu, beta, gamma);
    //interpolate(x, M, N, n_dim, t, n, d, kappa, kappa_den, mu, beta, gamma);
    for (int j = 0; j < N; j++)
    {
      x_temp[j][0] = 0.;
      x_temp[j][1] = 0.;
    }
    
    
      
    // Evolve patches
    dt = runge_kutta45(x, dxdt, dxdt_k1, dxdt_k2, dxdt_k3, dxdt_k4, dxdt_k5,\
                  dxdt_k6, dxdt_RK4, dxdt_RK5, tol_rk45_time, 2*dt, M, N,\
                  mu, beta, gamma, t, n, alpha, tol_rk45_space, h);
    time += dt;
    printf("time = %1.15lf\n", time);
    
    // Compute area
    area1 = compute_area(x, 0, M);
    //area2 = compute_area(x, M, N);
    printf("area1 = %lf\n", area1);
    //printf("area2 = %lf\n\n", area2);
    printf("--------------------------\n");
    
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

    N_old  = N;
   // points_reloc(px, t, n, pN, kappa, mu, gamma, beta);
    printf("N = %d\n", N);
    printf(" \n");
    
    
/*    d = (double*)realloc(d, N*sizeof(double));
    t = (double**)realloc(t, N*sizeof(double*));
   	n = (double**)realloc(n, N*sizeof(double*));
   	kappa = (double*)realloc(kappa, N*sizeof(double));
   	mu = (double*)realloc(mu, N*sizeof(double));
    beta = (double*)realloc(beta, N*sizeof(double));
   	gamma = (double*)realloc(gamma, N*sizeof(double));*/
   /*	for (int i = 0; i < N; i++)
   	{
   	  printf("x[%d][0] = %e \n", i, x[i][0]);	 		
   	}*/
   	
    for (int i = 0; i < N_old; i++)
   	{
   	  free(t[i]);
      free(n[i]);	 		
   	}  
    free(t);
    free(n);
    free(d);
    free(kappa);
    free(mu);
    free(beta);
    free(gamma); 

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

  
  return 0;
}
