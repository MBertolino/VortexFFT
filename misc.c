#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "orig_functions.h"
#include "fft_functions.h"
#include <unistd.h>


void print_to_file(double* x, int M, int N, int k)
{
  //Print to file
  char str[80] = "../results/circle_";
  char str2[80] = "";
  sprintf(str2, "%d", k);
  strcat(str, str2);
  strcat(str, ".txt");
  FILE* f = fopen(str, "wb");
  for (int i = 0; i < N; i++)
  {
    if (i == M)
      fprintf(f, "\n");
    fprintf(f, "%lf %lf\n", x[2*i], x[2*i+1]);
  }
  fclose(f);
}


void allocate(double** d, double** kappa, double** mu, double** beta, double** gamma, double** t, double** n, double** norm, double** k1, double** k2, double** k3, double** k4, double** k5, double** k6, int N)
{
  int size, zeros;
  size = 2*N;
  zeros = size*sizeof(double);
  
  // Interpolation
  *d = (double*)malloc(N*sizeof(double));
  *kappa = (double*)malloc(N*sizeof(double));
  *mu = (double*)malloc(N*sizeof(double));
  *beta = (double*)malloc(N*sizeof(double));
  *gamma = (double*)malloc(N*sizeof(double));
  *t = (double*)malloc(size*sizeof(double)); 
  *n = (double*)malloc(size*sizeof(double));
  *norm = (double*)malloc(size*sizeof(double));
  
  // Vector fields
  *k1 = (double*)malloc(size*sizeof(double));
  *k2 = (double*)malloc(size*sizeof(double));
  *k3 = (double*)malloc(size*sizeof(double));
  *k4 = (double*)malloc(size*sizeof(double));
  *k5 = (double*)malloc(size*sizeof(double));
  *k6 = (double*)malloc(size*sizeof(double));
  
  memset(*k1, 0, zeros);
  memset(*k2, 0, zeros);
  memset(*k3, 0, zeros);
  memset(*k4, 0, zeros);
  memset(*k5, 0, zeros);
  memset(*k6, 0, zeros);
  
  return;
}


void free_step(double* d, double* kappa, double* mu, double* beta, double* gamma, double* t, double* n, double* norm, double* k1, double* k2, double* k3, double* k4, double* k5, double* k6)
{
  // Free interpolation
  free(d);
  free(kappa);
  free(mu);
  free(beta);
  free(gamma); 
  free(t);
  free(n);
  free(norm);
  
  // Free derivative
  free(k1);
  free(k2);
  free(k3);
  free(k4);
  free(k5);
  free(k6);
  
  return;
}


double distance(double a_x, double a_y, double b_x, double b_y)
{
  return sqrt((a_x - b_x)*(a_x - b_x) + (a_y - b_y)*(a_y - b_y));
}


double scalar_prod(double a_x, double a_y, double b_x, double b_y)
{
  return a_x*b_x + a_y*b_y;
}


double runge_kutta45(double* x, double* k1, double* k2, double* k3, double* k4, double* k5, double* k6, double tol_rk45_time, double dt, int M, int N, double* mu, double* beta, double* gamma, double* t, double* n, double alpha, double tol_rk45_space, double h, double theta, double* norm)
{
  double dt_new, R, R_new;
    
  // Allocate temporary postions
  double *x_temp, *x_RK4, *x_RK5;
  x_temp = (double*)malloc(2*N*sizeof(double));
  x_RK4 = (double*)malloc(2*N*sizeof(double));
  x_RK5 = (double*)malloc(2*N*sizeof(double));
  memset(x_temp, 0, 2*N*sizeof(double));
  memset(x_RK4, 0, 2*N*sizeof(double));
  memset(x_RK5, 0, 2*N*sizeof(double));
  
  // Runge-Kutta 45
  do
  {
    // Step 1 in RK
    //vfield_fft(k1, x, M, N, alpha, theta);
    vfield_orig(k1, x_temp, mu, beta, gamma, t, n, M, N, alpha, h, tol_rk45_space, theta, norm);
    for (int j = 0; j < 2*N; j++)
      x_temp[j] = x[j] + dt*0.25*k1[j];
    
    // Step 2 in RK
    //vfield_fft(k2, x, M, N, alpha, theta);
    vfield_orig(k2, x_temp, mu, beta, gamma, t, n, M, N, alpha, h, tol_rk45_space, theta, norm);
    for (int j = 0; j < 2*N; j++)
      x_temp[j] = x[j] + dt*(3.0*k1[j] + 9.0*k2[j])/32.0;

    // Step 3 in RK
    //vfield_fft(k3, x, M, N, alpha, theta);
    vfield_orig(k3, x_temp, mu, beta, gamma, t, n, M, N, alpha, h, tol_rk45_space, theta, norm);
    for (int j = 0; j < 2*N; j++)
      x_temp[j] = x[j] + dt*(1932.*k1[j] - 7200.*k2[j] + 7296.*k3[j])/2197.;
      
    // Step 4 in RK
    //vfield_fft(k4, x, M, N, alpha, theta);
    vfield_orig(k4, x_temp, mu, beta, gamma, t, n, M, N, alpha, h, tol_rk45_space, theta, norm);
    for (int j = 0; j < 2*N; j++)
      x_temp[j] = x[j] + dt*(439./216.*k1[j] - 8.*k2[j] + 3680./513.*k3[j] - 845./4104.*k4[j]);
    
    // Step 5 in RK
    //vfield_fft(k5, x, M, N, alpha, theta);
    vfield_orig(k5, x_temp, mu, beta, gamma, t, n, M, N, alpha, h, tol_rk45_space, theta, norm);
    for (int j = 0; j < 2*N; j++)
      x_temp[j] = x[j] + dt*(-8./27.*k1[j] + 2.*k2[j] - 3544./2565.*k3[j]\
                          + 1859./4104.*k4[j] - 11./40.*k5[j]);
      
    // Step 6 in RK
    //vfield_fft(k6, x, M, N, alpha, theta);
    vfield_orig(k6, x_temp, mu, beta, gamma, t, n, M, N, alpha, h, tol_rk45_space, theta, norm);
    
    // RK4 and RK5 approx
    for (int j = 0; j < 2*N; j++)
    {
    	x_RK4[j] = x[j] + dt*(25./216.*k1[j] + 1408./2565.*k3[j] + 2197./4104.*k4[j] - 0.2*k5[j]);
	  	x_RK5[j] = x[j] + dt*(16./135.*k1[j] + 6656./12825.*k3[j] + 28561./56430.*k4[j] - 9./50.*k5[j]+2./55.*k6[j]);
    }

    // Compute error
    R = 0.;
    for (int i = 0; i < 2*N; i++)
    {
      R_new = fabs(x_RK5[i] - x_RK4[i]);
      if (R_new > R)
        R = R_new;
    }
    
    // Calculate update factor
    if (R > tol_rk45_time)
    {
      dt = 0.5*dt;
      continue;
    }
    
    // Predict next time step
    dt_new = 0.9*dt*sqrt(sqrt(tol_rk45_time/R));
    //dt_new = sqrt(sqrt((1.e-12*dt/(2.*R))))*dt;

    if (dt_new > 2.5e-3)
      dt_new = 2.5e-3;
    printf("dt_new = %e\n", dt_new);    
    printf("dt = %e\n\n", dt);
  } while (R > tol_rk45_time);
  
  // Update position
  for (int j = 0; j < 2*N; j++)
    x[j] = x_RK5[j];

  // Free temporary positions
  free(x_temp);
  free(x_RK4);
  free(x_RK5);
  
  //printf("X dot dxdt = %e \n", x[0]*dxdt_RK5[0]+x[1]*dxdt_RK5[1]);
  
  return dt_new;
}
