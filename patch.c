#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "functions.h"
#include <unistd.h>


#define TWOPI 6.2831853071795864769



int main(int argc, char **argv) {
  
  // Number of points
  int M = atoi(argv[1]); // Number of points in each circle
  int M2 = 128;
  int N = M + M2;
  int n_dim = 2;
  int size = N*n_dim;
  int T = 5000;
  double tol_rk45_time = 1.e-10;
  long double tol_rk45_space = 1.e-10;
  long double h = 1.e-4;
  double alpha = 0.7; // Interpolation between 2D Euler and Quasi-geostrophic
  double theta = -1.0;
  double dt = 1.e-3;//1.*h;
  double F, tpi, time;
  int *pN, *pM1, *pM2;
  double** px;
  int zeros;
  int N_old;
  time = 0;
  pM1 = &M;
  pM2 = &M2;

  // Allocate points
  double* x = (double*)malloc(size*sizeof(double));
  double kappa_den[2];  

  px = &x;
  pN = &N;
  
  // Generate circle
  for (int j = 0; j < M; j++) {
    x[2*j] =cos(TWOPI*j/(double)M) - 1.1;//cos(TWOPI*j/(double)M) + 0.45*sin(TWOPI*5*j/(double)M); // 
    x[2*j+1] = sin(TWOPI*j/(double)M);//sin(TWOPI*j/(double)M) + 0.3*cos(TWOPI*3*j/(double)M);
    
    x[2*j+2*M] =  cos(TWOPI*j/(double)M) + 1.1;
    x[2*j+1+2*M] = sin(TWOPI*j/(double)M);
  }
  double area1, area2;
 
  // Print to file  
  char str_start[80] = "../results/circle_start.txt";
  FILE* f = fopen(str_start, "wb");
  for (int i = 0; i < N; i++) {
    fprintf(f, "%lf,%lf\n", x[2*i], x[2*i+1]);
  }
  fclose(f);
  
  // Evolve
  for (int k = 0; k <= T; k++) {
  
  printf(" \n \n--------------------------\n \n");
    printf("k = %d\n", k);
    
    double* d = (double*)malloc(N*sizeof(double));
    double* kappa = (double*)malloc(N*sizeof(double));
    double* mu = (double*)malloc(N*sizeof(double));
    double* beta = (double*)malloc(N*sizeof(double));
    double* gamma = (double*)malloc(N*sizeof(double));
    double* t = (double*)malloc(size*sizeof(double)); 
    double* n = (double*)malloc(size*sizeof(double));
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
    
    zeros = size*sizeof(double);
    memset(dxdt, 0, zeros);
    memset(dxdt_fft, 0, zeros);
    memset(dxdt_k1, 0, zeros);
    memset(dxdt_k2, 0, zeros);
    memset(dxdt_k3, 0, zeros);
    memset(dxdt_k4, 0, zeros);
    memset(dxdt_k5, 0, zeros);
    memset(dxdt_k6, 0, zeros);
    memset(dxdt_RK4, 0, zeros);
    memset(dxdt_RK5, 0, zeros);
    
    // Interpolate
    interpolate(x, 0, M, n_dim, t, n, d, kappa, kappa_den, mu, beta, gamma);
    interpolate(x, M, N, n_dim, t, n, d, kappa, kappa_den, mu, beta, gamma);
    
  //  compute_fft(dxdt_fft, x, N, alpha);


    // Compare FFT and Mancho
   /* for (int j = 0; j < N; j++)
    {
      
      compute_derivative(dxdt, x, mu, beta, gamma, t, n, M, N, alpha, h, tol_rk45_space, j);
    }*/
    
    // Print to file
  /*  char strdx_fft[80] = "../results/dx_fft.txt";
    char strdx_ama[80] = "../results/dx_ama.txt";
    FILE* f_fft = fopen(strdx_fft, "wb");
    FILE* f_ama = fopen(strdx_ama, "wb");
    for (int j = 0; j < 8; j++)
    {
      printf("fft_x[%d, %d] = %e, \t %e\n", 2*j, 2*j+1, -dxdt_fft[2*j]/TWOPI, -dxdt_fft[2*j+1]/TWOPI);
      printf("ama_x[%d, %d] = %e, \t %e\n\n", 2*j, 2*j+1, -dxdt[2*j]/TWOPI, -dxdt[2*j+1]/TWOPI);
    }
    for (int j = 0; j < N; j++)
    {
      fprintf(f_fft, "%lf %lf\n", -dxdt_fft[2*j]/TWOPI, -dxdt_fft[2*j+1]/TWOPI);
      fprintf(f_ama, "%lf %lf\n", -dxdt[2*j]/TWOPI, -dxdt[2*j+1]/TWOPI);
    } 
    printf("Done\n");
    fclose(f_fft);
    fclose(f_ama);
    sleep(5);
    
    */
    // Evolve patches
    dt = runge_kutta45(x, dxdt_k1, dxdt_k2, dxdt_k3, dxdt_k4, dxdt_k5,\
                  dxdt_k6, dxdt_RK4, dxdt_RK5, tol_rk45_time, dt, M, N,\
                  mu, beta, gamma, t, n, alpha, tol_rk45_space, h, &time);
    printf("time = %1.15lf\n", time);
    area1 = compute_area(x, 0, M, t, n, mu, beta, gamma);
    area2 = compute_area(x, M, N, t, n, mu, beta, gamma);
    printf("area1 = %e, area2 = %e \n", area1, area2); 
    
    N_old  = N;
    interpolate(x, 0, M, n_dim, t, n, d, kappa, kappa_den, mu, beta, gamma);
    interpolate(x, M, N, n_dim, t, n, d, kappa, kappa_den, mu, beta, gamma);
    
    points_reloc(&x, t, n, pN, kappa, mu, beta, gamma, pM1, pM2, 2);
    
    //Print to file
    char str[80] = "../results/circle_";
    char str2[80] = "";
    sprintf(str2, "%d", k);
    strcat(str, str2);
    strcat(str, ".txt");
    FILE* f = fopen(str, "wb");
    for (int i = 0; i < N; i++) {
      fprintf(f, "%lf,%lf\n", x[2*i], x[2*i+1]);
    }
    fclose(f);

    // Redistribute the nodes
    printf("N = %d, N_old = %d \n", N, N_old);
    size = N*n_dim;
    
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
