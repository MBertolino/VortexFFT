#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "orig_functions.h"
#include "fft_functions.h"
#include "bh_functions.h"
#include <unistd.h>
#include <sys/time.h>

/* Here we define a routine called get_wall_seconds() the gets number
   of seconds and microseconds since the Epoch (1970-01-01 00:00:00
   +0000 (UTC)). The seconds and microseconds values are combined in a
   double number giving the number of seconds since the Epoch. */
double get_wall_seconds() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
  return seconds;
}


void print_to_file(double* x, int M, int N, int k)
{
  // Print to file
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
    fprintf(f, "%.12f %.12f\n", x[2*i], x[2*i+1]);
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
    vfield_fft(k1, x, M, N, alpha, theta);
    //vfield_orig(k1, x, mu, beta, gamma, t, n, M, N, alpha, h, tol_rk45_space, theta, norm);
    for (int j = 0; j < 2*N; j++)
      x_temp[j] = x[j] + dt*0.25*k1[j];
    
    // Step 2 in RK
    vfield_fft(k2, x_temp, M, N, alpha, theta);
    //vfield_orig(k2, x_temp, mu, beta, gamma, t, n, M, N, alpha, h, tol_rk45_space, theta, norm);
    for (int j = 0; j < 2*N; j++)
      x_temp[j] = x[j] + dt*(3.0*k1[j] + 9.0*k2[j])/32.0;

    // Step 3 in RK
    vfield_fft(k3, x_temp, M, N, alpha, theta);
    //vfield_orig(k3, x_temp, mu, beta, gamma, t, n, M, N, alpha, h, tol_rk45_space, theta, norm);
    for (int j = 0; j < 2*N; j++)
      x_temp[j] = x[j] + dt*(1932.*k1[j] - 7200.*k2[j] + 7296.*k3[j])/2197.;
      
    // Step 4 in RK
    vfield_fft(k4, x_temp, M, N, alpha, theta);
    //vfield_orig(k4, x_temp, mu, beta, gamma, t, n, M, N, alpha, h, tol_rk45_space, theta, norm);
    for (int j = 0; j < 2*N; j++)
      x_temp[j] = x[j] + dt*(439./216.*k1[j] - 8.*k2[j] + 3680./513.*k3[j] - 845./4104.*k4[j]);
    
    // Step 5 in RK
    vfield_fft(k5, x_temp, M, N, alpha, theta);
    //vfield_orig(k5, x_temp, mu, beta, gamma, t, n, M, N, alpha, h, tol_rk45_space, theta, norm);
    for (int j = 0; j < 2*N; j++)
      x_temp[j] = x[j] + dt*(-8./27.*k1[j] + 2.*k2[j] - 3544./2565.*k3[j]\
                          + 1859./4104.*k4[j] - 11./40.*k5[j]);
      
    // Step 6 in RK
    vfield_fft(k6, x_temp, M, N, alpha, theta);
    //vfield_orig(k6, x_temp, mu, beta, gamma, t, n, M, N, alpha, h, tol_rk45_space, theta, norm);
    
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


void compare_algo(double* x, double* mu, double* beta, double* gamma, double* t, double* n, int M, int N, double h, double tol_rk45_space, double alpha, double theta)
{
    // Allocate
    double* dxdt_ama = (double*)malloc(2*N*sizeof(double));
    double* dxdt_fft = (double*)malloc(2*N*sizeof(double));
    /*double* dxdt_bh = (double*)malloc(2*N*sizeof(double));
    node_t *tree = NULL;
    double* dx = (double*)malloc(2*N*sizeof(double));
    derivative_fft(x, dx, 0, M);
    */double t_drit, t_fft, t_bh;

    // Compare vector field calculations
    t_drit = get_wall_seconds();
    vfield_orig(dxdt_ama, x, mu, beta, gamma, t, n, M, N, alpha, h, tol_rk45_space, theta);
    t_drit = get_wall_seconds() - t_drit;

    sleep(0.1);

    t_fft = get_wall_seconds();
    vfield_fft(dxdt_fft, x, M, N, alpha, theta);
    t_fft = get_wall_seconds() - t_fft;

    sleep(0.1);

    /*t_bh = get_wall_seconds();
    for (int i = 0; i < N; i++)
      insert(&tree, 0.5, 0.5, 100, dxdt_bh[2*i], dxdt_bh[2*i+1], x[2*i], x[2*i+1], dx[2*i], dx[2*i+1], i);
    vfield_bh(tree, tree, theta, 0.001, N, alpha);
    t_bh = get_wall_seconds() - t_bh;
    write_tree(tree, &dxdt_bh);
   */ 
    // Print accuracy to file
    char strdx_ama[80] = "../results/dx_circle_alpha07_ama.txt";
    char strdx_fft[80] = "../results/dx_circle_alpha07_fft.txt";
    //char strdx_bh[80] = "../results/dx_other_alpha07_normal_bh.txt";
    FILE* fa_ama = fopen(strdx_ama, "a");
    FILE* fa_fft = fopen(strdx_fft, "a");
    //FILE* fa_bh = fopen(strdx_bh, "a");
    fprintf(fa_ama, "%.12f %.12f %d %lf\n", dxdt_ama[0], dxdt_ama[1], M, alpha);
    fprintf(fa_fft, "%.12f %.12f %d %lf\n", dxdt_fft[0], dxdt_fft[1], M, alpha);
    //fprintf(fa_bh, "%e %e %d %lf\n", dxdt_bh[0], dxdt_bh[1], M, alpha);
    fclose(fa_ama);
    fclose(fa_fft);
    //fclose(fa_bh);
    
    // Print time to file
    char time_ama[80] = "../results/time_drit.txt";
    char time_fft[80] = "../results/time_fft.txt";
    //char time_bh[80] = "../results/time_bh.txt";
    FILE* ft_ama = fopen(time_ama, "a");
    FILE* ft_fft = fopen(time_fft, "a");
    //FILE* ft_bh = fopen(time_bh, "a");
    fprintf(ft_ama, "%e %d %lf\n", t_drit, M, alpha);
    fprintf(ft_fft, "%e %d %lf\n", t_fft, M, alpha);
    //fprintf(ft_bh, "%e %d %lf\n", t_bh, M, alpha);
    fclose(ft_ama);
    fclose(ft_fft);
    //fclose(ft_bh);

    printf("time_ama = %e \n", t_drit);
    printf("time_fft = %e \n", t_fft);
    //printf("time__bh = %e \n\n", t_bh);

    printf("ama_x[%d, %d] = %e, \t %e\n", 0, 1, dxdt_ama[0], dxdt_ama[1]);
    printf("fft_x[%d, %d] = %e, \t %e\n", 0, 1, dxdt_fft[0], dxdt_fft[1]);
    //printf("bh__x[%d, %d] = %e, \t %e\n\n", 0, 1, dxdt_bh[0], dxdt_bh[1]);

    free(dxdt_ama);
    free(dxdt_fft);
    //free(dxdt_bh);

    return;

}
