#include "fft_functions.h"
#include "misc.h"
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>


void derivative_fft(double* x, double* dx, int start, int stop)
{ 
  int N = stop - start;
  int N_points = N/2;
  double k[N], filter[N];
  double a, m;
  a = 36;
  m = 36;
  
  // Setup FFT
  fftw_complex *in_x, *in_y, *out_x, *out_y;
  fftw_plan plan_for_x, plan_for_y, plan_back_x, plan_back_y;
  in_x = (fftw_complex*) fftw_malloc(N*sizeof(fftw_complex));
  in_y = (fftw_complex*) fftw_malloc(N*sizeof(fftw_complex));
  out_x = (fftw_complex*) fftw_malloc(N*sizeof(fftw_complex));
  out_y = (fftw_complex*) fftw_malloc(N*sizeof(fftw_complex));
  plan_for_x = fftw_plan_dft_1d(N, in_x, out_x, FFTW_FORWARD, FFTW_ESTIMATE);
  plan_for_y = fftw_plan_dft_1d(N, in_y, out_y, FFTW_FORWARD, FFTW_ESTIMATE);
  plan_back_x = fftw_plan_dft_1d(N, in_x, out_x, FFTW_BACKWARD, FFTW_ESTIMATE);
  plan_back_y = fftw_plan_dft_1d(N, in_y, out_y, FFTW_BACKWARD, FFTW_ESTIMATE);
  
  // Frequency
  for (int i = 0; i < N_points; i++)
  {
    k[i] = (double)i;
    k[N_points+i] = (double)(i-N_points); 
  }
  k[N_points] = 0;

  // Anti-aliasing filter
  for (int i = 0; i < N_points; i++)
    filter[i] = exp(-a*pow(k[i]/(double)N_points, m));
  for (int i = 1; i < N_points; i++)
    filter[i+N_points] = filter[N_points-i];  
  filter[N_points] = 0.;

  /*
  char str[80] = "../results/filter.txt";
  FILE* f = fopen(str, "wb");
  for (int i = 0; i < N; i++)
    fprintf(f, "%lf\n", filter[i]);
  fclose(f);
  */

  // Fourier transform x_i(p, t)
  for (int i = 0; i < N; i++)
  {
    in_x[i] = x[2*(start + i)];
    in_y[i] = x[2*(start + i)+1];
  }
  fftw_execute(plan_for_x); // Thread safe
  fftw_execute(plan_for_y);
  
  /*
  char str[80] = "../results/filter.txt";
  FILE* f = fopen(str, "wb");
  for (int i = 0; i < N; i++)
    fprintf(f, "%lf\n", creal(out_x[i]));
  fclose(f);
  */

  // Differentiate and transform back
  for (int i = 0; i < N; i++)
  {
    in_x[i] = filter[i]*I*k[i]*out_x[i];
    in_y[i] = filter[i]*I*k[i]*out_y[i];
  }
  fftw_execute(plan_back_x);
  fftw_execute(plan_back_y);
  
  for (int i = 0; i < N; i++)
  {
    dx[2*(start + i)] = creal(out_x[i])/((double)N);
    dx[2*(start + i)+1] = creal(out_y[i])/((double)N);
  }
  
  // Free
  fftw_destroy_plan(plan_for_x);
  fftw_destroy_plan(plan_for_y);
  fftw_destroy_plan(plan_back_x);
  fftw_destroy_plan(plan_back_y);
  fftw_free(in_x);
  fftw_free(in_y);
  fftw_free(out_x);
  fftw_free(out_y);
  
  return;
}


void vfield_fft(double* dxdt, double* x, int M, int N, double alpha, double theta)
{

  double alpha_d_x, aux_0, aux_1;
  
  double* dx = (double*)malloc(2*N*sizeof(double));
  
  // Compute derivative in nominator
  derivative_fft(x, dx, 0, M);
  derivative_fft(x, dx, M, N);

  // Estimate integral using Riemann sum 
  for (int j = 0; j < M; j++)
  {
  aux_0 = 0.;
  aux_1 = 0.;
  dxdt[2*j] = 0.;
  dxdt[2*j+1] = 0.;
    for (int i = 0; i < M; i++)
    {
      if (i == j)
        continue;

      // Denominator
      alpha_d_x = pow(distance(x[2*i], x[2*i+1], x[2*j], x[2*j+1]), alpha);
      
      // Minus tangent
      aux_0 += (dx[2*i] - dx[2*j])/alpha_d_x;
      aux_1 += (dx[2*i+1] - dx[2*j+1])/alpha_d_x;
      //aux_0 += dx[2*i]/alpha_d_x;
      //aux_1 += dx[2*i+1]/alpha_d_x;
    }
    
    // Normalize FFT and Riemann sum
    dxdt[2*j] = aux_0*theta/((double)M);
    dxdt[2*j+1] = aux_1*theta/((double)M);
    aux_0 = 0.;
    aux_1 = 0.;
    for (int i = M; i < N; i++)
    {
      // Denominator
      alpha_d_x = pow(distance(x[2*i], x[2*i+1], x[2*j], x[2*j+1]), alpha);
      
      // Normally 
      aux_0 += dx[2*i]/alpha_d_x;
      aux_1 += dx[2*i+1]/alpha_d_x;
    }
    dxdt[2*j] += aux_0*theta/((double)(N-M));
    dxdt[2*j+1] += aux_1*theta/((double)(N-M));
  }  
  
  // Estimate integral using Riemann sum 
  for (int j = M; j < N; j++)
  {
  aux_0 = 0.;
  aux_1 = 0.;
    for (int i = 0; i < M; i++)
    {
      // Denominator
      alpha_d_x = pow(distance(x[2*i], x[2*i+1], x[2*j], x[2*j+1]), alpha);
      
      // Normally 
      aux_0 += dx[2*i]/alpha_d_x;
      aux_1 += dx[2*i+1]/alpha_d_x;
    }
    
    // Normalize FFT and Riemann sum
    dxdt[2*j] += aux_0*theta/((double)M);
    dxdt[2*j+1] += aux_1*theta/((double)M);
    aux_0 = 0.;
    aux_1 = 0.;
    for (int i = M; i < N; i++)
    {
      if (i == j)
        continue;

      // Denominator
      alpha_d_x = pow(distance(x[2*i], x[2*i+1], x[2*j], x[2*j+1]), alpha);
      
      // Minus tangent
      aux_0 += (dx[2*i] - dx[2*j])/alpha_d_x;
      aux_1 += (dx[2*i+1] - dx[2*j+1])/alpha_d_x;
      //aux_0 += dx[2*i]/alpha_d_x;
      //aux_1 += dx[2*i+1]/alpha_d_x;
    }
    dxdt[2*j] += aux_0*theta/((double)(N-M));
    dxdt[2*j+1] += aux_1*theta/((double)(N-M));
  }  
  
  free(dx);
  
  return;
}


double area_fft(double* x, int start, int stop)
{
  int N, N_points;
  N = stop - start;
  double k[N], area;
  N_points = N/2;
  
  // Setup FFT
  fftw_complex *in_x, *out_x;
  fftw_plan plan_for_x, plan_back_x;
  in_x = (fftw_complex*) fftw_malloc(N*sizeof(fftw_complex));
  out_x = (fftw_complex*) fftw_malloc(N*sizeof(fftw_complex));
  plan_for_x = fftw_plan_dft_1d(N, in_x, out_x, FFTW_FORWARD, FFTW_ESTIMATE);
  plan_back_x = fftw_plan_dft_1d(N, in_x, out_x, FFTW_BACKWARD, FFTW_ESTIMATE);
  
  // Frequency
  for (int i = 0; i < N_points; i++)
  {
    k[i] = (double)i;
    k[N_points+i] = (double)(i-N_points); 
  }
  k[N_points] = 0;
  
  // Fourier transform x_i(p, t)
  for (int i = 0; i < N; i++)
    in_x[i] = x[2*(start + i)];
  fftw_execute(plan_for_x); // Thread safe
  
  // Differentiate and transform back
  for (int i = 0; i < N; i++)
    in_x[i] = I*k[i]*out_x[i];
  fftw_execute(plan_back_x);
  
  // Estimate integral using Riemann sum 
  area = 0;
  for (int i = 0; i < N; i++)
    area += creal(out_x[i])*(x[2*(start + i)+1]);
  area = -area/(2.*(double)(N*N));
  
  // Free
  fftw_destroy_plan(plan_for_x);
  fftw_destroy_plan(plan_back_x);
  fftw_free(in_x);
  fftw_free(out_x);
  fftw_cleanup();
  
  return area;
}

