#include "functions.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <complex.h>
#include <fftw3.h>

#define TWOPI 6.2831853071795864769
#define SQRTTWO 1.4142135623730950588
#define PRINT 0

void interpolate(double** x, int start, int N, int n_dim, double** t, double** n,\
                 double* d, double* kappa, double* kappa_den, double* mu,\
                 double* beta, double* gamma)
{
    //printf("N = %d \n", N);
  // Calculate t an n
  for (int j = start; j < N-1; j++) {
  // printf("j = %d \n", j);
    t[j][0] = x[j+1][0] - x[j][0];
    t[j][1] = x[j+1][1] - x[j][1];
    n[j][0] = -t[j][1];
    n[j][1] = t[j][0];
    d[j] = sqrt(t[j][0]*t[j][0] + t[j][1]*t[j][1]);
  }

  // Special case j = N-1
  t[N-1][0] = x[start][0] - x[N-1][0];
  t[N-1][1] = x[start][1] - x[N-1][1];
  n[N-1][0] = -t[N-1][1];
  n[N-1][1] = t[N-1][0];
  d[N-1] = sqrt(t[N-1][0]*t[N-1][0] + t[N-1][1]*t[N-1][1]);

  // kappa local curvature
  kappa_den[0] = d[N-1]*d[N-1]*t[start][0] + d[start]*d[start]*t[N-1][0];
  kappa_den[1] = d[N-1]*d[N-1]*t[start][1] + d[start]*d[start]*t[N-1][1];
  kappa[start] = 2.*(t[N-1][0]*t[start][1] - t[N-1][1]*t[start][0])\
    /sqrt(kappa_den[0]*kappa_den[0] + kappa_den[1]*kappa_den[1]);
    
  for (int j = start+1; j < N; j++) {
    // kappa local curvature
    kappa_den[0] = (d[j-1]*d[j-1]*t[j][0] + d[j]*d[j]*t[j-1][0]);
    kappa_den[1] = (d[j-1]*d[j-1]*t[j][1] + d[j]*d[j]*t[j-1][1]);
    kappa[j] = 2.*(t[j-1][0]*t[j][1] - t[j-1][1]*t[j][0])\
      /sqrt(kappa_den[0]*kappa_den[0] + kappa_den[1]*kappa_den[1]);
  }
  
  // Construct the cubic interpolation
  for (int j = start; j < N-1; j++) {
    // Cubic interpolation coefficients
    mu[j] = -1./3.*d[j]*kappa[j] - 1./6.*d[j]*kappa[j+1];
    beta[j] = 0.5*d[j]*kappa[j];
    gamma[j] = 1./6.*d[j]*(kappa[j+1] - kappa[j]);
  }
  mu[N-1] = -1./3.*d[N-1]*kappa[N-1] - 1./6.*d[N-1]*kappa[start];
  beta[N-1] = 0.5*d[N-1]*kappa[N-1];
  gamma[N-1] = 1./6.*d[N-1]*(kappa[start] - kappa[N-1]);

  return;
}

void autder(double* f, double* c_coeff, double alpha, int order)
{
  // Allocate memory for temporary coefficients
  double* a_ = (double*)malloc(order*sizeof(double));
  for (int n = 1; n < order; n++)
  {
    a_[n] = 0.;
    f[n] = 0.;
  }
    a_[0] = 1.;

  // Calculate temporary coefficients
  for (int n = 1; n < order; n++)
  {
    for (int j = 1; j <= n; j++)
    {
      a_[n] -= c_coeff[j]*a_[n-j];
    }
  }
  
  f[0] = 1.;
  
  // Calculate Taylor coefficients of order, "order" :D
  for (int n = 1; n < order; n++)
  {
    for (int j = 0; j < n; j++)
    {
      f[n] += (double)(n*alpha - j*(alpha + 1))*a_[n-j]*f[j];
    }
    f[n] /= (double)(n); 
  }
  free(a_);

  return;
}

void compute_fft(double* dxdt, double** x, int N, double alpha, int j)
{
  int N_points = 0.5*N;
  double k[N];
  double alpha_d_x;
  
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

  // Fourier transform x_i(p, t)
  for (int i = 0; i < N; i++)
  {
    in_x[i] = x[i][0];
    in_y[i] = x[i][1];
  }
  fftw_execute(plan_for_x); // Thread safe
  fftw_execute(plan_for_y);
  
  // Differentiate in time and transform back
  for (int i = 0; i < N; i++)
  {
    in_x[i] = I*k[i]*out_x[i];
    in_y[i] = I*k[i]*out_y[i];
  }
  fftw_execute(plan_back_x);
  fftw_execute(plan_back_y);
  
  // Estimate integral using Riemann sum 
  for (int i = 0; i < N; i++)
  {
    if (i == j)
      continue;
	  alpha_d_x = pow(sqrt((x[i][0] - x[j][0])*(x[i][0] - x[j][0])\
  	  + (x[i][1] - x[j][1])*(x[i][1] - x[j][1])), alpha);
    dxdt[0] += out_x[i]/alpha_d_x;
    dxdt[1] += out_y[i]/alpha_d_x;
  }
  
  // Normalize
  dxdt[0] = dxdt[0]*TWOPI/((double)(N*N));
  dxdt[1] = dxdt[1]*TWOPI/((double)(N*N));
  
  // Free
  fftw_destroy_plan(plan_for_x);
  fftw_destroy_plan(plan_for_y);
  fftw_destroy_plan(plan_back_x);
  fftw_destroy_plan(plan_back_y);
  fftw_free(in_x);
  fftw_free(in_y);
  fftw_free(out_x);
  fftw_free(out_y);
  fftw_cleanup();
  
  return;
}

void compute_derivative(double* dxdt, double** x, double* mu, double* beta, double* gamma, double** t, double** n, int M, int N, double alpha, double h, double tol_rk45_space, int j)
{
  //printf("Entering compute derivative\n");
  double d_x, d_ni, d_ti, d_xi;
  int order = 11;
  double Q = 0.01;
  double f = 1./sqrt(Q);
  
  // Generate coefficients
  double c[order];
  double g[order];
  double poly_coeff_c[order];
  double poly_coeff_g[order];
  double mu_2, beta_2;
  
  // Evolve the contour integrals
  for (int i = 0; i < N; i++) { // i < ???
    #if PRINT    
      printf("i = %d,  ", i);
      printf("j = %d,  ", j);
      //printf("i = %d, j = %d, \n", i, j);
    #endif
      
	  d_x = sqrt((x[i][0] - x[j][0])*(x[i][0] - x[j][0])\
				  + (x[i][1] - x[j][1])*(x[i][1] - x[j][1]));
    d_ti = -((x[j][0] - x[i][0])*t[j][0] + (x[j][1] - x[i][1])*t[j][1]);
    d_ni = -((x[j][0] - x[i][0])*n[j][0] + (x[j][1] - x[i][1])*n[j][1]);
   
    // Distance between x_i and x_{i+1}
    if (i+1 == N) {
      d_xi = sqrt((x[0][0] - x[i][0])*(x[0][0] - x[i][0])\
            + (x[0][1] - x[i][1])*(x[0][1] - x[i][1]));
    } else {
        d_xi = sqrt((x[i+1][0] - x[i][0])*(x[i+1][0] - x[i][0])\
              + (x[i+1][1] - x[i][1])*(x[i+1][1] - x[i][1])); 
    }
    
    // Initialize Taylor coefficients
    for (int n = 0; n < order; n++) {
      c[n] = 0.;
      g[n] = 0.;
      poly_coeff_c[n] = 0.;
      poly_coeff_g[n] = 0.;
    }
    poly_coeff_c[0] = 1.;
    poly_coeff_g[0] = 1.;
    
    // Evaluate integrals
    if ((i+1 == M && j == 0) || (i+1 == N && j == M)) {
      // Edge case
      // Case 2: Use formula (29) with shifted params
      #if PRINT
        printf("Edge Case 2\n"); 
      #endif

      // Update parameters
      mu_2 = mu[i] + 2.*beta[i] + 3.*gamma[i];
      beta_2 = -beta[i] - 3.*gamma[i];
      
      // Generate Taylor coefficients
      poly_coeff_c[1] = 2.*mu_2*beta_2/(1. + mu_2*mu_2);
      poly_coeff_c[2] = (beta_2*beta_2 + 2.*mu_2*gamma[i])/(1. + mu_2*mu_2);
      poly_coeff_c[3] = 2.*beta_2*gamma[i]/(1. + mu_2*mu_2);
      poly_coeff_c[4] = gamma[i]*gamma[i]/(1. + mu_2*mu_2);
      
      autder(c, poly_coeff_c, alpha/2., order);
      
      // Look at inputs in these functions
      evaluate_integral(dxdt, mu_2, beta_2, gamma[i], t[i], n[i], c, alpha); 
    
      } else {
      
      if ((x[i][0] == x[j][0]) && (x[i][0] == x[j][0])) {
        // Case 1: Use formula (29)
        #if PRINT
          printf("Case 1\n");
        #endif
        
        // Generate Taylor coefficients
        poly_coeff_c[1] = 2.*mu[i]*beta[i]/(1. + mu[i]*mu[i]);
        poly_coeff_c[2] = (beta[i]*beta[i] + 2.*mu[i]*gamma[i])/(1. + mu[i]*mu[i]);
        poly_coeff_c[3] = 2.*beta[i]*gamma[i]/(1. + mu[i]*mu[i]);
        poly_coeff_c[4] = gamma[i]*gamma[i]/(1. + mu[i]*mu[i]);

        autder(c, poly_coeff_c, alpha/2., order);
        
        evaluate_integral(dxdt, mu[i], beta[i], gamma[i], t[i], n[i], c, alpha);
      
      }
      else if ((i == j-1) && (i != M-1))
      {
        // Case 2: Use formula (29) with shifted params
        #if PRINT
          printf("Case 2\n");      
        #endif
        
		    // Update parameters
        mu_2 = mu[i] + 2.*beta[i] + 3.*gamma[i];
        beta_2 = -beta[i] - 3.*gamma[i];
      
        // Generate Taylor coefficients
        poly_coeff_c[1] = 2.*mu_2*beta_2/(1. + mu_2*mu_2);
        poly_coeff_c[2] = (beta_2*beta_2 + 2.*mu_2*gamma[i])/(1. + mu_2*mu_2);
        poly_coeff_c[3] = 2.*beta_2*gamma[i]/(1. + mu_2*mu_2);
        poly_coeff_c[4] = gamma[i]*gamma[i]/(1. + mu_2*mu_2);
      
        autder(c, poly_coeff_c, alpha/2., order);
        
        evaluate_integral(dxdt, mu_2, beta_2, gamma[i], t[i], n[i], c, alpha); // Look at inputs in these functions

      } else if (d_x > f*d_xi) {
        #if PRINT   
          printf("Case 3\n");
        #endif
        // Case 3: Use formula (31)
        
        // Generate Taylor coefficients
        poly_coeff_g[1] = (d_ti + d_ni*mu[i])/(d_x*d_x);
        
        poly_coeff_g[2] = ((t[i][0]*t[i][0] + t[i][1]*t[i][1])\
                        + mu[i]*mu[i]*(n[i][0]*n[i][0] + n[i][1]*n[i][1])\
                        + d_ni*beta[i])/(d_x*d_x);
        
        poly_coeff_g[3] = (2.*mu[i]*beta[i]*(n[i][0]*n[i][0]\
                        + n[i][1]*n[i][1]) + d_ni*gamma[i])/(d_x*d_x);
        
        poly_coeff_g[4] = ((beta[i]*beta[i] + 2.*mu[i]*gamma[i])\
                        *(n[i][0]*n[i][0] + n[i][1]*n[i][1]))/(d_x*d_x);
        
        poly_coeff_g[5] = 2.*beta[i]*gamma[i]*(n[i][0]*n[i][0]\
                        + n[i][1]*n[i][1])/(d_x*d_x);
        
        poly_coeff_g[6] = (gamma[i]*gamma[i]*(n[i][0]*n[i][0]\
                        + n[i][1]*n[i][1]))/(d_x*d_x);
        
        autder(g, poly_coeff_g, alpha/2., order);

        evaluate_integral_g(dxdt, mu[i], beta[i], gamma[i], d_x, d_ni, d_ti, t[i], n[i], g, alpha);

      } else {
        #if PRINT
          printf("Case 4\n");
        #endif
        // Case 4: Use Runge-Kutta 4-5
        evaluate_integral_RK(dxdt, x[i], x[j], mu[i], beta[i], gamma[i], tol_rk45_space, h, t[i], n[i], alpha);
      }
    }
    //if (i == 3)
    //  printf("\n");
    //sleep(1);
  }
  //printf("----------------------------\n\n");
  
  return;
}

void evaluate_integral(double* dxdt, double mu_i, double beta_i, double gamma_i, double* t_i, double* n_i, double* c, double alpha) {
  
  // Compute the integrals
  double first = 0.;
  double second = 0.;
  double p_coef, psq_coef;
  double t_abs = pow((t_i[0]*t_i[0] + t_i[1]*t_i[1]), alpha);
  
  double alpha_mu = pow((1 + mu_i*mu_i), 0.5*alpha);
  
  for (int n = 0; n < 11; n++) {
    first += c[n]/(double)(n - alpha + 1.);
    
    p_coef = 2.*beta_i/(double)(n - alpha + 2.);
    psq_coef = 3.*gamma_i/(double)(n - alpha + 3.);
    second += c[n]*(p_coef + psq_coef);
  }
  
  first = first/(t_abs*alpha_mu);
  second = second/(t_abs*alpha_mu);
  
  // Sum together
  dxdt[0] += first*(t_i[0] + mu_i*n_i[0]) + second*n_i[0];
  dxdt[1] += first*(t_i[1] + mu_i*n_i[1]) + second*n_i[1];

  return;
}

void evaluate_integral_g(double* dxdt, double mu_i, double beta_i, double gamma_i, double d_x, double d_ni, double d_ti, double* t_i, double* n_i, double* g, double alpha)
{
  // Compute the integrals
  double first = 0.;
  double second = 0.;
  double p_coef, psq_coef;

  for (int n = 0; n < 11; n++) {
    first += g[n];
    
    p_coef = 2.*beta_i/(double)(n - alpha + 2.);
    psq_coef = 3.*gamma_i/(double)(n - alpha + 3.);
    second += g[n]*(p_coef + psq_coef);
  }
  
  first = first/pow(d_x, alpha);
  second = second/pow(d_x, alpha);
  
  // Sum together
  dxdt[0] += first*(t_i[0] + mu_i*n_i[0]) + second*n_i[0];
  dxdt[1] += first*(t_i[1] + mu_i*n_i[1]) + second*n_i[1];
  
  return;
}

void evaluate_integral_RK(double* dxdt, double* x_i, double* x_j, double mu_i, double beta_i, double gamma_i,\
                          double tol_rk45_space, double h, double* t_i, double* n_i, double alpha)
{
  double first = evaluate_integral1_RK(x_i, x_j, tol_rk45_space, h, t_i, n_i, mu_i, beta_i, gamma_i, alpha);
  double second = evaluate_integral2_RK(x_i, x_j, tol_rk45_space, h, t_i, n_i, mu_i, beta_i, gamma_i, alpha);

  dxdt[0] += first*(t_i[0] + mu_i*n_i[0]) + second*n_i[0];
  dxdt[1] += first*(t_i[1] + mu_i*n_i[1]) + second*n_i[1];
  
  return;
}

double evaluate_integral1_RK(double* x_i, double* x_j, double tol_rk45_space, double h, double* t_i,\
									double* n_i, double mu_i, double beta_i, double gamma_i, double alpha)
{
  
	double p = 0.;
	double p_end = 1.; 
  double k1, k2, k3, k4, k5, k6;
	double Y = 0., Y1, Y2, delta, p_temp;
  long double R;
    while (p < p_end)
  {
    do
    {
        
		k1 = h*integrand1(x_i, x_j, p, t_i, n_i, mu_i, beta_i, gamma_i, alpha); 
		p_temp = p + 0.25*h;
		
		k2 = h*integrand1(x_i, x_j, p_temp, t_i, n_i, mu_i, beta_i, gamma_i, alpha);
		p_temp = p + 3.0*h/8.0;
		
		k3 = h*integrand1(x_i, x_j, p_temp, t_i, n_i, mu_i, beta_i, gamma_i, alpha);
		p_temp = p + 12.0*h/13.0;
		
		k4 = h*integrand1(x_i, x_j, p_temp, t_i, n_i, mu_i, beta_i, gamma_i, alpha);
		p_temp = p + h;
		
		k5 = h*integrand1(x_i, x_j, p_temp, t_i, n_i, mu_i, beta_i, gamma_i, alpha);
		p_temp = p + 0.5*h;
		
		k6 = h*integrand1(x_i, x_j, p_temp, t_i, n_i, mu_i, beta_i, gamma_i, alpha);
    
		//RK4 approx
		Y1 = Y + 25.0*k1/216.0 + 1408.0*k3/2565.0 + 2197.0*k4/4104.0 - 0.2*k5;

    //RK5 approx
		Y2 = Y + 16.0*k1/135.0 + 6656.0*k3/12825.0 + 28561.0*k4/56430.0\
				-9.0*k5/50.0 + 2.0*k6/55.0;
    
    // Compute error
    R = fabs(Y2 - Y1);
    if (R > tol_rk45_space)
      h = 0.5*h;
    } while (R > tol_rk45_space);
    
    p = p + h;
    Y = Y2;
    //printf("R = %e \n", R);
    if (R < 1.e-10*tol_rk45_space)
    {
      h = 1.5*h;
    }
    else 
    {
      h = 0.9*h*sqrt(sqrt((h*tol_rk45_space)/R));
    }
  }  
  return Y;
}

double integrand1(double* x_i, double* x_j, double p, double* t_i, double* n_i, double mu_i,\
                  double beta_i, double gamma_i, double alpha)
{
  double eta_i = (mu_i + (beta_i + gamma_i*p)*p)*p;
  
  // Integrand
  double x_part = (x_j[0] - x_i[0] - t_i[0]*p - eta_i*n_i[0]);
  double y_part = (x_j[1] - x_i[1] - t_i[1]*p - eta_i*n_i[1]);
  double w = 1.0/pow(x_part*x_part + y_part*y_part, 0.5*alpha);
  
  return w;
}

double evaluate_integral2_RK(double* x_i, double* x_j, double tol_rk45_space, double h,\
            double* t_i, double* n_i, double mu_i, double beta_i, double gamma_i, double alpha) 
{
  
	double p = 0.;
	double p_end = 1.; 
  double k1, k2, k3, k4, k5, k6;
	double Y = 0., Y1, Y2, R, delta, p_temp, h_new;
  h_new = 0;
  
  while (p < p_end)
  {
    do
    {
		k1 = h*integrand2(x_i, x_j, p, t_i, n_i, mu_i, beta_i, gamma_i, alpha); 
		p_temp = p + 0.25*h;
		
		k2 = h*integrand2(x_i, x_j, p_temp, t_i, n_i, mu_i, beta_i, gamma_i, alpha);
		p_temp = p + 3.0*h/8.0;
		
		k3 = h*integrand2(x_i, x_j, p_temp, t_i, n_i, mu_i, beta_i, gamma_i, alpha);
		p_temp = p + 12.0*h/13.0;
		
		k4 = h*integrand2(x_i, x_j, p_temp, t_i, n_i, mu_i, beta_i, gamma_i, alpha);
		p_temp = p + h;
		
		k5 = h*integrand2(x_i, x_j, p_temp, t_i, n_i, mu_i, beta_i, gamma_i, alpha);
		p_temp = p + 0.5*h;
		
		k6 = h*integrand2(x_i, x_j, p_temp, t_i, n_i, mu_i, beta_i, gamma_i, alpha);
    
		//RK4 approx
		Y1 = Y + 25.0*k1/216.0 + 1408.0*k3/2565.0 + 2197.0*k4/4104.0 - 0.2*k5;

    //RK5 approx
		Y2 = Y + 16.0*k1/135.0 + 6656.0*k3/12825.0 + 28561.0*k4/56430.0\
				-9.0*k5/50.0 + 2.0*k6/55.0;
    
    // Compute error
    R = fabs(Y2 - Y1);
    if (R > tol_rk45_space)
      h = 0.5*h;
    } while (R > tol_rk45_space);
    
    p = p + h;
    Y = Y2;
    //printf("R = %e \n", R);
    if (R < 1.e-10*tol_rk45_space)
    {
      h = 2*h;
    }
    else 
    {
      h = 0.9*h*sqrt(sqrt((h*tol_rk45_space)/R));
    }
  }

  return Y;
}

double integrand2(double* x_i, double* x_j, double p, double* t_i, double* n_i,\
                  double mu_i, double beta_i, double gamma_i, double alpha)
{
  double eta_i = (mu_i + (beta_i + gamma_i*p)*p)*p;
  
  // Integrand
  double x_part = (x_i[0] - x_j[0] - t_i[0]*p + eta_i*n_i[0]);
  double y_part = (x_i[1] - x_j[1] - t_i[1]*p + eta_i*n_i[1]);
  
  double w = p*(2.0*beta_i + 3.0*gamma_i*p)/pow(x_part*x_part + y_part*y_part, 0.5*alpha);
  
	return w;
}

// Relocating nodes
void points_reloc(double*** px, double** t, double** n, int* pN, double* kappa,\
									double* mu, double* beta, double* gamma) {
  
  double **x;
  x = *px;
  int N;
  N = *pN;
  double epsilon = 1.e-6;
  double L = 3.0;
  double a = 2.0/3.0;
  double v = 0.05; // 0.01, 0.03, 0.05 (Smaller value => Larger densities of points on the curve)
  
  // Either calculate all d[j]'s here or use input
  // Same goes with kappa and h
  double *d, *kappa_bar, *kappai_breve, *kappai_tilde, *sigmai_prim;
  double  *kappai_hat, *rho, *sigmai, *sigmai_tilde, **h, **x_temp;
  d = (double*)malloc(N*sizeof(double));
  kappa_bar = (double*)malloc(N*sizeof(double));
  kappai_breve = (double*)malloc(N*sizeof(double));
  kappai_tilde = (double*)malloc(N*sizeof(double));  
  kappai_hat = (double*)malloc(N*sizeof(double));
  rho = (double*)malloc(N*sizeof(double));
  sigmai = (double*)malloc(N*sizeof(double));
  sigmai_tilde = (double*)malloc(N*sizeof(double));
  sigmai_prim = (double*)malloc(N*sizeof(double));
  h = (double**)malloc(N*sizeof(double*));

  for (int k = 0; k < N; k++)
    h[k] = (double*)malloc(N*sizeof(double));
 
  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < N; j++)
    {
      if (j == N-1)
	    {
      h[i][j] = sqrt((x[i][0] - (x[0][0] + x[j][0])/2)\
                  * (x[i][0] - (x[0][0] + x[j][0])/2)\
                  + (x[i][1] - (x[0][1] + x[j][1])/2)\
                  * (x[i][1] - (x[0][1] + x[j][1])/2));
      }
      else
      {
        h[i][j] = sqrt((x[i][0] - (x[j+1][0] + x[j][0])/2)\
                  * (x[i][0] - (x[j+1][0] + x[j][0])/2)\
                  + (x[i][1] - (x[j+1][1] + x[j][1])/2)\
                  * (x[i][1] - (x[j+1][1] + x[j][1])/2));
      }
    //printf("h[%d][%d] = %e\n", i, j, h[i][j]);
    }
    if (i == N-1)
    {
      kappa_bar[i] = 0.5*(kappa[i] + kappa[0]);
      d[i] = sqrt((x[0][0] - x[i][0])*(x[0][0] - x[i][0])\
          + (x[0][1] - x[i][1])*(x[0][1] - x[i][1]));
    }
    else
    {
      kappa_bar[i] = 0.5*(kappa[i] + kappa[i+1]);
      d[i] = sqrt((x[i+1][0] - x[i][0])*(x[i+1][0] - x[i][0])\
           + (x[i+1][1] - x[i][1])*(x[i+1][1] - x[i][1]));
    }
  }
  
  double kappa_breve_temp;
  for (int i = 0; i < N; i++)
  {
    kappai_breve[i] = 0.;
  }

  for (int i = 0; i < N; i++)
  {
    kappa_breve_temp = 0.;

    for (int j = 0; j < N; j++)
      kappa_breve_temp += d[j]/(h[i][j]*h[i][j]);

    for (int j = 0; j < N; j++)
    {
      kappai_breve[i] += (d[j]*fabs(kappa_bar[j])/(h[i][j]*h[i][j]))*(1./kappa_breve_temp);

    }
    kappai_tilde[i] = pow((kappai_breve[i]*L), a)/(v*L)\
                    + SQRTTWO*kappai_breve[i];

    kappa_breve_temp = 0;
  }

  for (int i = 0; i < N-1; i++)
  {
    kappai_hat[i] = 0.5*(kappai_tilde[i] + kappai_tilde[i+1]);

    rho[i] = kappai_hat[i]/(1. + epsilon*kappai_hat[i]/SQRTTWO);

    sigmai[i] = rho[i]*d[i];
  } 
  
  kappai_hat[N-1] = 0.5*(kappai_tilde[N-1] + kappai_tilde[0]);
  rho[N-1] = kappai_hat[N-1]/(1. + epsilon*kappai_hat[N-1]/SQRTTWO);
  sigmai[N-1] = rho[N-1]*d[N-1];
  //*kappai_hat[N-1]
  double q, rest, p, dp, p_min, rest_dp;
  int i_min;
  q = 0.;
  dp = 0.005;
  p = 0.;
  
  int N_tilde;
  for (int i = 0; i < N; i++)
    q += sigmai[i];
  // printf("q = %e\n", q);
  N_tilde = round(q) + 2; 
  printf("N = %d         N_tilde = %d\n", N, N_tilde);
  for (int i = 0; i < N; i++)
  {
    sigmai_prim[i] = sigmai[i]*N_tilde/q;
    //printf("sigmai_prim[%d] = %e\n", i, sigmai_prim[i]);
  }
  
  // Copy x into temporary matrix
  x_temp = (double**)malloc(N*sizeof(double*));
  for (int i = 0; i < N; i++)
  {
  		x_temp[i] = (double*)malloc(2*sizeof(double));
  }
  
  for (int i = 0; i < N; i++)
  {
  		memcpy(&x_temp[i], &x[i], sizeof(x[i]));
  }
  
  printf("x_temp[0][0] = %lf       x_temp[0][1] = %lf\n", x_temp[0][0], x_temp[0][1]);
  printf("x_[0][0] = %lf       x_[0][1] = %lf\n", x[0][0], x[0][1]);
  
  printf("before realloc\n");
 	x = (double**)realloc(x, N_tilde*sizeof(double*));

 	for (int i = 0; i < N_tilde; i++)
 	{
 		x[i] = realloc(x[i], 2*sizeof(double));
 		 		
 	}
 	
 	printf("after realloc\n");
 	for (int i = 1; i < N_tilde; i++)
 	{
 		x[i][0] = 0;
 		x[i][1] = 0;
 	}
  
  printf("after zeroing \n");
  
 
  for (int j = 1; j < N_tilde; j++)
  {
  	// Assume minimum at i=0, p=0
  	i_min = 0;
  	p_min = 0.;
    for (int i = 0; i < N; i++)
    {
      while (p < 1-dp)
      {
      	rest = 0.;
      	rest_dp = 0.;
      	for (int l = 1; l < i; l++)
      	{
          rest += (sigmai_prim[l] + sigmai_prim[i]*p);
          rest_dp += (sigmai_prim[l] + sigmai_prim[i]*(p+dp));
        }
          
        if (j-1 > rest && j-1 < rest_dp)
				{
				//printf("test\n");
					i_min = i;
					p_min = p;
				}
        p += dp;
      }
	    p = 0;
    }
    x[j][0] = x_temp[i_min][0] + p_min*t[i_min][0] + (mu[i_min]*p_min + \
      					beta[i_min]*p_min*p_min + gamma[i_min]*p_min*p_min*p_min)*n[i_min][0];  
	  x[j][1] = x_temp[i_min][1] + p_min*t[i_min][1] + (mu[i_min]*p_min + \
      					beta[i_min]*p_min*p_min + gamma[i_min]*p_min*p_min*p_min)*n[i_min][1];
  }
   	for (int i = 0; i < N_tilde; i++)
      printf("x[%d][0] =  %e, x[%d][1] =  %e \n", i, x[i][0],\
  i, x[i][1]);
  // Free
  
  free(kappa_bar);
  free(d);
  free(kappai_breve);
  free(kappai_tilde);
  free(kappai_hat);
  free(rho);
  free(sigmai);
  free(sigmai_tilde);
  free(sigmai_prim);
  for (int k = 0; k < N; k++)
    free(h[k]);
  free(h); 
  for (int i = 0; i < N; i++)
  {
  	free(x_temp[i]);
  }
  *pN = N_tilde;
  *px = x;
  free(x_temp);
  printf("Exiting points_reloc\n");
  
  return;
}


double runge_kutta45(double** x, double** dxdt_k1, double** dxdt_k2, double** dxdt_k3, double** dxdt_k4, double** dxdt_k5, double** dxdt_k6, double** dxdt_RK4, double** dxdt_RK5, double tol_rk45_time, double dt, int M, int N, double* mu, double* beta, double* gamma, double** t, double** n, double alpha, double tol_rk45_space, double h)
{
  double F, delta, theta, tpi, dt_new;
  theta = -1.;
  F = dt*theta/(TWOPI);
  tpi = theta/(TWOPI);
  double R[2];
  double R_old[2];
  double R_max;
 // R[0] = 0.;
 // R[1] = 0.;
  R_max = 0.;
  double **x_temp, **x_RK4, **x_RK5;
  //x_temp = (double*)malloc(2*N*sizeof(double));
  x_temp = (double**)malloc(N*sizeof(double*));
  x_RK4 = (double**)malloc(N*sizeof(double*));
  x_RK5 = (double**)malloc(N*sizeof(double*));
  for (int j = 0; j < N; j++)
  {
    x_temp[j] = (double*)malloc(2*sizeof(double)); 
    x_RK4[j] = (double*)malloc(2*sizeof(double)); 
    x_RK5[j] = (double*)malloc(2*sizeof(double)); 
    x_temp[j][0] = 0.;
    x_temp[j][1] = 0.;
    x_RK4[j][0] = 0.;
    x_RK4[j][1] = 0.;
    x_RK5[j][0] = 0.;
    x_RK5[j][1] = 0.;
  } 
  // Runge-Kutta 45
  do
  {
    // Step 1 in RK
    for (int j = 0; j < N; j++)
    {
      compute_derivative(dxdt_k1[j], x, mu, beta, gamma, t, n, M, N, alpha, h, tol_rk45_space, j);
      dxdt_k1[j][0] = F*dxdt_k1[j][0];
      dxdt_k1[j][1] = F*dxdt_k1[j][1];
      //printf("dxdt dot x = %e\n", dxdt_k1[j][0]*x[j][0] + dxdt_k1[j][1]*x[j][1]);
    }
    for (int j = 0; j < N; j++)
    {
      x_temp[j][0] = x[j][0] + 0.25*dxdt_k1[j][0];
      x_temp[j][1] = x[j][1] + 0.25*dxdt_k1[j][1];
    }
    
    // Step 2 in RK
    for (int j = 0; j < N; j++)
    {
      compute_derivative(dxdt_k2[j], x_temp, mu, beta, gamma, t, n, M, N, alpha, h, tol_rk45_space, j);
      dxdt_k2[j][0] = F*dxdt_k2[j][0];
      dxdt_k2[j][1] = F*dxdt_k2[j][1];
    }
    for (int j = 0; j < N; j++)
    {
      x_temp[j][0] = x[j][0] + 3.0/32.0*dxdt_k1[j][0] + 9.0/32.0*dxdt_k2[j][0];
      x_temp[j][1] = x[j][1] + 3.0/32.0*dxdt_k1[j][1] + 9.0/32.0*dxdt_k2[j][1];
    }
		
    // Step 3 in RK
    for (int j = 0; j < N; j++)
    {
      compute_derivative(dxdt_k3[j], x_temp, mu, beta, gamma, t, n, M, N, alpha, h, tol_rk45_space, j);
      dxdt_k3[j][0] = F*dxdt_k3[j][0];
      dxdt_k3[j][1] = F*dxdt_k3[j][1];
    }
    for (int j = 0; j < N; j++)
    {
      x_temp[j][0] = x[j][0] + 1932.0/2197.0*dxdt_k1[j][0] - 7200.0/2197.0*dxdt_k2[j][0] + 7296.0/2197.0*dxdt_k3[j][0];
      x_temp[j][1] = x[j][1] + 1932.0/2197.0*dxdt_k1[j][1] - 7200.0/2197.0*dxdt_k2[j][1] + 7296.0/2197.0*dxdt_k3[j][1];
    }
    
    // Step 4 in RK
    for (int j = 0; j < N; j++)
    {
      compute_derivative(dxdt_k4[j], x_temp, mu, beta, gamma, t, n, M, N, alpha, h, tol_rk45_space, j);
      dxdt_k4[j][0] = F*dxdt_k4[j][0];
      dxdt_k4[j][1] = F*dxdt_k4[j][1];
    }
    for (int j = 0; j < N; j++)
    {
      x_temp[j][0] = x[j][0] + 439.0/216.0*dxdt_k1[j][0] - 8.0*dxdt_k2[j][0] + 3680.0/513.0*dxdt_k3[j][0] - 845.0/4104.0*dxdt_k4[j][0];
      x_temp[j][1] = x[j][1] + 439.0/216.0*dxdt_k1[j][1] - 8.0*dxdt_k2[j][1] + 3680.0/513.0*dxdt_k3[j][1] - 845.0/4104.0*dxdt_k4[j][1];
    }
    
    // Step 5 in RK
    for (int j = 0; j < N; j++)
    {
      compute_derivative(dxdt_k5[j], x_temp, mu, beta, gamma, t, n, M, N, alpha, h, tol_rk45_space, j);
      dxdt_k5[j][0] = F*dxdt_k5[j][0];
      dxdt_k5[j][1] = F*dxdt_k5[j][1];
    }
    for (int j = 0; j < N; j++)
    {
      x_temp[j][0] = x[j][0] - 8.0/27.0*dxdt_k1[j][0] + 2.0*dxdt_k2[j][0] - 3544.0/2565.0*dxdt_k3[j][0] + 1859.0/4104*dxdt_k4[j][0] - 11.0/40.0*dxdt_k5[j][0];
      x_temp[j][1] = x[j][1] - 8.0/27.0*dxdt_k1[j][1] + 2.0*dxdt_k2[j][1] - 3544.0/2565.0*dxdt_k3[j][1] + 1859.0/4104*dxdt_k4[j][1] - 11.0/40.0*dxdt_k5[j][1];
    }
    
    // Step 6 in RK
    for (int j = 0; j < N; j++)
    {
      compute_derivative(dxdt_k6[j], x_temp, mu, beta, gamma, t, n, M, N, alpha, h, tol_rk45_space, j);
      dxdt_k6[j][0] = F*dxdt_k6[j][0];
      dxdt_k6[j][1] = F*dxdt_k6[j][1];
    }
    
    for (int j = 0; j < N; j++)
    {
		// RK4 approx
  	dxdt_RK4[j][0] = 25.0*dxdt_k1[j][0]/216.0 + 1408.0*dxdt_k3[j][0]/2565.0 + 2197.0*dxdt_k4[j][0]/4104.0 - 0.2*dxdt_k5[j][0];
  	dxdt_RK4[j][1] = 25.0*dxdt_k1[j][1]/216.0 + 1408.0*dxdt_k3[j][1]/2565.0 + 2197.0*dxdt_k4[j][1]/4104.0 - 0.2*dxdt_k5[j][1];
    
    // RK5 approx
		dxdt_RK5[j][0] = 16.0*dxdt_k1[j][0]/135.0 + 6656.0*dxdt_k3[j][0]/12825.0 + 28561.0*dxdt_k4[j][0]/56430.0\
		 	-9.0*dxdt_k5[j][0]/50.0 + 2.0*dxdt_k6[j][0]/55.0;
		dxdt_RK5[j][1] = 16.0*dxdt_k1[j][1]/135.0 + 6656.0*dxdt_k3[j][1]/12825.0 + 28561.0*dxdt_k4[j][1]/56430.0\
		 	-9.0*dxdt_k5[j][1]/50.0 + 2.0*dxdt_k6[j][1]/55.0;

    // RK x approximations
    x_RK4[j][0] = x[j][0] + dxdt_RK4[j][0];
    x_RK4[j][1] = x[j][1] + dxdt_RK4[j][1];
    
    x_RK5[j][0] = x[j][0] + dxdt_RK5[j][0];
    x_RK5[j][1] = x[j][1] + dxdt_RK5[j][1];
    }
    
  	// Compute maximum error
    /*R_old[0] = fabs(x_RK5[0][0] - x_RK4[0][0]);
    R_old[1] = fabs(x_RK5[0][1] - x_RK4[0][1]);
    for (int j = 1; j < N; j++)
    {
	    R[0] = fabs(x_RK5[j][0] - x_RK4[j][0]);
      if (R[0] > R_old[0])
        R_old[0] = R[0];
	  	  
      R[1] = fabs(x_RK5[j][1] - x_RK4[j][1]);
      if (R[1] > R_old[1])
        R_old[1] = R[1];
    }
    R_max = R_old[0];
    if (R_old[1] > R_old[0])
      R_max = R_old[1];*/
     // Compute average error
    R[0] = 0.;
    R[1] = 0.;
    for (int i = 0; i < N; i++)
    {
      R[0] = R[0] + fabs(x_RK5[i][0] - x_RK4[i][0]);
      R[1] = R[1] + fabs(x_RK5[i][1] - x_RK4[i][1]);
    }
    
    R_max = 0.5*(R[0] + R[1])/N;
    
    
    printf("dt1 = %e\n", dt);
    // Calculate update factor
    if (R_max > tol_rk45_time)
      dt = 0.5*dt;
    
    // Make step smaller
    F = dt*tpi;
    printf("R_max = %e\n", R_max);
    printf("dt2 = %e\n\n", dt);
  } while (R_max > tol_rk45_time);
  dt_new = 0.9*dt*sqrt(sqrt(tol_rk45_time/R_max));
  printf("dt_new = %e\n", dt_new);
  
  // Update
  for (int j = 0; j < N; j++)
  {
    x[j][0] = x[j][0] + dxdt_RK5[j][0];
    x[j][1] = x[j][1] + dxdt_RK5[j][1];
  }
  
  for (int j = 0; j < N; j++)
  {
    free(x_temp[j]);
    free(x_RK4[j]);
    free(x_RK5[j]);
  }
  free(x_temp);
  free(x_RK4);
  free(x_RK5);
  
  printf("X dot dxdt = %e \n", x[0][0]*dxdt_RK5[0][0]+x[0][1]*dxdt_RK5[0][1]);

  return dt_new;
}

double compute_area(double** x, int start, int stop, double** t, double** n,\
                    double* mu, double* beta, double* gamma)
{
  double area, area2;
  area = 0;
  area2 = 0;
  for (int i = start; i < stop; i++)
  {
    area += (t[i][1] + (mu[i])*n[i][1])*x[i][0];
    
    //printf("%e, %e, %e, %e, %e\n",x[i][1], t[i][0]+mu[i]*n[i][0], t[i][1]+mu[i]*n[i][1], x[i+1][1] - x[i][1], fabs(x[i+1][1] - x[i][1] - t[i][1]+mu[i]*n[i][1]));
  }
  
  /*+ 2*beta[i] + 3*gamma[i]*/
  /*for (int j = start; j < stop-1; j++)
    area2 = area2 + x[j][0]*(x[j+1][1] - x[j][1]);
  area2 =  area2 + x[stop-1][0]*(x[start][1] - x[stop-1][1]);
  //area = area;
  printf("area2 = %e\n", area2);*/
  return area;
}
