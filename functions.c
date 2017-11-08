#include "functions.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#define TWOPI 6.2831853071795864769
#define SQRTTWO 1.4142135623730950588
#define PRINT 0


void interpolate(double** x, int N, int P, int n_dim, double** t, double** n, double* p,\
                 double* eta, double* d, double* kappa, double* kappa_den, double* mu,\
                 double* beta, double* gamma)
{

  for (int i = 0; i < P; i++)
    p[i] = (double)i/P;

  // Calculate t an n
  for (int j = 0; j < N-1; j++) {
    t[j][0] = x[(j+1)*P][0] - x[j*P][0];
    t[j][1] = x[(j+1)*P][1] - x[j*P][1];
    n[j][0] = -t[j][1];
    n[j][1] = t[j][0];
    d[j] = sqrt(t[j][0]*t[j][0] + t[j][1]*t[j][1]);
  }  
  // Special case j = N-1
  t[N-1][0] = x[0][0] - x[(N-1)*P][0];
  t[N-1][1] = x[0][1] - x[(N-1)*P][1];
  n[N-1][0] = -t[N-1][1];
  n[N-1][1] = t[N-1][0];
  d[N-1] = sqrt(t[N-1][0]*t[N-1][0] + t[N-1][1]*t[N-1][1]);
  
  // kappa local curvature
  kappa_den[0] = d[N-1]*d[N-1]*t[0][0] + d[0]*d[0]*t[N-1][0];
  kappa_den[1] = d[N-1]*d[N-1]*t[0][1] + d[0]*d[0]*t[N-1][1];
  kappa[0] = 2*(t[N-1][0]*t[0][1] - t[N-1][1]*t[0][0])\
    /sqrt(kappa_den[0]*kappa_den[0] + kappa_den[1]*kappa_den[1]);
  
  for (int j = 1; j < N; j++) {
    // kappa local curvature
    kappa_den[0] = (d[j-1]*d[j-1]*t[j][0] + d[j]*d[j]*t[j-1][0]);
    kappa_den[1] = (d[j-1]*d[j-1]*t[j][1] + d[j]*d[j]*t[j-1][1]);
    kappa[j] = 2*(t[j-1][0]*t[j][1] - t[j-1][1]*t[j][0])\
      /sqrt(kappa_den[0]*kappa_den[0] + kappa_den[1]*kappa_den[1]);
  }
  
  // Construct the cubic interpolation
  for (int j = 0; j < N; j++) {
    
    // Cubic interpolation coefficients
    mu[j] = -(double)1/3*d[j]*kappa[j] - (double)1/6*d[j]*kappa[j+1];
    beta[j] = 0.5*d[j]*kappa[j];
    gamma[j] = (double)1/6*d[j]*(kappa[j+1] - kappa[j]);
    
    for (int i = 0; i < P; i++) {
      eta[i] = mu[j]*p[i] + beta[j]*p[i]*p[i] + gamma[j]*p[i]*p[i]*p[i];
      x[j*P + i][0] = x[j*P][0] + p[i]*t[j][0] + eta[i]*n[j][0];
      x[j*P + i][1] = x[j*P][1] + p[i]*t[j][1] + eta[i]*n[j][1];
    }
  }

  return;
}

void local_coeffs(int NP, double** x, double** t_loc, double** n_loc, double* mu_loc, double* beta_loc, double* gamma_loc)
{
  double* d_loc = (double*)malloc(NP*sizeof(double));
  double* kappa_loc = (double*)malloc(NP*sizeof(double));
  double kappa_den[2];
  //	printf("Entering Local coeffs\n");
  // Calculate t an n
  for (int i = 0; i < NP-1; i++) {
    t_loc[i][0] = x[(i+1)][0] - x[i][0];
    t_loc[i][1] = x[(i+1)][1] - x[i][1];
    n_loc[i][0] = -t_loc[i][1];
    n_loc[i][1] = t_loc[i][0];
    d_loc[i] = sqrt(t_loc[i][0]*t_loc[i][0] + t_loc[i][1]*t_loc[i][1]);
  }  
  // Special case j = P-1
  t_loc[NP-1][0] = x[0][0] - x[NP-1][0];
  t_loc[NP-1][1] = x[0][1] - x[NP-1][1];
  n_loc[NP-1][0] = -t_loc[NP-1][1];
  n_loc[NP-1][1] = t_loc[NP-1][0];
  d_loc[NP-1] = sqrt(t_loc[NP-1][0]*t_loc[NP-1][0] + t_loc[NP-1][1]*t_loc[NP-1][1]);
  
  // kappa local curvature
  kappa_den[0] = d_loc[NP-1]*d_loc[NP-1]*t_loc[0][0] + d_loc[0]*d_loc[0]*t_loc[NP-1][0];
  kappa_den[1] = d_loc[NP-1]*d_loc[NP-1]*t_loc[0][1] + d_loc[0]*d_loc[0]*t_loc[NP-1][1];
  kappa_loc[0] = 2*(t_loc[NP-1][0]*t_loc[0][1] - t_loc[NP-1][1]*t_loc[0][0])\
    /sqrt(kappa_den[0]*kappa_den[0] + kappa_den[1]*kappa_den[1]);
  
  for (int i = 1; i < NP; i++) {
    // kappa local curvature
    kappa_den[0] = (d_loc[i-1]*d_loc[i-1]*t_loc[i][0] + d_loc[i]*d_loc[i]*t_loc[i-1][0]);
    kappa_den[1] = (d_loc[i-1]*d_loc[i-1]*t_loc[i][1] + d_loc[i]*d_loc[i]*t_loc[i-1][1]);
    kappa_loc[i] = 2*(t_loc[i-1][0]*t_loc[i][1] - t_loc[i-1][1]*t_loc[i][0])\
      /sqrt(kappa_den[0]*kappa_den[0] + kappa_den[1]*kappa_den[1]);
  }
  
  for (int i = 0; i < NP; i++) {
    // Cubic interpolation coefficients
    mu_loc[i] = -(double)1/3*d_loc[i]*kappa_loc[i] - (double)1/6*d_loc[i]*kappa_loc[i+1];
    beta_loc[i] = 0.5*d_loc[i]*kappa_loc[i];
    gamma_loc[i] = (double)1/6*d_loc[i]*(kappa_loc[i+1] - kappa_loc[i]);
  }
  
  free(d_loc);
  free(kappa_loc);
  //	printf("Exiting local coeffs\n");
  return;
}

void autder(double* f, double* c_coeff, double alpha, int order)
{
  //	printf("Entering autoder\n");
  // Allocate memory for temporary coefficients
  double* a_ = (double*)malloc(order*sizeof(double));
  for (int n = 1; n < order; n++)
    a_[n] = 0;
  a_[0] = 1;

  // Calculate temporary coefficients
  for (int n = 1; n < order; n++)
  {
    for (int j = 1; j <= n; j++)
    {
      a_[n] -= c_coeff[j]*a_[n-j];
    }
  }
  
  f[0] = 1;
  
  // Calculate Taylor coefficients of order, "order" :D
  for (int n = 1; n < order; n++)
  {
    for (int j = 0; j < n; j++)
    {
      f[n] += (n*alpha - j*(alpha + 1))*a_[n-j]*f[j];
    }
    f[n] /= (n+1); // What exactly should it be here?
  }
  
  free(a_);
  //	printf("Exiting autoder\n");
  return;
}

void compute_derivative(double* dxdt, double** x, double* mu, double* beta, double* gamma, double** t, double** n, int N, int P, double alpha, double h, double eps, int j)
{
  //	printf("Entering compute derivative\n");
  double d_x, d_ni, d_ti, d_xi;
  int order = 11;
  int ij, jp; // Index i + j
  double Q = 0.01;
  double f = 1/sqrt(Q);
  jp = j*P;
  // Generate coefficients
  double c[order];
  double g[order];
  double poly_coeff_c[order];
  double poly_coeff_g[order];
  
  // Evolve the contour integrals
  for (int i = 0; i < P; i++) { // i < ???
    ij = i+jp;
    #if PRINT    
      printf("ij = %d\n", ij);
      printf("jp = %d\n", jp);
    #endif    
    if (jp+P == N*P)
    {
      d_x = sqrt((x[0][0] - x[jp][0])*(x[0][0] - x[jp][0]) + (x[0][1] - x[jp][1])*(x[0][1] - x[jp][1]));
      d_ti = -((x[jp][0] - x[0][0])*t[j][0] + (x[jp][1] - x[0][1])*t[j][1]);
      d_ni = -((x[jp][0] - x[0][0])*n[j][0] + (x[jp][1] - x[0][1])*n[j][1]);
    } else {
      d_x = sqrt((x[ij][0] - x[jp][0])*(x[ij][0] - x[jp][0]) + (x[ij][1] - x[jp][1])*(x[ij][1] - x[jp][1]));
      d_ti = -((x[jp][0] - x[ij][0])*t[j][0] + (x[jp][1] - x[ij][1])*t[j][1]);
      d_ni = -((x[jp][0] - x[ij][0])*n[j][0] + (x[jp][1] - x[ij][1])*n[j][1]);
    } 
    
    // Distance between x_i and x_{i+1}
    if (ij+1 == N*P) {
      d_xi = sqrt((x[0][0] - x[ij][0])*(x[0][0] - x[ij][0]) + (x[0][1] - x[ij][1])*(x[0][1] - x[ij][1]));
    } else {
      d_xi = sqrt((x[ij+1][0] - x[ij][0])*(x[ij+1][0] - x[ij][0]) + (x[ij+1][1] - x[ij][1])*(x[ij+1][1] - x[ij][1]));
    }
    
    // Initialize Taylor coefficients
    for (int n = 0; n < order; n++) {
      c[n] = 0;
      g[n] = 0;
      poly_coeff_c[n] = 0;
      poly_coeff_g[n] = 0;
    }
    poly_coeff_c[0] = 1;
    poly_coeff_g[0] = 1;
    
    // Evaluate integrals
    if (ij+1 != N*P) {
      if (jp == ij) {
        //Case 1: Use formula (29)
        #if PRINT
          printf("Case 1\n");
        #endif
        
        // Generate Taylor coefficients
        poly_coeff_c[1] = 2*mu[j]*beta[j]/(1 + mu[j]*mu[j]);
        poly_coeff_c[2] = (beta[j]*beta[j] + 2*mu[j]*gamma[j])/(1 + mu[j]*mu[j]);
        poly_coeff_c[3] = 2*beta[j]*gamma[j]/(1 + mu[j]*mu[j]);
        poly_coeff_c[4] = gamma[j]*gamma[j]/(1 + mu[j]*mu[j]);

        autder(c, poly_coeff_c, alpha, order); // Should these be divided by two maybe?
      
        evaluate_integral(dxdt, mu[j], beta[j], gamma[j], t[j], n[j], c, alpha); // Look at inputs in these functions
      
      } else if (jp+P == ij) {
        // Case 2: Use formula (29) with shifted params
        #if PRINT
          printf("Case 2\n");      
        #endif
        
        mu[j] = mu[j] + 2*beta[j] + 3*gamma[j];
        beta[j] = -beta[j] - 3*gamma[j];
      
        // Generate Taylor coefficients
        poly_coeff_c[1] = 2*mu[j]*beta[j]/(1 + mu[j]*mu[j]);
        poly_coeff_c[2] = (beta[j]*beta[j] + 2*mu[j]*gamma[j])/(1 + mu[j]*mu[j]);
        poly_coeff_c[3] = 2*beta[j]*gamma[j]/(1 + mu[j]*mu[j]);
        poly_coeff_c[4] = gamma[j]*gamma[j]/(1 + mu[j]*mu[j]);
      
        autder(c, poly_coeff_c, alpha, order); // Should these be divided by two maybe?
        
        evaluate_integral(dxdt, mu[j], beta[j], gamma[j], t[j], n[j], c, alpha); // Look at inputs in these functions
        //printf("dxdt[%d] = %lf\n", j, dxdt[0]);
      } else if (sqrt((x[jp][0] - x[ij][0])*(x[jp][0] - x[ij][0]) + (x[jp][1] - x[ij][1])*(x[jp][1] - x[ij][1])) > f*d_xi) {
        #if PRINT   
          printf("Case 3\n");
        #endif
        // Case 3: Use formula (31)

        // Generate Taylor coefficients
        poly_coeff_g[1] = (d_ti + d_ni*mu[j+i])/(d_x*d_x);
        poly_coeff_g[2] = ((t[j][0]*t[j][0] + t[j][1]*t[j][1]) + mu[j]*mu[j]*(n[j][0]*n[j][0] + n[j][1]*n[j][1]) + d_x*beta[j])/(d_x*d_x);
        poly_coeff_g[3] = (2*mu[j]*beta[j]*(n[j][0]*n[j][0] + n[j][1]*n[j][1]) + d_ni*gamma[j])/(d_x*d_x);
        poly_coeff_g[4] = ((beta[j] + 2*mu[j]*gamma[j])*(n[j][0]*n[j][0] + n[j][1]*n[j][1]))/(d_x*d_x);
        poly_coeff_g[5] = 2*beta[j]*gamma[j]*(n[j][0]*n[j][0] + n[j][1]*n[j][1])/(d_x*d_x);
        poly_coeff_g[6] = gamma[j]*gamma[j]/6/(d_x*d_x);
    
        autder(g, poly_coeff_g, alpha, order);

        evaluate_integral_g(dxdt, mu[j], beta[j], gamma[j], d_x, d_ni, d_ti, t[j], n[j], g, alpha);

      } else {
        #if PRINT
          printf("Case 4\n");
        #endif
        // Case 4: Use Runge-Kutta 4-5
        evaluate_integral_RK(dxdt, x[ij], x[jp], mu[j], beta[j], gamma[j], eps, h, t[j], n[j], alpha);
      }
    } else {
      // Edge case
      // Case 2: Use formula (29) with shifted params
        #if PRINT
          printf("Edge Case 2\n"); 
          printf("jp = %d \n", jp);     
        #endif
        mu[j] = mu[j] + 2*beta[j] + 3*gamma[j];
        beta[j] = -beta[j] - 3*gamma[j];
      
        // Generate Taylor coefficients
        poly_coeff_c[1] = 2*mu[j]*beta[j]/(1 + mu[j]*mu[j]);
        poly_coeff_c[2] = (beta[j]*beta[j] + 2*mu[j]*gamma[j])/(1 + mu[j]*mu[j]);
        poly_coeff_c[3] = 2*beta[j]*gamma[j]/(1 + mu[j]*mu[j]);
        poly_coeff_c[4] = gamma[j]*gamma[j]/(1 + mu[j]*mu[j]);
      
        autder(c, poly_coeff_c, alpha, order); // Should these be divided by two maybe?
      
        evaluate_integral(dxdt, mu[j], beta[j], gamma[j], t[j], n[j], c, alpha); // Look at inputs in these functions
    }
  }

  return;
}

void evaluate_integral(double* dxdt, double mu_i, double beta_i, double gamma_i, double* t_i, double* n_i, double* c, double alpha) {
  
  // Compute the integrals
  double first = 0;
  double second = 0;
  double p_coef, psq_coef;
  double t_abs = pow((t_i[0]*t_i[0] + t_i[1]*t_i[1]), alpha);
  
  double alpha_mu = pow((1 + mu_i*mu_i), 0.5*alpha);
  
  for (int n = 0; n < 11; n++) {
    first += c[n]/(double)(n - alpha + 1);
    
    p_coef = 2*beta_i/(double)(n - alpha + 2);
    psq_coef = 3*gamma_i/(double)(n - alpha + 3);
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
  double first = 0;
  double second = 0;
  double p_coef, psq_coef;

  for (int n = 0; n < 11; n++) {
    first += g[n];
    
    p_coef = 2*beta_i/(double)(n - alpha + 2);
    psq_coef = 3*gamma_i/(double)(n - alpha + 3);
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
                          double eps, double h, double* t_i, double* n_i, double alpha)
{
  double first = evaluate_integral1_RK(x_i, x_j, eps, h, t_i, n_i, mu_i, beta_i, gamma_i, alpha);
  double second = evaluate_integral2_RK(x_i, x_j, eps, h, t_i, n_i, mu_i, beta_i, gamma_i, alpha);
  
  dxdt[0] += first*(t_i[0] + mu_i*n_i[0]) + second*n_i[0];
  dxdt[1] += first*(t_i[1] + mu_i*n_i[1]) + second*n_i[1];
  
  return;
}

double evaluate_integral1_RK(double* x_i, double* x_j, double eps, double h, double* t_i,\
									double* n_i, double mu_i, double beta_i, double gamma_i, double alpha)
{
  
	double p = 0;
	double p_end = 1; 
  double k1, k2, k3, k4, k5, k6;
	double Y = 0, Y1, Y2, R, delta, p_temp;
  
  /*
  // Print to file
  char str[80] = "../w_spike.csv";
  FILE* f = fopen(str, "wb");
    
  while (p < 1) {
	  printf("p = %lf\n", p);
    w = integrand1(x_i, x_j, p, t_i, n_i, mu_i, beta_i, gamma_i, alpha);
    fprintf(f, "%lf,%lf\n", p, Y);
	  printf("w = %lf\n\n", Y);
    p = p + 0.001;
  }
  fclose(f);
  sleep(5);
  */
  
	while (p < p_end) 
	{
		if ((p_end - p) < h)
		{
			h = p_end - p;
		}
    
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
		R = sqrt((Y2 - Y1)*(Y2 - Y1))/h;
    
		//Calculate update factor
		delta = 0.84*pow((eps/R), 0.25);
		
		// Check if to progress to next step or recalculate current step with
		// new step size. 		
		if (R <= eps)
		{
      // Update
			Y = Y1;
			p = p + h;
			h = delta*h;
		}
		else
		{
      // Make step smaller
			h = delta*h;
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

double evaluate_integral2_RK(double* x_i, double* x_j, double eps, double h,\
            double* t_i, double* n_i, double mu_i, double beta_i, double gamma_i, double alpha) 
{
  
	double p = 0;
	double p_end = 1; 
  double k1, k2, k3, k4, k5, k6;
	double Y = 0, Y1, Y2, R, delta, p_temp;
  
	while (p < p_end) 
	{
		if ((p_end - p) < h)
		{
			h = p_end - p;
		}
    
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
		R = sqrt((Y2 - Y1)*(Y2 - Y1))/h;
    
		//Calculate update factor
		delta = 0.84*pow((eps/R), 0.25);
		
		// Check if to progress to next step or recalculate current step with
		// new step size. 		
		if (R <= eps)
		{
      // Update
			Y = Y1;
			p = p + h;
			h = delta*h;
		}
		else
		{
      // Make step smaller
			h = delta*h;
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


/*
void runge_kutta_2D(double** dxdt, double** x)
{
  
  double k1, k2, k3, k4, k5, k6;
  double l1, l2, l3, l4, l5, l6;
	double w = int_IC, w1, w2, R, delta, w_temp, p_temp;
	int i = 0;
  
  k1 = dxdt[j][0];
  l1 = dxdt[j][1];
  
  		
  
  return;
}
*/

/* To be implemented later
void points_reloc(double** x, double NP) {
  
  double epsilon = 10^-6;
  double L = 3.0;
  double a = 2.0/3.0;
  double v = 0.05; // 0.01, 0.03, 0.05 (Smaller value => Larger densities of points on the curve)
  
  // Either calculate all d[j]'s here or use input
  // Same goes with kappa and h
  
  double kappa_bar = 0.5*(kappa[j] + kappa[j+1]);
  double kappai_breve = 0;
  double kappa_breve_temp = 0;
  for (int j = 0; j < NP; j++)
    kappa_breve_temp += d[j]/(h[i][j]*h[i][j]);
  for (int j = 0; j < NP; j++)
    kappai_breve += d[j]*abs(kappa_bar[j])/(h[i][j]*h[i][j])/kappa_breve_temp;
  double kappai_tilde = pow((kappai_breve*L), a)/(v*L) + SQRTTWO*kappai_breve;
  double kappai_hat = 0.5*(kappai_tilde + kappaip_tilde);
  
  double density = kappai_hat*kappa_hat/(1 + epsilon*kappai_hat/SQRTTWO);
  double sigmai = 
  
  return;
}*/
