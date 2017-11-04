#include "functions.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#define TWOPI 6.2831853071795864769
#define SQRTTWO 1.4142135623730950588

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
  double d_loc[NP], kappa_loc[NP], kappa_den[2];
  
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
  
  return;
}

void autder(double* f, double* c_coeff, double alpha, int order)
{
  // Allocate memory for temporary coefficients
  double* a_ = (double*)malloc(order*sizeof(double));
  a_[0] = 1/c_coeff[0];
  
  // calculate temporary coefficients
  for (int n = 1; n <= order; n++)
  {
    for (int j = 1; j <= n; j++)
    {
      a_[n] -= c_coeff[j]*a_[n-j];
    }
    a_[n] /= c_coeff[0];
  }
  
  f[0] = 1;
  
  // Calculate Taylor coefficients of order, "order" :D
  for (int n = 1; n <= order; n++)
  {
    for (int j = 0; j < n; j++)
    {
      f[n] += (n*alpha - j*(alpha + 1))*a_[n-j]*f[j];
    }
    f[n] /= (n*a_[0]);
  }
  
  free(a_);
  
  return;
}

void compute_derivative(double* dxdt, double** x, double* mu, double* beta, double* gamma, double** t, double** n, int P, double alpha, double h, double eps, int j)
{
  
  double d_x, d_ni, d_ti, d_xi;
  int int_IC = 0;      
  int order = 11;
  int ij; // Index i + j
  double Q = 0.025;
  double f = 1/sqrt(Q);

  // Generate coefficients
  double c[11];
  double g[11];
  double poly_coeff_c[11];
  double poly_coeff_g[11];
  for (int i = 0; i < order; i++) {
    c[i] = 0;
    g[i] = 0;
    poly_coeff_c[i] = 0;
    poly_coeff_g[i] = 0;
  }
  poly_coeff_c[0] = 1;
  poly_coeff_g[0] = 1;
  
  // Evolve the contour integrals
  for (int i = 0; i < P; i++) { // i < ???
    printf("i = %d\n", i);
    ij = i+j;
    d_x = sqrt((x[j][0] - x[ij][0])*(x[j][0] - x[ij][0]) + (x[j][1] - x[ij][1])*(x[j][1] - x[ij][1]));
    d_ti = -((x[j][0] - x[ij][0])*t[j][0] + (x[j][1] - x[ij][1])*t[j][1]);
    d_ni = -((x[j][0] - x[ij][0])*n[j][0] + (x[j][1] - x[ij][1])*n[j][1]);
    
    // Distance between x_i and x_{i+1}
    d_xi = sqrt((x[ij+1][0] - x[ij][0])*(x[ij+1][0] - x[ij][0]) + (x[ij+1][1] - x[ij][1])*(x[ij+1][1] - x[ij][1]));
    
    // Generate Taylor coefficients
    poly_coeff_c[1] = 2*mu[ij]*beta[ij]/(1 + mu[ij]*mu[ij]);
    poly_coeff_c[2] = (beta[ij]*beta[ij] + 2*mu[ij]*gamma[ij])/(1 + mu[ij]*mu[ij]);
    poly_coeff_c[3] = 2*beta[ij]*gamma[ij]/(1 + mu[ij]*mu[ij]);
    poly_coeff_c[4] = gamma[ij]*gamma[ij]/(1 + mu[ij]*mu[ij]);
      
    poly_coeff_g[1] = (d_ti + d_ni*mu[j+i])/(d_x*d_x);
    poly_coeff_g[2] = ((t[ij][0]*t[ij][0] + t[ij][1]*t[ij][1]) + mu[ij]*mu[ij]*(n[ij][0]*n[ij][0] + n[ij][1]*n[ij][1]) + d_x*beta[ij])/(d_x*d_x);
    poly_coeff_g[3] = (2*mu[ij]*beta[ij]*(n[ij][0]*n[ij][0] + n[ij][1]*n[ij][1]) + d_ni*gamma[ij])/(d_x*d_x);
    poly_coeff_g[4] = ((beta[ij] + 2*mu[ij]*gamma[ij])*(n[ij][0]*n[ij][0] + n[ij][1]*n[ij][1]))/(d_x*d_x);
    poly_coeff_g[5] = 2*beta[ij]*gamma[ij]*(n[ij][0]*n[ij][0] + n[ij][1]*n[j][1])/(d_x*d_x);
    poly_coeff_g[6] = gamma[ij]*gamma[ij]/6/(d_x*d_x);
    
    autder(c, poly_coeff_c, alpha, order); // Should these be divided by two maybe?
    autder(g, poly_coeff_g, 0.5*alpha, order);
    
    // Evaluate integrals
    if ((x[j][0] == x[ij][0]) && (x[j][1] == x[ij][1])) {
      printf("Case 1\n");
      // Case 1: Use formula (29)
      evaluate_integral(dxdt, mu[ij], beta[ij], gamma[ij], t[ij], n[ij], c, alpha); // Look at inputs in these functions
      
    } else if ((x[j][0] == x[ij-1][0]) && (x[j][1] == x[ij-1][1])) {
      printf("Case 2\n");
      // Case 2: Use formula (29) with shifted params
      evaluate_integral(dxdt, mu[ij] + 2*beta[ij] + 3*gamma[ij], - beta[ij] \
                                - 3*gamma[ij], gamma[ij], t[ij], n[ij], c, alpha);
    
    } else if (sqrt((x[j][0] - x[ij][0])*(x[j][0] - x[ij][0]) + (x[j][1] - x[ij][1])*(x[j][1] - x[ij][1])) > f*d_xi) {
      printf("Case 3\n");
      // Case 3: Use formula (31)      
      evaluate_integral_g(dxdt, mu[ij], beta[ij], gamma[ij], d_x, d_ni, d_ti, t[ij], n[ij], g, alpha);

    //} else if (sqrt((x[j][0] - x[i][0])*(x[j][0] - x[i][0]) + (x[j][1] - x[i][1])*(x[j][1] - x[i][1])) < 0.01) {
    } else {
      printf("Case 4\n");
      // Case 4: Use Runge-Kutta 4-5
      evaluate_integral_RK(dxdt, mu[ij], beta[ij], gamma[ij], eps, h, int_IC, t[ij], n[ij]);
    }
  }
   
  // What is theta?
  //derivative *= theta/(double)TWOPI;
  
  return;
}

void evaluate_integral(double* dxdt, double mu_i, double beta_i, double gamma_i, double* t_i, double* n_i, double* c, double alpha) {
  
  // Compute the integrals
  double first = 0;
  double second = 0;
  double p_coef, psq_coef;
  double t_abs = pow((t_i[0]*t_i[0] + t_i[1]*t_i[1]), 0.5*alpha);
  
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
  
  // COmpute the integrals
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

void evaluate_integral_RK(double* dxdt, double mu_i, double beta_i, double gamma_i,\
                          double eps, double h, double int_IC, double* t_i, double* n_i)
{

  printf("integral RK\n");
  double first = evaluate_integral1_RK(eps, h, int_IC, t_i, n_i, mu_i, beta_i, gamma_i);
  double second = evaluate_integral2_RK(eps, h, int_IC, t_i, n_i, mu_i, beta_i, gamma_i);
  printf("integral RK after\n");
  
  dxdt[0] += first*(t_i[0] + mu_i*n_i[0]) + second*n_i[0];
  dxdt[1] += first*(t_i[1] + mu_i*n_i[1]) + second*n_i[1];
  
  return;
}

double evaluate_integral1_RK(double eps, double h, double int_IC, double* t_i,\
									double* n_i, double mu_i, double beta_i, double gamma_i)
{
	double p = h;
	double p_end = 1; // Declare all variables.
  double k1, k2, k3, k4, k5, k6;
	double w = int_IC, w1, w2, R, delta, w_temp, p_temp;  
	int i = 0;
  
	
	while (p < p_end) 
	{
		if ((p_end - p) < h)
		{
			h = p_end - p;
		}
		
		k1 = h*integrand1(p, w, t_i, n_i, mu_i, beta_i, gamma_i);     //Func should be integrand1-function
		w_temp = w + 0.25*k1;
		p_temp = p + 0.25*h;
		
		k2 = h*integrand1(p_temp, w_temp, t_i, n_i, mu_i, beta_i, gamma_i);
		w_temp = w + 3.0*k1/32.0 + 9.0*k2/32.0;
		p_temp = p + 3.0*h/8.0;
		
		k3 = h*integrand1(p_temp, w_temp, t_i, n_i, mu_i, beta_i, gamma_i);
		w_temp = w + 1932.0*k1/2197.0 -7200.0*k2/2197.0 + 7296.0*k3/2197.0;
		p_temp = p + 12.0*h/13.0;
		
		k4 = h*integrand1(p_temp, w_temp, t_i, n_i, mu_i, beta_i, gamma_i);
		w_temp = w + 439.0*k1/216.0 - 8.0*k2 + 3680.0*k3/513.0 - 845.0*k4/4104.0;
		p_temp = p + h;
		
		k5 = h*integrand1(p_temp, w_temp, t_i, n_i, mu_i, beta_i, gamma_i);
		w_temp = w - 8.0*k1/27.0 + 2.0*k2 - 3544.0*k3/2565.0 + 1859.0*k4/4104.0 - 11.0*k5/40.0;
		p_temp = p + 0.5*h;
		
		k6 = h*integrand1(p_temp, w_temp, t_i, n_i, mu_i, beta_i, gamma_i);
						
		//RK4 approx
		w1 = w + 25.0*k1/216.0 + 1408.0*k3/2565.0 + 2197.0*k4/4104.0 - 0.2*k5;
		//RK5 approx
		w2 = w + 16.0*k1/135.0 + 6656.0*k3/12825.0 + 28561.0*k4/56430.0\
				-9.0*k5/50.0 + 2.0*k6/55.0;
		
		//Compute error
		R = sqrt(w2*w2 - w1*w1)/h;
		
		//Calculate update factor
		delta = 0.84*pow((eps/R), 0.25);
		
		//Check if to progress to next step or recalculate current step with
		// new step size. 		
		if (R <= eps)
		{
			w = w1;
			p = p + h;
			h = delta*h;
			i++;
		}
		else
		{
			h = delta*h;
		}
	}
	
	return w;
}

double integrand1(double p, double w, double* t_i, double* n_i, double mu_i,\
                  double beta_i, double gamma_i)
{
	double eta_i = mu_i*p + beta_i*p*p + gamma_i*p*p*p;
	double alpha = 0.5;
	double func;
	func = 1.0/pow(sqrt(p*p*t_i[0]*t_i[0] + 2*p*t_i[0]*eta_i*n_i[0] + eta_i*eta_i*n_i[0]*n_i[0] \
						+ p*p*t_i[1]*t_i[1] +2*p*t_i[1]*eta_i*n_i[1] + eta_i*eta_i*n_i[1]*n_i[1]), alpha);
						
	return func;
}

double evaluate_integral2_RK(double eps, double h, double int_IC,\
            double* t_i, double* n_i, double mu_i, double beta_i, double gamma_i) 
{
  
	double p = 0.1;
	double p_end = 1;
  double k1, k2, k3, k4, k5, k6;
	double w = int_IC, w1, w2, R, delta, w_temp, p_temp;
	int i = 0;
	
	while (p < p_end) 
	{
		if ((p_end - p) < h)
		{
			h = p_end - p;
		}
		k1 = h*integrand2(p, w, t_i, n_i, mu_i, beta_i, gamma_i); 
		w_temp = w + 0.25*k1;
		p_temp = p + 0.25*h;
		
		k2 = h*integrand2(p_temp, w_temp, t_i, n_i, mu_i, beta_i, gamma_i);
		w_temp = w + 3.0*k1/32.0 + 9.0*k2/32.0;
		p_temp = p + 3.0*h/8.0;
		
		k3 = h*integrand2(p_temp, w_temp, t_i, n_i, mu_i, beta_i, gamma_i);
		w_temp = w + 1932.0*k1/2197.0 -7200.0*k2/2197.0 + 7296.0*k3/2197.0;
		p_temp = p + 12.0*h/13.0;
		
		k4 = h*integrand2(p_temp, w_temp, t_i, n_i, mu_i, beta_i, gamma_i);
		w_temp = w + 439.0*k1/216.0 - 8.0*k2 + 3680.0*k3/513.0 - 845.0*k4/4104.0;
		p_temp = p + h;
		
		k5 = h*integrand2(p_temp, w_temp, t_i, n_i, mu_i, beta_i, gamma_i);
		w_temp = w - 8.0*k1/27.0 + 2.0*k2 - 3544.0*k3/2565.0 + 1859.0*k4/4104.0 - 11.0*k5/40.0;
		p_temp = p + 0.5*h;
		
		k6 = h*integrand2(p_temp, w_temp, t_i, n_i, mu_i, beta_i, gamma_i);
    
		//RK4 approx
		w1 = w + 25.0*k1/216.0 + 1408.0*k3/2565.0 + 2197.0*k4/4104.0 - 0.2*k5;
		//RK5 approx
		w2 = w + 16.0*k1/135.0 + 6656.0*k3/12825.0 + 28561.0*k4/56430.0\
				-9.0*k5/50.0 + 2.0*k6/55.0;
    printf("k1 = %lf\n", k1);
    printf("k3 = %lf\n", k3);
    printf("k4 = %lf\n", k4);
    printf("k5 = %lf\n", k5);
    printf("k6 = %lf\n", k6);
    printf("w1 = %lf\n", w1);
    printf("w2 = %lf\n", w2);
		
		//Compute error
		R = sqrt(w2*w2 - w1*w1)/h;
    printf("w2*w2 - w1*w1 = %lf\n", w2*w2 - w1*w1);
    printf("h = %lf\n", h);
    printf("R = %lf\n", R);
		
		//Calculate update factor
		delta = 0.84*pow((eps/R), 0.25);
		
		//Check if to progress to next step or recalculate current step with
		// new step size. 		
		if (R <= eps)
		{
			w = w1;
			p = p + h;
			h = delta*h;
			i++;
		}
		else
		{
			h = delta*h;
		}
	}
	
	return w;
}

double integrand2(double p, double w, double* t_i, double* n_i,\
                  double mu_i, double beta_i, double gamma_i)
{
	double func;
  //double func2;
	double alpha = 0.5;
  //double eta_i = (mu_i + (beta_i + gamma_i*p)*p)*p;
  
  // Tried using Horner's Scheme for better accuracy
  double a2[2], a3[2], a4[2], a5[2], a6[2];
  a2[0] = t_i[0]*t_i[0] + 2*t_i[0]*mu_i*n_i[0] + mu_i*mu_i*n_i[0]*n_i[0];
  a2[1] = t_i[1]*t_i[1] + 2*t_i[1]*mu_i*n_i[1] + mu_i*mu_i*n_i[1]*n_i[1];
  a3[0] = 2*t_i[0]*beta_i*n_i[0] + 2*mu_i*beta_i*n_i[0]*n_i[0];
  a3[1] = 2*t_i[1]*beta_i*n_i[1] + 2*mu_i*beta_i*n_i[1]*n_i[1];
  a4[0] = 2*t_i[0]*gamma_i*n_i[0] + beta_i*beta_i*n_i[0]*n_i[0] + 2*mu_i*gamma_i*n_i[0]*n_i[0];
  a4[1] = 2*t_i[1]*gamma_i*n_i[1] + beta_i*beta_i*n_i[1]*n_i[1] + 2*mu_i*gamma_i*n_i[1]*n_i[1];
  a5[0] = 2*beta_i*gamma_i*n_i[0]*n_i[0];
  a5[1] = 2*beta_i*gamma_i*n_i[1]*n_i[1];
  a6[0] = gamma_i*gamma_i*n_i[0];
  a6[1] = gamma_i*gamma_i*n_i[1];
  double x_sq = a2[0] + (a3[0] + (a4[0] + (a5[0] + a6[0]*p)*p)*p)*p;
  double y_sq = a2[1] + (a3[1] + (a4[1] + (a5[1] + a6[1]*p)*p)*p)*p;
  
  func = p*(2.0*beta_i + 3.0*gamma_i*p)/pow(p*sqrt(x_sq + y_sq), alpha);
  /*
  func2 = p*(2.0*beta_i + 3.0*gamma_i*p)/pow(sqrt(p*p*t_i[0]*t_i[0] + 2*p*t_i[0]*eta_i*n_i[0] + eta_i*eta_i*n_i[0]*n_i[0] \
						+ p*p*t_i[1]*t_i[1] +2*p*t_i[1]*eta_i*n_i[1] + eta_i*eta_i*n_i[1]*n_i[1]), alpha);
  */
  sleep(1);
  printf("func = %lf\n", func);
  
	return func;
}

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
