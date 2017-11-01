#include "functions.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TWOPI 6.2831853071795864769

void interpolate(int N, int P, double* x, double* y) {
	
  double* p = (double*)malloc(P*sizeof(double));
  double* eta = (double*)malloc(P*sizeof(double));
  for (int i = 0; i < P; i++)
    p[i] = (double)i/P;
  double* t_x = (double*)malloc(N*sizeof(double));
  double* t_y = (double*)malloc(N*sizeof(double));
  double* n_x = (double*)malloc(N*sizeof(double));
  double* n_y = (double*)malloc(N*sizeof(double));
  double* d = (double*)malloc(N*sizeof(double));
  double* kappa = (double*)malloc(N*sizeof(double));
  double kappa_den_x, kappa_den_y;
  double mu, beta, gamma;
  
  // Generate circle (j*P to avoid the interpolated nodes)
  for (int j = 0; j < N; j++) {
    x[j*P] = cos(TWOPI*j/(double)N);
    y[j*P] = sin(TWOPI*j/(double)N);
  }
  
  // Calculate t an n
  for (int j = 0; j < N-1; j++) {
    t_x[j] = x[(j+1)*P] - x[j*P];
    t_y[j] = y[(j+1)*P] - y[j*P];
    n_x[j] = -t_y[j];
    n_y[j] = t_x[j];
    d[j] = sqrt(t_x[j]*t_x[j] + t_y[j]*t_y[j]);
  }  
  // Special case j = N-1
  t_x[N-1] = x[0] - x[(N-1)*P];
  t_y[N-1] = y[0] - y[(N-1)*P];
  n_x[N-1] = -t_y[N-1];
  n_y[N-1] = t_x[N-1];
  d[N-1] = sqrt(t_x[N-1]*t_x[N-1] + t_y[N-1]*t_y[N-1]);
  
  // kappa local curvature
  kappa_den_x = (d[N-1]*d[N-1]*t_x[0] + d[0]*d[0]*t_x[N-1]);
  kappa_den_y = (d[N-1]*d[N-1]*t_y[0] + d[0]*d[0]*t_y[N-1]);
  kappa[0] = 2*(t_x[N-1]*t_y[0] - t_y[N-1]*t_x[0])\
    /sqrt(kappa_den_x*kappa_den_x + kappa_den_y*kappa_den_y);
  
  for (int j = 1; j < N; j++) {
    // kappa local curvature
    kappa_den_x = (d[j-1]*d[j-1]*t_x[j] + d[j]*d[j]*t_x[j-1]);
    kappa_den_y = (d[j-1]*d[j-1]*t_y[j] + d[j]*d[j]*t_y[j-1]);
    kappa[j] = 2*(t_x[j-1]*t_y[j] - t_y[j-1]*t_x[j])\
      /sqrt(kappa_den_x*kappa_den_x + kappa_den_y*kappa_den_y);
  }
  
  // Construct the cubic interpolation
  for (int j = 0; j < N; j++) {
    
    // Cubic interpolation coefficients
    mu = -(double)1/3*d[j]*kappa[j] - (double)1/6*d[j]*kappa[j+1];
    beta = 0.5*d[j]*kappa[j];
    gamma = (double)1/6*d[j]*(kappa[j+1] - kappa[j]);
    
    for (int i = 0; i < P; i++) {
      eta[i] = mu*p[i] + beta*p[i]*p[i] + gamma*p[i]*p[i]*p[i];
      x[j*P + i] = x[j*P] + p[i]*t_x[j] + eta[i]*n_x[j];
      y[j*P + i] = y[j*P] + p[i]*t_y[j] + eta[i]*n_y[j];
    }
  }
    
  // Free memory
  free(t_x);
  free(t_y);
  free(n_x);
  free(n_y);
  free(d);
  free(p);
  free(eta);
  free(kappa);
  
  return;
}

double compute_derivative(double* x, double* y, double* p, double* mu, \
  double* beta, double* gamma, double* t_x, double* t_y, double* n_x, double* n_y, double alpha) {
  
  double derivative = 0;
  double d_x, d_ni, d_ti;
  
  // Evolve the contour integrals
  for (int i = 0; i < N*P; i++) { // i < ???
    if ((x[j] == x[i]) && (y[j] == y[i])) {
      derivative += evaluate_integral(p, mu[i], beta[i], gamma[i], t_x[i], t_y[i], n_x[i], n_y[i], alpha);
      
    } else if ((x[j] == x[i+1]) && (y[j] == y[i+1])) {
    
      derivative += evaluate_integral(1-p, mu[i] + 2*beta[i] + 3*gamma[i], -beta[i] \
                                - 3*gamma[i], gamma[i], t_x[i], t_y[i], n_x[i], n_y[i], alpha);
    
    } else if (sqrt((x[j] - x[i])*(x[j] - x[i]) + (y[j] - y[i])*(y[j] - y[i])) > 0.5) {
    
      d_x = sqrt((x[j] - x[i])*(x[j] - x[i]) + (y[j] - y[i])*(y[j] - y[i]));
      d_ni = (x[j] - x[i])*t_y[j] - (y[j] - y[i])*t_x[j];
      d_ti = -(x[j] - x[i])*t_x[j] + (y[j] - y[i])*t_y[j];
      
      derivative += evaluate_integral_g(mu[i], d_x, d_ni, d_ti, alpha);
    
    } else if (sqrt((x[j] - x[i])*(x[j] - x[i]) + (y[j] - y[i])*(y[j] - y[i])) < 0.01) {
      
      // Compute integral with RK5.
      // derivative += evaluate_integral_RK();
      
      d_x = sqrt((x[j] - x[i])*(x[j] - x[i]) + (y[j] - y[i])*(y[j] - y[i]));
      d_ni = (x[j] - x[i])*t_y[j] - (y[j] - y[i])*t_x[j];
      d_ti = -(x[j] - x[i])*t_x[j] + (y[j] - y[i])*t_y[j];
      
      derivative += evaluate_integral_g(mu[i], d_x, d_ni, d_ti, alpha);
    }
  }
  
  derivative *= theta/(double)TWOPI;
  
  return derivative;
}

double evaluate_integral(double p, double mu_i, double beta_i, double gamma_i, \
  double t_xi, double t_yi, double* n_x, double* n_y, double alpha) {
  
  // Coefficients
  double* c = (double*)malloc(10*sizeof(double));
  
  // Compute the integrals
  double first = 0;
  double second = 0;
  double p_coef, psq_coef;
  double t_abs = pow(sqrt((t_xi*t_xi + t_yi*t_yi)), alpha);
  double alpha_mu = pow((1 + mu_i*mu_i), 0.5*alpha);
  
  for (int n = 0; n < 11; n++) {
    first += c[n]/(double)(n - alpha + 1);
    
    p_coef = 2*beta_i/(double)(n - alpha + 2);
    psq_coef = 3*gamma_i/(double)(n - alpha + 3);
    second += c[n]*(p_coef + psq_coef);
  }
  
  first = first/(t_abs*alpha_mu); // Multiply by (t_i + mu*n_i)
  second = second/(t_abs*alpha_mu); // Multiply by n_i
  
  free(c);
  
  return first + second;
}

double evaluate_integral_g(double mu_i, double d_x, double d_ni, double d_ti, double alpha) {
  
  // Coefficients
  double* g = (double*)malloc(10*sizeof(double));
  
  // COmpute the integrals
  double first = 0;
  double second = 0;
  for (int n = 0; n < 11; n++) {
    first += g[n];
    second += ...;
  }
  
  first = first/pow(d_x, alpha) // Multiply by (t_i + mu*n_i)
  second = second/... // Multiply by n_i
  
  
  free(g);
  
  return first + second;
}


double evaluate_integral1_RK(double eps, double h, double int_IC, double t_xi, \ 
									double t_yi, double n_xi, double n_yi, double eta_i) 
{
	double p0 = 0;
	double p_end = 1;
	double p;
	double w = int_IC, w1, w2, R, delta, w_temp, p_temp;
	int i = 0;
	
	while (p < p_end) 
	{
		if ((p_end - p) < h)
		{
			h = p_end - p;
		}
		
		k1 = h*integrand1(p, w, t_xi, t_yi, n_xi, n_yi, eta_i);     //Func should be integrand1-function
		w_temp = w + 0.25*k1;
		p_temp = p + 0.25*h;
		
		k2 = h*func(p_temp, w_temp, t_xi, t_yi, n_xi, n_yi, eta_i);
		w_temp = w + 3.0*k1/32.0 + 9.0*k2/32.0;
		p_temp = p + 3.0*h/8.0;
		
		k3 = h*func(p_temp, w_temp, t_xi, t_yi, n_xi, n_yi, eta_i);
		w_temp = w + 1932.0*k1/2197.0 -7200.0*k2/2197.0 + 7296.0*k3/2197.0;
		p_temp = p + 12.0*h/13.0;
		
		k4 = h*func(p_temp, w_temp, t_xi, t_yi, n_xi, n_yi, eta_i);
		w_temp = w + 439.0*k1/216.0 - 8.0*k2 + 3680.0*k3/513.0 - 845.0*k4/4104.0;
		t_temp = p + h;
		
		
		k5 = h*func(p_temp, w_temp, t_xi, t_yi, n_xi, n_yi, eta_i);
		w_temp = w - 8.0*k1/27.0 + 2.0*k2 - 3544.0*k3/2565.0 + 1859.0*k4/4104.0 - 11.0*k5/40.0;
		p_temp = p + 0.5*h;
		
		k6 = h*func(p_temp, w_temp, t_xi, t_yi, n_xi, n_yi, eta_i);
						
		//RK4 approx
		w1 = w + 25.0*k1/216.0 + 1408.0*k3/2565.0 + 2197.0*k4/4104.0 - 0.2*k5;
		//RK5 approx
		w2 = w + 16.0*k1/135.0 + 6656.0*k3/12825.0 + 28561.0*k4/56430.0 \
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

double integrand1(double p, double w, double t_xi, double t_yi, \
						double n_xi, double n_yi, double eta_i)
{
	double func;
	func = 1.0/sqrt((p*t_xi*t_xi + p*t_yi*t_yi) \ 
						+ (eta_i*n_xi*n_xi + eta_i*n_yi*n_yi));
						
	return func;
}

double evaluate_integral2_RK(double eps, double h, double int_IC, double t_xi, \ 
									double t_yi, double n_xi, double n_yi, double eta_i, \
									double beta_i, double gamma_i) 
{
	double p0 = 0;
	double p_end = 1;
	double p;
	double w = int_IC, w1, w2, R, delta, w_temp, p_temp;
	int i = 0;
	
	while (p < p_end) 
	{
		if ((p_end - p) < h)
		{
			h = p_end - p;
		}
		
		k1 = h*integrand2(p, w, t_xi, t_yi, n_xi, n_yi, eta_i, beta_i, gamma_i); 
		w_temp = w + 0.25*k1;
		p_temp = p + 0.25*h;
		
		k2 = h*func(p_temp, w_temp, t_xi, t_yi, n_xi, n_yi, eta_i, beta_i, gamma_i);
		w_temp = w + 3.0*k1/32.0 + 9.0*k2/32.0;
		p_temp = p + 3.0*h/8.0;
		
		k3 = h*func(p_temp, w_temp, t_xi, t_yi, n_xi, n_yi, eta_i, beta_i, gamma_i);
		w_temp = w + 1932.0*k1/2197.0 -7200.0*k2/2197.0 + 7296.0*k3/2197.0;
		p_temp = p + 12.0*h/13.0;
		
		k4 = h*func(p_temp, w_temp, t_xi, t_yi, n_xi, n_yi, eta_i, beta_i, gamma_i);
		w_temp = w + 439.0*k1/216.0 - 8.0*k2 + 3680.0*k3/513.0 - 845.0*k4/4104.0;
		t_temp = p + h;
		
		
		k5 = h*func(p_temp, w_temp, t_xi, t_yi, n_xi, n_yi, eta_i, beta_i, gamma_i);
		w_temp = w - 8.0*k1/27.0 + 2.0*k2 - 3544.0*k3/2565.0 + 1859.0*k4/4104.0 - 11.0*k5/40.0;
		p_temp = p + 0.5*h;
		
		k6 = h*func(p_temp, w_temp, t_xi, t_yi, n_xi, n_yi, eta_i, beta_i, gamma_i);
						
		//RK4 approx
		w1 = w + 25.0*k1/216.0 + 1408.0*k3/2565.0 + 2197.0*k4/4104.0 - 0.2*k5;
		//RK5 approx
		w2 = w + 16.0*k1/135.0 + 6656.0*k3/12825.0 + 28561.0*k4/56430.0 \
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

double integrand2(double p, double w, double t_xi, double t_yi, double n_xi, \
						double n_yi, double eta_i, double beta_i, double gamma_i)
{
	double func;
	func = (2.0*beta_i*p + 3.0*gamma_i*p*p)/sqrt((p*t_xi*t_xi + p*t_yi*t_yi) \ 
														+ (eta_i*n_xi*n_xi + eta_i*n_yi*n_yi));
	
	return func;
}
