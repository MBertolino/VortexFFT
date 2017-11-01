#include "functions.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TWOPI 6.2831853071795864769
#define SQRTTWO 1.4142135623730950588

void interpolate(int N, int P, double** x) {
	
  double* p = (double*)malloc(P*sizeof(double));
  double* eta = (double*)malloc(P*sizeof(double));
  for (int i = 0; i < P; i++)
    p[i] = (double)i/P;
  double** t = (double**)malloc(n_dim*sizeof(double*));
  double** n = (double**)malloc(n_dim*sizeof(double*));
  for (int i = 0; i < N; i++) {
    t[i] = (double*)malloc(N*sizeof(double));
    n[i] = (double*)malloc(N*sizeof(double));
  }
  double* d = (double*)malloc(N*sizeof(double));
  double* kappa = (double*)malloc(N*sizeof(double));
  double kappa_den[2];
  double mu, beta, gamma;
  
  // Generate circle (j*P to avoid the interpolated nodes)
  for (int j = 0; j < N; j++) {
    x[0][j*P] = cos(TWOPI*j/(double)N);
    x[1][j*P] = sin(TWOPI*j/(double)N);
  }
  
  // Calculate t an n
  for (int j = 0; j < N-1; j++) {
    t[0][j] = x[0][(j+1)*P] - x[0][j*P];
    t[1][j] = x[1][(j+1)*P] - x[1][j*P];
    n[0][j] = -t[1][j];
    n[1][j] = t[0][j];
    d[j] = sqrt(t[0][j]*t[0][j] + t[1][j]*t[1][j]);
  }  
  // Special case j = N-1
  t[0][N-1] = x[0][0] - x[0][(N-1)*P];
  t[1][N-1] = x[1][0] - x[1][(N-1)*P];
  n[0][N-1] = -t[1][N-1];
  n[1][N-1] = t[0][N-1];
  d[N-1] = sqrt(t[0][N-1]*t[0][N-1] + t[1][N-1]*t[1][N-1]);
  
  // kappa local curvature
  kappa_den[0] = d[N-1]*d[N-1]*t[0][0] + d[0]*d[0]*t[0][N-1];
  kappa_den[1] = d[N-1]*d[N-1]*t[1][0] + d[0]*d[0]*t[1][N-1];
  kappa[0] = 2*(t[0][N-1]*t[1][0] - t[1][N-1]*t[0][0])\
    /sqrt(kappa_den[0]*kappa_den[0] + kappa_den[1]*kappa_den[1]);
  
  for (int j = 1; j < N; j++) {
    // kappa local curvature
    kappa_den[0] = (d[j-1]*d[j-1]*t[0][j] + d[j]*d[j]*t[0][j-1]);
    kappa_den[1] = (d[j-1]*d[j-1]*t[1][j] + d[j]*d[j]*t[1][j-1]);
    kappa[j] = 2*(t[0][j-1]*t[1][j] - t[1][j-1]*t[0][j])\
      /sqrt(kappa_den[0]*kappa_den[0] + kappa_den[1]*kappa_den[1]);
  }
  
  // Construct the cubic interpolation
  for (int j = 0; j < N; j++) {
    
    // Cubic interpolation coefficients
    mu = -(double)1/3*d[j]*kappa[j] - (double)1/6*d[j]*kappa[j+1];
    beta = 0.5*d[j]*kappa[j];
    gamma = (double)1/6*d[j]*(kappa[j+1] - kappa[j]);
    
    for (int i = 0; i < P; i++) {
      eta[i] = mu*p[i] + beta*p[i]*p[i] + gamma*p[i]*p[i]*p[i];
      x[0][j*P + i] = x[0][j*P] + p[i]*t[0][j] + eta[i]*n[0][j];
      x[1][j*P + i] = y[1][j*P] + p[i]*t[1][j] + eta[i]*n[1][j];
    }
  }
    
  // Free memory
  for (int i = 0; i < n_dim; i++) {
    free(t[i]);
    free(n[i]);
  }
  free(t);
  free(n);
  free(d);
  free(p);
  free(eta);
  free(kappa);
  
  return;
}

double compute_derivative(double** x, double* p, double* mu, \
  double* beta, double* gamma, double* t, double* n, int NP, double alpha) {
  
  double* derivative;
  double d_x, d_ni, d_ti;
  
  // Evolve the contour integrals
  for (int i = 0; i < N*P; i++) { // i < ???
    if ((x[0][j] == x[0][i]) && (y[j] == y[i])) {
      // Case 1: Use formula (29)
      derivative += evaluate_integral(p, mu[i], beta[i], gamma[i], t[0][i], t[1][i], n[0][i], n[1][i], alpha);
      
    } else if ((x[0][j] == x[0][i+1]) && (x[1][j] == x[1][i+1])) {
      // Case 2: Use formula (29) with shifted params
      derivative += evaluate_integral(1-p, mu[i] + 2*beta[i] + 3*gamma[i], -beta[i] \
                                - 3*gamma[i], gamma[i], t[0][i], t[1][i], n[0][i], n[1][i], alpha);
    
    } else if (sqrt((x[0][j] - x[0][i])*(x[0][j] - x[0][i]) + (x[1][j] - x[1][i])*(x[1][j] - x[1][i])) > 0.5) {
      // Case 3: Use formula (31)
      d_x = sqrt((x[0][j] - x[0][i])*(x[0][j] - x[0][i]) + (x[1][j] - x[1][i])*(x[1][j] - x[1][i]));
      d_ni = (x[0][j] - x[0][i])*t[1][j] - (x[1][j] - x[1][i])*t[0][j];
      d_ti = -(x[0][j] - x[0][i])*t[0][j] + (x[1][j] - x[1][i])*t[1][j];
      
      derivative += evaluate_integral_g(mu[i], d_x, d_ni, d_ti, alpha);
    
    } else if (sqrt((x[0][j] - x[0][i])*(x[0][j] - x[0][i]) + (x[1][j] - x[1][i])*(x[1][j] - x[1][i])) < 0.01) {
      // Case 4: Use Runge-Kutta 4-5

      derivative += evaluate_integral_RK(mu[i], d_x, d_ni, d_ti, alpha);
    }
  }
  
  derivative *= theta/(double)TWOPI;
  
  return derivative;
}

double evaluate_integral(double p, double mu_i, double beta_i, double gamma_i, \
  double* t_i, double* n, double alpha) {
  
  // Coefficients
  double* c = (double*)malloc(10*sizeof(double));
  
  // Compute the integrals
  double first = 0;
  double second = 0;
  double eval[2];
  double p_coef, psq_coef;
  double t_abs = pow(sqrt((t_i[0]*t_i[0] + t_i[1]*t_i[1])), alpha);
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
  eval[0] = first*(t_i[0] + mu_i*n_i[0]) + second*n_i[0];
  eval[1] = first*(t_i[1] + mu_i*n_i[1]) + second*n_i[1];
  
  free(c);
  
  return eval;
}

double evaluate_integral_g(double mu_i, double d_x, double d_ni, double d_ti, double alpha) {
  
  // Coefficients
  double* g = (double*)malloc(10*sizeof(double));
  
  // COmpute the integrals
  double first = 0;
  double second = 0;
  double eval[2];
  double p_coef, psq_coef;

  for (int n = 0; n < 11; n++) {
    first += g[n];
    
    p_coef = 2*beta_i/(double)(n - alpha + 2);
    psq_coef = 3*gamma_i/(double)(n - alpha + 3);
    second += g[n]*(p_coef + psq_coef);
  }
  
  first = first/pow(d_x, alpha);
  second = second/...;
  
  // Sum together
  eval[0] = first*(t_i[0] + mu_i*n_i[0]) + second*n_i[0];
  eval[1] = first*(t_i[1] + mu_i*n_i[1]) + second*n_i[1];
  
  free(g);
  
  return eval;
}

<<<<<<< HEAD

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
=======
//double evaluate_integral_RK() {
//  
//  return ;
//}

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
}
>>>>>>> fd5bb46268edf7aa9225e769ce3e1e14205e76ba
