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
  
  double dxdt = 0;
  double d_xi, d_ni, d_ti;
  
  // Evolve the contour integrals
  for (int i = 0; i < N*P; i++) { // i < ???
    if ((x[j] == x[i]) && (y[j] == y[i])) {
      dxdt += evaluate_integral(p, mu[i], beta[i], gamma[i], t_x[i], t_y[i], n_x[i], n_y[i], alpha);
      
    } else if ((x[j] == x[i+1]) && (y[j] == y[i+1])) {
    
      dxdt += evaluate_integral(1-p, mu[i] + 2*beta[i] + 3*gamma[i], -beta[i] \
                                - 3*gamma[i], gamma[i], t_x[i], t_y[i], n_x[i], n_y[i], alpha);
    
    } else if (sqrt((x[j] - x[i])*(x[j] - x[i]) + (y[j] - y[i])*(y[j] - y[i])) > 0.5) {
    
      d_xi = sqrt((x[j] - x[i])*(x[j] - x[i]) + (y[j] - y[i])*(y[j] - y[i]));
      d_ni = (x[j] - x[i])*t_y[j] - (y[j] - y[i])*t_x[j];
      d_ti = -(x[j] - x[i])*t_x[j] + (y[j] - y[i])*t_y[j];
      
      dxdt += evaluate_integral_g(mu[i], d_xi, d_ni, d_ti, alpha);
    
    } else if (sqrt((x[j] - x[i])*(x[j] - x[i]) + (y[j] - y[i])*(y[j] - y[i])) < 0.01) {
    
      d_xi = sqrt((x[j] - x[i])*(x[j] - x[i]) + (y[j] - y[i])*(y[j] - y[i]));
      d_ni = (x[j] - x[i])*t_y[j] - (y[j] - y[i])*t_x[j];
      d_ti = -(x[j] - x[i])*t_x[j] + (y[j] - y[i])*t_y[j];
      
      dxdt += evaluate_integral_g(mu[i], d_xi, d_ni, d_ti, alpha);
    }
  }
  
  dxdt *= theta/(double)TWOPI;
  
  
  
  return dxdt;
}

double evaluate_integral(double p, double mu_i, double beta_i, double gamma_i, \
  double t_xi, double t_yi, double* n_x, double* n_y, double alpha) {
  
  // Coefficients
  double* c = (double*)malloc(10*sizeof(double));
  
  // Compute the first integral
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
  
  first = first/(t_abs*alpha_mu);
  second = second/(t_abs*alpha_mu);
  
  free(c);
  
  return first + second;
}

double evaluate_integral_g(double mu_i, double d_xi, double d_ni, double d_ti, double alpha) {
  
  // Coefficients
  double* g = (double*)malloc(10*sizeof(double));
  
  
  free(g);
  
  return integral;
}
