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