#include "functions.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#define TWOPI 6.2831853071795864769
#define SQRTTWO 1.4142135623730950588
#define PRINT 0


void interpolate(double** x, int N, int n_dim, double** t, double** n,\
                 double* d, double* kappa, double* kappa_den, double* mu,\
                 double* beta, double* gamma)
{

  // Calculate t an n
  for (int j = 0; j < N-1; j++) {
    t[j][0] = x[j+1][0] - x[j][0];
    t[j][1] = x[j+1][1] - x[j][1];
    n[j][0] = -t[j][1];
    n[j][1] = t[j][0];
    d[j] = sqrt(t[j][0]*t[j][0] + t[j][1]*t[j][1]);
  }  
  // Special case j = N-1
  t[N-1][0] = x[0][0] - x[N-1][0];
  t[N-1][1] = x[0][1] - x[N-1][1];
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
  }

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

void compute_derivative(double* dxdt, double** x, double* mu, double* beta, double* gamma, double** t, double** n, int N, double alpha, double h, double eps, int j)
{
  //printf("Entering compute derivative\n");
  double d_x, d_ni, d_ti, d_xi;
  int order = 11;
  double Q = 0.01;
  double f = 1/sqrt(Q);
  
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
      c[n] = 0;
      g[n] = 0;
      poly_coeff_c[n] = 0;
      poly_coeff_g[n] = 0;
    }
    poly_coeff_c[0] = 1;
    poly_coeff_g[0] = 1;
    // Evaluate integrals
    //printf("dx = %lf, d_xi = %lf\n", d_x, d_xi);
    if (i+1 == N && j == 0) {
          // Edge case
      // Case 2: Use formula (29) with shifted params
      #if PRINT
        printf("Edge Case 2\n"); 
      #endif
      //printf("test\n");   
		  // Update parameters
      mu_2 = mu[i] + 2*beta[i] + 3*gamma[i];
      beta_2 = -beta[i] - 3*gamma[i];
      
      // Generate Taylor coefficients
      poly_coeff_c[1] = 2*mu_2*beta_2/(1 + mu_2*mu_2);
      poly_coeff_c[2] = (beta_2*beta_2 + 2*mu_2*gamma[i])/(1 + mu_2*mu_2);
      poly_coeff_c[3] = 2*beta_2*gamma[i]/(1 + mu_2*mu_2);
      poly_coeff_c[4] = gamma[i]*gamma[i]/(1 + mu_2*mu_2);
      
      autder(c, poly_coeff_c, alpha, order); // Should these be divided by two maybe?
      
      evaluate_integral(dxdt, mu_2, beta_2, gamma[i], t[i], n[i], c, alpha); // Look at inputs in these functio
    
      } else {
      
      if (i == j) {
        //Case 1: Use formula (29)
        #if PRINT
          printf("Case 1\n");
        #endif
        
        // Generate Taylor coefficients
        poly_coeff_c[1] = 2*mu[i]*beta[i]/(1 + mu[i]*mu[i]);
        poly_coeff_c[2] = (beta[i]*beta[i] + 2*mu[i]*gamma[i])/(1 + mu[i]*mu[i]);
        poly_coeff_c[3] = 2*beta[i]*gamma[i]/(1 + mu[i]*mu[i]);
        poly_coeff_c[4] = gamma[i]*gamma[i]/(1 + mu[i]*mu[i]);

        autder(c, poly_coeff_c, alpha, order); // Should these be divided by two maybe?
      
        evaluate_integral(dxdt, mu[i], beta[i], gamma[i], t[i], n[i], c, alpha); // Look at inputs in these functions
      
      } else if (i == j - 1) {
        // Case 2: Use formula (29) with shifted params
        #if PRINT
          printf("Case 2\n");      
        #endif
        
		// Update parameters
        mu_2 = mu[i] + 2*beta[i] + 3*gamma[i];
        beta_2 = -beta[i] - 3*gamma[i];
      
        // Generate Taylor coefficients
        poly_coeff_c[1] = 2*mu_2*beta_2/(1 + mu_2*mu_2);
        poly_coeff_c[2] = (beta_2*beta_2 + 2*mu_2*gamma[i])/(1 + mu_2*mu_2);
        poly_coeff_c[3] = 2*beta_2*gamma[i]/(1 + mu_2*mu_2);
        poly_coeff_c[4] = gamma[i]*gamma[i]/(1 + mu_2*mu_2);
      
        autder(c, poly_coeff_c, alpha, order); // Should these be divided by two maybe?
        
        evaluate_integral(dxdt, mu_2, beta_2, gamma[i], t[i], n[i], c, alpha); // Look at inputs in these functions
        //printf("dxdt[%d] = %lf\n", j, dxdt[0]);
      } else if (d_x > f*d_xi) {
        #if PRINT   
          printf("Case 3\n");
        #endif
        // Case 3: Use formula (31)

        // Generate Taylor coefficients
        poly_coeff_g[1] = (d_ti + d_ni*mu[i])/(d_x*d_x);
        poly_coeff_g[2] = ((t[i][0]*t[i][0] + t[i][1]*t[i][1]) + mu[i]*mu[i]*(n[i][0]*n[i][0] + n[i][1]*n[i][1]) + d_x*beta[i])/(d_x*d_x);
        poly_coeff_g[3] = (2*mu[i]*beta[i]*(n[i][0]*n[i][0] + n[i][1]*n[i][1]) + d_ni*gamma[i])/(d_x*d_x);
        poly_coeff_g[4] = ((beta[i] + 2*mu[i]*gamma[i])*(n[i][0]*n[i][0] + n[i][1]*n[i][1]))/(d_x*d_x);
        poly_coeff_g[5] = 2*beta[i]*gamma[i]*(n[i][0]*n[i][0] + n[i][1]*n[i][1])/(d_x*d_x);
        poly_coeff_g[6] = gamma[i]*gamma[i]/6/(d_x*d_x);
    
        autder(g, poly_coeff_g, alpha, order);

        evaluate_integral_g(dxdt, mu[i], beta[i], gamma[i], d_x, d_ni, d_ti, t[i], n[i], g, alpha);

      } else {
        #if PRINT
          printf("Case 4\n");
        #endif
        // Case 4: Use Runge-Kutta 4-5
        evaluate_integral_RK(dxdt, x[i], x[j], mu[i], beta[i], gamma[i], eps, h, t[i], n[i], alpha);
      }
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
  double tol = 0.001;
  
  /*
  // Print to file
  char str[80] = "../w_spike.csv";
  FILE* f = fopen(str, "wb");
    
  while (p < 1) {
	  printf("p = %lf\n", p);
    Y = integrand1(x_i, x_j, p, t_i, n_i, mu_i, beta_i, gamma_i, alpha);
    fprintf(f, "%lf,%lf\n", p, Y);
	  printf("w = %lf\n\n", Y);
    p = p + 0.001;
  }
  fclose(f);
  sleep(5);
  */
  
	do
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
		
    //printf("Y1 = %lf,   Y2 = %lf\n", Y1, Y2);
    //printf("R = %lf\n", R);
		//printf("\n");
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
	} while (p_end - p > tol || p - p_end > tol);
  
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
  double tol = 0.001;

  do
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
  } while (p_end - p > tol || p - p_end > tol);
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
  double kappai_hat = 0.5*(kappai_tilde + kappai_tilde);
  
  double density = kappai_hat*kappa_hat/(1 + epsilon*kappai_hat/SQRTTWO);
  double sigmai = 
  
  return;
}*/
