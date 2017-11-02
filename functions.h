void interpolate(double** x, int N, int P, int n_dim, double** t, double** n, double* p,\
                 double* eta, double* d, double* kappa, double* kappa_den, double* mu,\
                 double* beta, double* gamma); 
// Return the derivative of x
double compute_derivative(double** x, double* mu, double* beta, double* gamma, double** t, double** n, int NP, double alpha, int j);

// Evaluate integral in case 1 and 2
double evaluate_integral(double mu_i, double beta_i, double gamma_i,\
  double* t_i, double* n_i, double alpha);

// Evaluate integral in case 3
double evaluate_integral_g(double mu_i, double d_x, double d_ni, double d_ti, double alpha);

// Evaluate integral in case 4
double evaluate_integral_RK(double mu_i, double eps, double h, double int_IC, double* t_i,\
									double* n_i, double eta_i);

double evaluate_integral1_RK(double eps, double h, double int_IC, double* t_i,\
									double* n_i, double eta_i);

double evaluate_integral2_RK(double eps, double h, double int_IC, double* t_i,\
									double* n_i, double eta_i, double beta_i, double gamma_i);

