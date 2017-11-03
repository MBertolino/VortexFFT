void interpolate(double** x, int N, int P, int n_dim, double** t, double** n, double* p,\
                 double* eta, double* d, double* kappa, double* kappa_den, double* mu,\
                 double* beta, double* gamma); 
// Return the derivative of x
void compute_derivative(double* dxdt, double** x, double* mu, double* beta, double* gamma, double** t, double** n, int NP, double alpha, double h, double eps, int j);

// Evaluate integral in case 1 and 2
void evaluate_integral(double* dxdt, double mu_i, double beta_i, double gamma_i,\
  double* t_i, double* n_i, double alpha);

// Evaluate integral in case 3
void evaluate_integral_g(double* dxdt, double mu_i, double beta_i, double gamma_i, double d_x, double d_ni, double d_ti, double* t_i, double* n_i, double alpha);

// Evaluate integral in case 4
void evaluate_integral_RK(double* dxdt, double mu_i, double beta_i, double gamma_i,\
                          double eps, double h, double int_IC, double* t_i, double* n_i);

double evaluate_integral1_RK(double eps, double h, double int_IC, double* t_i,\
									double* n_i, double mu_i, double beta_i, double gamma_i);

double evaluate_integral2_RK(double eps, double h, double int_IC, double* t_i,\
									double* n_i, double beta_i, double gamma_i);

double integrand1(double p, double w, double* t_i, double* n_i, double mu_i,\
                  double beta_i, double gamma_i);

double integrand2(double p, double w, double* t_i, double* n_i,\
                  double eta_i, double beta_i, double gamma_i);

