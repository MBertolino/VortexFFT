void interpolate(double** x, int N, int n_dim, double** t, double** n,\
                 double* d, double* kappa, double* kappa_den, double* mu,\
                 double* beta, double* gamma); 

void autder(double* f, double* c_coeff, double alpha, int order);

// Return the derivative of x
void compute_derivative(double* dxdt, double** x, double* mu, double* beta, double* gamma, double** t, double** n, int N, double alpha, double h, double eps, int j);

// Evaluate integral in case 1 and 2
void evaluate_integral(double* dxdt, double mu_i, double beta_i, double gamma_i,\
  double* t_i, double* n_i, double* c, double alpha);

// Evaluate integral in case 3
void evaluate_integral_g(double* dxdt, double mu_i, double beta_i, double gamma_i, double d_x, double d_ni, double d_ti, double* t_i, double* n_i, double* g, double alpha);

// Evaluate integral in case 4
void evaluate_integral_RK(double* dxdt, double* x_i, double* x_j, double mu_i, double beta_i, double gamma_i,\
                          double eps, double h, double* t_i, double* n_i, double alpha);

double evaluate_integral1_RK(double* x_i, double* x_j, double eps, double h, double* t_i,\
									double* n_i, double mu_i, double beta_i, double gamma_i, double alpha);

double evaluate_integral2_RK(double* x_i, double* x_j, double eps, double h, double* t_i,\
									double* n_i, double mu_i, double beta_i, double gamma_i, double alpha);

double integrand1(double* x_i, double* x_j, double p, double* t_i, double* n_i, double mu_i,\
                  double beta_i, double gamma_i, double alpha);

double integrand2(double* x_i, double* x_j, double p, double* t_i, double* n_i,\
                  double mu_i, double beta_i, double gamma_i, double alpha);

