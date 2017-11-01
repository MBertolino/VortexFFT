void interpolate(int N, int P, double** x);

// Return the derivative of x
double compute_derivative(double** x, double* p, double* mu, \
  double* beta, double* gamma, double* t_x, double* t_y, double* n_x, double* n_y, int NP, double alpha);

// Evaluate integral in case 1 and 2
double evaluate_integral(double p, double mu_i, double beta_i, double gamma_i, \
  double* t_i, double* n_i, double alpha);

// Evaluate integral in case 3
double evaluate_integral_g(double mu_i, double d_xi, double d_ni, double d_ti, double alpha);