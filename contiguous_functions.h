void interpolate(double* x, int start, int N, int n_dim, double* t, double* n,\
                 double* d, double* kappa, double* kappa_den, double* mu,\
                 double* beta, double* gamma); 


void autder(double* f, double* c_coeff, double alpha, int order);

// Return the derivative of x
void compute_derivative(double* dxdtx, double* dxdty, double* x, double* mu, double* beta, double* gamma, double* t, double* n, int M, int N, double alpha, double h, double eps, int j);

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

void points_reloc(double** px, double* t, double* n, int* pN, double* kappa,\
									double* mu, double* beta, double* gamma);

double runge_kutta45(double* x, double* dxdt_k1, double* dxdt_k2, double* dxdt_k3, double* dxdt_k4, double* dxdt_k5, double* dxdt_k6, double* dxdt_RK4, double* dxdt_RK5, double tol, double dt, int M, int N, double* mu, double* beta, double* gamma, double* t, double* n, double alpha, double eps, double h);

double compute_area(double* x, int start, int stop, double* t, double* n,\
                    double* mu, double* beta, double* gamma);
