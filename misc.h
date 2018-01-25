double get_wall_seconds();

void print_to_file(double* x, int M, int N, int k);

void allocate(double** d, double** kappa, double** mu, double** beta, double** gamma, double** t, double** n, double** norm, double** k1, double** k2, double** k3, double** k4, double** k5, double** k6, int N);

void free_step(double* d, double* kappa, double* mu, double* beta, double* gamma, double* t, double* n, double* norm, double* k1, double* k2, double* k3, double* k4, double* k5, double* k6);

double distance(double a_x, double a_y, double b_x, double b_y);

double scalar_prod(double a_x, double a_y, double b_x, double b_y);

double runge_kutta45(double* x, double* k1, double* k2, double* k3, double* k4, double* k5, double* k6, double tol_rk45_time, double dt, int M, int N, double* mu, double* beta, double* gamma, double* t, double* n, double alpha, double tol_rk45_space, double h, double theta, double* norm);

void compare_algo(double* x, double* mu, double* beta, double* gamma, double* t, double* n, int M, int N, double h, double tol_rk45_space, double alpha, double theta);

