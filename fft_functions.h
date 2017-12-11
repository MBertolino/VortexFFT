

// Return the derivative of x with FFT
void derivative_fft(double* x, double* dx, int start, int stop);

// Return the vector field with FFT
void vfield_fft(double* dxdt, double* x, int M, int N, double alpha, double theta);

// Compute the area with FFT
double area_fft(double* x, int start, int stop);