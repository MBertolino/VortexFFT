  
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "functions.h"

#define TWOPI 6.2831853071795864769

int main()
{
  // Number of points
  int N = 10;
  
  // Number of interpolation points
  int P = 100;
  // Coordinates
  double* x = (double*)malloc(P*N*sizeof(double));
  double* y = (double*)malloc(P*N*sizeof(double));
  for (int i = 0; i < N*P; i++) {
    x[i] = 0;
    y[i] = 0;
  }
  // Generate circle (j*P to avoid the interpolated nodes)
  for (int j = 0; j < N; j++) {
    x[j*P] = cos(TWOPI*j/(double)N);
    y[j*P] = sin(TWOPI*j/(double)N);
  }
  
  interpolate(N, P, x, y);
  
    // Print to file
  char str[80] = "../circle.csv";
  FILE* f = fopen(str, "wb");
  for (int i = 0; i < N*P; i++) {
    fprintf(f, "%lf,%lf\n", x[i], y[i]);
  }
  fclose(f);
  
}