#include <stdio.h>
#include <math.h>
#include <fftw3.h>

#define SWAP(a, b) temp = (a); (a) = (b); (b) = temp

// 256 bits
#define R2(n)     n,     n + 2*64,     n + 1*64,     n + 3*64
#define R4(n) R2(n), R2(n + 2*16), R2(n + 1*16), R2(n + 3*16)
#define R6(n) R4(n), R4(n + 2*4 ), R4(n + 1*4 ), R4(n + 3*4 )


/* #define REV(n) n,
0   1   2   3   4   5   6   7
000 001 010 011 100 101 110 111
000 100 010 110 001 101 011 111
0   4   2   6   1   5   3   7
// */ 


static const unsigned char BitReverseTable256[256] = {R6(0), R6(2), R6(1), R6(3)};

unsigned char ReverseBitsLookupTable(unsigned char v) {
  return BitReverseTable256[v];
}

unsigned char reverse_bit(unsigned char b, int N) {

  // Reversal 8-bit
  unsigned char size = N;
  unsigned char b_rev = 0;
  for (int position = size - 1; position > 0; position --) {
    b_rev += ((b&1) << position);
    b >>= 1;
  }

  return b_rev;  
}

void swap(float *a, float *b) {
  
  unsigned temp = *a;
  *a = *b;
  *b = temp;
}

int main() {
  
  int N = 8;
  float x[N], x_F[N];
  unsigned char i_rev; // Reverse bit
  
  // Real and imaginary twiddle factors
  float twiddle_R[N], twiddle_I[N];
  for (int k = 0; k < N; k++) {
    twiddle_R[k] = cos(2*M_PI*k/N);
    twiddle_I[k] = sin(2*M_PI*k/N);
  }
  
  // Input x
  for (int i = 0; i < N; i++)
    x[i] = 0;
  for (int i = 0; i < N; i+=2)
    x[i] = 1;
  for (int i = 0; i < N; i++)
    printf("%f ", x[i]);
  printf("\n");
  
  // Reverse input
  for (int i = 0; i < 0.5*N; i++) { 
    i_rev = reverse_bit(i, 4); // reverse this (8-bit) byte
    if (i < i_rev) {
      swap(&x[i], &x[i_rev]);
    }
  }
  
  // Print reversed input x
  printf("\n");
  for (int i = 0; i < N; i++)
    printf("%f ", x[i]);
  printf("\n");
  
  // Two-point DFT Butterfly
  /* int M = N >> 1; // N >> 1 = N/2
  int m = M
  for (int j = 0; j < N; j+=2) {
    x_F[j] = x[j] + twiddle_R[0]*x[j+1];
    x_F[j+1] = x[j] - twiddle_R[0]*x[j+1];
  }
  ///
  int M = 2;
  int m = M;
  while (N > M) {
    m << 1;
    theta = 2*M_PI/M;
    
    for (m = 1; m < M; m += 2) {
      
    }
  }
  
  /*
    // External loop
  while (n > M) {
      
      // Initi
    i_step = M << 1;
    theta = 2*M_PI/M;
    w_temp = sin(0.5*theta);
    wpr = -2.0*w_temp*w_temp;
    wpi = sin(theta);
    wr = 1.0;
    wi = 0.0;
    
    // Internal loops
    for (m = 1; m < M; m += 2) {
      for (i = m; i <= n; i += i_step) {
        j = i + M;
        tempr = wr*data[j-1] - wi*data[j];
        tempi = wr*data[j] + wi*data[j-1];
        data[j-1] = data[i-1] - tempr;
        data[j] = data[i] - tempi;
        data[i-1] += tempr;
        data[i] += tempi;
      }
      wr = (w_temp = wr)*wpr - wi*wpi + wr;
      wi = wi*wpr + wtemp*wpi + wi;
    }
    M = i_step;
  }
  //*/
  
  return(0);
}