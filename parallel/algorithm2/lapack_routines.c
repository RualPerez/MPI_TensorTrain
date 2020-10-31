#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "lapacke.h"
#include "alloc.h"

/* How to compile these functions in DTU hpc:

   module load mpi/3.1.3-gcc-8.2.0

*/ 

double* linear_system_solver_singular(int N, double* b, double** A) {
  /*
    Solve the equations A*x = b --> https://software.intel.com/en-us/mkl-developer-reference-c-gesv#90C462DB-A8BF-48A1-AE76-5E49D4EA04AF

    Input:
    int N --> int dimension of the squared matrix A and vector b
    double* b --> vector of size N
    double** A --> squared matrix (N x N), can be singular

    Output:
    double* x --> vector of size N

  */
  
  double rcond = 0.00001;
  lapack_int n, nrhs, lda, ldb, info;
  n = N; lda = n; ldb = 1; nrhs = 1;

  int rank;
  lapack_int* jpvt = (lapack_int *)malloc(n*sizeof(lapack_int)) ;
  for(int i=0; i<n; ++i)
    jpvt[0] = 0;

  // Solve the equations A*X = B 
  info = LAPACKE_dgelsy( LAPACK_ROW_MAJOR, n, n, nrhs, A[0], lda, b, ldb, jpvt, rcond, &rank );

  // Check for the exact singularity 
  if( info > 0 ) {
    printf( "The diagonal element of the triangular factor of A,\n" );
    printf( "U(%i,%i) is zero, so that A is singular;\n", info, info );
    printf( "the solution could not be computed.\n" );
    exit( 1 );
  }
  if (info < 0) exit( 1 );

  // Free jpvt
  free(jpvt);

  // Overwritten by the solution x
  return b;
}

double* linear_system_solver(int N, double* b, double** A) {
  /*
    Solve the equations A*x = b --> https://software.intel.com/en-us/mkl-developer-reference-c-gesv#90C462DB-A8BF-48A1-AE76-5E49D4EA04AF

    Input:
    int N --> int dimension of the squared matrix A and vector b
    double* b --> vector of size N
    double** A --> squared matrix (N x N)

    Output:
    double* x --> vector of size N

  */

  lapack_int n, nrhs, lda, ldb, info;
  n = N; lda = n; ldb = 1; nrhs = 1;

  lapack_int* ipiv = (lapack_int *)malloc(n*sizeof(lapack_int)) ;

  // Solve the equations A*X = B 
  info = LAPACKE_dgesv( LAPACK_ROW_MAJOR, n, nrhs, A[0], lda, ipiv, b, ldb );

  // Check for the exact singularity 
  if( info > 0 ) {
    printf( "The diagonal element of the triangular factor of A,\n" );
    printf( "U(%i,%i) is zero, so that A is singular;\n", info, info );
    printf( "the solution could not be computed.\n" );
    exit( 1 );
  }
  if (info < 0) exit( 1 );

  // Free ipiv
  free(ipiv);

  // Overwritten by the solution x
  return b;
}

void qr(int M, int N, double** A, double** Q, double** R) {
  /*
    Compute QR decomposition of A. It returns the same output as numpy.linalg.qr(A, mode='reduced')

    Input:
    int M --> int dimension, rows of A
    int N --> int dimension, columns of A
    double** A --> matrix input of size (M x N).
    double** Q --> ortoghonal matrix of size (M x min(M,N))
    double** R --> upper triangular matrix of size (min(M,N) x N)

    Output:
    double** Q --> Output of Q in 1D-format. Use afterwards d_malloc_2d to save the matrix in 2D format.
    double** R --> Output of R in 1D-format. Use afterwards d_malloc_2d to save the matrix in 2D format.

  */

  // Allocate Q and R
  // Q size (m x  rank)
  double* Q0;
  if ( M < N)
    Q0 = malloc(M*M * sizeof(double));
  else
    Q0 = malloc(M*N * sizeof(double));
  *Q = Q0;

  // R size (rank x n)
  double* R0;
  if ( M < N)
    R0 = malloc(M*N * sizeof(double));
  else
    R0 = malloc(M*M * sizeof(double));
  *R = R0;

  // Create argument variables
  lapack_int m, n, lda, info;
  m = M, n = N; lda = n; 

  // rank = min(m,n)
  int rank;
  if ( m < n)
    rank = m;
  else
    rank = n;

  double* tau = malloc(rank * sizeof(double)) ;

  // Computes the QR factorization of a general m-by-n matrix.
  info = LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, m, n, A[0], lda, tau);
  if (info <0) exit( 1 );

  // Copy the upper triangular Matrix R (rank x n) into position
  // Set the entire R to 0
  if (m < n)
    for ( int i = 0 ; i < M*N ; i++ ) 
      R0[i] = 0.;
  else
    for ( int i = 0 ; i < M*M ; i++ ) 
      R0[i] = 0.;
  
  for(int i=0; i < rank; ++i) 
    memcpy(R0+i*n+i, A[0]+n*i+i, (n-i)*sizeof(double)); // Copy upper triangular part from Lapack result.
  

  // Create orthogonal matrix Q (in tmpA)
  if (m < n)
    LAPACKE_dorgqr(LAPACK_ROW_MAJOR, m, m, rank, A[0], n, tau);
  else
    LAPACKE_dorgqr(LAPACK_ROW_MAJOR, m, n, rank, A[0], n, tau);

  // Copy results to Q
  if (m < n) {
    for(int i=0; i<m; ++i)
      memcpy(Q0+i*m, A[0]+i*n, sizeof(double)*m);
  }
  else if (m==n)
    memcpy(Q0, A[0], sizeof(double)*(m*m));
  else {
    for(int i=0; i<m; ++i)
      memcpy(Q0+i*n, A[0]+i*n, sizeof(double)*n);
  }  

  // free
  free(tau);

}

/*
int main() {
  
  /////////////////////////////
  //    TEST LINEAR SOLVE   //
  ///////////////////////////

  double** B = d_malloc_2d(3,3);
  double* b = (double *)malloc(3*sizeof(double)) ;

  B[0][0] = 0; B[0][1] =  5; B[0][2] = 2; 
  B[1][0] = 6; B[1][1] = 7; B[1][2] =  3; 
  B[2][0] = 4; B[2][1] = -2; B[2][2] =  -1;
  b[0] = 1; b[1] =  2; b[2] = 3; 

  b = linear_system_solver(3, b, B);

  for (int i = 0; i < 3; i++) 
    printf("%f ", b[i]);
  
  printf("\n \n \n");
  
  ///////////////////
  //    TEST QR   //
  /////////////////

  free(B);
  int m, n;
  m = 3;
  n = 5;
  B = d_malloc_2d(m,n);  
  // We have to reconstruct B in order to do QR decomposition
  B[0][0] = 0; B[0][1] =  5; B[0][2] = 2; B[0][3] = 1; B[0][4] = 0; 
  B[1][0] = 6; B[1][1] = 7; B[1][2] =  3; B[1][3] = 1; B[1][4] = -13; 
  B[2][0] = 4; B[2][1] = -2; B[2][2] = -1; B[2][3] = 1; B[2][4] = -21; 

  double *R1, *Q1;
  
  qr(m, n, B, &Q1, &R1);

  // m < n
  if (m<n) {
    double**Q = d_malloc_2d(m,m) ; memcpy(Q[0], Q1, sizeof(double)*m*m);
    double**R = d_malloc_2d(m,n) ; memcpy(R[0], R1, sizeof(double)*m*n);

    int ii = 0;
    for ( int i = 0 ; i < m ; i++ ) {
      for ( int j = 0 ; j < m ; j++ ) {
        printf("%f ", Q[0][ii++]);
      }
      printf("\n");
    }

    printf("\n \n");

    ii = 0;
    for ( int i = 0 ; i < m ; i++ ) {
      for ( int j = 0 ; j < n ; j++ ) {
        printf("%f ", R[0][ii++]);
      }
      printf("\n");
    }
  }
  // m >= n
  else {
    double**Q = d_malloc_2d(m,n) ; memcpy(Q[0], Q1, sizeof(double)*m*n);
    double**R = d_malloc_2d(n,n) ; memcpy(R[0], R1, sizeof(double)*n*n);

   int ii = 0;
    for ( int i = 0 ; i < m ; i++ ) {
      for ( int j = 0 ; j < n ; j++ ) {
        printf("%f ", Q[0][ii++]);
      }
      printf("\n");
    }

    printf("\n \n");

    ii = 0;
    for ( int i = 0 ; i < n ; i++ ) {
      for ( int j = 0 ; j < n ; j++ ) {
        printf("%f ", R[0][ii++]);
      }
      printf("\n");
    } 
  }

}

*/