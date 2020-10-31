#include <stdio.h>
#include <stdlib.h>
#include "alloc.h"

double** matmult_nat(int m, int n, int k, int free_memory, double **A, double **B) {
    
    int i, j, h;
    double **C = d_malloc_2d(m, n);

    // Initializing elements of matrix mult to 0.
    for(i = 0; i < m; ++i) 
        for(j = 0; j < n; ++j) 
            C[i][j] = 0;

    // Multiplying matrix firstMatrix and secondMatrix and storing in array mult.
    for(i = 0; i < m; ++i) 
        for(j = 0; j < n; ++j) 
            for(h=0; h < k; ++h) 
                C[i][j] += A[i][h] * B[h][j];

    if (free_memory==1) {
        free(A); free(B);
    }

    return C;
}

double* matmult_vec(int m, int n, double* v, double** A) {
    /*
    Matrix vector multiplication: u = A * v

    Input:
        int m --> row dimension of A
        int n --> column dimension of A and dimension of v
        double* v --> vec of size n
        double** A --> matrix of size (m x n)

    Output:
        double* u --> vector of size m

    */

    double* u = malloc(m*sizeof(double));

    for(int i = 0; i < m; ++i) 
        u[i] = 0;

    for(int i = 0; i < m; ++i) 
        for(int j = 0; j < n; ++j) 
            u[i] += A[i][j] * v[j];

    return u;
}

double** transpose(int m, int n, int free_memory, double** A) {
    
    double **B = d_malloc_2d(n, m);

    // Transpose matrix
    for(int i = 0; i < m; ++i) 
        for(int j = 0; j < n; ++j) 
            B[j][i] = A[i][j];
    
    if (free_memory)
        free(A);
    
    return B;
}

double mse(int m, double* predictions, double* y) {
    /*
    Compute the Mean Squared Error between the predictions and the y vectors.

    Input:
        int m --> dimension of predictions and y, i.e., number of instances
        double* predictions --> vec of size m
        double* y --> vec of size m (target values)

    Output:
        double error --> MSE

    */
    double error = 0;
    double l1;

    for(int i=0; i<m; ++i) {
        l1 = predictions[i]-y[i];
        error += l1*l1;
    }

    error = error / ((double) m);

    return error;
}