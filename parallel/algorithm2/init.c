#include <stdlib.h>
#include <stdio.h>
#include <math.h>

double gaussrand() {
    /*
     Generate a random number with a normal or Gaussian distribution.
     Method discussed in Knuth and due originally to Marsaglia.
    */

    static double V1, V2, S;
    static int phase = 0;
    double X;

    if(phase == 0) {
        do {
            double U1 = (double)rand() / RAND_MAX;
            double U2 = (double)rand() / RAND_MAX;

            V1 = 2 * U1 - 1;
            V2 = 2 * U2 - 1;
            S = V1 * V1 + V2 * V2;
            } while(S >= 1 || S == 0);

        X = V1 * sqrt(-2 * log(S) / S);
    } else
        X = V2 * sqrt(-2 * log(S) / S);

    phase = 1 - phase;

    return X;
}

void zeros_3D_tensor(int rank_k_1, int n_k, int rank_k, double ***tn) {

    // Initialise the values of tn to one

    int i,j,k;
    
    for (i = 0; i < rank_k_1; ++i)
        for (j = 0; j < n_k; ++j)
            for (k = 0; k < rank_k; ++k)
                tn[i][j][k] = 0;

}

void ones_3D_tensor(int rank_k_1, int n_k, int rank_k, double ***tn) {

    // Initialise the values of tn to one

    int i,j,k;
    
    for (i = 0; i < rank_k_1; ++i)
        for (j = 0; j < n_k; ++j)
            for (k = 0; k < rank_k; ++k)
                tn[i][j][k] = 1;

}


void rdm_3D_tensor(int rank_k_1, int n_k, int rank_k, double ***tn) {

    /* 
    Initialise the values of tn based on a normal distribution.
    Moreover, afterwards it is divided by its norm.
    */
    int i,j,k;
    double  norm = 0;
    
    for (i = 0; i < rank_k_1; ++i)
        for (j = 0; j < n_k; ++j)
            for (k = 0; k < rank_k; ++k) {
                tn[i][j][k] = gaussrand();
                norm += tn[i][j][k]*tn[i][j][k];
            }

    norm = sqrt(norm);

    // Divide each element by the norm of the tensor
    for (i = 0; i < rank_k_1; ++i)
        for (j = 0; j < n_k; ++j)
            for (k = 0; k < rank_k; ++k) 
                tn[i][j][k] /= norm;

}

void ones_matrix_void(int m, int n, double **A) {

    // Initialise the values of A to one

    int i,j;

    for (i = 0; i < m; ++i)
        for (j = 0; j < n; ++j)
            A[i][j] = 1;

}


void ones_vector(int m, double *v) {

    // Initialise the values of v to one
    
    int i;

    for (i = 0; i < m; ++i)
        v[i] = 1;

}
