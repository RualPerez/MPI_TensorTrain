#include <stdio.h>
#include <stdlib.h>
#include "alloc.h"

double** reshape_tensor_to_matrix(int rank_k_1, int poly_order, int rank_k, int type, double*** tn) {
    /*
    Reshape tensor tn of size (rank_k_1 x poly_order x rank_k) to matrix A of size:
        - type = 0: ( (rank_k_1 * poly_order) x rank_k )
        - type = 1: ( rank_k_1 x (poly_order* rank_k) )

    Input:
        int rank_k_1 --> first integer dimension of tn
        int poly_order --> second interger dimension of tn
        int rank_k --> third integer dimension of tn
        int type --> type 0 compresses the first and second dimension
                     type 1 compresses the second and third dimension
        double*** tn --> tensor of size (rank_k_1 x poly_order x rank_k)

    Output:
        double** V --> vector of size either rank_k_1 or rank_k

    */

    double** A;

    // reshape(tn, rank_k_1 * poly_order, rank_k) in Matlab
    if (type==0) {
        A = d_malloc_2d(rank_k_1 * poly_order, rank_k);
        
        for (int i = 0; i < rank_k_1; i++)
            for (int j = 0; j < poly_order; j++) 
                for (int k = 0; k < rank_k; k++)  
                    A[j*rank_k_1+i][k] = tn[i][j][k];
    }

    // reshape(tn, rank_k_1, poly_order * rank_k) in Matlab
    else {
        A = d_malloc_2d(rank_k_1, poly_order * rank_k);

        for (int i = 0; i < rank_k_1; i++)
            for (int j = 0; j < poly_order; j++) 
                for (int k = 0; k < rank_k; k++)  
                    A[i][k*poly_order+j] = tn[i][j][k];
    }

    return A;
}

double** reshape_matrix_to_matrix(int m, int n, int m_f, int n_f, double** A) {
    /*
    Reshape matrix A of size (m x n) to matrix of size (m_f x n_f).
    It reshapes as Matlab, in a columnwise way.

    Input:
        int m --> row integer dimension of A
        int n --> column interger dimension of A
        int m_f --> row integer dimension of A reshaped
        int n_f --> column integer dimension of A reshaped
        double** A --> matrix of size (m x n)

    Output:
        double** B --> matrix reshaped to (m_f x n_f)

    */

    double** B = d_malloc_2d(m_f, n_f);
    int global_pos;

    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++) {
            global_pos = j * m + i;
            B[global_pos%m_f][global_pos/m_f] = A[i][j];
        }
    
    free(A);

    return B;
}

double*** reshape_matrix_to_tensor(int m, int n, int m_f, int n_f, int k_f, double** A) {
    /*
    Reshape matrix A of size (m x n) to matrix of size (m_f x n_f).
    It reshapes as Matlab, in a columnwise way.

    Input:
        int m --> row integer dimension of A
        int n --> column interger dimension of A
        int m_f --> row integer dimension of A reshaped
        int n_f --> column integer dimension of A reshaped
        double** A --> matrix of size (m x n)

    Output:
        double** tn --> matrix reshaped to (m_f x n_f)

    */

    double*** tn = d_malloc_3d(m_f, n_f, k_f);
    int global_pos, global_pos_in_slice;

    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++) {
            // global position flattening the matrix in a column-wise way
            global_pos = j * m + i;
            // global position flattening the slice tn[:,:,global_pos/(m_f*n_f)]
            global_pos_in_slice = global_pos % (m_f*n_f);

            tn[global_pos_in_slice%m_f][global_pos_in_slice/m_f][global_pos/(m_f*n_f)] = A[i][j];
        }

    return tn;
}

double** reshape_vector_to_matrix(int m, int n, double* v) {
    /*
    Reshape vector v size m*n to matrix A of size (m x n).
    It reshapes as Matlab, in a columnwise way.

    Input:
        int m --> row integer dimension of A
        int n --> column interger dimension of A

    Output:
        double** A --> matrix of size (m x n)

    */

    double** A = d_malloc_2d(m, n);

    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++) 
            A[i][j] = v[j*m+i];
        
    free(v);

    return A;
}