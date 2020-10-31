#include <stdlib.h>

void ones_3D_tensor(int rank_k_1, int n_k, int rank_k, double ***tn) {

    // Initialise the values of tn to one

    int i,j,k;
    
    for (i = 0; i < rank_k_1; ++i)
        for (j = 0; j < n_k; ++j)
            for (k = 0; k < rank_k; ++k)
                tn[i][j][k] = 1;

}


void ones_matrix(int m, int n, double **A) {

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