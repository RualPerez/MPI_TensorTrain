#include <stdlib.h>
#include <stdio.h>

double*** d_malloc_3d(int m, int n, int k) {

    if (m <= 0 || n <= 0 || k <= 0)
        return NULL;

    double ***array3D = malloc(m * sizeof(double **) +
                               m * n * sizeof(double *) +
                               m * n * k * sizeof(double));
    if (array3D == NULL) {
        return NULL;
    }

    for(int i = 0; i < m; i++) {
        array3D[i] = (double **) array3D + m + i * n ;
    }

    for(int i = 0; i < m; i++) {
        for(int j = 0; j < n; j++) {
            array3D[i][j] = (double *) array3D + m + m * n + i * n * k + j * k;
        }
    }

    return array3D;
}


// allocate a double-prec m x n matrix
double** d_malloc_2d(int m, int n) {
    if (m <= 0 || n <= 0) return NULL;
    double **A = malloc(m * sizeof(double *) +
                        m*n*sizeof(double) );

    if (A == NULL) {
        printf("Fail to allocate matrix \n");
        return NULL;
    }

    for (int i = 0; i < m; i++)
        A[i] = (double *) A + m + i * n;

    return A;
}


// allocate a int m x n matrix
int** d_malloc_2d_int(int m, int n) {
    if (m <= 0 || n <= 0) return NULL;
    int **A = malloc(m * sizeof(int*) +
                     m*n*sizeof(int) );

    if (A == NULL) {
        printf("Fail to allocate matrix \n");
        return NULL;
    }

    for (int i = 0; i < m; i++)
        A[i] = (int *) A + m + i * n;

    return A;
}
