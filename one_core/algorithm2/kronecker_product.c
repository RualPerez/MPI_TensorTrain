#include <stdio.h>
#include <stdlib.h>
#include "alloc.h"

double** kroneckerproduct(int m, int n, int p, int q, double** A, double** B) {
	/*
	Kronecker product between two matrices: https://ide.geeksforgeeks.org/YnaiKYwcwR

	Input:
		int m --> row dimension of A
		int n --> column dimension of A
		int p --> row dimension of B
		int  q--> column dimension of B
		double** A --> left matrix in the kronecker product
		double** B --> right matrix in the kronecker product

	Output:
		double** C --> matrix of size (p*m x q*n)

    */
	
	double** C = d_malloc_2d(p*m, q*n);

	for (int i = 0; i < p*m; i++) 
		for (int j = 0; j < q*n; j++) 
			C[i][j] = A[i/p][j/q]*B[i%p][j%q];
	
	return C;
}

double** repmat(int repetition_row, int repetition_col, int m, int n, double** A) {
	/*
	creates a large matrix B consisting of an repetition_row-by-repetition_col tiling of copies of A. 
	If A is a matrix, the size of B is [m*repetition_row, n*repetition_col].

	Input:
		int repetition_row --> times to repeat A along vertical/row dimension
		int repetition_row --> times to repeat A along horizontal/column dimension
		int m --> row dimension of A
		int  n--> column dimension of A
		double** A --> input matrix

	Output:
		double** B --> matrix of size (m*repetition_row x n*repetition_col)

    */

	double** B = d_malloc_2d(m*repetition_row, n*repetition_col);

	for (int i = 0; i < m*repetition_row; i++) 
		for (int j = 0; j < n*repetition_col; j++) 
			B[i][j] = A[i%m][j%n];

	return B;

}

double** element_wise_mult(int m, int n, double** A, double** B) {
	/*
	Element-wise multiplication / Hadamard product between two matrices 

	Input:
		int m --> row dimension of A and B
		int n --> column dimension of A and B
		double** A --> input matrix
		double** B --> input matrix

	Output:
		double** C --> output matrix 

    */

	double** C = d_malloc_2d(m, n);

	for (int i = 0; i < m; i++) 
		for (int j = 0; j < n; j++) 
			C[i][j] = A[i][j]*B[i][j];

	return C;

}

double** ones_matrix(int m, int n) {

	// Initialise the values of A to one

	double **A = d_malloc_2d(m,n);

    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
            A[i][j] = 1;

    return A;

}


double** dotkron(int m, int n, int p, int q, int u, int v, double** A, double** B, double** C) {
	/*
	The Kroneck product of same rows for given 2 or 3 matrices 

	Input:
		int m --> row dimension of A
		int n --> column dimension of A
		int p --> row dimension of B
		int q--> column dimension of B
		double** A --> left matrix in the kronecker product of same rows
		double** B --> middle/right matrix in the kronecker product of same rows (2/3 matrices)
		double** B --> right matrix in the kronecker product of same rows (3 matrices)

	Output:
		double** D --> matrix of size (m x n*q) for 2 matrices
					   matrix of size (m x n*q*v) for 3 matrices
	*/

	double** y;

	// 2-matrix case
	if ((u==0) && (v==0)) {
		if (m != p)
			printf("Matrices should have equal rows in dotkron. \n" );
		else {
			double** repmatrix = repmat(1, q, m, n, A);
			double** kron = kroneckerproduct(p, q, 1, n, B, ones_matrix(1, n)); 
			y = element_wise_mult(p, q*n, repmatrix, kron);
		}
	}

	// 3-matrix case
	else {
		if ((m != p) || (p != u))
			printf("Matrices should have equal rows in dotkron. \n" );
		else
			y= dotkron(m, n, m, q*v, 0, 0, A, dotkron(p, q, u, v, 0, 0, B, C, NULL), NULL);
	}

	return y;
}



/*
int main() {

	double** A = d_malloc_2d(3,2);
	double** B = d_malloc_2d(3,3); 

    A[0][0] = 1; A[0][1] =  2; A[1][0] = 3; A[1][1] = 4; A[2][0] = 1; A[2][1] =  0;
    B[0][0] = 0; B[0][1] =  5; B[0][2] = 2; B[1][0] = 6; B[1][1] = 7; B[1][2] =  3; B[2][0] = 4; B[2][1] = -2; B[2][2] =  -1;
  
    double**C = dotkron(3, 2, 3, 3, 3, 3, A, B, B); 

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 18; j++)  
			printf("%f ", C[i][j]);
	printf("\n");
	}

}
*/