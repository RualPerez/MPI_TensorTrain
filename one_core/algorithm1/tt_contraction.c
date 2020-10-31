#include <stdio.h>
#include <stdlib.h>

void tt_core_contraction(int rank_k_1, int n_k, int rank_k, double* v, double** V, double*** G) {
	/*
	2-mode multiplication between 3-tensor and vector: V = G x_2 vec

	Input:
		int rank_k_1 --> integer dimension
		int n_k --> integer dimension
		int rank_k --> integer dimension
		double* vec --> vandermonde vector of size n_k
		double** V --> empty matrix where to save the output
		double*** G --> 3-tensor of size (rank_k_1 x n_k x rank_k)

	Output:
		double** V --> matrix of size (rank_k_1 x rank_k)

    */
	
	// Init with zeros
	for (int i = 0; i < rank_k_1; i++)
    	for (int j = 0; j < rank_k; j++) 
    		V[i][j] = 0;
    	
    // G(j) is a rank_k_1 x rank_k matrix for each index j
	for (int i = 0; i < rank_k_1; i++)
    	for (int j = 0; j < n_k; j++) 
    		for (int k = 0; k < rank_k; k++)  
    			V[i][k] += G[i][j][k] * v[j];
    		
}


void reshape_matrix_to_vector(int rank_k_1, int rank_k, double* v, double** V) {
	/*
	Reshape matrix V of size (rank_k_1 x rank_k) to vector v
		note that either rank_k_1 = 1  or  rank_k = 1

	Input:
		int rank_k_1 --> integer dimension
		int rank_k --> integer dimension
		double* v --> empty vector where to save the output
		double** V --> matrix of size (rank_k_1 x rank_k)

	Output:
		double** V --> vector of size either rank_k_1 or rank_k

    */

	// Check correct input dimensions
    if (rank_k_1 != 1 && rank_k != 1) {
    	printf("Error in reshaping, either rank_k_1 or rank_k must be equal 1 \n");
        return;
    }

    // Reshape for last V obtained in the tt_core_contraction
    if (rank_k_1 > 1) {
    	for (int i = 0; i < rank_k_1; i++)
    		v[i] = V[i][0];
    }

    // Reshape for first V obtained in the tt_core_contraction
    else {
    	for (int j = 0; j < rank_k; j++) 
    		v[j] = V[0][j];
    	
    }

}

void vector_matrix(int rank_k, int rank_k1, double** f, double** V) {
	/*
	Vector matrix multiplication: f = f * V

	Input:
		int rank_k --> integer dimension
		int rank_k1 --> integer dimension, note that indicates rank_(k+1)
		double** f --> pointer to vector of size rank_k + pointer to vector where to save the output
		double** V --> matrix of size (rank_k x rank_k1)

	Output:
		double* f --> vector of size either rank_k1

    */

	// Auxiliar vector where to save the result
	double* output = malloc( rank_k1 * sizeof(double));

	// Init with zeros
	for (int j = 0; j < rank_k1; j++) 
			output[j] = 0;
	
	// Vector matrix multiplication
	for (int i = 0; i < rank_k; i++) 
		for (int j = 0; j < rank_k1; j++) 
			output[j] += (*f)[i] * V[i][j];

	// Free and save result to f
	free(*f);
	*f = output;
}
 

void vector_vector(int rank_k, double** f, double* V) {
	/*
	Vector vector multiplication: f = f * V

	Input:
		int rank_k --> integer dimension (of last V)
		double* f --> vector of size rank_k + vector where to save the output
		double* V --> vector of size rank_k (last V)

	Output:
		double* f --> regression prediction, a double

    */

	// Auxiliar number where to save the result
	double* output = malloc(1 * sizeof(double));
	*output = 0;

	// Vector vector multiplication
	for (int i = 0; i < rank_k; i++) 
		*output += (*f)[i] * V[i];

	// Free and save result to f
	free(*f);
	*f = output;
}