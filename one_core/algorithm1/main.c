#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "alloc.h"
#include "init.h"
#include "tt_contraction.h"
#include "vandermonde.h"
#include "read_write.h"
#include "read_write_int.h"

int main(int argc, char *argv []) {

	int num_features = read_int(argv[1]);
	int* poly_order = read_vector_int(argv[1]);
	int* rank = read_vector_int(argv[2]);
	
	/**************************************/
	// Allocate and initialize TT-cores (G)
	/**************************************/

	double **** G_list = malloc(num_features * sizeof(double ***) );
	char filename[256];
	for(int i=0; i<num_features; i++) {
		G_list[i] = d_malloc_3d(rank[i], poly_order[i], rank[i+1]);
		//ones_3D_tensor(rank[i], poly_order[i], rank[i+1], G_list[i]);
		
		// Create filename: ("G_i/G_%d.txt", i)
		snprintf(filename, sizeof(filename), "../../Data/algorithm1/G_%d.txt", i);
		G_list[i] = read_3Dtensor(filename);
	}


	/*********************************************/
	// Read dataset and create vandermonde vectors
	/*********************************************/

	int vandermonde_size = 0;
	for(int i=0; i<num_features; i++) 
		vandermonde_size += poly_order[i];
	
	double** X = read_matrix("Data/X.txt");
	int num_instances = read_int("Data/X.txt");

	double** X_vandermonde = d_malloc_2d(num_instances, vandermonde_size);
	vandermonde_vec(poly_order, num_instances, num_features, X, X_vandermonde);
	free(X); 
	
	// Only one instace in the dataset X
	/*double* X_vandermonde = malloc(num_features*vandermonde_size * sizeof(double ***) );
	X_vandermonde = read_vector("Data/X_vandermonde.txt");
	//ones_vector(num_features*poly_order, X_vandermonde);*/


	
	/******************************/
	// Algorithm 1: TT contraction
	/******************************/

	double *** V_list = malloc(num_features * sizeof(double **) );
	double* vandermonde;

	int column_pos;
	double * first_V;
	double * last_V = malloc(rank[num_features-1] * sizeof(double) );
	double * f;

	// Repeat algorithm 1 for each instance
	for (int ins=0; ins<num_instances; ins++) {

		column_pos = 0;
		// Multiply vandermonde vector with its corresponding TT_core (G)
		for (int i=0; i<num_features; i++) {
			vandermonde = malloc(poly_order[i] * sizeof(double) );
			memcpy(vandermonde, &(X_vandermonde[ins][column_pos]), poly_order[i]  * sizeof(double));
			column_pos += poly_order[i];
			
			V_list[i] = d_malloc_2d(rank[i], rank[i+1]); 
			tt_core_contraction(rank[i], poly_order[i], rank[i+1], vandermonde, V_list[i], G_list[i]);
		
			free(vandermonde);
		}

		// Reshape output from first and last TT-core
		first_V = malloc(rank[1] * sizeof(double) );
		reshape_matrix_to_vector(rank[0], rank[1], first_V, V_list[0]);
		reshape_matrix_to_vector(rank[num_features-1], rank[num_features], last_V, V_list[num_features-1]);
		free(V_list[0]); free(V_list[num_features-1]);

		// Compress multiplications of V
		f = first_V;
		for (int i=1; i<num_features-1; i++) {
			vector_matrix(rank[i], rank[i+1], &f, V_list[i]);
			free(V_list[i]);
		}
		
		vector_vector(rank[num_features-1], &f, last_V);

		// Print result
		printf("Regression prediction %f \n", *f);
		free(f);
	}

	// Free
	free(last_V);
	free(X_vandermonde); free(poly_order); free(rank);
	for (int i=0; i<num_features; i++) 
		free(G_list[i]);
	free(G_list);
	
    
    return 0;

}
