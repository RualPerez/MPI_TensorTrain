#include <stdlib.h>
#include <stdio.h>

void vandermonde_vec(int* poly_order, int num_instances, int num_features, double** X, double** X_vandermonde) {
    /*
	Creates the vandermonde vectors.

	Input:
		int* poly_order --> polynomial order for the vandermonde vectors / maxim power of the features
        int num_instances --> number of rows in X 
        int num_features --> number of columns in X
        double** X --> dataset
		double** X_vandermonde --> empty matrix where to save the output
		
	Output:
		double** X_vandermonde --> dataset with vandermonde vectors

    */

    // Get powers
    int column_pos;
    double init_feature;

    for (int i = 0; i < num_instances; i++) {
        column_pos = 0;
    	for (int j = 0; j < num_features; j++) {
    		init_feature = X[i][j];
    		X_vandermonde[i][column_pos] = 1; // power of 0
            column_pos += 1;
    		for (int k = 1; k < poly_order[j]; k++) { // power from 1 to poly_order[j]
                X_vandermonde[i][column_pos] = X_vandermonde[i][column_pos-1] * init_feature;
                column_pos += 1;
            }
    	}
    }
    return; // X_vandermonde;
}