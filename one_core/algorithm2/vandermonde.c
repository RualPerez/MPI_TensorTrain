#include <stdlib.h>
#include <stdio.h>
#include "alloc.h"

double** vandermonde_vec(int poly_order, int num_instances, int feature, double** X) {
    /*
    Create the vandermonde vector given a column of X.

    Input:
        int poly_order --> polynomial order for the vandermonde vectors / maxim power of the features
        int num_instances --> number of rows in X 
        int feature --> column of X that we will take to create the vander. vector
        double** X --> dataset
        
    Output:
        double** X_vandermonde --> matrix with vandermonde vectors

    */

    // Get powers
    int column_pos;
    double init_feature;

    double** X_vandermonde = d_malloc_2d(num_instances, poly_order);

    for (int i = 0; i < num_instances; i++) {
        column_pos = 0;
        init_feature = X[i][feature];
        X_vandermonde[i][column_pos] = 1; // power of 0
        column_pos += 1;
        for (int k = 1; k < poly_order; k++) { // power from 1 to poly_order[feature]
            X_vandermonde[i][column_pos] = X_vandermonde[i][column_pos-1] * init_feature;
            column_pos += 1;
        }
    }
    return X_vandermonde;
}