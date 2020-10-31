#include <stdio.h>
#include <stdlib.h>
#include "alloc.h"

double* read_vector(char* filename) {
    /*
    Read from a file a vector. 
    The filename must contain in the first line an integer containing the size of the vector.
    Elements of the vector must be separated by an space.
    
    Code inspired from https://cboard.cprogramming.com/c-programming/153674-reading-file-storing-values-array-help-please-o.html

    Input:
        char* filename --> name of the file

    Output:
        double* vector --> read vector

    */
    
    int n;
    FILE * fp;
 
    fp = fopen(filename, "r");
    
    // read the size of the vector
    fscanf(fp, "%d", &n);
    double* vector = malloc(n * sizeof(double) );

    // read rest of lines: the vector
    for (int i = 0; i < n; i++)
        fscanf(fp, "%lf", &vector[i]);

    // Same using a while:
    /*
    int i = 0;
    while (fscanf(fp, "%lf", &vector[i]) != EOF) 
        ++i;
    */
    
    /*
    // Display vector
    for (int i = 0; i < n; i++)
        printf("vector[%d] = %f \n",i, vector[i]);
    */

    fclose(fp);

    return vector;
}   

void write_vector(char* filename, int n, double* vector) {
    /*
    Write a vector to a file.
    The filename will contain in the first line an integer containing the size of the vector.
    Elements of the vector will be separated by an space.
    
    Input:
        char* filename --> name of the file
        int n --> vector size
        double* vector --> vector to be written

    Output:
        

    */
    
    FILE * fp;
    fp = fopen(filename, "w+");

    // write the size of the vector
    fprintf(fp, "%d \n", n); 

    // write the vector
    for (int i = 0; i < n; i++)
        fprintf(fp, "%lf ", vector[i]); 

    fclose(fp);
}   

double** read_matrix(char* filename) {
    /*
    Read from a file a matrix. 
    The filename must contain in the first line two integers containing the size of the matrix.
    Elements of the matrix must be separated by an space.

    Input:
        char* filename --> name of the file

    Output:
        double** matrix --> read matrix

    */
    
    int m, n;
    FILE * fp;
 
    fp = fopen(filename, "r");
    
    // read the size of the matrix
    fscanf(fp, "%d", &m); fscanf(fp, "%d", &n);
    double** matrix = d_malloc_2d(m, n);

    // read rest of lines: the matrix
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            fscanf(fp, "%lf", &matrix[i][j]);
    
    fclose(fp);
    
    /*
    // Display matrix
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            printf("matrix[%d][%d] = %f \n",i,j, matrix[i][j]);
    */ 

    return matrix;
}   

void write_matrix(char* filename, int m, int n, double** matrix) {
    /*
    Write a matrix to a file.
    The filename will contain in the first line two integers containing the size of the matrix.
    Elements of the matrix will be separated by an space and an end line after each row.
    
    Input:
        char* filename --> name of the file
        int m --> row size
        int n --> column size
        double** matrix --> matrix to be written

    Output:
        
    */
    
    FILE * fp;
    fp = fopen(filename, "w+");

    // write the size of the matrix
    fprintf(fp, "%d ", m);  fprintf(fp, "%d \n", n); 

    // write the matrix
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++)
            fprintf(fp, "%lf ", matrix[i][j]);
        fprintf(fp, "\n");
    }
        

    fclose(fp);
}   

double*** read_3Dtensor(char* filename) {
    /*
    Read from a file a tensor. 
    The filename must contain in the first line three integers containing the size of the tensor.
    Elements of the tensor must be separated by an space.

    Input:
        char* filename --> name of the file

    Output:
        double*** tensor --> read tensor

    */
    
    int m, n, k;
    FILE * fp;
 
    fp = fopen(filename, "r");
    
    // read the size of the tensor
    fscanf(fp, "%d", &m); fscanf(fp, "%d", &n); fscanf(fp, "%d", &k);
    double*** tensor = d_malloc_3d(m, n, k);

    // read rest of lines: the tensor
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            for (int h = 0; h < k; h++)
                fscanf(fp, "%lf", &tensor[i][j][h]);
    
    fclose(fp);
    
    /*
    // Display tensor
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            for (int h = 0; h < k; h++)
                printf("tensor[%d][%d][%d] = %f \n",i,j,h, tensor[i][j][h]);
    */ 

    return tensor;
}   

void write_3Dtensor(char* filename, int m, int n, int k, double*** tensor) {
    /*
    Write a tensor to a file.
    The filename will contain in the first line three integers containing the size of the tensor.
    Elements of the tensor will be separated by an space and an end line after each row and column.
    
    Input:
        char* filename --> name of the file
        int m --> row size
        int n --> column size
        int k --> third order size
        double*** tensor --> tensor to be written

    Output:
        
    */
    
    FILE * fp;
    fp = fopen(filename, "w+");

    // write the size of the tensor
    fprintf(fp, "%d ", m); fprintf(fp, "%d ", n); fprintf(fp, "%d \n", k); 

    // write the tensor
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            for (int h = 0; h < k; h++)
                fprintf(fp, "%lf ", tensor[i][j][h]);
            fprintf(fp, "\n");    
        }
        fprintf(fp, "\n");
    }        

    fclose(fp);
}   


// How to call the functions
/*
int main(int argc, char *argv []) {

    char* filename = argv[1];
    double* vector = read_vector(filename);
    write_vector("cp_vector.txt", 10, vector);

    filename = argv[2];
    double** matrix = read_matrix(filename);
    write_matrix("cp_matrix.txt", 2, 5, matrix);

    filename = argv[3];
    double*** tensor = read_3Dtensor(filename);
    write_3Dtensor("cp_tensor.txt", 2, 5, 1, tensor);
}
*/