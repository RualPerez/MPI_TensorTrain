#include <stdio.h>
#include <stdlib.h>
#include "alloc.h"

int read_int(char* filename) {
    /*
    Read from a file the first integer. 
    
    Input:
        char* filename --> name of the file

    Output:
        int integer --> read integer

    */
    
    int integer;
    FILE * fp;
 
    fp = fopen(filename, "r");
    fscanf(fp, "%d", &integer);
    fclose(fp);

    return integer;
}   

int* read_vector_int(char* filename) {
    /*
    Read from a file a vector. 
    The filename must contain in the first line an integer containing the size of the vector.
    Elements of the vector must be separated by an space.
    
    Code inspired from https://cboard.cprogramming.com/c-programming/153674-reading-file-storing-values-array-help-please-o.html

    Input:
        char* filename --> name of the file

    Output:
        int* vector --> read vector

    */
    
    int n;
    FILE * fp;
 
    fp = fopen(filename, "r");
    
    // read the size of the vector
    fscanf(fp, "%d", &n);
    int* vector = malloc(n * sizeof(int) );

    // read rest of lines: the vector
    for (int i = 0; i < n; i++)
        fscanf(fp, "%d", &vector[i]);

    // Same using a while:
    /*
    int i = 0;
    while (fscanf(fp, "%d", &vector[i]) != EOF) 
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

void write_vector_int(char* filename, int n, int* vector) {
    /*
    Write a vector to a file.
    The filename will contain in the first line an integer containing the size of the vector.
    Elements of the vector will be separated by an space.
    
    Input:
        char* filename --> name of the file
        int n --> vector size
        int* vector --> vector to be written

    Output:
        

    */
    
    FILE * fp;
    fp = fopen(filename, "w+");

    // write the size of the vector
    fprintf(fp, "%d \n", n); 

    // write the vector
    for (int i = 0; i < n; i++)
        fprintf(fp, "%d ", vector[i]); 

    fclose(fp);
}   

int** read_matrix_int(char* filename) {
    /*
    Read from a file a matrix. 
    The filename must contain in the first line two integers containing the size of the matrix.
    Elements of the matrix must be separated by an space.

    Input:
        char* filename --> name of the file

    Output:
        int** matrix --> read matrix

    */
    
    int m, n;
    FILE * fp;
 
    fp = fopen(filename, "r");
    
    // read the size of the matrix
    fscanf(fp, "%d", &m); fscanf(fp, "%d", &n);
    int** matrix = d_malloc_2d_int(m, n);

    // read rest of lines: the matrix
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            fscanf(fp, "%d", &matrix[i][j]);
    
    fclose(fp);
    
    /*
    // Display matrix
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            printf("matrix[%d][%d] = %f \n",i,j, matrix[i][j]);
    */ 

    return matrix;
}   

void write_matrix_int(char* filename, int m, int n, int** matrix) {
    /*
    Write a matrix to a file.
    The filename will contain in the first line two integers containing the size of the matrix.
    Elements of the matrix will be separated by an space and an end line after each row.
    
    Input:
        char* filename --> name of the file
        int m --> row size
        int n --> column size
        int** matrix --> matrix to be written

    Output:
        
    */
    
    FILE * fp;
    fp = fopen(filename, "w+");

    // write the size of the matrix
    fprintf(fp, "%d ", m);  fprintf(fp, "%d \n", n); 

    // write the matrix
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++)
            fprintf(fp, "%d ", matrix[i][j]);
        fprintf(fp, "\n");
    }
        

    fclose(fp);
}   

// How to call the functions
/*
int main(int argc, char *argv []) {

    char* filename = argv[1];
    int* vector = read_vector(filename);
    write_vector("cp_vector.txt", 10, vector);

    filename = argv[2];
    int** matrix = read_matrix(filename);
    write_matrix("cp_matrix.txt", 2, 5, matrix);

}
*/