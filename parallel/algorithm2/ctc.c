#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include "alloc.h"
#include "index_size_functions.h"
#include "ctc.h"

double** ctc(double** C, int N, int M, int rank, int size, int blocksize) {
    // Computes the matrix-matrix product L = transpose(C)*C.
    
    // need to distribute the rows in C, and to start
    // to compute the product right away
    
    // store part results in 1d-array, reduce this afterwards
    int COUNT = M*(M+1)/2;
    double *L_elements_local = malloc(COUNT*sizeof(double));
    for (int i = 0; i < COUNT; i++)
        L_elements_local[i] = 0.0;
    
    // Local C-matrix on rank
    int N_local = local_size(N, blocksize, rank, size, 1, 'r');
    double** C_local = d_malloc_2d(N_local, M);
    
    
    // create data in arrays to be used in send-receive loop
    
    int full_blocks = N / blocksize;
    int elements_in_block;
    
    int current_row;
    int current_blocksize;
    int dest;
    int K = full_blocks;
    int non_zero_dests = 0;
    
    int *dests = malloc((full_blocks+1)*sizeof(int));
    int *current_rows = malloc((full_blocks+1)*sizeof(int));
    int *current_blocksizes = malloc((full_blocks+1)*sizeof(int));
    
    for (int i = 0; i <= full_blocks; i++) {
        current_row = i*blocksize;
        if (i == full_blocks) {
            if (current_row < N) {
                K = full_blocks + 1;
                dests[i] = global2rank(current_row, blocksize, size);
                if (dests[i] != 0) {
                    non_zero_dests++;
                }
                current_rows[i] = current_row;
                current_blocksizes[i] = N-current_row;
            } else {
                dests[i] = -1; // if used then error in code
                current_rows[i] = -1; // if used then error in code
                current_blocksizes[i] = -1; // if used then error in code
            }
        } else {
            dests[i] = global2rank(current_row, blocksize, size);
            if (dests[i] != 0) {
                non_zero_dests++;
            }
            current_rows[i] = current_row;
            current_blocksizes[i] = blocksize;
        }
    }
    
    // do send-receive loop
    
    MPI_Request requests[non_zero_dests];
    
    if (rank == 0) {
        int request_counter = 0;
        for (int i = 0; i < K; i++) {
            
            // fetch variables
            current_row = current_rows[i];
            dest = dests[i];
            current_blocksize = current_blocksizes[i];
            elements_in_block = M*current_blocksize;
            //printf("rank %d, dest %d, currentblocksize %d, elements in block %d\n", rank, dest, current_blocksize, elements_in_block);
            
            // compute for rank 0 or send to other ranks
            if (dest == 0) {
                int counter = 0;
                for (int h = 0; h < M; h++) {
                    for (int j = h; j < M; j++) {
                        // For the chosen combinations of columns, do a
                        // dot product of columns i and j in C.
                        for (int i_loc = current_row; i_loc < current_row+current_blocksize; i_loc++) {
                            L_elements_local[counter] += C[i_loc][h]*C[i_loc][j];
                        }
                        counter++;
                    }
                }
                //printf("L_elements_local[0] %f\n", L_elements_local[0]);
            } else {
                //printf("SEND PARAMS dest %d, tag %d, request_counter %d, elements_in_block %d\n", dest, i, request_counter, elements_in_block);
                MPI_Isend(&(C[current_row][0]), elements_in_block, MPI_DOUBLE, dest, i, MPI_COMM_WORLD, &requests[request_counter]);
                request_counter++;
                //printf("SENDING DONE\n");
                // int MPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request * request)
            }
        }
        MPI_Waitall(non_zero_dests, requests, MPI_STATUSES_IGNORE);
    } else {
        int current_local_row;
        // receive the rows
        for (int k = 0; k < K; k++) {
            dest = dests[k];
            if (dest == rank) {
                
                current_row = current_rows[k]; // global
                current_blocksize = current_blocksizes[k];
                elements_in_block = M*current_blocksize;
                current_local_row = global2local(current_row, blocksize, size);
                
                // receive row(s) of C-matrix
                MPI_Recv(&(C_local[current_local_row][0]), elements_in_block, MPI_DOUBLE, 0, k, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                
                // compute the dot products of the columens in the local C
                int counter = 0;
                for (int h = 0; h < M; h++) {
                    for (int j = h; j < M; j++) {
                        // For the chosen combinations of columns, do a
                        // dot product of columns i and j in C.
                        for (int i_loc = current_local_row; i_loc < current_local_row+current_blocksize; i_loc++) {
                            L_elements_local[counter] += C_local[i_loc][h]*C_local[i_loc][j];
                        }
                        counter++;
                    }
                }
            }
        }
    }
    
    // free temp arrays
    free(C_local);
    free(dests);
    free(current_rows);
    free(current_blocksizes);

    // sum these array element by element on rank 0
    double *L_elements = malloc(COUNT*sizeof(double));
    
    MPI_Reduce(L_elements_local,    // void* send_data
               L_elements,         // void* recv_data
               COUNT,               // int count
               MPI_DOUBLE,          // MPI_Datatype datatype
               MPI_SUM,             // MPI_Op op
               0,                   // int root
               MPI_COMM_WORLD);     // MPI_Comm communicator

    // free the used arrays
    free(L_elements_local);
    
    // L = transpose(C)*C
    double **L = d_malloc_2d(M, M);
    
    if (rank == 0) {
        // place the summed array in the correct places in L
        int counter = 0;
        for (int i = 0; i < M; i++) {
            L[i][i] = L_elements[counter];
            counter++;
            for (int j = i+1; j < M; j++) {
                L[i][j] = L_elements[counter];
                L[j][i] = L_elements[counter];
                counter++;
            }
        }
    } else {
        L = NULL;
    }
    // free the used arrays
    free(L_elements);
    return L;
}
