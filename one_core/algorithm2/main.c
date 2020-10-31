#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "alloc.h"
#include "lapack_routines.h"
#include "kronecker_product.h"
#include "read_write.h"
#include "read_write_int.h"
#include "reshape.h"
#include "operations.h"
#include "init.h"
#include "vandermonde.h"

int main(int argc, char *argv []) {
    /*
    The executeable takes the following input arguments:
    Example:
        Data/poly_order_2.txt   --> File with number vandermonde vectors and dimensions
        Data/rank_2.txt         --> File with number of TT-ranks and each Grank
        Data/X_20640.txt        --> File with number of instances, dimension and all data variables
        Data/y_20640.txt        --> File with number of instancesand all data outputs variables
        48                      --> Blocksize
    */

    // Checking number of input arguments
    if (argc < 6) printf("\nERROR! Too few input arguments!\n");

    // Read hyperparameters (command-line inputs)
    int num_features = read_int(argv[1]);
    int* poly_order = read_vector_int(argv[1]);
    int* Grank = read_vector_int(argv[2]);

    // Read dataset
    int num_instances = read_int(argv[3]);
    double** X = read_matrix(argv[3]);
    double* y = read_vector(argv[4]);
     
    // Read blocksize
    int block_size = atoi(argv[5]);

    /**************************************/
    // Allocate and initialize TT-cores (G)
    /**************************************/

    srand48(time(NULL)); srand(time(NULL));

    double **** G_list = malloc(num_features * sizeof(double ***) );
    for(int i=0; i<num_features; i++) {
        G_list[i] = d_malloc_3d(Grank[i], poly_order[i], Grank[i+1]);
        rdm_3D_tensor(Grank[i], poly_order[i], Grank[i+1], G_list[i]);
    }

    /**************************************/
    //        Left-to-right sweep: 
    //      Orthogonalize the G_list
    /**************************************/

    double** Q, ** R, ** temp;
    double *R1, *Q1;
    int m_r, n_r, m, n;

    R = d_malloc_2d(1,1); R[0][0] = 1;
    m_r = 1; n_r = 1;

    for(int i=0; i<num_features; ++i) {
        temp = reshape_tensor_to_matrix(Grank[i], poly_order[i], Grank[i+1], 1, G_list[i]); free(G_list[i]);
        temp = matmult_nat(m_r,poly_order[i]*Grank[i+1],n_r,1,R,temp); // frees R and temp
        Grank[i] = m_r;
        temp = reshape_matrix_to_matrix(m_r, poly_order[i]*Grank[i+1], m_r*poly_order[i], Grank[i+1], temp); // frees temp
        qr(m_r*poly_order[i], Grank[i+1], temp, &Q1, &R1); free(temp);
        // set size of R
        n_r = Grank[i+1];

        m = m_r*poly_order[i]; n = Grank[i+1];
        // Assign Q to the i-th TT-core: G_list[i]
        if (m < n) {// QR we have m < n --> Q is (m x m) and R is (m x n)
            // convert Q1 and R1 in a 2d-allocation
            Q = d_malloc_2d(m,m); memcpy(Q[0], Q1, sizeof(double)*m*m);
            R = d_malloc_2d(m,n); memcpy(R[0], R1, sizeof(double)*m*n);
            
            G_list[i] = reshape_matrix_to_tensor(m_r*poly_order[i], m_r*poly_order[i], m_r, poly_order[i], m_r*poly_order[i], Q); free(Q);
            
            // set size of R
            m_r = m_r*poly_order[i];
        }
        else {// QR we have m >= n --> Q is (m x n) and R is (n x n)
            // convert Q1 and R1 in a 2d-allocation
            Q = d_malloc_2d(m,n); memcpy(Q[0], Q1, sizeof(double)*m*n);
            R = d_malloc_2d(n,n); memcpy(R[0], R1, sizeof(double)*n*n);

            G_list[i] = reshape_matrix_to_tensor(m_r*poly_order[i], Grank[i+1], m_r, poly_order[i], Grank[i+1], Q); free(Q);

            //set size of R
            m_r = Grank[i+1]; 
        }
        free(Q1); free(R1);
    }


    /**********************************************/
    // Allocate and initialize matrices Matd,
    // which is defined as p_k and q_k in the paper
    /**********************************************/

    double *** Matd = malloc( (num_features+1) * sizeof(double ***) );
    Matd[num_features] = ones_matrix(num_instances, Grank[num_features]);

    double** X_vandermonde;

    /**************************************/
    //        Left-to-right sweep: 
    //      Orthogonalize the G_list
    //                &
    //           Obtain Matd
    /**************************************/

    double** tmp;

    for(int i=num_features-1; i>0; --i) {
        temp = reshape_tensor_to_matrix(Grank[i], poly_order[i], Grank[i+1], 0, G_list[i]); free(G_list[i]);
        temp = matmult_nat(Grank[i]*poly_order[i], n_r, Grank[i+1], 1, temp, R); // frees R and temp
        Grank[i+1] = n_r;
        temp = reshape_matrix_to_matrix(Grank[i]*poly_order[i], Grank[i+1], Grank[i], poly_order[i]*Grank[i+1], temp); // frees temp
        temp = transpose(Grank[i], poly_order[i]*Grank[i+1], 1, temp); // frees temp
        qr(poly_order[i]*Grank[i+1], Grank[i], temp, &Q1, &R1); free(temp);
        // set size of R
        m_r = Grank[i];

        m = poly_order[i]*Grank[i+1]; n = Grank[i];
        // Assign Q to the i-th TT-core: G_list[i]
        if (m < n) {// QR we have m < n --> Q is (m x m) and R is (m x n)
            // convert Q1 and R1 in a 2d-allocation
            Q = d_malloc_2d(m,m); memcpy(Q[0], Q1, sizeof(double)*m*m);
            R = d_malloc_2d(m,n); memcpy(R[0], R1, sizeof(double)*m*n);

            // transpose Q and R
            Q = transpose(m,m,1,Q); // frees Q
            R = transpose(m,n,1,R); // frees R

            G_list[i] = reshape_matrix_to_tensor(poly_order[i]*Grank[i+1], poly_order[i]*Grank[i+1], poly_order[i]*Grank[i+1], poly_order[i], Grank[i+1], Q); 
            
            // set size of R
            n_r = poly_order[i]*Grank[i+1];

            // Compute Matd[i]
            X_vandermonde = vandermonde_vec(poly_order[i], num_instances, i, X);
            tmp = dotkron(num_instances, poly_order[i], num_instances, Grank[i+1], 0, 0, X_vandermonde, Matd[i+1], NULL); // tmp = (num_instances x (poly_order[i]*Grank[i+1]) )
            Matd[i] = matmult_nat(num_instances, m, m, 1, tmp, transpose(m,m,1,Q)); // frees tmp and Q
        }
        else {// QR we have m >= n --> Q is (m x n) and R is (n x n)
            // convert Q1 and R1 in a 2d-allocation
            Q = d_malloc_2d(m,n); memcpy(Q[0], Q1, sizeof(double)*m*n);
            R = d_malloc_2d(n,n); memcpy(R[0], R1, sizeof(double)*n*n);

            // transpose Q and R
            Q = transpose(m,n,1,Q); // frees Q
            R = transpose(n,n,1,R); // frees R

            G_list[i] = reshape_matrix_to_tensor(Grank[i], poly_order[i]*Grank[i+1], Grank[i], poly_order[i], Grank[i+1], Q); 

            //set size of R
            n_r = Grank[i];

            // Compute Matd[i]
            X_vandermonde = vandermonde_vec(poly_order[i], num_instances, i, X);
            tmp = dotkron(num_instances, poly_order[i], num_instances, Grank[i+1], 0, 0, X_vandermonde, Matd[i+1], NULL); // tmp = (num_instances x (poly_order[i]*Grank[i+1]) )
            Matd[i] = matmult_nat(num_instances, n, m, 1, tmp, transpose(n,m,1,Q)); // frees tmp and Q
        }
        free(X_vandermonde); free(Q1); free(R1);
    }

    temp = reshape_tensor_to_matrix(Grank[0], poly_order[0], Grank[1], 0, G_list[0]); free(G_list[0]);
    temp = matmult_nat(Grank[0]*poly_order[0], n_r, Grank[1], 1, temp, R); // frees R and temp
    Grank[1] = n_r;
    G_list[0] = reshape_matrix_to_tensor(Grank[0]*poly_order[0], n_r, Grank[0], poly_order[0], Grank[1], temp); 
    Matd[0] = ones_matrix(num_instances, Grank[0]);
    

    /****************************/
    //         TRAINING
    /****************************/

    int ite=0; int itemax=(num_features-1)*4;
    int loopind, dir, ind;

    double** C, **C_transpose, **transpose_by_matrix;
    double* vec_G, *predictions, *transpose_by_labels;

    double* error = malloc(itemax*sizeof(double));

    while (ite < itemax ) {
        ite = ite + 1; 
        loopind = (ite-1) % (2*(num_features-1)); ++loopind;

        // dir indicates if we are in a right-to-left or left-to-right sweep
        if (loopind < num_features-0.5)
            dir = 1;
        else
            dir = 0;
        
        // ind indicates which TT-core index are we gonna optimize
        if (loopind <= num_features)
            ind = loopind;
        else
            ind = 2*num_features-loopind;
        ind -= 1;

        // Define and solve the linear system C_k vec(G_k) = y (equation 19 in the paper)
        X_vandermonde = vandermonde_vec(poly_order[ind], num_instances, ind, X);
        C = dotkron(num_instances, Grank[ind], num_instances, poly_order[ind], num_instances, Grank[ind+1], Matd[ind+1], X_vandermonde, Matd[ind+1]); // tmp = (num_instances x (Grank[ind]*poly_order[ind]*Grank[ind+1]) )
        C_transpose = transpose(num_instances, Grank[ind]*poly_order[ind]*Grank[ind+1], 0, C);
        
        transpose_by_matrix = matmult_nat(Grank[ind]*poly_order[ind]*Grank[ind+1],Grank[ind]*poly_order[ind]*Grank[ind+1],num_instances,0, C_transpose, C); 
        transpose_by_labels = matmult_vec(Grank[ind]*poly_order[ind]*Grank[ind+1], num_instances, y, C_transpose);  
        
        vec_G = (double *) malloc(Grank[ind]*poly_order[ind]*Grank[ind+1] * sizeof(double));
        vec_G = linear_system_solver_singular(Grank[ind]*poly_order[ind]*Grank[ind+1], transpose_by_labels, transpose_by_matrix);
        free(transpose_by_matrix); free(C_transpose); // freeing transpose_by_labels is freeing vec_G, as it is overwritten in linear_system_solver

        // Error
        predictions = matmult_vec(num_instances, Grank[ind]*poly_order[ind]*Grank[ind+1], vec_G, C); 
        error[ite] = mse(num_instances, predictions, y); free(predictions);
        printf("Error iteration %d: %f \n", ite, error[ite]);

        if (dir==1) {
            temp = reshape_vector_to_matrix(Grank[ind]*poly_order[ind], Grank[ind+1], vec_G); // frees vec_G
            qr(Grank[ind]*poly_order[ind], Grank[ind+1], temp, &Q1, &R1); free(temp);

            m = Grank[ind]*poly_order[ind]; n = Grank[ind+1];
            
            // Assign Q to the i-th TT-core: G_list[i]
            if (m < n) {// QR we have m < n --> Q is (m x m) and R is (m x n)
                // convert Q1 and R1 in a 2d-allocation
                Q = d_malloc_2d(m,m); memcpy(Q[0], Q1, sizeof(double)*m*m);
                R = d_malloc_2d(m,n); memcpy(R[0], R1, sizeof(double)*m*n);

                // include R to the next TT-core: G_list[ind+1]
                temp = reshape_tensor_to_matrix(Grank[ind+1], poly_order[ind+1], Grank[ind+2], 1, G_list[ind+1]); free(G_list[ind+1]);
                temp = matmult_nat(Grank[ind]*poly_order[ind], poly_order[ind+1]*Grank[ind+2], Grank[ind+1], 1, R, temp); // frees R and temp

                Grank[ind+1] = Grank[ind]*poly_order[ind]; //size(Q,2)
                G_list[ind+1] = reshape_matrix_to_tensor(Grank[ind+1], poly_order[ind+1]*Grank[ind+2], Grank[ind+1], poly_order[ind+1], Grank[ind+2], temp); free(temp);

                // assign Q to the current TT-core: G_list[ind]
                G_list[ind] = reshape_matrix_to_tensor(Grank[ind]*poly_order[ind], Grank[ind]*poly_order[ind], Grank[ind], poly_order[ind], Grank[ind+1], Q); 

                // upload Matd
                Matd[ind+1] = dotkron(num_instances, Grank[ind], num_instances, poly_order[ind], 0, 0, Matd[ind], X_vandermonde, NULL); // Matd[ind+1] = (num_instances x (Grank[ind]*poly_order[ind]) );
                Matd[ind+1] = matmult_nat(num_instances, Grank[ind+1], Grank[ind]*poly_order[ind], 1, Matd[ind+1], Q); // frees Matd[ind+1] and Q
            }

            else {// QR we have m >= n --> Q is (m x n) and R is (n x n)
                // convert Q1 and R1 in a 2d-allocation
                Q = d_malloc_2d(m,n); memcpy(Q[0], Q1, sizeof(double)*m*n);
                R = d_malloc_2d(n,n); memcpy(R[0], R1, sizeof(double)*n*n);

                // include R to the next TT-core: G_list[ind+1]
                temp = reshape_tensor_to_matrix(Grank[ind+1], poly_order[ind+1], Grank[ind+2], 1, G_list[ind+1]); free(G_list[ind+1]);
                temp = matmult_nat(Grank[ind+1], poly_order[ind+1]*Grank[ind+2], Grank[ind+1], 1, R, temp); // frees R and temp

                Grank[ind+1] = Grank[ind+1]; //size(Q,2)
                G_list[ind+1] = reshape_matrix_to_tensor(Grank[ind+1], poly_order[ind+1]*Grank[ind+2], Grank[ind+1], poly_order[ind+1], Grank[ind+2], temp); free(temp);

                // assign Q to the current TT-core: G_list[ind]
                G_list[ind] = reshape_matrix_to_tensor(Grank[ind]*poly_order[ind], Grank[ind+1], Grank[ind], poly_order[ind], Grank[ind+1], Q); 

                // upload Matd
                Matd[ind+1] = dotkron(num_instances, Grank[ind], num_instances, poly_order[ind], 0, 0, Matd[ind], X_vandermonde, NULL); // Matd[ind+1] = (num_instances x (Grank[ind]*poly_order[ind]) );
                Matd[ind+1] = matmult_nat(num_instances, Grank[ind+1], Grank[ind]*poly_order[ind], 1, Matd[ind+1], Q); // frees Matd[ind+1] and Q

            }

            
        }

        if (dir==0) {
            temp = reshape_vector_to_matrix(Grank[ind], poly_order[ind]*Grank[ind+1], vec_G); // frees vec_G
            temp = transpose(Grank[ind], poly_order[ind]*Grank[ind+1],1,temp); // frees temp
            qr(poly_order[ind]*Grank[ind+1], Grank[ind], temp, &Q1, &R1); free(temp);

            m = poly_order[ind]*Grank[ind+1]; n = Grank[ind];
            
            // Assign Q to the i-th TT-core: G_list[i]
            if (m < n) {// QR we have m < n --> Q is (m x m) and R is (m x n)
                // convert Q1 and R1 in a 2d-allocation
                Q = d_malloc_2d(m,m); memcpy(Q[0], Q1, sizeof(double)*m*m);
                R = d_malloc_2d(m,n); memcpy(R[0], R1, sizeof(double)*m*n);

                // transpose Q and R
                Q = transpose(m,m,1,Q); // frees Q
                R = transpose(m,n,1,R); // frees R

                // include R to the previous TT-core: G_list[ind-1]
                temp = reshape_tensor_to_matrix(Grank[ind-1], poly_order[ind-1], Grank[ind], 0, G_list[ind-1]); free(G_list[ind-1]);
                temp = matmult_nat(Grank[ind-1]*poly_order[ind-1], poly_order[ind]*Grank[ind+1], Grank[ind], 1, temp, R); // frees R and temp

                Grank[ind] = poly_order[ind]*Grank[ind+1]; //size(Q,1)
                G_list[ind-1] = reshape_matrix_to_tensor(Grank[ind-1]*poly_order[ind-1], poly_order[ind]*Grank[ind+1], Grank[ind-1], poly_order[ind-1], Grank[ind], temp); free(temp);

                // assign Q to the current TT-core: G_list[ind]
                G_list[ind] = reshape_matrix_to_tensor(Grank[ind], Grank[ind], Grank[ind], poly_order[ind], Grank[ind+1], Q); 

                // upload Matd
                Q = transpose(m,m,1,Q); // frees Q
                Matd[ind] = dotkron(num_instances, poly_order[ind], num_instances, Grank[ind+1], 0, 0, X_vandermonde, Matd[ind+1], NULL); // Matd[ind] = (num_instances x (poly_order[ind]*Grank[ind+1]) );
                Matd[ind] = matmult_nat(num_instances, Grank[ind], poly_order[ind]*Grank[ind+1], 1, Matd[ind], Q); // frees Matd[ind] and Q            
            }
            else {// QR we have m >= n --> Q is (m x n) and R is (n x n)
                // convert Q1 and R1 in a 2d-allocation
                Q = d_malloc_2d(m,n); memcpy(Q[0], Q1, sizeof(double)*m*n);
                R = d_malloc_2d(n,n); memcpy(R[0], R1, sizeof(double)*n*n);

                // transpose Q and R
                Q = transpose(m,n,1,Q); // frees Q
                R = transpose(n,n,1,R); // frees R

                // include R to the previous TT-core: G_list[ind-1]
                temp = reshape_tensor_to_matrix(Grank[ind-1], poly_order[ind-1], Grank[ind], 0, G_list[ind-1]); free(G_list[ind-1]);
                temp = matmult_nat(Grank[ind-1]*poly_order[ind-1], Grank[ind], Grank[ind], 1, temp, R); // frees R and temp

                Grank[ind] = Grank[ind]; //size(Q,1)
                G_list[ind-1] = reshape_matrix_to_tensor(Grank[ind-1]*poly_order[ind-1], Grank[ind], Grank[ind-1], poly_order[ind-1], Grank[ind], temp); free(temp);

                // assign Q to the current TT-core: G_list[ind]
                G_list[ind] = reshape_matrix_to_tensor(Grank[ind], poly_order[ind]*Grank[ind+1], Grank[ind], poly_order[ind], Grank[ind+1], Q); 

                // upload Matd
                Q = transpose(n,m,1,Q); // frees Q
                Matd[ind] = dotkron(num_instances, poly_order[ind], num_instances, Grank[ind+1], 0, 0, X_vandermonde, Matd[ind+1], NULL); // Matd[ind] = (num_instances x (poly_order[ind]*Grank[ind+1]) );
                Matd[ind] = matmult_nat(num_instances, Grank[ind], poly_order[ind]*Grank[ind+1], 1, Matd[ind], Q); // frees Matd[ind] and Q 
            }

        }
    }

    // Save weights
    char filename[256];
    for(int i=0; i<num_features; i++) {
        snprintf(filename, sizeof(filename), "Data/trained_G_p%d_r%d/G_%d.txt", poly_order[0], Grank[1], i);
        write_3Dtensor(filename, Grank[i], poly_order[i], Grank[i+1], G_list[i]);
    }
    // Save rank
    snprintf(filename, sizeof(filename), "Data/trained_G_p%d_r%d/rank.txt", poly_order[0], Grank[1]);
    write_vector_int(filename, num_features + 1, Grank);
    // Save poly_order
    snprintf(filename, sizeof(filename), "Data/trained_G_p%d_r%d/poly_order.txt", poly_order[0], Grank[1]);
    write_vector_int(filename, num_features, poly_order);




}
