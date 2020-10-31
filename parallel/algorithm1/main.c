#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "alloc.h"
#include "init.h"
#include "tt_contraction.h"
#include "vandermonde.h"
#include "read_write.h"
#include "read_write_int.h"
#include "operations.h"


int main(int argc, char *argv []) {
    /*
    The executeables takes the following input arguments:
    Example:
        Data/feat_200/poly_order_2.txt --> File with number vandermonde vectors and dimensions
        Data/feat_200/rank_2.txt       --> File with number of TT-ranks and each rank
        Data/X_20000_200.txt           --> File with number of instances, dimension and all data variables
        Data/y_20000.txt              --> File with number of instancesand all data outputs variables - ONLY USED WHEN CHECKING RESULTS WITH EXECUTABLE test_mse
    */

    // Initialize MPI
    MPI_Init(&argc,&argv);

    // Query size and rank
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // No buffering of stdout
    setbuf(stdout, NULL);

    // Checking number of input arguments
	#ifndef MSE
    	if (argc < 4 && rank == 0) printf("\nERROR! Too few input arguments!\n\n");
    #endif
    // MSE makes the 'test_mse' executable for testing correctness of the code
    #ifdef MSE
    	if (argc < 5 && rank == 0) printf("\nERROR! Too few input arguments!\n\n");
    #endif

    // Read hyperparameters (command-line inputs)
	int num_features = read_int(argv[1]);
	int* poly_order = read_vector_int(argv[1]);
	int* Grank = read_vector_int(argv[2]);

	#ifndef DATA
		// ALG1 makes the 'algorithm1' executable for timing the code
		#ifdef ALG1
    		// total timing
    		double start_t = MPI_Wtime();
		#endif
		// Read dataset
		int num_instances = read_int(argv[3]);
		double** X = d_malloc_2d(num_instances, num_features);
		X = d_malloc_2d(num_instances, num_features);
		#ifdef ALG1
    		// data load timing
    		double data_t0 = MPI_Wtime();
		#endif
		X = read_matrix(argv[3]);
		#ifdef ALG1
    		// Final data timing
    		double data_t1 = MPI_Wtime();
		#endif
	#endif

	// DATA makes the 'load_data' executable for timing the different data loading
	#ifdef DATA

		int num_instances = read_int(argv[3]);
		double** X = d_malloc_2d(num_instances, num_features);
	
    	// data load timing
    	double data1_t0 = MPI_Wtime();
		
		// In case data should only be loaded at process 0 and then distributed
		
		//int num_instances = read_int(argv[3]);
		//double** X = d_malloc_2d(num_instances, num_features);
		if (rank == 0){
			X = read_matrix(argv[3]);
    	}
    	MPI_Bcast(*X, num_instances*num_features, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    	
    	// Final data timing
    	double data1_t1 = MPI_Wtime();
	
    	free(X);
		
		// data load timing
    	double data2_t0 = MPI_Wtime();
	
    	// Read dataset
		//int num_instances = read_int(argv[3]);
		//double** X = d_malloc_2d(num_instances, num_features);
		X = d_malloc_2d(num_instances, num_features);
		X = read_matrix(argv[3]);
	
    	// Final data timing
    	double data2_t1 = MPI_Wtime();

    #endif

	
	/**************************************/
	// Allocate and initialize TT-cores (G)
	/**************************************/

	double **** G_list = malloc(num_features * sizeof(double ***) );
	char filename[256];
	for(int i=0; i<num_features; i++) {
		G_list[i] = d_malloc_3d(Grank[i], poly_order[i], Grank[i+1]);
		//ones_3D_tensor(general_rank, poly_order, general_rank, G_list[i]);
		
		// Create filename: ("G_i/G_%d.txt", i)
		snprintf(filename, sizeof(filename), "../../Data/algorithm1/G_%d.txt", num_features, poly_order[0], Grank[1], i);
		G_list[i] = read_3Dtensor(filename);
	}


	/*********************************************/
	// Read dataset and create vandermonde vectors
	/*********************************************/

	int vandermonde_size = 0;
	for(int i=0; i<num_features; i++) 
		vandermonde_size += poly_order[i];

    // convert to local X and local number of instances
    
    int ins_per_rank = num_instances / size;
/*    int num_lost_obs = num_instances % size;
    if (rank==0) {
        printf("Instances per rank: %d\n", ins_per_rank);
        printf("Instances not included: %d\n", num_lost_obs);
    }
*/

    int ins_start = rank*ins_per_rank;
    int ins_stop = (rank+1)*ins_per_rank;
    
    double** X_local = d_malloc_2d(ins_per_rank, num_features);
    
    int k = 0;
    for (int j = ins_start; j < ins_stop; j++) {
        for (int feat = 0; feat < num_features; feat++) {
            X_local[k][feat] = X[j][feat];
        }
        k++;
    }
    
    
    // create vandermonde vectors
	double** X_vandermonde = d_malloc_2d(ins_per_rank, vandermonde_size);
	vandermonde_vec(poly_order, ins_per_rank, num_features, X_local, X_vandermonde);
    
	free(X);
    free(X_local);

	
	/******************************/
	// Algorithm 1: TT contraction
	/******************************/

	double *** V_list = malloc(num_features * sizeof(double **) );
	double* vandermonde;

	int column_pos;
	double * first_V;
	double * last_V = malloc(Grank[num_features-1] * sizeof(double) );
	double * f;
	#ifdef MSE
		double * f_local = malloc(ins_per_rank*sizeof(double));
	#endif

	// ALG1 makes the 'algorithm1' executable for timing the code
	#ifdef ALG1
    	// start alg1 timing
    	double alg1_t0 = MPI_Wtime();
    #endif
    

	// Repeat algorithm 1 for each instance
	for (int ins=0; ins<ins_per_rank; ins++) {

		column_pos = 0;
		// Multiply vandermonde vector with its corresponding TT_core (G)
		for (int i=0; i<num_features; i++) {
			vandermonde = malloc(poly_order[i] * sizeof(double) );
			memcpy(vandermonde, &(X_vandermonde[ins][column_pos]), poly_order[i]  * sizeof(double));
			column_pos += poly_order[i];
			
			V_list[i] = d_malloc_2d(Grank[i], Grank[i+1]);
			tt_core_contraction(Grank[i], poly_order[i], Grank[i+1], vandermonde, V_list[i], G_list[i]);
		
			free(vandermonde);
		}

		// Reshape output from first and last TT-core
		first_V = malloc(Grank[1] * sizeof(double) );
		reshape_matrix_to_vector(Grank[0], Grank[1], first_V, V_list[0]);
		reshape_matrix_to_vector(Grank[num_features-1], Grank[num_features], last_V, V_list[num_features-1]);
		free(V_list[0]); free(V_list[num_features-1]);

		// Compress multiplications of V
		f = first_V;
		for (int i=1; i<num_features-1; i++) {
			vector_matrix(Grank[i], Grank[i+1], &f, V_list[i]);
			free(V_list[i]);
		}
		
		vector_vector(Grank[num_features-1], &f, last_V);

		#ifdef MSE
			// Print result
/*			if (ins == 0)
				printf("For 1st instance of rank %2d, regression prediction %.1f\n", rank, 	*f);
*/			f_local[ins] = *f;
		#endif
		free(f);
	}

	#ifdef ALG1
		// Final alg1 timing
    	double alg1_t1 = MPI_Wtime();
    	double end_t = MPI_Wtime();
    #endif

    // Makes sure all ranks has reached this point
    MPI_Barrier(MPI_COMM_WORLD);
    
    // Timings, min/avg/max
    #ifdef DATA
    	double data1_t = data1_t1 - data1_t0;
    	double data1_t_avg;
    	MPI_Reduce(&data1_t, &data1_t0, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    	MPI_Reduce(&data1_t, &data1_t1, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    	MPI_Reduce(&data1_t, &data1_t_avg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    	data1_t_avg /= (double) size;
    	
    	double data2_t = data2_t1 - data2_t0;
    	double data2_t_avg;
    	MPI_Reduce(&data2_t, &data2_t0, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    	MPI_Reduce(&data2_t, &data2_t1, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    	MPI_Reduce(&data2_t, &data2_t_avg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    	data2_t_avg /= (double) size;
    #endif

    #ifdef ALG1
    	double data_t = data_t1 - data_t0;
    	double data_t_avg;
    	MPI_Reduce(&data_t, &data_t0, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    	MPI_Reduce(&data_t, &data_t1, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    	MPI_Reduce(&data_t, &data_t_avg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    	data_t_avg /= (double) size;
	
    	double alg1_t = alg1_t1 - alg1_t0;
    	double alg1_t_avg;
    	MPI_Reduce(&alg1_t, &alg1_t0, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    	MPI_Reduce(&alg1_t, &alg1_t1, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    	MPI_Reduce(&alg1_t, &alg1_t_avg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    	alg1_t_avg /= (double) size;
	
    	double total_t = end_t - start_t;
    	double avg_t;
    	MPI_Reduce(&total_t, &start_t, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    	MPI_Reduce(&total_t, &end_t, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    	MPI_Reduce(&total_t, &avg_t, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    	avg_t /= (double) size;
    #endif

	// Print inputs
	if ( rank == 0 ) {
    	printf("Arguments:\n");
    	printf("  Number of instances = %d\n", num_instances);
    	printf("  Number of features = %d\n", num_features);
    	printf("  Poly orders =");
    	for (int i=0; i<num_features; i++) printf(" %d", poly_order[i]);
    	printf("\n  TT-ranks =");
    	for (int i=0; i<num_features+1; i++) printf(" %d", Grank[i]);
    	printf("\n  MPI processors: %d\n", size);
	}

	// Print timings
    if (rank == 0) {
    	#ifdef DATA
        	printf("\nTime b-cast loading data min/avg/max  %12.8f / %12.8f / %12.8f\n\n", 	data1_t0, data1_t_avg, data1_t1);
        	printf("Time all loading data min/avg/max  %12.8f / %12.8f / %12.8f\n\n", data2_t0	, data2_t_avg, data2_t1);
        	fflush(stdout);
        #endif
        #ifdef ALG1
        	printf("\nTime to load data min/avg/max  %12.8f / %12.8f / %12.8f\n\n", data_t0	, data_t_avg, data_t1);
        	printf("Time algorithm 1 min/avg/max  %12.8f / %12.8f / %12.8f\n\n", alg1_t0, 	alg1_t_avg, alg1_t1);
        	printf("Total time min/avg/max  %12.8f / %12.8f / %12.8f\n\n", start_t, avg_t, 	end_t);
        	fflush(stdout);
        #endif
    }

	#ifdef MSE
    	double *y = malloc(ins_per_rank*size*sizeof(double));
    	MPI_Gather(f_local, ins_per_rank, MPI_DOUBLE, y, ins_per_rank, MPI_DOUBLE, 0, 	MPI_COMM_WORLD);
    	if (rank == 0){
    		double *y_test = read_vector(argv[4]);
    		double error = mse(ins_per_rank*size, y, y_test);
    		printf("\nMSE = %lf\n\n", error);
		}
		free(y);
		free(f_local);
	#endif
    
	// Free
	free(last_V);
	free(X_vandermonde); free(poly_order); free(Grank);
	for (int i=0; i<num_features; i++) 
		free(G_list[i]);
	free(G_list);
	
    MPI_Finalize();
    return 0;

}
