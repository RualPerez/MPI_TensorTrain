void tt_core_contraction(int rank_k_1, int n_k, int rank_k, double* v, double** V, double*** G);
void reshape_matrix_to_vector(int rank_k_1, int rank_k, double* v, double** V);
void vector_matrix(int rank_k, int rank_k1, double** f, double** V);
void vector_vector(int rank_k, double** f, double* V);