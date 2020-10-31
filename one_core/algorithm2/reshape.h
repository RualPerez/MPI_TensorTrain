double** reshape_tensor_to_matrix(int rank_k_1, int poly_order, int rank_k, int type, double*** tn);
double** reshape_matrix_to_matrix(int m, int n, int m_f, int n_f, double** A);
double*** reshape_matrix_to_tensor(int m, int n, int m_f, int n_f, int k_f, double** A);
double** reshape_vector_to_matrix(int m, int n, double* v);