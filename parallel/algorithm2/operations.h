double** matmult_nat(int m, int n, int k, int free_memory, double **A, double **B);
double* matmult_vec(int m, int n, double* v, double** A);
double** transpose(int m, int n, int free_memory, double** A);
double mse(int m, double* predictions, double* y);