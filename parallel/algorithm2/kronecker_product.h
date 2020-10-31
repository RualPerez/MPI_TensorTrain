double** kroneckerproduct(int m, int n, int p, int q, double** A, double** B);
double** repmat(int repetition_row, int repetition_col, int m, int n, double** A);
double** element_wise_mult(int m, int n, double** A, double** B);
double** ones_matrix(int m, int n);
double** dotkron(int m, int n, int p, int q, int u, int v, double** A, double** B, double** C);