double* read_vector(char* filename);
void write_vector(char* filename, int n, double* vector);
double** read_matrix(char* filename);
void write_matrix(char* filename, int m, int n, double** matrix);
double*** read_3Dtensor(char* filename);
void write_3Dtensor(char* filename, int m, int n, int k, double*** tensor);