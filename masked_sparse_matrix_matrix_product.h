#ifndef MASKED_SPARSE_MATRIX_MATRIX_PRODUCT
#define MASKED_SPARSE_MATRIX_MATRIX_PRODUCT

void sequential_masked_sparse_matrix_matrix_product(
    int * a_lower_row, int * a_lower_col,
    int * a_csc_row, int * a_csc_col,
    int   a_nnz,     int   a_N,
    int * c_csc_row, int * c_csc_col, int * c_csc_val
);



void openCilk_masked_sparse_matrix_matrix_product(
    int * a_lower_row, int * a_lower_col,
    int * a_csc_row, int * a_csc_col,
    int   a_nnz,     int   a_N,
    int * c_csc_row, int * c_csc_col, int * c_csc_val,
    int   num_of_threads
);

void openMP_masked_sparse_matrix_matrix_product(
    int * a_lower_row, int * a_lower_col,
    int * a_csc_row, int * a_csc_col,
    int   a_nnz,     int   a_N,
    int * c_csc_row, int * c_csc_col, int * c_csc_val
);

void boundary_sequential_masked_sparse_matrix_matrix_product(void *arg);

void pThreads_masked_sparse_matrix_matrix_product(
    int * a_lower_row,  int * a_lower_col,
    int * a_csc_row,    int * a_csc_col,
    int   a_nnz,        int   a_N,
    int * c_csc_row,    int * c_csc_col,    int * c_csc_val,
    int num_of_threads
);

#endif
