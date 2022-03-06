sequential:
$ gcc -O3 -pthread convert_coo_to_csc.c masked_sparse_matrix_matrix_product.c mmio.c main.c
$ ./a.out
pthreads:
$ gcc -O4 -pthread -D NUM_OF_THREADS=16 convert_coo_to_csc.c masked_sparse_matrix_matrix_product.c mmio.c main.c
$ ./a.out
openMP:
$ gcc -O4 -fopenmp -D NUM_OF_THREADS=16 convert_coo_to_csc.c masked_sparse_matrix_matrix_product.c mmio.c main.c
$ ./a.out
openCilk:
$ clang -O3 -fcilkplus -pthread -D NUM_OF_THREADS=16 convert_coo_to_csc.c masked_sparse_matrix_matrix_product.c mmio.c main.c
$ ./a.out