#include <stdio.h>
#include <stdlib.h>
#include "convert_coo_to_csc.h"
#include <cilk/cilk.h>
#include <pthread.h>
#include <sys/types.h>
#include <omp.h>

// IMPORTANT

/* if you want to run the open cilk function then comment the #include <opp.h> and the whole openMP function
    if you want to run any other function then comment the #include <cilk/cilk.h> and the whole openCilk function */


pthread_mutex_t lock; //define the lock and use when needed in the parallel functions

void sequential_masked_sparse_matrix_matrix_product(
    int * a_lower_row, int * a_lower_col,
    int * a_csc_row, int * a_csc_col,
    int   a_nnz,     int   a_N,
    int * c_csc_row, int * c_csc_col, int * c_csc_val
) {
    /* A is symmetric so C is symmetric, so we only need to compute the lower or upper triangular
        matrix of C and then convert it to the whole matrix like we did with our A matrix in the main.c */

    int * c_temp_lower_row = calloc(a_nnz/2, sizeof(int));
    int * c_temp_lower_col = calloc(a_nnz/2, sizeof(int));
    int * c_temp_lower_val = calloc(a_nnz/2, sizeof(int));

    int index = 0;
    int c_lower_nnz = 0;

    for (int i=0; i < a_N; ++i) {

        /* start by taking each column of A matrix*/
        int a_col_start = a_lower_col[i]; // this is where it starts
        int a_col_end   = a_lower_col[i+1]; // and this where it ends

        for (int j=a_col_start; j < a_col_end; ++j) { // now take every non zero element of each column

            /* we call a1 the first A matrix of the A*A product and a2 the second A*/

            int a1_col = i; // so this is each column of the a1 matrix because i goes from 0 to N
            int a2_col = a_lower_row[j]; /* a2_col is every row (we call it col in the name because
                                            A is symmetric so rows and columns are the same thing )
                                            of the a2 matrix THAT corresponds to the
                                            non zero elements of the a1_col column. So now we are computing
                                            only the indices of the C matrix that have non zero values in
                                            the A matrix. In other words we accomplice Hadamard product*/

            int a1_col_start = a_csc_col[a1_col];
            int a1_col_end = a_csc_col[a1_col+1];

            int a2_col_start = a_csc_col[a2_col];
            int a2_col_end = a_csc_col[a2_col+1];

            int value = 0;

            int k = a1_col_start;
            int l = a2_col_start;

            /*  see how many times our a1_col and a2_col have a nonzero value in the same place
                the number of times is the value of the element (a2_col, a1_col) of the C matrix */

            while ( (k<a1_col_end) && (l<a2_col_end) ) {

                if (a_csc_row[k] == a_csc_row[l]) {
                    ++value;
                    ++k;
                    ++l;
                    continue;
                }

                else if (a_csc_row[k] > a_csc_row[l]) {
                    ++l;
                    continue;
                }

                else {
                    ++k;
                    continue;
                }

            }

            if (value > 0) {
                c_temp_lower_row[index] = a1_col;
                c_temp_lower_col[index] = a2_col;
                c_temp_lower_val[index] = value;
                ++index;
                ++c_lower_nnz;
            }

        }
    }

    /* now we have the lower triangular of the C matrix in COO
        we convert it to the CSC form of the whole C matrix as we did for A in main.c
        so there is nothing new here */

    int c_nnz = 2 * c_lower_nnz;

    printf("There are %d non zero values. \n", c_nnz);

    int * c_lower_row = calloc(c_lower_nnz, sizeof(int));
    int * c_lower_col = calloc(c_lower_nnz, sizeof(int));
    int * c_lower_val = calloc(c_lower_nnz, sizeof(int));

    for (int i=0; i < c_lower_nnz; ++i) {
        c_lower_row[i] = c_temp_lower_row[i];
        c_lower_col[i] = c_temp_lower_col[i];
        c_lower_val[i] = c_temp_lower_val[i];
    }

    free(c_temp_lower_row);
    free(c_temp_lower_col);
    free(c_temp_lower_val);

    int * c_upper_row = calloc(c_lower_nnz, sizeof(int));
    int * c_upper_col = calloc(c_lower_nnz, sizeof(int));
    int * c_upper_val = calloc(c_lower_nnz, sizeof(int));

    for (int i=0; i < c_lower_nnz; ++i) {
        c_upper_row[i] = c_lower_col[i];
        c_upper_col[i] = c_lower_row[i];
        c_upper_val[i] = c_lower_val[i];
    }



    int * c_coo_row = calloc(c_nnz, sizeof(int));
    int * c_coo_col = calloc(c_nnz, sizeof(int));
    int * c_coo_val = calloc(c_nnz, sizeof(int));

    for (int i=0; i < c_nnz; ++i) {

        if (i < c_lower_nnz) {
            c_coo_row[i] = c_lower_row[i];
            c_coo_col[i] = c_lower_col[i];
            c_coo_val[i] = c_lower_val[i];
        }

        else {
            c_coo_row[i] = c_lower_row[i-c_lower_nnz];
            c_coo_col[i] = c_lower_col[i-c_lower_nnz];
            c_coo_val[i] = c_lower_val[i-c_lower_nnz];
        }
    }



    c_csc_row = calloc(c_nnz, sizeof(int));
    c_csc_col = calloc( (a_N + 1), sizeof(int) );
    c_csc_val = calloc(c_nnz, sizeof(int));

    coo2csc_col(c_csc_col, c_coo_col, c_nnz, a_N);

    int * c_upper_csc_row = calloc(c_lower_nnz, sizeof(int));
    int * c_upper_csc_col = calloc(c_lower_nnz, sizeof(int));
    int * c_upper_csc_val = calloc(c_lower_nnz, sizeof(int));

    coo2csc(c_upper_csc_row, c_upper_csc_col, c_upper_csc_val,
            c_upper_row, c_upper_col, c_upper_val,
            c_lower_nnz, a_N);

    int * c_lower_csc_row = calloc(c_lower_nnz, sizeof(int));
    int * c_lower_csc_col = calloc(c_lower_nnz, sizeof(int));
    int * c_lower_csc_val = calloc(c_lower_nnz, sizeof(int));

    coo2csc(c_lower_csc_row, c_lower_csc_col, c_lower_csc_val,
            c_lower_row, c_lower_col, c_lower_val,
            c_lower_nnz, a_N);

    coo2csc_row_val(c_upper_csc_row, c_upper_csc_col, c_upper_csc_val,
                   c_lower_csc_row, c_lower_csc_col, c_lower_csc_val,
                   c_csc_row, c_csc_val, a_N);


    /* compute the number of triangles
        just get the sum of the values of all elements of the C matrix
        and then divide them by 6 */

    int counter = 0;

    for (int i=0; i < c_nnz; ++i)
        counter = counter + c_csc_val[i];

    int num_of_triangles = counter / 6;

    printf("there are %d triangles.\n", num_of_triangles);


    free(c_lower_row);
    free(c_lower_col);
    free(c_lower_val);

    free(c_upper_row);
    free(c_upper_col);
    free(c_upper_val);

    free(c_coo_row);
    free(c_coo_col);
    free(c_coo_val);

    free(c_upper_csc_row);
    free(c_upper_csc_col);
    free(c_upper_csc_val);

    free(c_lower_csc_row);
    free(c_lower_csc_col);
    free(c_lower_csc_val);

}

/* same thing with sequential, but instead of for loops we have cilk_for loops
    and we add pragma grainsize to change the number of threads */

void openCilk_masked_sparse_matrix_matrix_product(
    int * a_lower_row, int * a_lower_col,
    int * a_csc_row, int * a_csc_col,
    int   a_nnz,     int   a_N,
    int * c_csc_row, int * c_csc_col, int * c_csc_val,
    int   num_of_threads
) {
    int * c_temp_lower_row = calloc(a_nnz/2, sizeof(int));
    int * c_temp_lower_col = calloc(a_nnz/2, sizeof(int));
    int * c_temp_lower_val = calloc(a_nnz/2, sizeof(int));

    int index = 0;
    int c_lower_nnz = 0;



    pthread_mutex_init(&lock, NULL); //initialize the lock


    const int gs_1 = (2048 < a_N / (8 * num_of_threads)) ? 2048 : a_N / (8 * num_of_threads);

    #pragma grainsize gs_1
    cilk_for (int i=0; i < a_N; ++i) {

        int a_col_start = a_lower_col[i];
        int a_col_end   = a_lower_col[i+1];


        const int gs_2 =  (2048 <(a_col_end - a_col_start) / (8 * num_of_threads)) ? 2048 : (a_col_end - a_col_start) / (8 * num_of_threads);
        #pragma grainsize gs_2
        cilk_for (int j=a_col_start; j < a_col_end; ++j) {


            int a1_col = i;
            int a2_col = a_lower_row[j];

            int a1_col_start = a_csc_col[a1_col];
            int a1_col_end = a_csc_col[a1_col+1];

            int a2_col_start = a_csc_col[a2_col];
            int a2_col_end = a_csc_col[a2_col+1];

            int value = 0;

            int k = a1_col_start;
            int l = a2_col_start;

            while ( (k<a1_col_end) && (l<a2_col_end) ) {

                if (a_csc_row[k] == a_csc_row[l]) {
                    ++value;
                    ++k;
                    ++l;
                    continue;
                }

                else if (a_csc_row[k] > a_csc_row[l]) {
                    ++l;
                    continue;
                }

                else {
                    ++k;
                    continue;
                }

            }

            if (value > 0) {

                pthread_mutex_lock(&lock);

                c_temp_lower_row[index] = a1_col;
                c_temp_lower_col[index] = a2_col;
                c_temp_lower_val[index] = value;
                ++index;
                ++c_lower_nnz;

                pthread_mutex_unlock(&lock);
            }

        }
    }

    pthread_mutex_destroy (&lock);

    int c_nnz = 2 * c_lower_nnz;

    //printf("There are %d non zero values. \n", c_nnz);

    int * c_lower_row = calloc(c_lower_nnz, sizeof(int));
    int * c_lower_col = calloc(c_lower_nnz, sizeof(int));
    int * c_lower_val = calloc(c_lower_nnz, sizeof(int));

    for (int i=0; i < c_lower_nnz; ++i) {
        c_lower_row[i] = c_temp_lower_row[i];
        c_lower_col[i] = c_temp_lower_col[i];
        c_lower_val[i] = c_temp_lower_val[i];
    }

    free(c_temp_lower_row);
    free(c_temp_lower_col);
    free(c_temp_lower_val);

    int * c_upper_row = calloc(c_lower_nnz, sizeof(int));
    int * c_upper_col = calloc(c_lower_nnz, sizeof(int));
    int * c_upper_val = calloc(c_lower_nnz, sizeof(int));

    for (int i=0; i < c_lower_nnz; ++i) {
        c_upper_row[i] = c_lower_col[i];
        c_upper_col[i] = c_lower_row[i];
        c_upper_val[i] = c_lower_val[i];
    }



    int * c_coo_row = calloc(c_nnz, sizeof(int));
    int * c_coo_col = calloc(c_nnz, sizeof(int));
    int * c_coo_val = calloc(c_nnz, sizeof(int));

    for (int i=0; i < c_nnz; ++i) {

        if (i < c_lower_nnz) {
            c_coo_row[i] = c_lower_row[i];
            c_coo_col[i] = c_lower_col[i];
            c_coo_val[i] = c_lower_val[i];
        }

        else {
            c_coo_row[i] = c_lower_row[i-c_lower_nnz];
            c_coo_col[i] = c_lower_col[i-c_lower_nnz];
            c_coo_val[i] = c_lower_val[i-c_lower_nnz];
        }
    }



    c_csc_row = calloc(c_nnz, sizeof(int));
    c_csc_col = calloc( (a_N + 1), sizeof(int) );
    c_csc_val = calloc(c_nnz, sizeof(int));

    coo2csc_col(c_csc_col, c_coo_col, c_nnz, a_N);

    int * c_upper_csc_row = calloc(c_lower_nnz, sizeof(int));
    int * c_upper_csc_col = calloc(c_lower_nnz, sizeof(int));
    int * c_upper_csc_val = calloc(c_lower_nnz, sizeof(int));

    coo2csc(c_upper_csc_row, c_upper_csc_col, c_upper_csc_val,
            c_upper_row, c_upper_col, c_upper_val,
            c_lower_nnz, a_N);

    int * c_lower_csc_row = calloc(c_lower_nnz, sizeof(int));
    int * c_lower_csc_col = calloc(c_lower_nnz, sizeof(int));
    int * c_lower_csc_val = calloc(c_lower_nnz, sizeof(int));

    coo2csc(c_lower_csc_row, c_lower_csc_col, c_lower_csc_val,
            c_lower_row, c_lower_col, c_lower_val,
            c_lower_nnz, a_N);

    coo2csc_row_val(c_upper_csc_row, c_upper_csc_col, c_upper_csc_val,
                   c_lower_csc_row, c_lower_csc_col, c_lower_csc_val,
                   c_csc_row, c_csc_val, a_N);

    int counter = 0;

    for (int i=0; i < c_nnz; ++i)
        counter = counter + c_csc_val[i];

    int num_of_triangles = counter / 6;

    //printf("there are %d triangles.\n", num_of_triangles);


    free(c_lower_row);
    free(c_lower_col);
    free(c_lower_val);

    free(c_upper_row);
    free(c_upper_col);
    free(c_upper_val);

    free(c_coo_row);
    free(c_coo_col);
    free(c_coo_val);

    free(c_upper_csc_row);
    free(c_upper_csc_col);
    free(c_upper_csc_val);

    free(c_lower_csc_row);
    free(c_lower_csc_col);
    free(c_lower_csc_val);

}


/* we just use #pragma omp parallel for ( and omp_set_num_threads(NUM_OF_THREADS) in main.c )*/

void openMP_masked_sparse_matrix_matrix_product(
    int * a_lower_row, int * a_lower_col,
    int * a_csc_row, int * a_csc_col,
    int   a_nnz,     int   a_N,
    int * c_csc_row, int * c_csc_col, int * c_csc_val
) {
    int * c_temp_lower_row = calloc(a_nnz/2, sizeof(int));
    int * c_temp_lower_col = calloc(a_nnz/2, sizeof(int));
    int * c_temp_lower_val = calloc(a_nnz/2, sizeof(int));

    int index = 0;
    int c_lower_nnz = 0;


    pthread_mutex_init(&lock, NULL); //initialize the lock

    #pragma omp parallel for
    for (int i=0; i < a_N; ++i) {

        int a_col_start = a_lower_col[i];
        int a_col_end   = a_lower_col[i+1];

        //#pragma omp parallel for
        for (int j=a_col_start; j < a_col_end; ++j) {


            int a1_col = i;
            int a2_col = a_lower_row[j];

            int a1_col_start = a_csc_col[a1_col];
            int a1_col_end = a_csc_col[a1_col+1];

            int a2_col_start = a_csc_col[a2_col];
            int a2_col_end = a_csc_col[a2_col+1];

            int value = 0;

            int k = a1_col_start;
            int l = a2_col_start;

            while ( (k<a1_col_end) && (l<a2_col_end) ) {

                if (a_csc_row[k] == a_csc_row[l]) {
                    ++value;
                    ++k;
                    ++l;
                    continue;
                }

                else if (a_csc_row[k] > a_csc_row[l]) {
                    ++l;
                    continue;
                }

                else {
                    ++k;
                    continue;
                }

            }

            if (value > 0) {

                pthread_mutex_lock(&lock);

                c_temp_lower_row[index] = a1_col;
                c_temp_lower_col[index] = a2_col;
                c_temp_lower_val[index] = value;
                ++index;
                ++c_lower_nnz;

                pthread_mutex_unlock(&lock);
            }

        }
    }

    pthread_mutex_destroy (&lock);

    int c_nnz = 2 * c_lower_nnz;

    //printf("There are %d non zero values. \n", c_nnz);

    int * c_lower_row = calloc(c_lower_nnz, sizeof(int));
    int * c_lower_col = calloc(c_lower_nnz, sizeof(int));
    int * c_lower_val = calloc(c_lower_nnz, sizeof(int));

    for (int i=0; i < c_lower_nnz; ++i) {
        c_lower_row[i] = c_temp_lower_row[i];
        c_lower_col[i] = c_temp_lower_col[i];
        c_lower_val[i] = c_temp_lower_val[i];
    }

    free(c_temp_lower_row);
    free(c_temp_lower_col);
    free(c_temp_lower_val);

    int * c_upper_row = calloc(c_lower_nnz, sizeof(int));
    int * c_upper_col = calloc(c_lower_nnz, sizeof(int));
    int * c_upper_val = calloc(c_lower_nnz, sizeof(int));

    for (int i=0; i < c_lower_nnz; ++i) {
        c_upper_row[i] = c_lower_col[i];
        c_upper_col[i] = c_lower_row[i];
        c_upper_val[i] = c_lower_val[i];
    }



    int * c_coo_row = calloc(c_nnz, sizeof(int));
    int * c_coo_col = calloc(c_nnz, sizeof(int));
    int * c_coo_val = calloc(c_nnz, sizeof(int));

    for (int i=0; i < c_nnz; ++i) {

        if (i < c_lower_nnz) {
            c_coo_row[i] = c_lower_row[i];
            c_coo_col[i] = c_lower_col[i];
            c_coo_val[i] = c_lower_val[i];
        }

        else {
            c_coo_row[i] = c_lower_row[i-c_lower_nnz];
            c_coo_col[i] = c_lower_col[i-c_lower_nnz];
            c_coo_val[i] = c_lower_val[i-c_lower_nnz];
        }
    }



    c_csc_row = calloc(c_nnz, sizeof(int));
    c_csc_col = calloc( (a_N + 1), sizeof(int) );
    c_csc_val = calloc(c_nnz, sizeof(int));

    coo2csc_col(c_csc_col, c_coo_col, c_nnz, a_N);

    int * c_upper_csc_row = calloc(c_lower_nnz, sizeof(int));
    int * c_upper_csc_col = calloc(c_lower_nnz, sizeof(int));
    int * c_upper_csc_val = calloc(c_lower_nnz, sizeof(int));

    coo2csc(c_upper_csc_row, c_upper_csc_col, c_upper_csc_val,
            c_upper_row, c_upper_col, c_upper_val,
            c_lower_nnz, a_N);

    int * c_lower_csc_row = calloc(c_lower_nnz, sizeof(int));
    int * c_lower_csc_col = calloc(c_lower_nnz, sizeof(int));
    int * c_lower_csc_val = calloc(c_lower_nnz, sizeof(int));

    coo2csc(c_lower_csc_row, c_lower_csc_col, c_lower_csc_val,
            c_lower_row, c_lower_col, c_lower_val,
            c_lower_nnz, a_N);

    coo2csc_row_val(c_upper_csc_row, c_upper_csc_col, c_upper_csc_val,
                   c_lower_csc_row, c_lower_csc_col, c_lower_csc_val,
                   c_csc_row, c_csc_val, a_N);

    int counter = 0;

    for (int i=0; i < c_nnz; ++i)
        counter = counter + c_csc_val[i];

    int num_of_triangles = counter / 6;

    //printf("there are %d triangles.\n", num_of_triangles);


    free(c_lower_row);
    free(c_lower_col);
    free(c_lower_val);

    free(c_upper_row);
    free(c_upper_col);
    free(c_upper_val);

    free(c_coo_row);
    free(c_coo_col);
    free(c_coo_val);

    free(c_upper_csc_row);
    free(c_upper_csc_col);
    free(c_upper_csc_val);

    free(c_lower_csc_row);
    free(c_lower_csc_col);
    free(c_lower_csc_val);

}


/* this struct has the input arguments for the next function */

typedef struct {
    int * a_lower_row,      * a_lower_col;
    int * a_csc_row,        * a_csc_col;
    int * c_temp_lower_row, * c_temp_lower_col, * c_temp_lower_val;
    int * index,            * c_lower_nnz;
    int   c_col_start,        c_col_end;
} parm;

/* this is like the sequential but it computes the lower triangular C matrix
    from the column c_col_start to the c_col_end */

void *boundary_sequential_masked_sparse_matrix_matrix_product(void *arg) {

    parm *p = (parm *)arg;




    for (int i=p->c_col_start; i < p->c_col_end; ++i) {

        int a_col_start = p->a_lower_col[i];
        int a_col_end   = p->a_lower_col[i+1];

        for (int j=a_col_start; j < a_col_end; ++j) {


            int a1_col = i;
            int a2_col = p->a_lower_row[j];

            int a1_col_start = p->a_csc_col[a1_col];
            int a1_col_end = p->a_csc_col[a1_col+1];

            int a2_col_start = p->a_csc_col[a2_col];
            int a2_col_end = p->a_csc_col[a2_col+1];

            int value = 0;

            int k = a1_col_start;
            int l = a2_col_start;

            while ( (k<a1_col_end) && (l<a2_col_end) ) {

                if (p->a_csc_row[k] == p->a_csc_row[l]) {
                    ++value;
                    ++k;
                    ++l;
                    continue;
                }

                else if (p->a_csc_row[k] > p->a_csc_row[l]) {
                    ++l;
                    continue;
                }

                else {
                    ++k;
                    continue;
                }

            }

            if (value > 0) {

                pthread_mutex_lock(&lock);

                p->c_temp_lower_row[p->index[0]] = a1_col;
                p->c_temp_lower_col[p->index[0]] = a2_col;
                p->c_temp_lower_val[p->index[0]] = value;
                p->index[0] = p->index[0] + 1 ;
                p->c_lower_nnz[0] = p->c_lower_nnz[0] + 1 ;

                pthread_mutex_unlock(&lock);
            }

        }
    }


    pthread_exit(NULL);
}


/* this function just sets every thread to compute the same portion of the C matrix
    using the previews function. So if we have num_of_threads threads, then each thread
    computes 1/num_of_threads part of the C matrix.
    this function works exactly like the hello.c example we had for the pthreads */

void pThreads_masked_sparse_matrix_matrix_product(
    int * a_lower_row,  int * a_lower_col,
    int * a_csc_row,    int * a_csc_col,
    int   a_nnz,        int   a_N,
    int * c_csc_row,    int * c_csc_col,    int * c_csc_val,
    int num_of_threads
) {

    int * c_temp_lower_row = calloc(a_nnz/2, sizeof(int));
    int * c_temp_lower_col = calloc(a_nnz/2, sizeof(int));
    int * c_temp_lower_val = calloc(a_nnz/2, sizeof(int));

    int index = 0;
    int c_lower_nnz = 0;

    pthread_t *threads;
    pthread_attr_t pthread_custom_attr;
    parm *p;

    threads = (pthread_t *) malloc(num_of_threads * sizeof(pthread_t));
    pthread_attr_init( &pthread_custom_attr );

    p = (parm *) malloc(num_of_threads * sizeof(parm));

    int iterations = a_N /num_of_threads;

    for (int i=0; i < num_of_threads; ++i) {
        p[i].a_lower_row = a_lower_row;
        p[i].a_lower_col = a_lower_col;
        p[i].a_csc_row = a_csc_row;
        p[i].a_csc_col = a_csc_col;
        p[i].c_temp_lower_row = c_temp_lower_row;
        p[i].c_temp_lower_col = c_temp_lower_col;
        p[i].c_temp_lower_val = c_temp_lower_val;
        p[i].index = &index;
        p[i].c_lower_nnz = &c_lower_nnz;

        if (i == num_of_threads-1) {
            p[i].c_col_start = i * iterations;
            p[i].c_col_end = a_N;
        }

        else {
            p[i].c_col_start = i * iterations;
            p[i].c_col_end = (i * iterations) + iterations;
        }
    }


    pthread_mutex_init(&lock, NULL); //initialize the lock

    for (int i=0; i < num_of_threads; ++i) {
        pthread_create(&threads[i], &pthread_custom_attr,
                       boundary_sequential_masked_sparse_matrix_matrix_product, (void *)(p+i));
    }

    for (int i=0; i < num_of_threads; ++i) {
        pthread_join(threads[i],NULL);
    }

    free(p);

    pthread_mutex_destroy (&lock);


    // now we have the lower triangular C matrix and the rest is not new

    int c_nnz = 2 * c_lower_nnz;

    //printf("There are %d non zero values. \n", c_nnz);

    int * c_lower_row = calloc(c_lower_nnz, sizeof(int));
    int * c_lower_col = calloc(c_lower_nnz, sizeof(int));
    int * c_lower_val = calloc(c_lower_nnz, sizeof(int));

    for (int i=0; i < c_lower_nnz; ++i) {
        c_lower_row[i] = c_temp_lower_row[i];
        c_lower_col[i] = c_temp_lower_col[i];
        c_lower_val[i] = c_temp_lower_val[i];
    }

    free(c_temp_lower_row);
    free(c_temp_lower_col);
    free(c_temp_lower_val);

    int * c_upper_row = calloc(c_lower_nnz, sizeof(int));
    int * c_upper_col = calloc(c_lower_nnz, sizeof(int));
    int * c_upper_val = calloc(c_lower_nnz, sizeof(int));

    for (int i=0; i < c_lower_nnz; ++i) {
        c_upper_row[i] = c_lower_col[i];
        c_upper_col[i] = c_lower_row[i];
        c_upper_val[i] = c_lower_val[i];
    }



    int * c_coo_row = calloc(c_nnz, sizeof(int));
    int * c_coo_col = calloc(c_nnz, sizeof(int));
    int * c_coo_val = calloc(c_nnz, sizeof(int));

    for (int i=0; i < c_nnz; ++i) {

        if (i < c_lower_nnz) {
            c_coo_row[i] = c_lower_row[i];
            c_coo_col[i] = c_lower_col[i];
            c_coo_val[i] = c_lower_val[i];
        }

        else {
            c_coo_row[i] = c_lower_row[i-c_lower_nnz];
            c_coo_col[i] = c_lower_col[i-c_lower_nnz];
            c_coo_val[i] = c_lower_val[i-c_lower_nnz];
        }
    }



    c_csc_row = calloc(c_nnz, sizeof(int));
    c_csc_col = calloc( (a_N + 1), sizeof(int) );
    c_csc_val = calloc(c_nnz, sizeof(int));

    coo2csc_col(c_csc_col, c_coo_col, c_nnz, a_N);

    int * c_upper_csc_row = calloc(c_lower_nnz, sizeof(int));
    int * c_upper_csc_col = calloc(c_lower_nnz, sizeof(int));
    int * c_upper_csc_val = calloc(c_lower_nnz, sizeof(int));

    coo2csc(c_upper_csc_row, c_upper_csc_col, c_upper_csc_val,
            c_upper_row, c_upper_col, c_upper_val,
            c_lower_nnz, a_N);

    int * c_lower_csc_row = calloc(c_lower_nnz, sizeof(int));
    int * c_lower_csc_col = calloc(c_lower_nnz, sizeof(int));
    int * c_lower_csc_val = calloc(c_lower_nnz, sizeof(int));

    coo2csc(c_lower_csc_row, c_lower_csc_col, c_lower_csc_val,
            c_lower_row, c_lower_col, c_lower_val,
            c_lower_nnz, a_N);

    coo2csc_row_val(c_upper_csc_row, c_upper_csc_col, c_upper_csc_val,
                   c_lower_csc_row, c_lower_csc_col, c_lower_csc_val,
                   c_csc_row, c_csc_val, a_N);

    int counter = 0;

    for (int i=0; i < c_nnz; ++i)
        counter = counter + c_csc_val[i];

    int num_of_triangles = counter / 6;

    //printf("there are %d triangles.\n", num_of_triangles);


    free(c_lower_row);
    free(c_lower_col);
    free(c_lower_val);

    free(c_upper_row);
    free(c_upper_col);
    free(c_upper_val);

    free(c_coo_row);
    free(c_coo_col);
    free(c_coo_val);

    free(c_upper_csc_row);
    free(c_upper_csc_col);
    free(c_upper_csc_val);

    free(c_lower_csc_row);
    free(c_lower_csc_col);
    free(c_lower_csc_val);



}
