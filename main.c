#include <stdio.h>
#include <stdlib.h>
#include "mmio.h"
#include "convert_coo_to_csc.h"
#include "masked_sparse_matrix_matrix_product.h"
#include <sys/time.h>
#include <math.h>
#include <string.h>
#include <omp.h>

// IMPORTANT

/* if you want to run the open cilk function then comment the #include <opp.h> and the whole openMP function
    if you want to run any other function then comment the #include <cilk/cilk.h> and the whole openCilk function */


int main(int argc, char *argv[])
{
    int ret_code;
    MM_typecode matcode;
    FILE *f = fopen("/home/csal/pds/pds-codebase/matrix_market/dblp-2010.mtx", "r"); /* you might want to change this
                                                                                        if you want to run this program */
    int M, N, nz;


    int *I,             *J;
    int *lower_row,     *lower_col;
    int *upper_row,     *upper_col;
    int *csc_row,       *csc_col;
    int *upper_csc_row, *upper_csc_col;
    int *lower_csc_row, *lower_csc_col;
    int * c_csc_row = NULL;
    int * c_csc_col = NULL;
    int * c_csc_val = NULL;
    //double *val;



    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }


    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */

    if ( !mm_is_pattern(matcode) )
    {
        printf("Sorry, this application supports only ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }

    /* find out size of sparse matrix .... */

    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
        exit(1);

    int nnz = 2 * nz; // number of non zero values for the whole matrix


    /* reseve memory for matrices */

    I =             (int *) malloc(nnz * sizeof(int));
    J =             (int *) malloc(nnz * sizeof(int));

    lower_row =     (int *) malloc(nz * sizeof(int));
    lower_col =     (int *) malloc(nz * sizeof(int));

    upper_row =     (int *) malloc(nz * sizeof(int));
    upper_col =     (int *) malloc(nz * sizeof(int));

    csc_row =       (int *) malloc(nnz * sizeof(int));
    csc_col =       (int *) malloc((N + 1) * sizeof(int));

    upper_csc_row = (int *) malloc(nz * sizeof(int));
    upper_csc_col = (int *) malloc((N + 1) * sizeof(int));

    lower_csc_row = (int *) malloc(nz * sizeof(int));
    lower_csc_col = (int *) malloc((N + 1) * sizeof(int));
    //val = (double *) malloc(nz * sizeof(double));


    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

    /* MM has only the lower triangular matrix */

    for (int i=0; i<nz; i++)
    {
        if( fscanf(f, "%d %d \n", &lower_row[i], &lower_col[i]) == 0)
            return -1;
        lower_row[i]--;  /* adjust from 1-based to 0-based */
        lower_col[i]--;
    }

    if (f !=stdin) fclose(f);

    /* save here the lower so we can make the whole matrix in coo form*/

    for (int i=0; i<nz; i++)
    {
        I[i] = lower_row[i];
        J[i] = lower_col[i];
    }



    for (int i=0; i<nz; i++)
    {
        I[nz+i] = J[i];
        J[nz+i] = I[i];
    }



    /*make the upper*/

    for ( int i=0; i < nz; ++i){
        upper_row[i] = lower_col[i];
        upper_col[i] = lower_row[i];
    }




    /************************/
    /* now write out basic characteristics of matrix */
    /************************/

    mm_write_banner(stdout, matcode);
    mm_write_mtx_crd_size(stdout, M, N, nnz);


    /* Convert COO to CSC */

    coo2csc_col(csc_col, J, nnz, N); // create only the col array of csc form for our matrix

    coo2csc_pattern(upper_csc_row, upper_csc_col,   // create csc form for the upper...
            upper_row, upper_col,
            nz, N);

    coo2csc_pattern(lower_csc_row, lower_csc_col,   // ...and the lower
            lower_row, lower_col,
            nz, N);

    coo2csc_row(upper_csc_row, upper_csc_col,       // create the row array of csc form for our matrix
                   lower_csc_row, lower_csc_col,    // by adding the csc forms of the upper and the lower matrices
                   csc_row, N);                     // now we have the whole csc form of our matrix.




    // make experiments manually to "benchmark" our functions

   

    double myTime = 0.0;

    struct timeval start,end;

    printf("Start of calculation \n\n");
    
    /*
    int num_of_threads = (int) pow(2,10);

    FILE *file_ptr = fopen("/home/csal/pds/pds-codebase/times/dblp-2010/Pthreads/10.txt", "w");



    for (int j=0; j < N_SAMPLES; ++j) {



        gettimeofday(&start,NULL); //Start timing the computation

        pThreads_masked_sparse_matrix_matrix_product(lower_csc_row, lower_csc_col,
                                                    csc_row, csc_col,
                                                    nnz, N,
                                                    c_csc_row, c_csc_col, c_csc_val,
                                                    num_of_threads);

        gettimeofday(&end,NULL); //Stop timing the computation

        myTime = (end.tv_sec+(double)end.tv_usec/1000000) - (start.tv_sec+(double)start.tv_usec/1000000);

        //printf("time: %lf seconds \n",myTime);

        fprintf(file_ptr, "%lf\n", myTime);
    }


    fclose(file_ptr);


    */




    // save here all the function calls and uncomment if you need to test them


    /************************/
        /* SEQUENTIAL */
    /************************/

    /*
    gettimeofday(&start,NULL); //Start timing the computation

    sequential_masked_sparse_matrix_matrix_product(lower_csc_row, lower_csc_col,
                                                    csc_row, csc_col,
                                                    nnz, N,
                                                    c_csc_row, c_csc_col, c_csc_val);

    gettimeofday(&end,NULL); //Stop timing the computation

    myTime = (end.tv_sec+(double)end.tv_usec/1000000) - (start.tv_sec+(double)start.tv_usec/1000000);

    printf("time: %lf seconds \n",myTime);


    printf("\n");

    */



    /************************/
        /* OPENCILK */
    /************************/

    /*
    gettimeofday(&start,NULL); //Start timing the computation

    openCilk_masked_sparse_matrix_matrix_product(lower_csc_row, lower_csc_col,
                                                    csc_row, csc_col,
                                                    nnz, N,
                                                    c_csc_row, c_csc_col, c_csc_val,
                                                    NUM_OF_THREADS);

    gettimeofday(&end,NULL); //Stop timing the computation

    myTime = (end.tv_sec+(double)end.tv_usec/1000000) - (start.tv_sec+(double)start.tv_usec/1000000);

    printf("time: %lf seconds \n",myTime);

    */




    /************************/
        /* OPENMP */
    /************************/

    /*
    gettimeofday(&start,NULL); //Start timing the computation

    omp_set_num_threads(NUM_OF_THREADS);
    openMP_masked_sparse_matrix_matrix_product(lower_csc_row, lower_csc_col,
                                                    csc_row, csc_col,
                                                    nnz, N,
                                                    c_csc_row, c_csc_col, c_csc_val);

    gettimeofday(&end,NULL); //Stop timing the computation

    myTime = (end.tv_sec+(double)end.tv_usec/1000000) - (start.tv_sec+(double)start.tv_usec/1000000);

    printf("time: %lf seconds \n",myTime);

    */



    /************************/
        /* PTHREADS */
    /************************/

    /*
    gettimeofday(&start,NULL); //Start timing the computation


    pThreads_masked_sparse_matrix_matrix_product(lower_csc_row, lower_csc_col,
                                                    csc_row, csc_col,
                                                    nnz, N,
                                                    c_csc_row, c_csc_col, c_csc_val,
                                                    NUM_OF_THREADS);

    gettimeofday(&end,NULL); //Stop timing the computation

    myTime = (end.tv_sec+(double)end.tv_usec/1000000) - (start.tv_sec+(double)start.tv_usec/1000000);

    printf("time: %lf seconds \n",myTime);

    */




    printf("\nEnd of calculation \a \n");





    free(lower_col);
    free(lower_row);
    free(upper_col);
    free(upper_row);
    free(lower_csc_col);
    free(lower_csc_row);
    free(upper_csc_col);
    free(upper_csc_row);
    free(I);
    free(J);
    free(csc_row);
    free(csc_col);

	return 0;
}

