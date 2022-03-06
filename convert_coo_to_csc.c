/**
 *   \file coo2csc.c
 *   \brief An example of COO to CSC matrix conversion
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

/*****************************************************************************/
/*                             routine definition                            */
/*****************************************************************************/

/**
 *  \brief COO to CSC conversion
 *
 *  Converts a square matrix from COO to CSC format.
 *
 *  Note: The routine assumes the input COO and the output CSC matrix
 *  to be square.
 *
 */
void coo2csc_pattern(
  int       * row,       /*!< CSC row start indices */
  int       * col,       /*!< CSC column indices */
  int       * row_coo,   /*!< COO row indices */
  int       * col_coo,   /*!< COO column indices */
  int             nnz,       /*!< Number of nonzero elements */
  int               n         /*!< Number of rows/columns */
) {


  // ----- cannot assume that input is already 0!
  for (int l = 0; l < n+1; l++) col[l] = 0;


  // ----- find the correct column sizes
  for (int l = 0; l < nnz; l++)
    col[col_coo[l]]++;

  // ----- cumulative sum
  for (int i = 0, cumsum = 0; i < n; i++) {
    int temp = col[i];
    col[i] = cumsum;
    cumsum += temp;
  }
  col[n] = nnz;
  // ----- copy the row indices to the correct place
  for (int l = 0; l < nnz; l++) {
    int col_l;
    col_l = col_coo[l];

    int dst = col[col_l];
    row[dst] = row_coo[l];

    col[col_l]++;
  }
  // ----- revert the column pointers
  for (int i = 0, last = 0; i < n; i++) {
    int temp = col[i];
    col[i] = last;
    last = temp;
  }

}

/* it's the same with coo2csc_pattern */

void coo2csc(
  int       * row,       /*!< CSC row start indices */
  int       * col,       /*!< CSC column indices */
  int       * val,       /*!< CSC val indices*/
  int       * row_coo,   /*!< COO row indices */
  int       * col_coo,   /*!< COO column indices */
  int       * val_coo,   /*!< COO val indices*/
  int             nnz,       /*!< Number of nonzero elements */
  int               n         /*!< Number of rows/columns */
) {


  // ----- cannot assume that input is already 0!
  for (int l = 0; l < n+1; l++) col[l] = 0;


  // ----- find the correct column sizes
  for (int l = 0; l < nnz; l++)
    col[col_coo[l]]++;

  // ----- cumulative sum
  for (int i = 0, cumsum = 0; i < n; i++) {
    int temp = col[i];
    col[i] = cumsum;
    cumsum += temp;
  }
  col[n] = nnz;
  // ----- copy the row indices to the correct place
  for (int l = 0; l < nnz; l++) {
    int col_l;
    col_l = col_coo[l];

    int dst = col[col_l];
    row[dst] = row_coo[l];
    val[dst] = val_coo[l];

    col[col_l]++;
  }
  // ----- revert the column pointers
  for (int i = 0, last = 0; i < n; i++) {
    int temp = col[i];
    col[i] = last;
    last = temp;
  }

}

/*it only creates the row indices of csc form */

void coo2csc_row(
     int *upper_csc_row,
     int *upper_csc_col,
     int *lower_csc_row,
     int *lower_csc_col,
     int *csc_row,
     int N
) {

    int index = 0;

    for( int i=0; i < N; ++i){

        int upper_col_start = upper_csc_col[i];
        int upper_col_end   = upper_csc_col[i+1];

        for( int k=upper_col_start; k < upper_col_end; ++k){
            csc_row[index] = upper_csc_row[k];
            ++index;
        }

        int lower_col_start = lower_csc_col[i];
        int lower_col_end   = lower_csc_col[i+1];

        for ( int k=lower_col_start; k < lower_col_end; ++k){
            csc_row[index] = lower_csc_row[k];
            ++index;
        }
    }

}

/* it creates the row and value indices of csc form*/

void coo2csc_row_val(
     int *upper_csc_row,
     int *upper_csc_col,
     int *upper_csc_val,
     int *lower_csc_row,
     int *lower_csc_col,
     int *lower_csc_val,
     int *csc_row,
     int *csc_val,
     int N
) {

    int index = 0;

    for( int i=0; i < N; ++i){

        int upper_col_start = upper_csc_col[i];
        int upper_col_end   = upper_csc_col[i+1];

        for( int k=upper_col_start; k < upper_col_end; ++k){
            csc_row[index] = upper_csc_row[k];
            csc_val[index] = upper_csc_val[k];
            ++index;
        }

        int lower_col_start = lower_csc_col[i];
        int lower_col_end   = lower_csc_col[i+1];

        for ( int k=lower_col_start; k < lower_col_end; ++k){
            csc_row[index] = lower_csc_row[k];
            csc_val[index] = lower_csc_val[k];
            ++index;
        }
    }

}

/*it only creates the col indices of csc form */

void coo2csc_col(
  int       * col,       /*!< CSC column indices */
  int       * col_coo,   /*!< COO column indices */
  int             nnz,       /*!< Number of nonzero elements */
  int               n         /*!< Number of rows/columns */
) {


  // ----- cannot assume that input is already 0!
  for (int l = 0; l < n+1; l++) col[l] = 0;


  // ----- find the correct column sizes
  for (int l = 0; l < nnz; l++)
    col[col_coo[l]]++;

  // ----- cumulative sum
  for (int i = 0, cumsum = 0; i < n; i++) {
    int temp = col[i];
    col[i] = cumsum;
    cumsum += temp;
  }
  col[n] = nnz;


}

