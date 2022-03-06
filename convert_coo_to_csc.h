#ifndef CONVERT_COO_TO_CSC
#define CONVERT_COO_TO_CSC

void coo2csc_pattern(
  int       * row,       /*!< CSC row start indices */
  int       * col,       /*!< CSC column indices */
  int       * row_coo,   /*!< COO row indices */
  int       * col_coo,   /*!< COO column indices */
  int             nnz,       /*!< Number of nonzero elements */
  int               n         /*!< Number of rows/columns */
);

void coo2csc(
  int       * row,       /*!< CSC row start indices */
  int       * col,       /*!< CSC column indices */
  int       * val,
  int       * row_coo,   /*!< COO row indices */
  int       * col_coo,   /*!< COO column indices */
  int       * val_coo,
  int             nnz,       /*!< Number of nonzero elements */
  int               n         /*!< Number of rows/columns */
);

void coo2csc_row(
     int *upper_csc_row,
     int *upper_csc_col,
     int *lower_csc_row,
     int *lower_csc_col,
     int *csc_row,
     int N
);

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
);

void coo2csc_col(
  int       * col,       /*!< CSC column indices */
  int       * col_coo,   /*!< COO column indices */
  int             nnz,       /*!< Number of nonzero elements */
  int               n         /*!< Number of rows/columns */
);

#endif
