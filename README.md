# Parallel and distribution systems: 1st project code
\
\
If we want to run the open cilk function then we comment the `#include <opp.h>` at _main_ and the whole openMP function at the _maskedSparseMatrixMatrixProduct_
\
\
If we want to run any other function then we comment the `#include <cilk/cilk.h>` and the whole openCilk function at the _maskedSparseMatrixMatrixProduct_
