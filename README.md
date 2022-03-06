# Parallel and distribution systems: 1st project code

Με δεδομένο το adjacency matrix ενός γράφου μπορούμε να υπολογίσουμε τον αριθμό των τριγώνων
που σχηματίζονται ακολουθώντας τα μονοπάτια μήκους 3 που καταλήγουν στον αρχικό κόμβο.
Μαθηματικά περιγράφεται από τον τύπο:
\
`C = A .* (AA)`
\
Όπου το (AA) μας δίνει τον αριθμό των μονοπατιών μηκους 2 από έναν κόμβο σε έναν άλλο και το
Hadamard Product το συνολικό τρίγωνο.
\
Όλα τα adjacency matrices αποτελούν τετράγωνους,δυαδικούς και συμμετρικούς πίνακες και εύκολα
αποδεικνήουμε πώς τοσο το inner product αλλα και το outer product ,το row-by-row και το column-by-
column product μας δίνουν το τετράγωνο του Α.
\
Στην προκειμένη περίπτωση χρησιμοποιούμε CSC μορφή για την αναπρασταση των αραίων adjacency
matrices και έτσι ευνοείται ο υπολογισμός του col-by-col γινομένου καθώς παράγονται μη μηδενικά
στοιχεία μόνο μεταξύ στηλών που παρουσιάζουν μη-μηδενικα στοιχεία (1) σε κοινές στήλες. Τέλος το
Hadamard product λειτουργεί ως μάσκα και μας δίνει την δυνατότητα να αναζητήσουμε τα non zero
στοιχεία μόνο στις θέσεις που υποδεικνύει ο πίνακας Α.
\
\
If we want to run the open cilk function then we comment the `#include <opp.h>` at _main_ and the whole openMP function at the _maskedSparseMatrixMatrixProduct_
\
\
If we want to run any other function then we comment the `#include <cilk/cilk.h>` and the whole openCilk function at the _maskedSparseMatrixMatrixProduct_
