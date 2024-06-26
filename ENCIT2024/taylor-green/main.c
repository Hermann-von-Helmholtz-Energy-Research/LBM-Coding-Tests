/* Includes */
#include <stdlib.h>

/* Lattice */
const unsigned int scale = 1; /* multiples of 2 up to 64 */
const unsigned int NX = 32 * scale;
const unsigned int NY = NX;
const unsigned int ndir = 9; /* D2Q9 */
const size_t mem_size_ndir = sizeof(double) * NX * NY * ndir;

/* Time Steps */
const unsigned int NSTEPS = 204800 / scale / scale;

// Index Functions
inline size_t scalar_index(unsigned int x, unsigned int y) {
    return (size_t)(NX) * (size_t)(y) + (size_t)(x);
}

/* Main */
int main (int argc, char* argv[]){
    /* allocate memory */
    double *f1 = (double *) malloc(mem_size_ndir);

    /* execute */

    /* deallocate memory */
    free(f1);

    /* return */
    return 0;
}

