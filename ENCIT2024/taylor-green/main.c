
/* Constants */
const unsigned int scale = 1;
const unsigned int NX = 32 * scale;
const unsigned int NY = NX;

/* Lattice */
const unsigned int ndir = 9; /* D2Q9 */
const size_t mem_size_ndir = sizeof(double) * NX * NY * ndir;

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

