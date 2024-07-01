/* Includes */
#include <stdlib.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Lattice */
const unsigned int scale = 1; /* multiples of 2 up to 64 */
const unsigned int NX = 32 * scale;
const unsigned int NY = NX;
const unsigned int ndir = 9; /* D2Q9 */
const size_t mem_size_ndir = sizeof(double) * NX * NY * ndir;
const size_t mem_size_scalar = sizeof(double) * NX * NY;
// Weights
const double w0 = 4/9;    // zero weight
const double ws = 1/9;    // adjacent weight
const double wd = 1/36;   // diagonal weight
// Array of lattice weights
const double wi[] = {w0,ws,ws,ws,ws,wd,wd,wd,wd};
const int dirx[] = {0,1,0,-1,0,1,-1,-1,1};
const int diry[] = {0,0,1,0,-1,1,1,-1,-1};
//Constants for physical and simulation parameters
const double nu = 1/6;
const double tau = 3 * nu + 0.5;
const double u_max = 0.04/scale;
const double rho0 = 1;
/* Time Steps */
const unsigned int NSTEPS = 204800 / scale / scale;

// Index Functions
extern inline size_t scalar_index(unsigned int x, unsigned int y) {
    return (size_t)(NX) * (size_t)(y) + (size_t)(x);
}
extern inline size_t field_index(unsigned int x, unsigned int y, unsigned int d) {
    return (size_t)(ndir) * ((size_t)(NX) * (size_t)(y) + (size_t)(x)) + (size_t)(d);
}

// Taylor-Green Functions
void taylor_green(unsigned int t, unsigned int x, unsigned int y, double *r, double *u, double *v) {
    double kx = 2 * M_PI/NX;
    double ky = 2 * M_PI/NY;
    double td = 1/(nu * (kx * kx + ky * ky));

    double X = x + 0.5;
    double Y = y + 0.5;
    double ux = -u_max * sqrt(ky/kx) * cos(kx*X) * sin(ky*Y) * exp(-1*t/td);
    double uy = u_max * sqrt(kx/ky) * sin(kx*X) * cos(ky*Y) * exp(-1*t/td);
    double P = -0.25*rho0*u_max*u_max*((ky/kx)*cos(2*kx*X)+(kx/ky)*cos(2*ky*Y))*exp(-2*t/td);
    double rho = rho0 + 3 * P;

    *r = rho;
    *u = ux;
    *v = uy;
}

void taylor_green2(unsigned int t, double *r, double *u, double *v) {
    for(unsigned int y = 0; y < NY; ++y)
    for(unsigned int x = 0; x < NX; ++x) {
        size_t sidx = scalar_index(x,y);
	taylor_green(t,x,y,&r[sidx],&u[sidx],&v[sidx]);
    }
}

// Particle Equilibrium
void init_equilibrium(double *f, double *r, double *u, double *v) {
    for(unsigned int y = 0; y < NY; ++y) {
	for(unsigned int x = 0; x < NX; ++x) {
      	    double rho = r[scalar_index(x,y)];
	    double ux = u[scalar_index(x,y)];
	    double uy = v[scalar_index(x,y)];

	    for(unsigned int i = 0; i < ndir; ++i) {
	        double cidotu = dirx[i] * ux + diry[i] * uy;
	        f[field_index(x,y,i)] = wi[i] * rho * (1 + 3 * cidotu + 4.5 * cidotu * cidotu - 1.5 * (ux*ux+uy*uy));
            }
        }
    }
}

// Stream function
void stream(double *f_src, double *f_dst) {
    for(unsigned int y = 0; y < NY; ++y) {
	for(unsigned int x = 0; x < NX; ++x) {
	    for(unsigned int i = 0; i < ndir; ++i) {
		// enforce periodicity, add NX to ensure that value is positive
		unsigned int xmd = (NX + x - dirx[i]) %NX;
		unsigned int ymd = (NY + y - diry[i]) %NY;

		f_dst[field_index(x,y,i)] = f_src[field_index(xmd,ymd,i)];
            }
	}
    }
}

// Compute density and velocity
void compute_rho_u(double *f, double *r, double *u, double *v) {
    for(unsigned int y = 0; y < NY; ++y) {
	for(unsigned int x = 0; x < NX; ++x) {
	    double rho = 0;
	    double ux = 0;
	    double uy = 0;

	    for(unsigned int i = 0; i < ndir; ++i) {
		rho += f[field_index(x,y,i)];
		ux += dirx[i] * f[field_index(x,y,i)];
		uy += diry[i] * f[field_index(x,y,i)];
	    }
	    r[scalar_index(x,y)] = rho;
	    u[scalar_index(x,y)] = ux/rho;
	    v[scalar_index(x,y)] = uy/rho;
	}
    }
}

// Collision operation
void collide(double *f, double *r, double *u, double *v) {
    // useful constance
    const double tauinv = 2/(6 * nu * 1);    // 1/tau
    const double omtauinv = 1 - tauinv;      // 1 - 1/tau

    for(unsigned int y = 0; y < NY; ++y) {
	for(unsigned int x = 0; x < NX; ++x) {
	    double rho = r[scalar_index(x,y)];
	    double ux = u[scalar_index(x,y)];
	    double uy = v[scalar_index(x,y)];

	    for(unsigned int i = 0; i < ndir; ++i) {
		// calculate dot product
		double cidotu = dirx[i] * ux + diry[i] * uy;
		// calculate equilibrium
		double feq = wi[i] * rho * (1 + 3*cidotu + 4.5*cidotu*cidotu - 1.5*(ux*ux + uy*uy));

		// relax to equilibrium
		f[field_index(x,y,i)] = omtauinv*f[field_index(x,y,i)] + tauinv * feq;
	    }
	}
    }
}


/* Main */
int main (int argc, char* argv[]){
    /* allocate memory */
    double *f1 = (double*) malloc(mem_size_ndir);
    double *f2 = (double*) malloc(mem_size_ndir);
    double *rho = (double*) malloc(mem_size_scalar);
    double *ux = (double*) malloc(mem_size_scalar);
    double *uy = (double*) malloc(mem_size_scalar);

    // Compute Taylor-Green flow at t = 0 to initialise rho, ux, uy fields
    taylor_green2(0,rho,ux,uy);

    // Initialise f1 as equilibrium for rho, ux, uy
    init_equilibrium(f1,rho,ux,uy);

    // Main simuation loop
    for(unsigned int n = 0; n < NSTEPS; ++n)  {
    // Stream from f1 storing to f2
    stream(f1,f2);

    // Calculate post-streaming density and velocity
    compute_rho_u(f2,rho,ux,uy);

    // Perform collision on f2
    collide(f2,rho,ux,uy);

    // Swap pointers
    double *temp = f1;
    f1 = f2;
    f2 = temp;
    }

    // Deallocate memory
    free(f1); free(f2);
    free(rho); free(ux); free(uy);

    /* return */
    return 0;
}

