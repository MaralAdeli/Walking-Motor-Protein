# define Nx 7 /* number of particles in chain along x */
# define Ny 1 /* number of particles in chain along y */
# define Nz 1 /* number of particles in chain along z */
# define N (Nx*Ny*Nz) /* Number of chain Particles*/
# define Mx 2 /* number of particles in motor along x */
# define My 1 /* number of particles in motor along y */
# define Mz 1 /* number of particles in motor along z */
# define M (Mx*My*Mz) /* Number of motor Particles*/
# define DIM 3 /* Coordinates*/
# define NORM 2147483647.0
# define rnd() ((double)rand()/NORM)
# define NSTEP 1000
# define NTRAJ 100
# define NGAP (NSTEP/NTRAJ)

# define rcut 3  /* cut off length for lj potential */
# define LJR0  (pow(2.0,(1.0/6.0))) 
# define F0   0  /* external force magnitude */
# define chsig 1.0 /* chain sigma*/
# define cheps 1.0 /*chain epsilon*/ 
# define sb  (LJR0*chsig) /* chain spring equilibrium */

double delt=5.0e-6;
double sc=40; /* chain spring constant  */
double psc =80;/*protein sprig constant  */
double psb=0.9*chsig;/*spring equillibrium */
double h=0.93*chsig; /* min height */
double De = 4*cheps;/* morse */
double kappa = 4.0/chsig; /*morse*/
double chm=1.0; /*chain mass */
double pm=1.0; /*protein mass */
double bm;
