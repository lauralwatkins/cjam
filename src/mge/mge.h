/* -----------------------------------------------------------------------------
  MGE PROGRAMS
    
    mge_addbh     : add a black hole component to an MGE
    mge_dens      : MGE volume density at a given position
    mge_deproject : MGE deprojection for a given inclination angle
    mge_qmed      : MGE median flattening
    mge_read      : read MGE from file into structure
    mge_surf      : MGE surface density at a given position
    multigaussexp : MGE structure
  
  Laura L Watkins [lauralwatkins@gmail.com]
----------------------------------------------------------------------------- */

struct multigaussexp {
    double *area;
    double *sigma;
    double *q;
    int ntotal;
};

struct multigaussexp mge_addbh( struct multigaussexp *, double, double );

double mge_dens( struct multigaussexp *, double, double );

struct multigaussexp mge_deproject( struct multigaussexp *, double );

double mge_qmed( struct multigaussexp *, double );

void mge_read( char *, int, struct multigaussexp * );

double* mge_surf( struct multigaussexp *, double *, double *, int );
