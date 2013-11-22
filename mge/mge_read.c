/* ----------------------------------------------------------------------------
  MGE_READ
    
    Reads in the components of an MGE from a file.
    
    INPUTS
      file   : name of file containing MGE
      length : number of MGE components
      mge    : variable to contain MGE information
  
  Mark den Brok
  Laura L Watkins [lauralwatkins@gmail.com]
  
  This code is released under a BSD 2-clause license.
  If you use this code for your research, please cite:
  Watkins et al. 2013, MNRAS
  "Discrete dynamical models of omega Centauri"
  http://adsabs.harvard.edu/abs/2013MNRAS.tmp.2480W
---------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mge.h"
#include "../tools/readcol.h"


void mge_read( char *file, int length, struct multigaussexp *mge ) {
    
    int ntot, *intarr;
    char *flags;
    char *fmt = "%i %lf %lf %lf";
    
    flags = (char *) malloc( ( 39 + strlen( file ) ) * sizeof( char ) );
    sprintf( flags, "quiet bufsize=1024 filelen=%i skipsym=# ", length );
    
    intarr = (int *) malloc( length * sizeof( int ) );
    mge->area = (double *) malloc( length * sizeof( double ) );
    mge->sigma = (double *) malloc( length * sizeof( double ) );
    mge->q = (double *) malloc( length * sizeof( double ) );
    
    ntot = readcol( file, flags, fmt, intarr, mge->area, mge->sigma, mge->q );
    mge->ntotal = ntot;
    
    return;
    
}
