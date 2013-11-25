/* ----------------------------------------------------------------------------
  MGE_SURF
    
    Calculates the surface density of an MGE at given projected (x',y').
    
    INPUTS
      mge : projected MGE
      xp  : projected x'
      yp  : projected y'
      nxy : number of (x',y') points
  
  Laura L Watkins [lauralwatkins@gmail.com]
  
  This code is released under a BSD 2-clause license.
  If you use this code for your research, please cite:
  Watkins et al. 2013, MNRAS, 436, 2598
  "Discrete dynamical models of omega Centauri"
  http://adsabs.harvard.edu/abs/2013MNRAS.436.2598W
---------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mge.h"


double* mge_surf( struct multigaussexp *mge, double *xp, double *yp, int nxy ) {
    
    int i, j;
    double *surf;
    
    surf = (double *) malloc( nxy * sizeof( double ) );
    
    for ( j = 0; j < nxy ; j++ ) { // positions
        
        surf[j] = 0.;
        for ( i = 0; i < mge->ntotal; i++ ) { // mges
            surf[j] += mge->area[i] * exp( -0.5 / pow( mge->sigma[i], 2 ) \
                * ( xp[j] * xp[j] + pow( yp[j] / mge->q[i], 2 ) ) );
        }
        
    }
    
    return surf;

}
