/* ----------------------------------------------------------------------------
  MGE_DENS
    
    Calculates the volume density of an MGE at given (R,z).
    
    INPUTS
      mge : intrinsic MGE
      r   : intrinsic R
      z   : intrinsic z
  
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


double mge_dens( struct multigaussexp *mge, double r, double z ) {
    
    int i;
    double dens;
    
    dens = 0.;
    for ( i = 0; i < mge->ntotal; i++ ) {
        dens += mge->area[i] * exp( -0.5 / pow( mge->sigma[i], 2 ) \
            * ( r * r + pow( z / mge->q[i], 2 ) ) );
    }
    
    return dens;
    
}
