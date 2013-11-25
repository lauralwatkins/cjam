/* ----------------------------------------------------------------------------
  MGE_DEPROJECT
    
    Deprojects an MGE given an inclination value.
    
    INPUTS
      pmge : projected MGE
      incl : inclination [radians]
  
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


struct multigaussexp mge_deproject( struct multigaussexp *pmge, double incl ) {
    
    struct multigaussexp imge;
    double si, ci;
    int i;
    
    si = sin( incl );
    ci = cos( incl );
    
    imge.ntotal = pmge->ntotal;
    imge.sigma = (double *) malloc( pmge->ntotal * sizeof( double ) );
    imge.area = (double *) malloc( pmge->ntotal * sizeof( double ) );
    imge.q = (double *) malloc( pmge->ntotal * sizeof( double ) );
    
    for ( i = 0; i < pmge->ntotal ; i++ ) {
        
        // sigmas stay the same
        imge.sigma[i] = pmge->sigma[i];
        
        // convert flattening values (and exit if q is too small)
        imge.q[i] = pmge->q[i] * pmge->q[i] - ci * ci;
        if ( imge.q[i] < 0. ) {
            printf( "Inclination %lf is too low: q<0\n", incl );
            exit(0);
        }
        imge.q[i] = sqrt( imge.q[i] ) / si;
        if ( imge.q[i] < 0.05 ) {
            printf( "There are components with q<0.05\n" );
            exit(0);
        }
        
        // convert surface density to volume density
        imge.area[i] = pmge->area[i] * pmge->q[i] / pmge->sigma[i] \
            / imge.q[i] / sqrt( 2. * M_PI );
    }
    
    return imge;
    
}
