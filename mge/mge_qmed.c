/* ----------------------------------------------------------------------------
  MGE_QMED
    
    Calculates median flattening of projected MGE.
    
    INPUTS
      pmge  : projected MGE
      lim   : upper limit on size of MGE dispersion
  
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
#include "../tools/tools.h"


double mge_qmed( struct multigaussexp *pmge, double lim ) {
    
    int i, m, *goodvals;
    double *median_arr, qmed;
    
    // keep only components with sigma < lim
    goodvals = (int *) malloc( pmge->ntotal * sizeof( int ) );
    m = 0;
    for ( i = 0; i < pmge->ntotal; i++ ) {
        if ( pmge->sigma[i] < lim ) {
            goodvals[m] = i;
            m++;
        }
    }
    
    // calculate median
    if ( m < 3 ) {
        median_arr = (double *) malloc( pmge->ntotal * sizeof( double ) );
        for ( i = 0; i < pmge->ntotal; i++ ) median_arr[i] = pmge->q[i];
        qmed = median( median_arr, pmge->ntotal );
    } else {
        median_arr = where( pmge->q, goodvals, m );
        qmed = median( median_arr, m );
    }
    
    free( median_arr );
    free( goodvals );
    
    return qmed;
    
}
