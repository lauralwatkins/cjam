/* ----------------------------------------------------------------------------
  WHERE
    
    Selects a given subset of an array.
    
    INPUTS
      in     : input array
      select : indices of selected elements
      n      : number of elements in select
    
  Laura L Watkins [lauralwatkins@gmail.com]
  
  This code is released under a BSD 2-clause license.
  If you use this code for your research, please cite:
  Watkins et al. 2013, MNRAS, 436, 2598
  "Discrete dynamical models of omega Centauri"
  http://adsabs.harvard.edu/abs/2013MNRAS.436.2598W
---------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>


double* where( double *in, int *select, int n ) {
    
    int i;
    double *out;
    
    // allocate memory
    out = (double *) malloc( n * sizeof( double ) );
    if ( out == NULL ) {
        printf( "Error: failed in memory declaration\n" );
        exit(1);
    }
    
    // only take selected elements
    for ( i = 0; i < n; i++ ) {
        out[i] = in[ select[i] ];
    }
    
    return out;
    
}
