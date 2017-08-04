/* ----------------------------------------------------------------------------
  RANGE
    
    Creates an array of n numbers between given limits.
    
    INPUTS
      lo   : lower limit
      hi   : upper limit
      n    : number of elements in array
      open : if set, end limits are not included in the array
    
  Laura L Watkins [lauralwatkins@gmail.com]
  
  This code is released under a BSD 2-clause license.
  If you use this code for your research, please cite:
  Watkins et al. 2013, MNRAS, 436, 2598
  "Discrete dynamical models of omega Centauri"
  http://adsabs.harvard.edu/abs/2013MNRAS.436.2598W
---------------------------------------------------------------------------- */

#include <stdlib.h>


double* range( double lo, double hi, int n, int open ) {
    
    int i;
    double *returnval;
    
    returnval = (double *) malloc( n * sizeof( double ) );
    
    if ( open ) {
        
        // create n bins between lo and hi, return the centre of the bins
        for ( i = 0; i < n; i++ ) {
            returnval[i] = lo + ( hi - lo ) * ( 0.5 + i ) / n;
        }
    
    } else {
        
        // for n < 0, create an integer array
        if ( n < 0 ) n = (int) ( hi - lo );
        
        // return n element array with limits included
        for ( i = 0; i < n; i++ ) {
            returnval[i] = lo + ( hi - lo ) * i / ( n - 1 );
        }
        
    }
    
    return (double *) returnval;
    
}
