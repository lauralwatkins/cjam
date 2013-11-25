/* ----------------------------------------------------------------------------
  MEDIAN
    
    Calculates the median of an array of values.  For an even number of
    elements, the mean of the middle two elements of the sorted array is
    returned.  For an odd number of elements, the middle element of the sorted
    array is returned.
    
    INPUTS
      x : input array
      n : number of elements in array
    
  Laura L Watkins [lauralwatkins@gmail.com]
  
  This code is released under a BSD 2-clause license.
  If you use this code for your research, please cite:
  Watkins et al. 2013, MNRAS, 436, 2598
  "Discrete dynamical models of omega Centauri"
  http://adsabs.harvard.edu/abs/2013MNRAS.436.2598W
---------------------------------------------------------------------------- */

#include "tools.h"


double median( double *x, int n ) {
    
    sort_dbl( x, n );
    
    if ( n % 2 == 0 ) return ( x[n/2] + x[n/2 - 1] ) / 2.;
    else return x[n/2];
    
}
