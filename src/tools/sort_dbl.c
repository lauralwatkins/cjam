/* ----------------------------------------------------------------------------
  SORT_DBL
    
    Sorts an array of doubles and returns the sorted array
    
    INPUTS
      x : input array to be sorted
      n : number of elements in array
    
  Laura L Watkins [lauralwatkins@gmail.com]
  
  This code is released under a BSD 2-clause license.
  If you use this code for your research, please cite:
  Watkins et al. 2013, MNRAS, 436, 2598
  "Discrete dynamical models of omega Centauri"
  http://adsabs.harvard.edu/abs/2013MNRAS.436.2598W
---------------------------------------------------------------------------- */

#include <stdlib.h>


// comparison function for doubles
int compare_dbl( const void *xx, const void *yy ) {
    
    double x = *( (double *) xx );
    double y = *( (double *) yy );
    
    if ( x > y ) return 1;
    else {
        if ( x < y ) return -1;
        else return 0;
    }
    
};


// sorting function for doubles using quicksort
void sort_dbl( double *array, int n ) {
    
    qsort( array, n, sizeof( double ), compare_dbl );
    
}
