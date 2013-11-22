/* ----------------------------------------------------------------------------
  MINIMUM
    
    Returns the minimum of an array of values.
    
    INPUTS
      x : input array
      n : number of elements in array
    
  Laura L Watkins [lauralwatkins@gmail.com]
  
  This code is released under a BSD 2-clause license.
  If you use this code for your research, please cite:
  Watkins et al. 2013, MNRAS
  "Discrete dynamical models of omega Centauri"
  http://adsabs.harvard.edu/abs/2013MNRAS.tmp.2480W
---------------------------------------------------------------------------- */


double minimum( double *x, int n ) {
    
    int i;
    double mx = x[0];
    
    for ( i = 0; i < n; i++ ) {
        if ( x[i] < mx ) mx = x[i];
    }
    
    return mx;
    
}
