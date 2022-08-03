/* ----------------------------------------------------------------------------
  JAM_AXI_RMS
    
    Wrapper for second moment calculator designed to interface with Python.
    
    INPUTS
      xp : projected x' [pc]
      yp : projected y' [pc]
      nxy : number of x' and y' values given
      incl : inclination [radians]
      lum_area : projected luminous MGE area
      lum_sigma : projected luminous MGE width
      lum_q : projected luminous MGE flattening
      pot_area : projected potential MGE area
      pot_sigma : projected potential MGE sigma
      pot_q : projected potential MGE flattening
      beta : velocity anisotropy (1-vz^2/vr^2)
      nrad : number of radial bins in interpolation grid
      nang : number of angular bins in interpolation grid
      rxx : array to hold the xx second moments calculated
      ryy : array to hold the yy second moments calculated
      rzz : array to hold the zz second moments calculated
      rxy : array to hold the xy second moments calculated
      rxz : array to hold the xz second moments calculated
      ryz : array to hold the yz second moments calculated
---------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "jam.h"
#include "../mge/mge.h"


void jam_axi_rms(double *xp, double *yp, int nxy, double incl, \
double *lum_area, double *lum_sigma, double *lum_q, int lum_total, \
double *pot_area, double *pot_sigma, double *pot_q, int pot_total, \
double *beta, int nrad, int nang, int* integrationFlag, \
double *rxx, double *ryy, double *rzz, double *rxy, double *rxz, double *ryz) {
    
    jam_axi_rms_axes(xp, yp, nxy, incl,
        lum_area, lum_sigma, lum_q, lum_total,
        pot_area, pot_sigma, pot_q, pot_total,
        beta, nrad, nang, integrationFlag,
        rxx, ryy, rzz, rxy, rxz, ryz,
        1, 1, 1);
    
    return;
}
