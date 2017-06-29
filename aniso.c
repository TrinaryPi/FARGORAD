/* C Header
	* @filename        : aniso.c
	* @author          : Clement Baruteau
	* @last_modified_by: trinarypi
	* @last_modified   : 2017/06/26 16:58
	* @description     :
*/
#include "mp.h"
/* In this routine we calculate the anisotropy coefficient alpha that
   models the non-axisymmetric component of the disk self-gravity, by
   means of an anisotropic pressure tensor (see Baruteau & Masset
   07b). Alpha is defined here as alpha = 1-beta/Q, where beta depends
   only on the smoothing length to disk thickness ratio
   (eta=eps/H): */
/*           eta           beta       */
/*           0.1          0.32(4)     */
/*           0.3          0.61(4)     */
/*           0.6          0.94(1)     */

void compute_anisotropic_pressurecoeff (sys)
     PlanetarySystem *sys;
{
  real xpl, ypl, rpl;
  real Q, beta, alpha;
  xpl = sys->x[0];
  ypl = sys->y[0];
  rpl = sqrt(xpl*xpl + ypl*ypl);
  //Q = pow(rpl,-2.0+SIGMASLOPE)*ASPECTRATIO/M_PI/SIGMA0;
  Q = pow(rpl,-2.0+SIGMASLOPE)*compute_aspectratio(xpl, ypl)/M_PI/SIGMA0;
  beta = 0.0;
  if ( fabs(THICKNESSSMOOTHING-0.1) < 1e-2 ) beta=0.324;
  if ( fabs(THICKNESSSMOOTHING-0.3) < 1e-2 ) beta=0.614;
  if ( fabs(THICKNESSSMOOTHING-0.6) < 1e-2 ) beta=0.941;
  alpha = 1.0 - beta/Q;
  if (SG_initcounter == 1) {
    printf ("Q = %lg, beta = %lg et alpha = %lg\n", Q, beta, alpha);
  }
  SG_aniso_coeff = alpha;
}
