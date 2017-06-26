/* C Header
	* @filename        : sgproto.h
	* @author          : Matthew Mutter (trinarypi)
	* @last_modified_by: trinarypi
	* @last_modified   : 2017/06/26 16:59
	* @description     :
*/
void compute_selfgravity ();
void init_sg ();
void compute_fftdensity ();
void compute_fftkernel ();
void compute_sgacc ();
void update_sgvelocity ();
void init_azimutalvelocity_withSG ();
void init_planetarysys_withSG ();
void compute_SGZeroMode ();
void compute_anisotropic_pressurecoeff ();
