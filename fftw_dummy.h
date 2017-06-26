/* C Header
	* @filename        : fftw_dummy.h
	* @author          : Matthew Mutter (trinarypi)
	* @last_modified_by: trinarypi
	* @last_modified   : 2017/06/26 16:58
	* @description     :
*/
/***********************************************/
/*                                             */
/*                                             */
/*  Fake functions library for non fftw built  */
/*                                             */
/*                                             */
/***********************************************/

#define fftw_real 2
#define fftw_complex 3

typedef int rfftwnd_mpi_plan;

void fftw_malloc();
void rfftw2d_mpi_create_plan();
void rfftwnd_mpi();
