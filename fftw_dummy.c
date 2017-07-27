/* C Header
	* @filename        : fftw_dummy.c
	* @author          : Frederic Masset
	* @last_modified_by: trinarypi
	* @last_modified   : 2017/06/26 16:59
	* @description     :
*/
/****************************************************/
/*                                                  */
/*                                                  */
/*  Fake fftw functions library for non-fftw build  */
/*                                                  */
/*                                                  */
/****************************************************/

#include <stdio.h>
#include "fftw_dummy.h"

int *fftw_malloc ()
{
}

int *rfftw2d_mpi_create_plan ()
{
}

int *rfftwnd_mpi_destroy_plan ()
{
}

int *rfftwnd_mpi ()
{
}

int *rfftwnd_mpi_local_sizes ()
{
}

int *rfftw_create_plan ()
{
}

int *rfftw_one ()
{
}