/* C Header
	* @filename        : fftw_dummy.h
	* @author          : Clement Baruteau
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
#define FFTW_REAL_TO_COMPLEX 0 
#define FFTW_COMPLEX_TO_REAL 0
#define FFTW_MEASURE 0
#define FFTW_IN_PLACE 0
#define FFTW_TRANSPOSED_ORDER 0

typedef double fftw_real;
typedef struct {
     fftw_real re, im;
} fftw_complex;
typedef int rfftwnd_mpi_plan;
typedef int rfftw_plan;

int *fftw_malloc();
int *rfftw2d_mpi_create_plan();
int *rfftwnd_mpi_destroy_plan();
int *rfftwnd_mpi();
int *rfftwnd_mpi_local_sizes();
int *rfftw_create_plan();
int *rfftw_one();
