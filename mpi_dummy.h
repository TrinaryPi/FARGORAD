/* C Header
	* @filename        : mpi_dummy.h
	* @author          : Frederic Masset
	* @last_modified_by: trinarypi
	* @last_modified   : 2017/06/26 16:59
	* @description     :
*/
/****************************************************/
/*                                                  */
/*                                                  */
/*  Fake MPI functions library for sequential built */
/*                                                  */
/*                                                  */
/****************************************************/

#define MPI_COMM_WORLD 0
#define MPI_DOUBLE 2
#define MPI_CHAR 1
#define MPI_LONG 3
#define MPI_INT 0
#define MPI_MIN 0
#define MPI_MAX 0
#define MPI_SUM 0

typedef int MPI_Request;
typedef int MPI_Status;

void MPI_Comm_rank ();
void MPI_Barrier ();
void MPI_Comm_size ();
void MPI_Init ();
void MPI_Finalize ();
void MPI_Bcast ();
void MPI_ISend ();
void MPI_IRecv ();
void MPI_Allreduce ();
void MPI_Gather ();
void MPI_Gatherv ();
void MPI_Scatterv ();
