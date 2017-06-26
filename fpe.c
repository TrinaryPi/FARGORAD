/* C Header
	* @filename        : fpe.c
	* @author          : Matthew Mutter (trinarypi)
	* @last_modified_by: trinarypi
	* @last_modified   : 2017/06/26 16:58
	* @description     :
*/
#include "mp.h"

void handfpe() {
  fprintf (stdout, "CPU #%d has been signaled with a floating point exception.\n",
	   CPU_Rank);
  fprintf (stdout, "Run is aborted.\n");
  prs_exit (1);
}

void setfpe () {
#ifdef _TRAP_FPE
  feenableexcept (FE_DIVBYZERO | FE_INVALID);
  /* Can also be any of the value below:
   FE_INEXACT           inexact result
   FE_DIVBYZERO         division by zero
   FE_UNDERFLOW         result not representable due to underflow
   FE_OVERFLOW          result not representable due to overflow
   FE_INVALID           invalid operation
   FE_ALL_EXCEPT        bitwise OR of all supported exceptions
  */
  signal (SIGFPE, handfpe);
#endif
}

