/* C Header
	* @filename        : merge.c
	* @author          : Frederic Masset
	* @last_modified_by: trinarypi
	* @last_modified   : 2017/06/26 16:59
	* @description     :
*/
#include "mp.h"

void merge (nb)
     int nb;
{
  extern boolean  Write_Density, Write_Velocity, Write_Energy, Write_Eccentricity;
  extern boolean  Write_Temperature, Write_DivV, Write_Qplus;
  extern boolean SelfGravity, SGZeroMode;
  extern boolean  Write_DiscHeight, Write_Kappa, Write_Qminus, Write_Qirr, Write_Residual, Write_Coeffs, Write_Rfld, Write_QirrRT, Write_Taus;
  extern boolean  Write_OpticalDepths;
  boolean bool = NO;
  int i, j, l;
  int one_if_odd;
  char radix[512];
  char command[1024];
  if (!CPU_Master) return;
  message ("Merging output files...");
  for (j = 0; j < 26+(AdvecteLabel == YES); j++) {
    switch (j) {
    case 0:
      strcpy (radix, "dens");
      bool = Write_Density;
      break;
    case 1: 
      strcpy (radix, "vrad");
      bool = Write_Velocity;
      break;
    case 2: 
      strcpy (radix, "vtheta");
      bool = Write_Velocity;
      break;
    case 3: 
      strcpy (radix, "energy");
      bool = Write_Energy;
      break;
    case 4: 
      strcpy (radix, "Temperature");
      bool = Write_Temperature;
      break;
    case 5: 
      strcpy (radix, "DivV");
      bool = Write_DivV;
      break;
    case 6: 
      strcpy (radix, "Qplus");
      bool = Write_Qplus;
      break;
    case 7: 
      strcpy (radix, "eccen");
      bool = Write_Eccentricity;
      break;
    case 8: 
      strcpy (radix, "DiscHeight");
      bool = Write_DiscHeight;
      break;
    case 9: 
      strcpy (radix, "RKappaval");
      bool = Write_Kappa;
      break;
    case 10: 
      strcpy (radix, "Qminus");
      bool = Write_Qminus;
      break;
    case 11: 
      strcpy (radix, "Qirr");
      bool = Write_Qirr;
      break;
    case 12: 
      strcpy (radix, "Residual");
      bool = Write_Residual;
      break;
    case 13: 
      strcpy (radix, "D");
      bool = Write_Coeffs;
      break;
    case 14: 
      strcpy (radix, "B");
      bool = Write_Coeffs;
      break;
    case 15: 
      strcpy (radix, "U1_");
      bool = Write_Coeffs;
      break;
    case 16: 
      strcpy (radix, "U2_");
      bool = Write_Coeffs;
      break;
    case 17: 
      strcpy (radix, "U3_");
      bool = Write_Coeffs;
      break;
    case 18: 
      strcpy (radix, "U4_");
      bool = Write_Coeffs;
      break;
    case 19: 
      strcpy (radix, "R");
      bool = Write_Rfld;
      break;
    case 20: 
      strcpy (radix, "lambda");
      bool = Write_Rfld;
      break;
    case 21: 
      strcpy (radix, "QirrRT");
      bool = Write_QirrRT;
      break;
    case 22: 
      strcpy (radix, "Tau_cell");
      bool = Write_Taus;
      break;
    case 23: 
      strcpy (radix, "Tau_grid");
      bool = Write_Taus;
      break;
    case 24:
      strcpy (radix, "tau");
      bool = Write_OpticalDepths;
      break;
    case 25:
      strcpy (radix, "taueff");
      bool = Write_OpticalDepths;
      break;
    case 26: 
      strcpy (radix, "label");
      bool = YES;
      break;
    }
    one_if_odd = (CPU_Number%2 == 0 ? 0 : 1);
    if (bool == YES) {
      if ( SelfGravity && !SGZeroMode ) {
	for (i = 0; i < (CPU_Number+one_if_odd)/2; i++) {
	  if ( i != 0 ) {
	    sprintf (command, "cd %s; cat gas%s%d.dat.%05d >> gas%s%d.dat", \
		     OUTPUTDIR, radix, nb, i, radix, nb);
	    system (command);
	  }
	  l = (i + (CPU_Number + one_if_odd)/2)%(CPU_Number + one_if_odd);
	  if ( i != CPU_Highest ) {
	    sprintf (command, "cd %s; cat gas%s%d.dat.%05d >> gas%s%d.dat", \
		     OUTPUTDIR, radix, nb, l, radix, nb);
	    system (command);
	  }
	}
	sprintf (command, "cd %s; rm -f gas%s%d.dat.0*",	\
		 OUTPUTDIR, radix, nb);
	system (command);
      }
      else {
	for (i = 1; i < CPU_Number; i++) {
	  sprintf (command, "cd %s; cat gas%s%d.dat.%05d >> gas%s%d.dat", \
		   OUTPUTDIR, radix, nb, i, radix, nb);
	  system (command);
	}
	sprintf (command, "cd %s; rm -f gas%s%d.dat.0*",	\
		 OUTPUTDIR, radix, nb);
	system (command);
      }
    }
  }
  message ("done\n");
  fflush (stdout);
}
