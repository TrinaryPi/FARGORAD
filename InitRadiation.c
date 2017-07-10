/* C Header
	* @filename        : InitRadiation.c
	* @author          : Matthew Mutter (trinarypi)
	* @last_modified_by: trinarypi
	* @last_modified   : 2017/06/26 13:32
	* @description     :
*/
#include "mp.h"
#include "radiation.h"

extern boolean BinaryOn;
extern boolean Adiabatic;
extern boolean HydroOn;
extern boolean RadiativeOnly;
extern boolean Cooling;
extern boolean CustomCooling;

void InitialiseRadiationModule(density, bsys)
	// Input
	PolarGrid *density;
	BinarySystem *bsys;
{
	// Function
	InitRadiationFields ();
	/* Calculate the initial disc height throughout the disc */
	ComputeDiscHeight (bsys);
	InitRosselandOpacity();
	ComputeRKappa (density);
	if (( Irradiation ) || ( RayTracingHeating )) {
		InitIrradiation();
	}
	if ( RadTransport ) {
		InitRadTransport();
	}
	if ( RayTracingHeating ) {
		InitRayTracing();
	}
	if (( RadCooling ) || ( Irradiation ) || ( RadTransport ) || ( RayTracingHeating )) {
		PrintGlobalConstants();
	}
	if ( RadCooling ) {
		ComputeQminus(density, bsys);
	}
	if ( Irradiation ) {
		ComputeQirr(density, bsys);
	}
	if ( RayTracingHeating ) {
		ComputeRayTracingHeating(density, bsys);
	}
	if ( RadTransport ) {
		ComputeRadTransCoeffs(density, DT);
	}
}

void InitRadiationFields ()
	// Input N/A
{
	//Function
	DiscHeight					= CreatePolarGrid(NRAD, NSEC, "DiscHeight");
	RKappaval 					= CreatePolarGrid(NRAD, NSEC, "RKappaval");
	if ( RadCooling ) {
		Qminus						= CreatePolarGrid(NRAD, NSEC, "Qminus");
	}
	if ( Irradiation ) {
		Qirr							= CreatePolarGrid(NRAD, NSEC, "Qirr");
	}
	if ( RadTransport ) {
		Darr							= CreatePolarGrid(NRAD, NSEC, "D");
		Barr							= CreatePolarGrid(NRAD, NSEC, "B");
		U1arr							= CreatePolarGrid(NRAD, NSEC, "U1_");
		U2arr							= CreatePolarGrid(NRAD, NSEC, "U2_");
		U3arr							= CreatePolarGrid(NRAD, NSEC, "U3_");
		U4arr							= CreatePolarGrid(NRAD, NSEC, "U4_");
		Residual					= CreatePolarGrid(NRAD, NSEC, "Residual");
		Rfld							= CreatePolarGrid(NRAD, NSEC, "R");
		lambdafld					= CreatePolarGrid(NRAD, NSEC, "lambda");
		TempSourcesSinks 	= CreatePolarGrid(NRAD, NSEC, "TSS");
		TempGuess 				= CreatePolarGrid(NRAD, NSEC, "TempGuess");
		TempGuess_old 		= CreatePolarGrid(NRAD, NSEC, "TempGuessOld");
	}
	if (( Irradiation ) || ( RadCooling )) {
		OpticalDepth 			= CreatePolarGrid(NRAD, NSEC, "tau");
		OpticalDepthEff 	= CreatePolarGrid(NRAD, NSEC, "taueff");
	}
}

void InitIrradiation()
	// Input N/A
{	
	// Declaration
	int ns;
	char *filename;

	// Function/Output
	filename = IRRADIATIONCONFIG;
	masterprint("Initialising the Stellar parameters for radiation sources...\n");
	masterprint("############################################################################\n");
	ns = FindNumberOfPlanets(filename);
	IrrSources = AllocIrradiationSources(ns);
	IrrSources = InitIrradiationSources(filename, ns);
	ListIrradiationSources();
	masterprint("############################################################################\n");
	masterprint("Done.\n\n");
}

IrradiationSource *InitIrradiationSources(filename, ns)
	// Input
	char *filename;
	int ns;
{	
	// Declaration
	int i = 0;
  char s[512], nm[512], test1[512], *s1;
  float rstar, tstar;
  real *rs, *ts, constant;
	boolean centralsource;
	FILE *input;

	// Constants
	constant = TEMPCGS;

	// Assignment
	rs = IrrSources->Rstar;
	ts = IrrSources->Tstar;

	// Function
	input = fopen (filename, "r");
	IrrSources->nb = ns;
	while ( fgets(s, 510, input) != NULL ) {
		sscanf(s, "%s ", nm);
  	if ( isalpha(s[0]) ) {
  		s1 = s + strlen(nm);
  		sscanf(s1 + strspn(s1, "\t :=>_"), "%f %f %s", &rstar, &tstar, test1);
    	centralsource = YES;
    	if ( tolower(*test1) == 'n' ) {
    		centralsource = NO;  
    	}
    	rs[i] = (real)rstar;
    	ts[i] = (real)tstar/constant;
    	IrrSources->CentralSource[i] = centralsource;
    	i++;
		}	
  }

  // Output
  return IrrSources;
}

IrradiationSource *AllocIrradiationSources(nsources)
	// Input
	int nsources;
{
	// Declaration
	int i;
  IrradiationSource *array;
  real *rstar, *tstar;
  boolean *centralsource;
  	
  // Function
  array  = (IrradiationSource *)malloc (sizeof(IrradiationSource));
  if ( array == NULL ) {
    fprintf (stderr, "Not enough memory.\n");
    prs_exit (1);
  }
  rstar = (real *)malloc (sizeof(real)*(nsources+1));
  tstar = (real *)malloc (sizeof(real)*(nsources+1));
  centralsource = (boolean *)malloc (sizeof(real)*(nsources+1));
  if (( rstar == NULL ) || ( tstar == NULL ) || ( centralsource == NULL )) {
    fprintf (stderr, "Not enough memory.\n");
    prs_exit (1);
  }
  array->Rstar = rstar;
  array->Tstar = tstar;
  array->CentralSource = centralsource;
  for (i = 0; i < nsources; i++) {
    rstar[i] = tstar[i] = 0.0;
    centralsource[i] = YES;
  }

  // Output
  return array;
}

void ListIrradiationSources ()
	// Input N/A
{
	// Declaration
  int nb, s;
  real constant, *rs, *ts;

  // Assignment
  nb = IrrSources->nb;
  rs = IrrSources->Rstar;
  ts = IrrSources->Tstar;

  // Constants
  constant = TEMPCGS;
 
  // Function/Output
  if ( !CPU_Master ) {
  	return;
  }
  for (s = 0; s < nb; s++) {
    printf ("# Source number %d:\n", s);
    printf ("# ---------------\n");
    printf ("# Rstar = %.3f AU   = %.3f Rsolar\n", rs[s], rs[s]/0.0046491);
    printf ("# Tstar = %.3f code = %.3f K\n", ts[s], ts[s]*constant);
    if ( IrrSources->CentralSource[s] ) {
      printf ("# Central Source.\n");
    } else {
      printf ("# Non-central Source.\n");
    }
  }
}

void InitRadTransport ()
	// Input N/A
{	
	// Declaration
	int i, J;
	real dx, dy, dxdy2, min, max, MINom, MAXom;

	// Constants
	max = 1.0;
	min = 2.0;
	J = GLOBALNRAD;
	dy = 2.0*PI/(real)NSEC;
	
	// Function/Output
	masterprint("Initialising the Radiation Transport Module...\n");
	masterprint("############################################################################\n");
	masterprint("# Options for SOR solver:\n");
	if ( OptimalOmega ) {
# pragma omp parallel for private (dx, dxdy2, rhoJac)
		for (i = 0; i < NRAD; i++) {
			dx = Rsup[i] - Rinf[i];
			dxdy2 = dx*dx/dy/dy;
			rhoJac[i] = (cos(2.0*(real)PI/(real)J) + dxdy2*cos(dy))/(1.0+dxdy2);
			omegaOpt[i] = 2.0/(1.0+sqrt(1.0-rhoJac[i]*rhoJac[i]));
			if ( omegaOpt[i] > max ) {
				max = omegaOpt[i];
			}
			if ( omegaOpt[i] < min ) {
				min = omegaOpt[i];
			}
		}
		MPI_Allreduce(&min, &MINom, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
		MPI_Allreduce(&max, &MAXom, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
		masterprint("# Optimal Omega       = YES\n");
		masterprint("# Omega               = %f->%f\n", MINom, MAXom);
	} else if ( ChebyshevOmega ) {
# pragma omp parallel for private (dx, dxdy2, rhoJac)
		for (i = 0; i < NRAD; i++) {
			dx = Rsup[i] - Rinf[i];
			dxdy2 = dx*dx/dy/dy;
			rhoJac[i] = (cos((real)PI/(real)J) + dxdy2*cos(dy))/(1.0+dxdy2);
		}
		masterprint("# Optimal Omega       = NO\n");
		masterprint("# Chebyshev Accl.     = YES\n");
	} else {
# pragma omp parallel for
		for (i = 0; i < NRAD; i++) {
			omegaOpt[i] = SOROMEGA;
		}
		masterprint("# Optimal Omega       = NO\n");
		masterprint("# Omega               = %f\n", SOROMEGA);
	}
	masterprint("# Max Iterations      = %d\n", MAXITERATIONS);
	masterprint("# Relative Tolerance  = %g\n", TOLERANCE);
	if ( Relative_Diff ) {
		masterprint("# Tolerance Calc.     = Relative Difference\n");
	}
	if ( Relative_Max ) {
		masterprint("# Tolerance Calc.     = Relative Maximum\n");
	}
	if ( Residual_Diff ) {
		masterprint("# Tolerance Calc.     = Residual Difference\n");
	}
	if ( Residual_Max ) {
		masterprint("# Tolerance Calc.     = Residual Maximum\n");
	}
	masterprint("############################################################################\n");
	masterprint("Done.\n\n");
	ResetTempSourcesSinks();
}

void PrintGlobalConstants()
	// Input N/A
{	
	// Function/Output
	masterprint("For future reference or for post processing\n");
	masterprint("############################################################################\n");
	masterprint("# GLOBAL CONSTANTS \n");
	masterprint("# R_MUcgs    = %g\n", R_MUCGS);
	masterprint("# E2cgs      = %g\n", E2CGS);
	masterprint("# STEFANKcgs = %g\n", STEFANKCGS);
	masterprint("# STEFANK    = %g\n", STEFANK);
	masterprint("# TEMPcgs    = %g\n", TEMPCGS);
	masterprint("# Ccgs       = %g\n", CCGS);
	masterprint("# C          = %g\n", Ccode);
	masterprint("# Acgs       = %g\n", (4.0*STEFANKCGS)/CCGS);
	masterprint("# A          = %g\n", Acode);
	masterprint("############################################################################\n\n");
}

void ListEOSParams()
	// Input N/A
{	
	// Function/Output
	masterprint("\n");
	if ( Adiabatic == NO ) {
		masterprint("Equation of State is Isothermal\n");
	} else {
		masterprint("Equation of State is Adiabatic.\nThe following Heating/Cooling/Transport processes are active:\n");
		masterprint("############################################################################\n");
		if ( HydroOn ) {
			masterprint("# HydroOn            = YES. Viscous heating effects are on.\n");
		} else {
	  	masterprint("# HydroOn            = NO.  Viscous heating effects are off.\n");
	  }

		if ( RadiativeOnly ) {
			masterprint("# RadiativeOnly      = YES. Heating/Cooling terms are radiative only (no Pressure heating).\n");
		} else {
			masterprint("# RadiativeOnly      = NO.  Pressure heating term is included.\n");
		}
		
		if ( Cooling ) {
			if ( CustomCooling ) {
				masterprint("# CustomCooling      = YES. Using a Custom cooling prescription, Cooling time = %g P(r=1).\n", COOLINGTIME0);
			} else {
				masterprint("# Cooling            = YES. Cooling time = %g.\n", COOLINGTIME0);
			}
		} else if ( RadCooling ) {
			masterprint("# RadCooling         = YES. Radiative Cooling from top and bottom surfaces of disc.\n");
		} else {
			masterprint("# Cooling            = NO.  No Cooling terms.\n");
		}

		if ( Irradiation ) {
			if ( STARTAPER > 1 ) {
				masterprint("# Irradiation        = YES. Heating of disc surfaces by stellar irradiation. Stellar luminosity tapered on over %g P(r=1).\n", STARTAPER);
			} else {
				masterprint("# Irradiation        = YES. Heating of disc surfaces by stellar irradiation.\n");
			}
		} else {
			masterprint("# Irradiation        = NO.  No Irradation of disc surfaces.\n");
		}

		if ( RayTracingHeating ) {
			masterprint("# RayTracingHeating  = YES.");
			if ( ExplicitRayTracingHeating ) {
				masterprint(" Disc midplane stellar irradiative heating solved explicitly.");
			} else {
				masterprint(" Disc midplane stellar irradiative heating solved implicitly in SOR solver.");
			}
			if ( STARTAPER > 1 ) {
				masterprint(" Stellar luminosity tapered on over %g P(r=1).", STARTAPER);
			}
			masterprint("\n");
		} else {
			masterprint("# RayTracingHeating  = NO.  No RT irradiative heating.\n");
		}

		if ( RadTransport ) {
			masterprint("# RadTransport       = YES. FLD of radiation in disc midplane by implicit SOR solver.\n");
		} else {
			masterprint("# RadTransport       = NO.  No FLD.\n");
		}
	
		masterprint("############################################################################\n\n");
	}
}

void SetTemperatureBoundaries()
	// Input n/a
{
	// Declaration
	int i, j, l, ns, nr;
	int lim, lip;
	real *T;

	nr = Temperature->Nrad;
	ns = Temperature->Nsec;
	T = Temperature->Field;

	if ( CPU_Rank == 0 ) {
		i = 1;
		for (j = 0; j < ns; j++) {
			l = j+i*ns;
	 		lim = l-ns;
	 		lip = l+ns;
			T[lim] = ApplyInnerBC(T[l], T[lip]);
		}
	}
	if ( CPU_Rank == CPU_Highest ) {
	 	i = nr-2;
	 	for (j = 0; j < ns; j++) {
	 		l = j+i*ns;
	 		lim = l-ns;
	 		lip = l+ns;
	 		T[lip] = ApplyOuterBC(T[l], T[lim]);
	 	}
	}
}


void BoundaryConditionsFLD(out_field, ref_field)
	// Input
	PolarGrid *out_field;
	PolarGrid *ref_field;
{
	// Declaration
	int i, j, l, nr, ns;
	int lim, lip;
	real *out_data, *ref_data;

	// Assignment
	nr = ref_field->Nrad;
	ns = ref_field->Nsec;
	out_data = out_field->Field;
	ref_data = ref_field->Field;

	if ( CPU_Rank == 0 ) {
		i = 1;
		for (j = 0; j < ns; j++) {
			l = j+i*ns;
	 		lim = l-ns;
	 		lip = l+ns;
			out_data[lim] = ApplyInnerBC(ref_data[l], ref_data[lip]);
		}
	}
	if ( CPU_Rank == CPU_Highest ) {
	 	i = nr-2;
	 	for (j = 0; j < ns; j++) {
	 		l = j+i*ns;
	 		lim = l-ns;
	 		lip = l+ns;
	 		out_data[lip] = ApplyOuterBC(ref_data[l], ref_data[lim]);
	 	}
	}
}
