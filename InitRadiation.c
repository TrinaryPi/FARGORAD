#include "mp.h"

extern boolean RadCooling, Irradiation, RadTransport, VarDiscHeight, BinaryOn, OptimalOmega, ChebyshevOmega;
extern boolean Relative_Diff, Residual_Diff, Relative_Max, Residual_Max;
extern boolean LinPap1985_Opacity, BellLin1994_Opacity, PowerLaw_Opacity;
extern boolean RayTracingHeating, ExplicitRayTracingHeating, Adiabatic, HydroOn, RadiativeOnly, Cooling, CustomCooling;
extern boolean	ExplicitRadTransport;

void InitialiseRadiationModule(density, bsys)
	// Input
	PolarGrid *density;
	BinarySystem *bsys;
{
	// Function
	InitRadiationFields ();
	/* Calculate the initial disc height throughout the disc */
	ComputeDiscHeight (bsys);
	if (( RadCooling ) || ( Irradiation ) || ( RadTransport ) || ( RayTracingHeating )) {
		InitRosselandOpacity();
		ComputeRKappa (density);
	}
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
		ComputeRayTracingHeating (density, bsys);
	}
	if ( RadTransport ) {
		ComputeRadTransCoeffs(density, DT);
	}
}

void InitRadiationFields ()
	// Input N/A
{
	//Function
	DiscHeight		= CreatePolarGrid(NRAD, NSEC, "DiscHeight");
	if ( RadCooling )
		Qminus		= CreatePolarGrid(NRAD, NSEC, "Qminus");
	if ( Irradiation )
		Qirr		= CreatePolarGrid(NRAD, NSEC, "Qirr");
	if ( RadTransport ) {
		Darr		= CreatePolarGrid(NRAD, NSEC, "D");
		Barr		= CreatePolarGrid(NRAD, NSEC, "B");
		U1arr		= CreatePolarGrid(NRAD, NSEC, "U1_");
		U2arr		= CreatePolarGrid(NRAD, NSEC, "U2_");
		U3arr		= CreatePolarGrid(NRAD, NSEC, "U3_");
		U4arr		= CreatePolarGrid(NRAD, NSEC, "U4_");
		Residual	= CreatePolarGrid(NRAD, NSEC, "Residual");
		Rfld		= CreatePolarGrid(NRAD, NSEC, "R");
		lambdafld	= CreatePolarGrid(NRAD, NSEC, "lambda");
		TempSourcesSinks = CreatePolarGrid(NRAD, NSEC, "TSS");
		TempGuess = CreatePolarGrid(NRAD, NSEC, "TempGuess");
		TempGuess_old = CreatePolarGrid(NRAD, NSEC, "TempGuessOld");
	}
	if (( RadCooling ) || ( Irradiation ) || ( RadTransport ) || ( RayTracingHeating ))
		RKappaval 	= CreatePolarGrid(NRAD, NSEC, "RKappaval");
	if (( Irradiation ) || ( RadCooling )) {
		OpticalDepth = CreatePolarGrid(NRAD, NSEC, "tau");
		OpticalDepthEff = CreatePolarGrid(NRAD, NSEC, "taueff");
	}


}

void InitRosselandOpacity()
	// Input N/A
{
	//Function
	if ( LinPap1985_Opacity ) {
		masterprint("Initialising the Lin & Papaloizou (1985) Rosseland Mean Opacity table...");
	}
	if ( BellLin1994_Opacity ) {
		masterprint("Initialising the Bell & Lin (1994) Rosseland Mean Opacity table...");
	}
	if ( PowerLaw_Opacity ) {
		masterprint("Initialising a Power Law Rosseland Mean Opacity...");
	}

	RKappa = CreateOpacTable();
	RKappa = InitOpacTable();
	masterprint("Done.\n\n");
}

OpacTable *CreateOpacTable ()
// Input N/A
{
	// Declaration
	OpacTable 	*array;
	int 				ni, i;
	real 				*kappa0, *alpha, *beta, *tempmin, *tempmax, *constant, *exponent;

	// Constants
	if ( LinPap1985_Opacity ) {
		ni = 6;
	}
	if ( BellLin1994_Opacity ) {
		ni = 7;
	}
	if ( PowerLaw_Opacity ) {
		ni = 0;
	}

	//Function
	array = (OpacTable *) malloc(sizeof(OpacTable));
	if ( array == NULL )
  	erreur("Insufficient memory for Opacity Table creation");

	kappa0 		= (real *) malloc(sizeof(real) * (ni+1));
	alpha 		= (real *) malloc(sizeof(real) * (ni+1));
	beta 		= (real *) malloc(sizeof(real) * (ni+1));
	tempmin 	= (real *) malloc(sizeof(real) * (ni+1));
	tempmax 	= (real *) malloc(sizeof(real) * (ni+1));
	constant 	= (real *) malloc(sizeof(real) * (ni+1));
	exponent 	= (real *) malloc(sizeof(real) * (ni+1));

	if (( kappa0 == NULL ) || ( alpha == NULL ) ||\
		( beta == NULL ) || ( tempmin == NULL ) ||\
		( tempmax == NULL ) || ( constant == NULL ) ||\
		( exponent == NULL )) {
		fprintf (stderr, "Not enough memory.\n");
		fprintf (stderr, "Turning off Radiative Cooling and Radiative Transport - \n");
		RadCooling = NO;
		fprintf (stderr, "Please provide more memory if this is what you require!\n");
	}
	array->Ni 			= ni;
	array->Kappa0 		= kappa0;
	array->Alpha 		= alpha;
	array->Beta 		= beta;
	array->Tmin 		= tempmin;
	array->Tmax 		= tempmax;
	array->Constant 	= constant;
	array->Exponent 	= exponent;
	for (i = 0; i < ni; i++)
 		kappa0[i] = alpha[i] = beta[i] = tempmin[i] = tempmax[i] = constant[i] = exponent[i] = 0.0;

 	// Output
	return array;
}

OpacTable *InitOpacTable ()
	// Input N/A
{
	// Declaration
	int 	i, ni;
	real 	ap, bp;
	real 	*kappa0, *alpha, *beta, *tmin, *tmax, *constant, *exponent;

	// Assignment
	ni 			= RKappa->Ni;
	kappa0 		= RKappa->Kappa0;
	alpha 		= RKappa->Alpha;
	beta 		= RKappa->Beta;
	tmin 		= RKappa->Tmin;
	tmax 		= RKappa->Tmax;
	constant 	= RKappa->Constant;
	exponent 	= RKappa->Exponent;

	// Function
	if ( LinPap1985_Opacity ) {
		kappa0[0] = 2.0E-4; 	//	#1
		kappa0[1] = 2.0E16; 	//	#2
		kappa0[2] = 5.0E-3; 	//	#3
		kappa0[3] = 2.0E34; 	//	#4
		kappa0[4] = 2.0E-8; 	//	#5
		kappa0[5] = 1.0E-36; 	//	#6
		kappa0[6] = 1.5E20;		//	#7
		
		alpha[0] = 0.0; 				//	#1
		alpha[1] = 0.0; 				//	#2
		alpha[2] = 0.0; 				//	#3
		alpha[3] = 0.6666667; 			//	#4
		alpha[4] = 0.6666667; 			//	#5
		alpha[5] = 0.3333333; 			//	#6
		alpha[6] = 1.0;					//	#7
		
		beta[0] = 2.0; 					//	#1
		beta[1] = -7.0; 				//	#2
		beta[2] = 1.0; 					//	#3
		beta[3] = -9.0; 				//	#4
		beta[4] = 3.0; 					//	#5
		beta[5] = 10.0; 				//	#6
		beta[6] = -2.5;				//	#7
		
		tmax[0] = 170.0; 				//	#1
		tmax[1] = 210.0; 				//	#2
		tmax[2] = 0.0; 					//	#3
		tmax[3] = 3000.0; 			//	#4
		tmax[4] = 0.0; 					//	#5
		tmax[5] = 0.0; 					//	#6
		tmax[6] = 0.0;					//	#7
		
		constant[2] = 4.6E3;
		constant[4] = 1.1E4;
		constant[5] = 3.0E4;
		
		exponent[2] = 0.0666667;
		exponent[4] = 0.0476190;
		exponent[5] = 0.0535354;
	}

	if ( BellLin1994_Opacity) {
		kappa0[0] = 2.0E-4;		//	#1
		kappa0[1] = 2.0E16;		//	#2
		kappa0[2] = 0.1;		//	#3
		kappa0[3] = 2.0E81;		//	#4
		kappa0[4] = 1.0E-8;		//	#5
		kappa0[5] = 1.0E-36;	//	#6
		kappa0[6] = 1.5E20;		//	#7
		kappa0[7] = 0.348;		//	#8

		alpha[0] = 0.0;			//	#1
		alpha[1] = 0.0;			//	#2
		alpha[2] = 0.0;			//	#3
		alpha[3] = 1.0;			//	#4
		alpha[4] = 0.6666667;	//	#5
		alpha[5] = 0.3333333;	//	#6
		alpha[6] = 1.0;			//	#7
		alpha[7] = 0.0;			//	#8

		beta[0] = 2.0;			//	#1
		beta[1] = -7.0;			//	#2
		beta[2] = 0.5;			//	#3
		beta[3] = -24;			//	#4
		beta[4] = 3.0;			//	#5
		beta[5] = 10;			//	#6
		beta[6] = -2.5;			//	#7
		beta[7] = 0.0;			//	#8

		tmax[0] = 170.0;		//	#1
		tmax[1] = 210.0;		//	#2
		tmax[2] = 0.0;			//	#3
		tmax[3] = 0.0;			//	#4
		tmax[4] = 0.0;			//	#5
		tmax[5] = 0.0;			//	#6
		tmax[6] = 0.0;			//	#7
		tmax[7] = 0.0;			//	#8

		constant[4] = 1.1E4;
		constant[5] = 3.0E4;

		exponent[4] = 0.04762;
		exponent[5] = 0.0533;
    
    	for (i = 0; i < ni; i++) {
	    	if (( tmax[i] == 0 ) && ( constant[i] == 0 )) {
				ap = alpha[i] - alpha[i+1];
				bp = beta[i+1] - beta[i];
				constant[i] = pow(kappa0[i]/kappa0[i+1], (1.0/bp));
				exponent[i] = ap/bp;
			}
    	}
	}
	// } else {
	// 	kappa0[0]  = 1.0E-4;
	// 	kappa0[1]  = 3.0E0;
	// 	kappa0[2]  = 1.0E-2;
	// 	kappa0[3]  = 5.0E4;
	// 	kappa0[4]  = 1.0E-1;
	// 	kappa0[5]  = 2.0E15;
	// 	kappa0[6]  = 2.0E-2;
	// 	kappa0[7]  = 2.0E81;
	// 	kappa0[8]  = 1.0E-8;
	// 	kappa0[9]  = 1.0E-36;
	// 	kappa0[10] = 1.5E20;

	// 	alpha[0] = 0.0;
	// 	alpha[1] = 0.0;
	// 	alpha[2] = 0.0;
	// 	alpha[3] = 0.0;
	// 	alpha[4] = 0.0;
	// 	alpha[5] = 0.0;
	// 	alpha[6] = 0.0;
	// 	alpha[7] = 1.0;
	// 	alpha[8] = 0.6666667;
	// 	alpha[9] = 0.3333333;
	// 	alpha[10] = 1.0;

	// 	beta[0] = 2.1;
	// 	beta[1] = -0.01;
	// 	beta[2] = 1.1;
	// 	beta[3] = -1.5;
	// 	beta[4] = 0.7;
	// 	beta[5] = -5.2;
	// 	beta[6] = 0.8;
	// 	beta[7] = -24.0;
	// 	beta[8] = 3.0;
	// 	beta[9] = 10.0;
	// 	beta[10] = -2.5;

	// 	tmax[0] = 132.0;
	// 	tmax[1] = 170.0;
	// 	tmax[2] = 375.0;
	// 	tmax[3] = 390.0;
	// 	tmax[4] = 580.0;
	// 	tmax[5] = 680.0;

	// 	for (i = 6; i < 10; i++) {
	// 		ap = alpha[i] - alpha[i+1];
	// 		bp = beta[i+1] - beta[i];
	// 		constant[i] = pow(kappa0[i]/kappa0[i+1], (1.0/bp));
	// 		exponent[i] = ap/bp;
	// 	}
	// 	tmin[0] = 0.0;
	// 	for (i = 1; i < 7; i++) {
	// 		tmin[i] = tmax[i-1];
	// 	}
	// }

	// Output
	return RKappa;
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
	    	if ( tolower(*test1) == 'n' )
	    		centralsource = NO;  
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
  // if (( nsources == 1 ) && ( centralsource[0] == NO ))
  // 	centralsource[0] = YES;

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
	max = 1;
	min = 2;
	J = GLOBALNRAD;
	dy = 2*PI/(real)NSEC;
	

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
			if ( omegaOpt[i] > max )
				max = omegaOpt[i];
			if ( omegaOpt[i] < min )
				min = omegaOpt[i];
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
		for (i = 0; i < NRAD; i++)
			omegaOpt[i] = SOROMEGA;
		masterprint("# Optimal Omega       = NO\n");
		masterprint("# Omega               = %f\n", SOROMEGA);
	}
	masterprint("# Max Iterations      = %d\n", MAXITERATIONS);
	masterprint("# Relative Tolerance  = %g\n", TOLERANCE);
	if ( Relative_Diff )
		masterprint("# Tolerance Calc.     = Relative Difference\n");
	if ( Relative_Max )
		masterprint("# Tolerance Calc.     = Relative Maximum\n");
	if ( Residual_Diff )
		masterprint("# Tolerance Calc.     = Residual Difference\n");
	if ( Residual_Max )
		masterprint("# Tolerance Calc.     = Residual Maximum\n");
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
	if ( Adiabatic == NO )
		masterprint("Equation of State is Isothermal\n");
	else {
		masterprint("Equation of State is Adiabatic.\nThe following Heating/Cooling/Transport processes are active:\n");
		masterprint("############################################################################\n");
		if ( HydroOn )
			masterprint("# HydroOn            = YES. Viscous heating effects are on.\n");
	  else
	  	masterprint("# HydroOn            = NO.  Viscous heating effects are off.\n");

		if ( RadiativeOnly )
			masterprint("# RadiativeOnly      = YES. Heating/Cooling terms are radiative only (no Pressure heating).\n");
		else
			masterprint("# RadiativeOnly      = NO.  Pressure heating term is included.\n");
		
		if ( Cooling ) {
			if ( CustomCooling )
				masterprint("# CustomCooling      = YES. Using a Custom cooling prescription, Cooling time = %g P(r=1).\n", COOLINGTIME0);
			else
				masterprint("# Cooling            = YES. Cooling time = %g.\n", COOLINGTIME0);

		} else if ( RadCooling )
			masterprint("# RadCooling         = YES. Radiative Cooling from top and bottom surfaces of disc.\n");
		else
			masterprint("# Cooling            = NO.  No Cooling terms.\n");

		if ( Irradiation ) {
			if ( STARTAPER > 1 )
				masterprint("# Irradiation        = YES. Heating of disc surfaces by stellar irradiation. Stellar luminosity tapered on over %g P(r=1).\n", STARTAPER);
			else
				masterprint("# Irradiation        = YES. Heating of disc surfaces by stellar irradiation.\n");
			
		} else
		masterprint("# Irradiation        = NO.  No Irradation of disc surfaces.\n");

		if ( RayTracingHeating ) {
		masterprint("# RayTracingHeating  = YES.");
			if ( ExplicitRayTracingHeating )
				masterprint(" Disc midplane stellar irradiative heating solved explicitly.");
			else
				masterprint(" Disc midplane stellar irradiative heating solved implicitly in SOR solver.");

			if ( STARTAPER > 1 )
				masterprint(" Stellar luminosity tapered on over %g P(r=1).", STARTAPER);

			masterprint("\n");
		} else
			masterprint("# RayTracingHeating  = NO.  No RT irradiative heating.\n");

		if ( RadTransport )
			masterprint("# RadTransport       = YES. FLD of radiation in disc midplane by implicit SOR solver.\n");
		else
			masterprint("# RadTransport       = NO.  No FLD.\n");
	
		masterprint("############################################################################\n\n");
	}
}

void SetTemperatureBoundaries()
	// Input n/a
{
	// Declaration
	int i, j, l, ns, nr;
	real *T;

	nr = Temperature->Nrad;
	ns = Temperature->Nsec;
	T = Temperature->Field;

	if ( CPU_Rank == 0 ) {
		i = 0;
		for (j = 0; j < ns; j++) {
			l = j+i*ns;
			T[l] = TINNER;
		}
	}
	if ( CPU_Rank == CPU_Highest ) {
		i = nr-1;
		for (j = 0; j < ns; j++) {
			l = j+i*ns;
			T[l] = TOUTER;
		}
	}
}
