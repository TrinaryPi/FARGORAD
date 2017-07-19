/* C Header
	* @filename        : Opacity.c
	* @author          : Matthew Mutter (trinarypi)
	* @last_modified_by: trinarypi
	* @last_modified   : 2017/06/26 16:57
	* @description     :
*/
#include "mp.h"
#include "radiation.h"

void InitRosselandOpacity()
	// Input N/A
{
	//Function
	if ( LinPap1985_Opacity ) {
		masterprint("Initialising the Lin & Papaloizou (1985) Rosseland Mean Opacity table...");
	} else if ( BellLin1994_Opacity ) {
		masterprint("Initialising the Bell & Lin (1994) Rosseland Mean Opacity table...");
	} else if ( PowerLaw_Opacity ) {
		masterprint("Initialising a Power Law Rosseland Mean Opacity...");
	} else {
		masterprint("Initialising a Constant Rosseland Mean Opacity...");
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
	real 				*kappa0, *alpha, *beta, *tempmax, *constant, *exponent;

	// Constants
	if ( LinPap1985_Opacity ) {
		ni = 6;
	} else if ( BellLin1994_Opacity ) {
		ni = 7;
	} else {
		ni = 0;
	}

	//Function
	array = (OpacTable *) malloc(sizeof(OpacTable));
	if ( array == NULL ) {
  	erreur("Insufficient memory for Opacity Table creation");
	}

	kappa0 		= (real *) malloc(sizeof(real) * (ni+1));
	alpha 		= (real *) malloc(sizeof(real) * (ni+1));
	beta 			= (real *) malloc(sizeof(real) * (ni+1));
	tempmax 	= (real *) malloc(sizeof(real) * (ni+1));
	constant 	= (real *) malloc(sizeof(real) * (ni+1));
	exponent 	= (real *) malloc(sizeof(real) * (ni+1));

	if (( kappa0 == NULL ) || ( alpha == NULL ) ||\
		( beta == NULL ) || ( tempmax == NULL ) ||\
		( constant == NULL ) || ( exponent == NULL )) {
		fprintf (stderr, "Not enough memory.\n");
		fprintf (stderr, "Please provide more memory if this is what you require!\n");
	}
	array->Ni				= ni;
	array->Kappa0		= kappa0;
	array->Alpha		= alpha;
	array->Beta			= beta;
	array->Tmax			= tempmax;
	array->Constant	= constant;
	array->Exponent	= exponent;
	for (i = 0; i < ni; i++) {
 		kappa0[i] = alpha[i] = beta[i] = tempmax[i] = constant[i] = exponent[i] = 0.0;
	}

 	// Output
	return array;
}

OpacTable *InitOpacTable ()
	// Input N/A
{
	// Declaration
	int 	i, ni;
	real 	ap, bp;
	real 	*kappa0, *alpha, *beta, *tmax, *constant, *exponent;

	// Assignment
	ni 			= RKappa->Ni;
	kappa0 		= RKappa->Kappa0;
	alpha 		= RKappa->Alpha;
	beta 		= RKappa->Beta;
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
	if ( BellLin1994_Opacity ) {
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

	// Output
	return RKappa;
}

void ComputeRKappa (Sigma)
	// Input
	PolarGrid *Sigma;
{
	// Declaration
	int i, j, l, ns, nr;
	real *dens, *temp, *H, *rkappa;
	real rhocgs, tempcgs;
	real BETA_POWER, GAMMA_POWER, c1, cgsconstant1, cgsconstant2, kap0;

	// Assignment
	ns = Sigma->Nsec;
	nr = Sigma->Nrad;
	dens = Sigma->Field;
	temp = Temperature->Field;
	H = DiscHeight->Field;
	rkappa = RKappaval->Field;

	// Constants
	BETA_POWER = -0.5;
	GAMMA_POWER = 0.0;
	kap0 = 3.0E1;
	c1 = 0.5;
	cgsconstant1 = c1*MCGS/pow(LCGS, 3.0);
	cgsconstant2 = TEMPCGS;

	// Function
	if ( Constant_Opacity ) {
		for (i = 0; i < nr; i++) {
			for (j = 0; j < ns; j++) {
				l = j+i*ns;
				rkappa[l] = KAPPAR_CONSTANT;
			}
		}
	} else if ( PowerLaw_Opacity ) {
		for (i = 0; i < nr; i++) {
			for (j = 0; j < ns; j++) {
				l = j+i*ns;
				tempcgs = cgsconstant2*temp[l];
				rkappa[l] = kap0*pow(tempcgs, BETA_POWER)*pow(Rmed[i], GAMMA_POWER);
			}
		}
	} else {
		for (i = 0; i < nr; i++) {
			for (j = 0; j < ns; j++) {
				l = j+i*ns;
				tempcgs = cgsconstant2*temp[l];
				rhocgs = cgsconstant1*dens[l]/H[l];
				rkappa[l] = CoolingRegime(rhocgs, tempcgs);
			}
		}
	}

	// Debug
	if ( RadiationDebug ) {
		int check_neg = 1;
    int check_zero = 1;
		CheckField(RKappaval, check_neg, check_zero, "ComputeRKappa");
	}
}

real CoolingRegime (rho, T)
	// Input
	real rho, T;
{
	// Declaration
	int i, max_i;
	real *kappa0, *alpha, *beta, *cons, *expn;
	real  kappacgs, tmax_max, tmax;

	// Assignment
	max_i = RKappa->Ni;
	kappa0 = RKappa->Kappa0;
	alpha = RKappa->Alpha;
	beta = RKappa->Beta;
	cons = RKappa->Constant;
	expn = RKappa->Exponent;

	// Function
	if ( OpacitySmoothing ) {
		if ( LinPap1985_Opacity ) {
			const real power1 = 4.444E-2;
			const real power2 = 2.381E-2;
			const real power3 = 2.267E-1;

			const real t234 = 1.6E3;
			const real t456 = 5.7E3;
			const real t678 = 2.28E6;

			/* coefficients for opacity laws 1, 2, and 3 in cgs units */
			const real ak1 = 2.0E-4;
			const real ak2 = 2.0E16;
			const real ak3 = 5.0E-3;

			/* coefficients for opacity laws 3, 4, 5, 6, 7, and 8 in T_4 units */
			const real bk3 = 50.0;
			const real bk4 = 2.0E-2;
			const real bk5 = 2.0E4;
			const real bk6 = 1.0E4;
			const real bk7 = 1.5E10;
			const real bk8 = 0.348;

			/* test T against (T_23 * T_34 * T_34)**0.333333333 */
			if ( T > t234*pow(rho, power1) ) {
				/* to avoid overflow */
				real ts4 = 1.0E-4*T;
				real density13 = pow(rho, 1.0/3.0);
				real density23 = density13*density13;
				real ts42 = ts4*ts4;
				real ts44 = ts42*ts42;
				real ts48 = ts44*ts44;

				/* test T against (T_45 * T_56)**0.5 */
				if ( T > t456 * pow(rho, power2) ) {
					if (( T < t678 * pow(rho, power3) ) || ( rho <= 1.0E-10 )) {
						/* disjoint opacity laws for 5, 6, and 7 */
						real o5 = bk5*density23*ts42*ts4;
						real o6 = bk6*density13*ts48*ts42;
						real o7 = bk7*rho/(ts42*sqrt(ts4));

						/* parameters used for smoothing */
						real o6an = o6*o6;
						real o7an = o7*o7;

						/* smoothed and continuous opacity law for regions 5, 6, and 7 */
						kappacgs = pow(pow(o6an*o7an/(o6an + o7an), 2.0) + pow(o5/(1.0 + pow(ts4/(1.1*pow(rho, 0.04762)), 10.0)), 4.0), 0.25);
					} else {
						/* disjoint opacity laws for 7 and 8 */
						real o7 = bk7*rho/(ts42*sqrt(ts4));
						real o8 = bk8;

						/* parameters used for smoothing */
						real o7an = o7*o7;
						real o8an = o8*o8;

						/* smoothed and continuous opacity law for regions 7 and 8 */
						kappacgs = pow(o7an*o7an + o8an*o8an, 0.25);
						/* no scattering */
						//return bk7*density/(ts42*sqrt(ts4));
					}
				} else {
					/*  disjoint opacity laws for 3, 4, and 5 */
					real o3 = bk3 * ts4;
					real o4 = bk4 * density23 / ( ts48 * ts4 );
					real o5 = bk5 * density23 * ts42 * ts4;
					/* parameters used for smoothing */
					real o4an = o4 * o4 * o4 * o4;
					real o3an = o3 * o3 * o3;

					/* smoothed and continuous opacity law for regions 3, 4, and 5 */
					kappacgs =  pow((o4an*o3an/(o4an + o3an)) + pow(o5/(1.0 + 6.561E-5/ts48), 4.0), 0.25);
				}
			} else {
				/* different powers of temperature */
				real t2 = T*T;
				real t4 = t2*t2;
				real t8 = t4*t4;
				real t10 = t8*t2;

				/* disjoint opacity laws */
				real o1 = ak1*t2;
				real o2 = ak2*T/t8;
				real o3 = ak3*T;

				/* parameters used for smoothing */
				real o1an = o1*o1;
				real o2an = o2*o2;

				/* smoothed and continuous opacity law for regions 1, 2, and 3 */
				kappacgs = pow(pow(o1an*o2an/(o1an + o2an), 2.0) + pow(o3/(1.0 + 1.0E22/t10), 4.0), 0.25);
			}
		}
		if ( BellLin1994_Opacity) {
			const real power1 = 2.8369E-2;
			const real power2 = 1.1464E-2;
			const real power3 = 2.2667E-1;

			const real t234 = 1.46E3;
			const real t456 = 4.51E3;
			const real t678 = 2.37E6;

			/* coefficients for opacity laws 1, 2, and 3 in cgs units */
			const real ak1 = 2.0E-4;
			const real ak2 = 2.0E16;
			const real ak3 = 0.1E0;

			/* coefficients for opacity laws 3, 4, 5, 6, 7, and 8 in T_4 units */
			const real bk3 = 10.0;
			const real bk4 = 2.0E-15;
			const real bk5 = 1.0E4;
			const real bk6 = 1.0E4;
			const real bk7 = 1.5E10;
			const real bk8 = 0.348;

			if( T < 1.0 ) {
				T = 10.0;
			}

			if ( T > t234*pow(rho, power1) ) {
				/* to avoid overflow */
				real ts4 = 1.0E-4*T;
				real density13 = pow(rho, 1.0/3.0);
				real density23 = density13*density13;
				real ts42 = ts4*ts4;
				real ts44 = ts42*ts42;
				real ts48 = ts44*ts44;

				/* test T against (T_45 * T_56)**0.5 */
				if ( T > t456*pow(rho, power2) ) {
					/* test T against (T67 * T78)**.5 */
					if (( T < t678*pow(rho, power3) ) || ((( rho <= 1.0E10 ) && ( T < 1.0E4 )))) {
						/* disjoint opacity laws for 5, 6, and 7 */
						real o5 = bk5*density23*ts42*ts4;
						real o6 = bk6*density13*ts48*ts42;
						real o7 = bk7*rho/(ts42*sqrt(ts4));

						/* parameters used for smoothing */
						real o6an = o6*o6;
						real o7an = o7*o7;

						/* smoothed and continuous opacity law for regions 5, 6, and 7 */
						kappacgs = pow(pow((o6an*o7an/(o6an + o7an)), 2.0) + pow((o5/(1.0 + pow((ts4/(1.1*pow(rho, 0.04762))), 10.0))), 4.0) , 0.25);
					} else {
						/* disjoint opacity laws for 7 and 8 */
						real o7 = bk7*rho/(ts42*sqrt(ts4));
						real o8 = bk8;

						/* parameters used for smoothing */
						real o7an = o7*o7;
						real o8an = o8*o8;

						/* smoothed and continuous opacity law for regions 7 and 8 */
						kappacgs = pow(o7an*o7an + o8an*o8an, 0.25);
					}
				} else {
					/* disjoint opacity laws for 3, 4, and 5 */
					real o3 = bk3*sqrt(ts4);
					real o4 = bk4*rho/(ts48*ts48*ts48);
					real o5 = bk5*density23*ts42*ts4;

					/* parameters used for smoothing */
					real o4an = pow(o4, 4.0);
					real o3an = pow(o3, 4.0);

					/* smoothed and continuous opacity law for regions 3, 4, and 5 */
					kappacgs = pow((o4an*o3an/(o4an + o3an)) + pow(o5/(1.0 + 6.561E-5/ts48*1.0E2*density23), 4.0), 0.25);
				}
			} else {
				/* different powers of temperature */
				real t2 = T*T;
				real t4 = t2*t2;
				real t8 = t4*t4;
				real t10 = t8*t2;

				/* disjoint opacity laws */
				real o1 = ak1*t2;
				real o2 = ak2*T/t8;
				real o3 = ak3*sqrt(T);

				/* parameters used for smoothing */
				real o1an = o1*o1;
				real o2an = o2*o2;

				/* smoothed and continuous opacity law for regions 1, 2, and 3 */
				kappacgs = pow(pow(o1an*o2an/(o1an + o2an), 2.0) + pow(o3/(1.0 + 1.0E22/t10), 4.0), 0.25);
			}
		}
	} else {
		if ( RKappa->Tmax[max_i-1] == 0 ) {
			tmax_max = cons[max_i-1]*pow(rho, expn[max_i-1]);
		}

		for (i = 0; i < max_i; i++) {
			if (T > tmax_max) {
				i = max_i;
				break;
			}
			if (RKappa->Tmax[i] == 0) {
				tmax = cons[i]*pow(rho, expn[i]);
			} else {
				tmax = RKappa->Tmax[i];
			}
			if (T <= tmax) {
				break;
			}
		}
		kappacgs = kappa0[i]*pow(rho, alpha[i])*pow(T, beta[i]);
	}

	// Output
	return kappacgs*MCGS/LCGS/LCGS;
}
