#include "mp.h"


#define LITTLE_NUMBER_FIX  1.0E-15

extern boolean RadCooling, Irradiation, RadTransport, VarDiscHeight, BinaryOn;
extern boolean InnerBCCons, InnerBCGrad, InnerBCExtrap, OuterBCCons, OuterBCGrad, OuterBCExtrap;
extern boolean OptimalOmega, Relative_Diff, Relative_Max, Residual_Diff, Residual_Max, Write_Coeffs;
extern boolean NoCFL, RadiationDebug;
extern boolean LinPap1985_Opacity, BellLin1994_Opacity;
extern boolean RayTracingHeating, ExplicitRayTracingHeating, ImplicitRayTracingHeating;
extern real StarTaper;
extern boolean	BitschSKappa, ImplicitRadiative;
extern boolean	ChebyshevOmega;
extern boolean	PowerLaw_Opacity;
extern boolean	OpacitySmoothing;
extern boolean	ExplicitRadTransport;

static real NormAv = 0;
static int IterAv = 0, counter = 1, counter2 = 1;;

// static real *SendInnerBoundary;
// static real *SendOuterBoundary;
// static real *RecvInnerBoundary;
// static real *RecvOuterBoundary;
// static int allocated_radcom = 0;
static int toggle_qirrrt = 0;

void ComputeDiscHeight (bsys)
	// Input
	BinarySystem *bsys;
{
	// Declaration
	int i, j, l, s, nb, ns, nr;
	real *xstar, *ystar, *mb, xc, yc, dist, angle, smooth;
	real *height, *cs;
	real average, csiso, drx, dry, A, q[2], dist_smooth;

	// Assignment of Fields etc.
	height = DiscHeight->Field;
	cs = SoundSpeed->Field;
	ns = DiscHeight->Nsec;
	nr = DiscHeight->Nrad;

	// Function
	if ( VarDiscHeight ) {
		if ( BinaryOn ) {
			xstar = bsys->x;
			ystar = bsys->y;
			nb = bsys->nb;
			mb = bsys->mass;
			drx = xstar[0] - xstar[1];
			dry = ystar[0] - ystar[1];
			q[0] = mb[0]/mb[1];
			q[1] = 1.0/q[0];
			A = sqrt(drx*drx + dry*dry);
#pragma omp parallel for private(average, j, l, angle, xc, yc, smooth, dist, MeanDiscHeight, csiso)
			for (i = 0; i < nr; i++) {
				average = 0.0;
				for (j = 0; j < ns; j++) {
					l = j + i*ns;
					angle = (real)j/(real)ns*2.0*PI;
					xc = Rmed[i]*cos(angle);
					yc = Rmed[i]*sin(angle);
					height[l] = 0.0;
					for (s = 0; s < nb; s++) {
						smooth = (0.49*pow(q[s], 2.0/3.0))/(0.6*pow(q[s], 2.0/3.0) + log(1.0 + pow(q[s], 1.0/3.0)))*A;
						smooth *= smooth;
						// smooth = 0.0;
						dist = (xc - xstar[s])*(xc -xstar[s]) + (yc - ystar[s])*(yc - ystar[s]);
						dist_smooth = sqrt(dist + smooth);
						csiso = cs[l]/sqrt(ADIABATICINDEX);
						height[l] += G*mb[s]/csiso/csiso/dist_smooth/dist_smooth/dist_smooth;
					}
					height[l] = pow(height[l], -0.5);
					average += height[l];
				}
				MeanDiscHeight[i] = average/(real)ns;
			}
		} else {
#pragma omp parallel for private(average, j, l, MeanDiscHeight)
			for (i = 0; i < nr; i++) {
				average = 0.0;
				for (j = 0; j < ns; j++) {
					l = j+i*ns;
					height[l] = cs[l]*sqrt(pow(Rmed[i], 3.0)/G/1.0); // M = mass of stellar system = 1.0
					average += height[l];
				}
				MeanDiscHeight[i] = average/(real)ns;
			}
		}
	} else {
#pragma omp parallel for private(j, l, average, MeanDiscHeight)
		for (i = 0; i < nr; i++) {
			MeanDiscHeight[i] = ASPECTRATIO*Rmed[i];
			for (j = 0; j < ns; j++) {
				l = j+i*ns;
				height[l] = MeanDiscHeight[i];
			}
		}
	}

	// Debug
	if ( RadiationDebug ) {
		int check_neg = 1;
    	int check_zero = 1;
		CheckField(DiscHeight, check_neg, check_zero, "ComputeDiscHeight");
	}
}

real compute_varheight_smoothing (xp, yp)
	// Input
	real xp, yp;
{
	// Declaration
	real dist, smoothing;

	// Function
	dist = sqrt(xp*xp + yp*yp);
	if (( dist > Rmed[NRAD-1] ) || ( dist < Rmed[0] ))
		return 0.0001*THICKNESSSMOOTHING;

	smoothing = compute_aspectratio(xp, yp);
	smoothing = smoothing * THICKNESSSMOOTHING * pow(dist, 1.0+FLARINGINDEX);

	// Debug
	if ( RadiationDebug ) {
		int check_neg = 1;
    	int check_zero = 1;
		int flag = 0, globalFlag = 0;
		int foo = CheckValue(smoothing, check_neg, check_zero);
		if (foo > flag) {
			flag = foo;
		}

		if ( globalFlag != 0 ) {
			printf("Error: Non-normal value in variable disc height smoothing. Exiting.\n");
			MPI_Finalize();
			exit(globalFlag);
		}
	}
	// Output
	return smoothing;
}

real compute_aspectratio (x, y)
	// Input
	real x, y;
{	
	// Declaration
	int i, i2, j, j2, l1, l2, l3, l4, ns;
	real angle, dangle, ang1, ang2,  dist;
	real fr1, fr2, H, HoverR;

	// Assignment
	ns = DiscHeight->Nsec;

	// Constants
	dangle = 2.0*PI/(real)ns;

	// Function
	dist = sqrt(x*x + y*y);
	i = 0;
	while ( Rmed[i] <= dist ) {
		i++;
	}
	i2 = i;
	i = i-1;
	
	j = 0;
	angle = atan2(y,x) + PI;
	if ( angle <= 0.5*dangle ) {
		j = ns-1;
		j2 = 0;
		ang1 = -0.5*dangle;
		ang2 = 0.5*dangle;
	} else if ( angle >= (2.0*PI)-(0.5*dangle) ) {
		j = ns-1;
		j2 = 0;
		ang1 = (2.0*PI)-(0.5*dangle);
		ang2 = (2.0*PI)+(0.5*dangle);
	} else {
		while ( dangle*((real)j + 0.5) < angle )
			j++;
		j2 = j;
		j = j2 - 1;
		ang1 = dangle*((real)j + 0.5);
		ang2 = dangle*((real)j2 + 0.5);
	}

	l1 = j+i*ns;
	l2 = j2+i*ns;
	l3 = j+i2*ns;
	l4 = j2+i2*ns;

	fr1 = (DiscHeight->Field[l1])*(ang2 - angle)/(dangle) + (angle - ang1)/(dangle)*(DiscHeight->Field[l2]);
	fr2 = (DiscHeight->Field[l3])*(ang2 - angle)/(dangle) + (angle - ang1)/(dangle)*(DiscHeight->Field[l4]);
	H = fr1*(Rmed[i2] - dist)/(Rmed[i2] - Rmed[i]) + fr2*(dist - Rmed[i])/(Rmed[i2] - Rmed[i]);
	HoverR = H/dist;

	// Debug
	if ( RadiationDebug ) {
		int flag = 0, globalFlag = 0;
		int check_neg = 1;
    	int check_zero = 1;
		int foo = CheckValue(HoverR, check_neg, check_zero);
		if (foo > flag) {
			flag = foo;
		}

		if ( globalFlag != 0 ) {
			printf("Error: Non-normal value in compute_aspectratio (H/r). Exiting.\n");
			MPI_Finalize();
			exit(globalFlag);
		}
	}

	// Output
	return HoverR;
}

// void ComputeQterms (Sigma, bsys)
// 	// Input
// 	PolarGrid *Sigma;
// 	BinarySystem *bsys;
// {
// 	// Declaration
// 	int i, j, l, s, ns, nr, nb;
// 	real *dens, *T, *qminus, *qirr, *H, *rkappa, *Rs, *Tstar, Ts;
// 	real *xs, *ys;
// 	real angle, x, y, dist, qirrterm;
// 	real 	*tau, *taueff, term1, term2, Wg, Tdisc4_eff;
// 	real disc_albedo, c2, constant1, constant2, tau_min;
// 	boolean *central_source;

// 	// Constants
// 	disc_albedo = 0.5;
// 	c2 = 0.5;
// 	constant1 = 2.0*STEFANK;
// 	tau_min = 0;

// 	// Assignment
// 	ns = Sigma->Nsec;
// 	nr = Sigma->Nrad;
// 	dens = Sigma->Field;
// 	T = Temperature->Field;
// 	H = DiscHeight->Field;
// 	rkappa = RKappaval->Field;
// 	if ( RadCooling ) {
// 		qminus = Qminus->Field;
// 	}
// 	if ( Irradiation ) {
// 		nb = IrrSources->nb;
// 		Rs = IrrSources->Rstar;
// 		Tstar = IrrSources->Tstar;
// 		central_source = IrrSources->CentralSource;
// 		qirr = Qirr->Field;
// 		if (nb > 1) {
// 			xs = bsys->x;
// 			ys = bsys->y;
// 		}
// 	}
// 	tau = OpticalDepth->Field;
// 	taueff = OpticalDepthEff->Field;

// 	// Function
// 	#pragma omp parallel for private(j, l, angle, Ts, constant2, dist, angle, x, y, term1, term2, Wg, qirrterm, Tdisc4_eff)
// 	for (i = 0; i < nr; i++) {
// 		for (j = 0; j < ns; j++) {
// 			l = j+i*ns;
// 			tau[l] = c2*rkappa[l]*dens[l];
// 			taueff[l] = (0.375*tau[l]) + (1.0/(4.0*tau[l] + tau_min)) + 0.866;

// 			if ( Irradiation ) {
// 				qirr[l] = 0.0;
// 				for (s = 0; s < nb; s++) {
// 					Ts = StarTaper*Tstar[s];
// 					constant2 = constant1*(1.0 - disc_albedo)*pow(Ts*Ts*Rs[s], 2.0);
// 					if (central_source[s] == YES) {
// 						dist = Rmed[i];	
// 					} else {
// 						angle = GlobalTheta[j];
// 						x = Rmed[i]*cos(angle);
// 						y = Rmed[i]*sin(angle);
// 						dist = sqrt((x - xs[s])*(x - xs[s]) + (y - ys[s])*(y - ys[s]));
// 					}
// 					term1 = 0.4*Rs[s]/dist;
// 					term2 = 0.286*H[l]/dist;
// 					term2 = 0.0;
// 					Wg = term1 + term2;
// 					if ( dist == 0.0 ) {
// 						qirrterm = 0.0;
// 					} else {
// 						qirrterm = constant2*pow(dist, -2.0)*Wg/taueff[l];
// 					}
// 					qirr[l] += qirrterm;
// 					if (isnan(qirr[l]) != 0) {
// 						printf("s = %d, nb = %d, StarTaper = %f, Tstar[s] = %f", s, nb, StarTaper, Tstar[s]);
// 					}
// 				}
// 			}
// 			if ( RadCooling ) {
// 				Tdisc4_eff = pow(T[l], 4.0)/taueff[l];
// 				//constant1 = 2*STEFANK
// 				qminus[l] = constant1*Tdisc4_eff;
// 			}
// 		}
// 	}

// 	// Debug
// 	if ( RadiationDebug ) {
// 		if ( Irradiation ) {
// 			int check_neg = 1;
//     		int check_zero = 1;
// 			CheckField(Qirr, check_neg, check_zero, "ComputeQterms");
// 		}
// 		if ( RadCooling ) {
// 			int check_neg = 1;
//     		int check_zero = 1;
// 			CheckField(Qminus, check_neg, check_zero, "ComputeQterms");
// 		}
// 	}
// }

void ComputeQirr (Sigma, bsys)
	// Input
	PolarGrid *Sigma;
	BinarySystem *bsys;
{
	// Declaration
	int i, j, l, s, ns, nr, nb;
	real *dens, *qirr, *H, *rkappa, *Rs, *Tstar, Ts;
	real *xs, *ys;
	real angle, x, y, dist, qirrterm;
	real 	*tau, *taueff, term1, term2, Wg;
	real disc_albedo, c2, constant1, constant2, tau_min;
	boolean *central_source;

	// Constants
	disc_albedo = 0.5;
	c2 = 0.5;
	constant1 = 2.0*STEFANK;
	tau_min = 0;

	// Assignment
	ns = Sigma->Nsec;
	nr = Sigma->Nrad;
	dens = Sigma->Field;
	H = DiscHeight->Field;
	rkappa = RKappaval->Field;
	nb = IrrSources->nb;
	Rs = IrrSources->Rstar;
	Tstar = IrrSources->Tstar;
	central_source = IrrSources->CentralSource;
	qirr = Qirr->Field;
	if ( nb > 1 ) {
		xs = bsys->x;
		ys = bsys->y;
	}
	tau = OpticalDepth->Field;
	taueff = OpticalDepthEff->Field;

	// Function
	#pragma omp parallel for private(j, l, angle, Ts, constant2, dist, angle, x, y, term1, term2, Wg, qirrterm)
	for (i = 0; i < nr; i++) {
		for (j = 0; j < ns; j++) {
			l = j+i*ns;
			tau[l] = c2*rkappa[l]*dens[l];
			taueff[l] = (0.375*tau[l]) + (1.0/(4.0*tau[l] + tau_min)) + 0.866;
			qirr[l] = 0.0;
			for (s = 0; s < nb; s++) {
				Ts = StarTaper*Tstar[s];
				constant2 = constant1*(1.0 - disc_albedo)*pow(Ts*Ts*Rs[s], 2.0);
				if (central_source[s] == YES) {
					dist = Rmed[i];	
				} else {
					angle = GlobalTheta[j];
					x = Rmed[i]*cos(angle);
					y = Rmed[i]*sin(angle);
					dist = sqrt((x - xs[s])*(x - xs[s]) + (y - ys[s])*(y - ys[s]));
				}
				term1 = 0.4*Rs[s]/dist;
				term2 = 0.286*H[l]/dist;
				term2 = 0.0;
				Wg = term1 + term2;
				if ( dist == 0.0 ) {
					qirrterm = 0.0;
				} else {
					qirrterm = constant2*pow(dist, -2.0)*Wg/taueff[l];
				}
				qirr[l] += qirrterm;
			}
		}
	}

	// Debug
	if ( RadiationDebug ) {
		int check_neg = 1;
    int check_zero = 1;
		CheckField(Qirr, check_neg, check_zero, "ComputeQplus");
	}
}

void ComputeQminus (Sigma)
	// Input
	PolarGrid *Sigma;
{
	// Declaration
	int i, j, l, ns, nr;
	real *dens, *T, *qminus, *rkappa;
	real 	*tau, *taueff, Tdisc4_eff;
	real c2, constant1, tau_min;

	// Assignment
	ns = Sigma->Nsec;
	nr = Sigma->Nrad;
	dens = Sigma->Field;
	T = Temperature->Field;
	rkappa = RKappaval->Field;
	qminus = Qminus->Field;
	tau = OpticalDepth->Field;
	taueff = OpticalDepthEff->Field;

	// Constants
	c2 = 0.5;
	constant1 = 2.0*STEFANK;
	tau_min = 0;

	// Function
	#pragma omp parallel for private(j, l, Tdisc4_eff)
	for (i = 0; i < nr; i++) {
		for (j = 0; j < ns; j++) {
			l = j+i*ns;
			tau[l] = c2*rkappa[l]*dens[l];
			taueff[l] = (0.375*tau[l]) + (1.0/(4.0*tau[l] + tau_min)) + 0.866;
			Tdisc4_eff = pow(T[l], 4.0)/taueff[l];
			qminus[l] = constant1*Tdisc4_eff;
		}
	}

	// Debug
	if ( RadiationDebug ) {
		int check_neg = 1;
    int check_zero = 1;
		CheckField(Qminus, check_neg, check_zero, "ComputeQminus");
	}
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

	// Constants
	BETA_POWER = -0.5;
	GAMMA_POWER = 0.0;
	kap0 = 3.0E1;
	c1 = 0.5;
	cgsconstant1 = c1*MCGS/pow(LCGS, 3.0);
	cgsconstant2 = TEMPCGS;
	
	// Assignment
	ns = Sigma->Nsec;
	nr = Sigma->Nrad;
	dens = Sigma->Field;
	temp = Temperature->Field;
	H = DiscHeight->Field;
	rkappa = RKappaval->Field;

	// Function
#pragma omp parallel for private(j, l, rhocgs, tempcgs)
	for (i = 0; i < nr; i++) {
		for (j = 0; j < ns; j++) {
			l = j+i*ns;
			
			rhocgs = cgsconstant1*dens[l]/H[l];
			tempcgs = cgsconstant2*temp[l];
			if ( PowerLaw_Opacity ) {
				// Power-Law Opacity 
				rkappa[l] = kap0 * pow(tempcgs, BETA_POWER) * pow(Rmed[i], GAMMA_POWER);
			} else {
				// Convert density and Temperature from code units to cgs units
				rhocgs = cgsconstant1*dens[l]/H[l];
				// Find Rosseland opacity from opacity table look-up
				rkappa[l] = CoolingRegime(rhocgs, tempcgs);
			}
			// rkappa[l] = 1.5E7;
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
	real *kappa0, *alpha, *beta, *cons, *expn, *tmin;
	real  kappacgs, tmax_max, tmax;

	// Assignment
	max_i = RKappa->Ni;
	kappa0 = RKappa->Kappa0;
	alpha = RKappa->Alpha;
	beta = RKappa->Beta;
	cons = RKappa->Constant;
	expn = RKappa->Exponent;
	tmin = RKappa->Tmin;

	// Function
	if ( OpacitySmoothing ) {
		if ( LinPap1985_Opacity ) {
	    	/*if (density < 1e-10)
	        	density = 1e-10;*/
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
		if (RKappa->Tmax[max_i-1] == 0) {
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

void SubStep4 (gas_density, gas_energynew, timestep)
	// Input
  PolarGrid *gas_density;
  PolarGrid *gas_energynew;
  real timestep;
{
	// Declaration
	int i, ii, j, l, nr, ns;
	int lim, lip, ljm, ljp;
	int iteration;
	real *T, *Tguess, *Tguess_old;
	real *B, *U1, *U2, *U3, *U4, *residual, *Rij;
	real qirrrt = 0.0, qirr = 0.0, qminus = 0.0;
	real norm_tmp[2], norm, cvfac, tol;

	// Assignment
	ns = gas_density->Nsec;
	nr = gas_density->Nrad;
	T = Temperature->Field;
	B = Barr->Field;
	U1 = U1arr->Field;
	U2 = U2arr->Field;
	U3 = U3arr->Field;
	U4 = U4arr->Field;
	residual = Residual->Field;
	Rij = TempSourcesSinks->Field;
	Tguess = TempGuess->Field;
	Tguess_old = TempGuess_old->Field;

	// Constants
	cvfac = timestep*MU*(ADIABATICINDEX-1.0)/R;

	if ( CPU_Rank == 0 ) {
		i = 0;
		for (j = 0; j < ns; j++) {
			l = j+i*ns;
			Tguess[l] = TINNER;
			T[l] = TINNER;
		}
	}
	if ( CPU_Rank == CPU_Highest ) {
		i = nr-1;
		for (j = 0; j < ns; j++) {
			l = j+i*ns;
			Tguess[l] = TOUTER;
			T[l] = TOUTER;
		}
	}

	// Function
	ComputeRadTransCoeffs(gas_density, timestep);

	if ( ExplicitRadTransport ) {
		SubStep4_Explicit(gas_density, gas_energynew, timestep);
	} else {
		norm_tmp[0] = 0.0;
	 	norm_tmp[1] = 0.0;

	# pragma omp parallel for private(j,l, qirrrt, qirr, qminus)
	 	for (i = 0; i < nr; i++) {
	 		if ( ChebyshevOmega ) {
	 			omegaOpt[i] = 1.0;
	 		}
	 		for (j = 0; j < ns; j++) {
	 			l = j+i*ns;
	 			if ( ImplicitRadiative ) {
	 				if  (RayTracingHeating ) {
		  			qirrrt = QirrRT->Field[l];
		  		}
		  		if ( Irradiation ) {
		  			qirr = Qirr->Field[l];
		  		}
		  		if ( RadCooling ) {
		  			qminus = Qminus->Field[l];
		  		}
		  	}
		  	Rij[l] = T[l] + qirrrt + cvfac*(qirr - qminus)/gas_density->Field[l];
	 			Tguess[l] = T[l];
	 			Tguess_old[l] = 0.0;
	 			norm_tmp[1] += fabs(Rij[l]);
	 		}
	 	}

	 	MPI_Allreduce(&norm_tmp[1], &norm_tmp[1], 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	 	tol = TOLERANCE*norm_tmp[1];
	 	norm = 2.0*tol;
	 	iteration = 0;
	 	// norm = 2.0*TOLERANCE;

	 	// while (( norm > TOLERANCE ) && ( iteration < MAXITERATIONS )) {
	 	while (( norm > tol ) && ( iteration < MAXITERATIONS )) {
	 		norm = 0.0;
	 		norm_tmp[0] = 0.0;
	 		// if ( Relative_Max )
	 		// 	norm_tmp[1] = 0.0;

	 		//Apply BC
	 		if ( CPU_Rank == 0 ) {
	 			i = 1;
	#pragma omp parallel for private(l, lim, lip)
	 			for (j = 0; j < ns; j++) {
	 				l = j+i*ns;
	 				lim = l-ns;
	 				lip = l+ns;
	 				Tguess[lim] = ApplyInnerBC(Tguess[l], Tguess[lip]);
	 			}
	 		}
	 		if ( CPU_Rank == CPU_Number-1 ) {
	 			i = nr-2;
	#pragma omp parallel for private(l, lim, lip)
	 			for (j = 0; j < ns; j++) {
	 				l = j+i*ns;
	 				lim = l-ns;
	 				lip = l+ns;
	 				Tguess[lip] = ApplyOuterBC(Tguess[l], Tguess[lim]);
	 			}
	 		}
	 		
			// Black-loop over all even numbered cells
	// Commented out .# pragma omp parallel for private (j, l, lim, lip, ljm, ljp, Qirrrt, Tguess_old)
	# pragma omp parallel for private (ii, j, l, lim, lip, ljm, ljp, Qirrrt, Tguess_old, norm_tmp)
	 		for (i = One_or_active; i < MaxMO_or_active; i++) {
	 			ii = i + IMIN;
	 			for (j = 0; j < ns; j++) {
	 				l = j+(ii)*ns;
	 				if ( (ii % 2) == 0 ) {
	 					l = l+1;
	 				}
	 				if ( (l % 2) == 0 ) {
	 					l = j+i*ns;
	 					lim = l-ns;
	 					lip = l+ns;
	 					if (j == 0) {
	  						ljm = ns*(i+1)-1;
	 					} else {
	  						ljm = l-1;
	 					}
	  					if (j == ns-1) {
	  						ljp = i*ns;
	  					} else {
	  						ljp = l+1;
	  					}
	  				
	  					Tguess_old[l] = Tguess[l];
	  					residual[l] = U1[l]*Tguess[lim] + U2[l]*Tguess[lip] + U3[l]*Tguess[ljm] + U4[l]*Tguess[ljp] + B[l]*Tguess[l] - Rij[l];
	  					Tguess[l] = Tguess_old[l] - omegaOpt[i]*residual[l]/B[l];
	  					norm_tmp[0] += fabs(residual[l]);
	  					// Internal Debug
	  					// if ( isnan(Tguess[l]) ) {
	  					// 	fprintf (stderr, "Tguess at cell %dx%d is NaN after iteration %d.\n", (i+IMIN), j, iteration);
	  					// }

		  				/* norm  = |x^{k} - x^{k-1}| / |x^{k}| --> norm_relative_l2 */
		  		// 		if ( Relative_Diff ) {
		  		// 			norm_tmp[0] += pow(Tguess_old[l] - Tguess[l], 2.0)/pow(Tguess_old[l], 2.0);
		  		// 		}
		  		// 		/* norm = max(|x^{k} - x^{k-1}| / |x^{k}|) --> norm_relative_max */
		  		// 		if ( Relative_Max ) {
		  		// 			norm_tmp[0] = fmax(norm_tmp[0], fabs(Tguess_old[l] - Tguess[l]));
		  		// 			norm_tmp[1] = fmax(norm_tmp[1], fabs(Tguess_old[l]));
		  		// 		}
			  	// 		/* norm = |r^{k}| / |b{k}| --> norm_residual_l2 */
			  	// 		if ( Residual_Diff ) {
		  		// 			norm_tmp[0] += pow(residual[l], 2.0);
		  		// 			norm_tmp[1] += pow(Rij[l], 2.0);
						// }
		  		// 		/* norm = max(|r^{k}| / |b{k}|) --> norm_residual_max */
		  		// 		if ( Residual_Max ) {
		  		// 			norm_tmp[0] = fmax(norm_tmp[0], fabs(residual[l]));
		  		// 			norm_tmp[1] = fmax(norm_tmp[1], fabs(Rij[l]));
		  		// 		}
	 				}
	 			}
	 		}

	 		if ( ChebyshevOmega ) {
	#pragma omp parallel for
	 			for (i = 0; i < nr; i++) {
	 				if ( iteration == 0 ) {
	 					omegaOpt[i] = 1.0/(1.0 - 0.5*rhoJac[i]*rhoJac[i]);
	 				} else {
	 					omegaOpt[i] = 1.0/(1.0 - 0.25*rhoJac[i]*rhoJac[i]*omegaOpt[i]);
	 				}
	 			}
	 		}

			//Overlap zones communication between neighbouring processors
	 		CommunicateFieldBoundaries(TempGuess);

			// Red-loop over all odd numbered cells
	// Commented out .# pragma omp parallel for private (j, l, lim, lip, ljm, ljp, Qirrrt, Tguess_old)
	# pragma omp parallel for private (j, l, lim, lip, ljm, ljp, Qirrrt, Tguess_old, norm_tmp)

	 		for (i = One_or_active; i < MaxMO_or_active; i++) {
	 			ii = i + IMIN;
	 			for (j = 0; j < ns; j++) {
	 				l = j+(ii)*ns;
	 				if ( (ii % 2) == 0 ) {
	 					l = l+1;
	 				}
	 				if ( (l % 2) != 0 ) {
	 					l = j+i*ns;
	 					lim = l-ns;
	 					lip = l+ns;
	 					if ( j == 0 ) {
	  						ljm = ns*(i+1)-1;
	 					} else {
	  						ljm = l-1;
	  					}
	  
		  				if ( j == ns-1 ) {
		  					ljp = i*ns;
		  				} else {
		  					ljp = l+1;
		  				}

		  				Tguess_old[l] = Tguess[l];
		  				residual[l] = U1[l]*Tguess[lim] + U2[l]*Tguess[lip] + U3[l]*Tguess[ljm] + U4[l]*Tguess[ljp] + B[l]*Tguess[l] - Rij[l];
		  				Tguess[l] = Tguess_old[l] - omegaOpt[i]*residual[l]/B[l];
		  				norm_tmp[0] += fabs(residual[l]);
		  				
		  				// Internal Debug
		  				// if ( isnan(Tguess[l]) ) {
		  				// 	fprintf (stderr, "Tguess at cell %dx%d is NaN after iteration %d.\n", (i+IMIN), j, iteration);
		  				// }

		  				/* norm  = |x^{k} - x^{k-1}| / |x^{k}| --> norm_relative_l2 */
		  		// 		if ( Relative_Diff ){
		  		// 			norm_tmp[0] += pow(Tguess_old[l] - Tguess[l], 2.0)/pow(Tguess_old[l], 2.0);
		  		// 		}
		  		// 		/* norm = max(|x^{k} - x^{k-1}| / |x^{k}|) --> norm_relative_max */
		  		// 		if ( Relative_Max ) {
		  		// 			norm_tmp[0] = fmax(norm_tmp[0], fabs(Tguess_old[l] - Tguess[l]));
		  		// 			norm_tmp[1] = fmax(norm_tmp[1], fabs(Tguess_old[l]));
		  		// 		}
			  	// 		/* norm = |r^{k}| / |b{k}| - norm_residual_l2 */
			  	// 		if ( Residual_Diff ) {
		  		// 			norm_tmp[0] += pow(residual[l], 2.0);
		  		// 			norm_tmp[1] += pow(Rij[l], 2.0);
						// }
		  		// 		/* norm = max(|r^{k}| / |b{k}|) - norm_residual_max */
		  		// 		if ( Residual_Max ) {
		  		// 			norm_tmp[0] = fmax(norm_tmp[0], fabs(residual[l]));
		  		// 			norm_tmp[1] = fmax(norm_tmp[1], fabs(Rij[l]));
		  		// 		}
	 				}
	 			}
	 		}

	 		if ( ChebyshevOmega ) {
	#pragma omp parallel for
	 			for ( i = 0; i < nr; i++) {
	 				omegaOpt[i] = 1.0/(1.0 - 0.25*rhoJac[i]*rhoJac[i]*omegaOpt[i]);
	 			}
	 		}

	 		MPI_Allreduce(&norm_tmp[0], &norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	 		// /* norm  = |x^{k} - x^{k-1}| / |x^{k}| - norm_relative_l2 */
	 		// if ( Relative_Diff ) { 
		 	// 	MPI_Allreduce(&norm_tmp[0], &norm_tmp[0], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		 	// 	norm = sqrt(norm_tmp[0]);
		 	// }
	 		// /* norm = max(|x^{k} - x^{k-1}| / |x^{k}|) - norm_relative_max */
	 		// if ( Relative_Max ) {
		 	// 	MPI_Allreduce(&norm_tmp[0], &norm_tmp[0], 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		 	// 	MPI_Allreduce(&norm_tmp[1], &norm_tmp[1], 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		 	// 	norm = norm_tmp[0]/norm_tmp[1];
		 	// }
	 		// /* norm = |r^{k}| / |b{k}| - norm_residual_l2 */
	 		// if ( Residual_Diff ) {
		 	// 	MPI_Allreduce(&norm_tmp[0], &norm_tmp[0], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		 	// 	norm_tmp[0] = sqrt(norm_tmp[0]);
		 	// 	MPI_Allreduce(&norm_tmp[1], &norm_tmp[1], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	 		// 	norm_tmp[1] = sqrt(norm_tmp[1]);
		 	// 	norm = norm_tmp[0]/norm_tmp[1];
		 	// }
	 		// /* norm = max(|r^{k}| / |b{k}|) - norm_residual_max */
	 		// if ( Residual_Max ) {
		 	// 	MPI_Allreduce(&norm_tmp[0], &norm_tmp[0], 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		 	// 	MPI_Allreduce(&norm_tmp[1], &norm_tmp[1], 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
		 	// 	norm = norm_tmp[0]/norm_tmp[1];
		 	// }
		 	// masterprint("Norm = %g\n", norm_tmp[0]);
			//Overlap zones communication between neighbouring processors
	 		CommunicateFieldBoundaries(TempGuess);

	  		iteration++;
	  		if ( RadiationDebug ) {
	 			int check_neg = 0;
	    		int check_zero = 0;
	    		char str[256];
	    		sprintf(str, "SubStep4 - Iteration %d", iteration);
				CheckField(TempGuess, check_neg, check_zero, str);
	 		}
	 	}

	 	NormAv += norm;
	 	IterAv += iteration;
	 	counter++;

	 	if ( RadiationDebug ) {
	 		int check_neg = 1;
	    int check_zero = 1;
			CheckField(TempGuess, check_neg, check_zero, "SubStep4");
	 	}

	# pragma omp parallel for private(j, l)
	 	for (i = 0; i < nr; i++) {
	 		for (j = 0; j < ns; j++) {
	 			l = j+i*ns;
	 			T[l] = Tguess[l];
	 			if (( ExplicitRayTracingHeating ) && (( i >= One_or_active ) || ( i < MaxMO_or_active ))) {
	 				T[l] += QirrRT->Field[l];
	 			}
	 		}
	 	}

	 	// if ( RayTracingHeating ) {
	 	// 	ToggleQirrRT(gas_density, timestep);
	 	// }

	 	// Debug
	 	if (( iteration >= MAXITERATIONS ) && ( norm > TOLERANCE )) {
	 		fprintf (stderr, "Max number of iterations has been reached without convergence. Exiting.\n");
			WriteDiskPolar (gas_density, 9999);    /* We write the HD arrays */
			WriteDiskPolar (Temperature, 9999);
			merge(9999);
	    prs_exit(1);
	 	}

	 	if ( RadiationDebug ) {
	 		int check_neg = 1;
	    	int check_zero = 1;
			CheckField(Temperature, check_neg, check_zero, "SubStep4");
	 	}

	 	// Output
	 	ComputeNewEnergyField(gas_density, gas_energynew);
 	}
}

void SubStep4_Explicit_Irr (gas_density, gas_energy, timestep)
	// Input
	PolarGrid *gas_density;
	PolarGrid *gas_energy;
  real timestep;
{
	// Declaration
	int i, j, l, ns, nr;
	real *temperature, *Qirrrt, *energy, *density;

	// Assignment
	nr = QirrRT->Nrad;
	ns = QirrRT->Nsec;
	Qirrrt = QirrRT->Field;
	temperature = Temperature->Field;
	energy = gas_energy->Field;
	density = gas_density->Field;


	// Function
# pragma omp parallel for private (j, l)
 	for (i = 0; i < nr; i++) {
 		for (j = 0; j < ns; j++) {
 			l = j+i*ns;
 			temperature[l] += timestep*Qirrrt[l];
 			energy[l] = CV*density[l]*temperature[l];
 		}
 	}
}

void ComputeRadTransCoeffs(gas_density, dt)
	// Input
	PolarGrid *gas_density;
	real dt;
{	
	// Declaration
	int nr, ns, i, im, ip, jm, jp;
	int j, l, lim, lip, ljm, ljp;
	real *sigma, *kappa, *height, *T, rho;
	real *D, *B, *U1, *U2, *U3, *U4;
	real Dip, Dim, Djm, Djp;
	real gij, gj, gip, gim;
	real cvfac, c1;
	real dr1, dr2;

  // Assignment
	nr = gas_density->Nrad;
	ns = gas_density->Nsec;
	sigma = gas_density->Field;
	kappa = RKappaval->Field;
	height = DiscHeight->Field;
	D = Darr->Field;
	B = Barr->Field;
	U1 = U1arr->Field;
	U2 = U2arr->Field;
	U3 = U3arr->Field;
	U4 = U4arr->Field;
	T = Temperature->Field;
	real elim, elip, eljm, eljp, el;
	real grad1, grad2, grad;
	real *rfld, *lambda;

	// Constants
	c1 = 0.5;
	if ( ExplicitRadTransport ) {
		cvfac = 1.0/CV;
	} else {
		cvfac = dt/CV;
	}

	// Function
	for (i = (One_or_active - 1); i <= MaxMO_or_active; i++) {
		if ( i == 0 ) {
			im = 0;
			ip = 1;
		} else if ( i == nr-1 ) {
			im = nr-2;
			ip = nr-1;
		} else {
			im = i-1;
			ip = i+1;
		}
		for (j = 0; j < ns; j++) {
			l = j+i*ns;
			if ( j == 0 ) {
				jm = ns - 1;
				jp = j + 1;
			} else if ( j == ns-1 ) {
				jm = j - 1;
				jp = 0;
			} else {
				jm = j - 1;
				jp = j + 1;
			}
			lim = j+im*ns;
			lip = j+ip*ns;
			ljm = jm+i*ns;
			ljp = jp+i*ns;

			// if ( CPU_Rank == CPU_Highest ) {
			// 	if ( j == 0 ) {
			// 		if (i == MaxMO_or_active) {
			// 			printf("i = %d, im = %d, ip = %d, j = %d, jm = %d, jp = %d\n", i, im, ip, j, jm, jp);
			// 			printf("l = %d, lim = %d, lip = %d, ljm = %d, ljp = %d\n", l , lim, lip, ljm, ljp);
			// 			printf("Tinner = %g, Touter = %g\n", TINNER, TOUTER);
			// 			printf("Tlip = %g\n", T[lip]);
			// 		}
			// 	}
			// }

			/* Commented out.
			if ( NoCFL )
				height[l] =  ASPECTRATIO*Rmed[i]; */
			rho = c1*sigma[l]/height[l];
			
			rfld = Rfld->Field;
			lambda = lambdafld->Field;

			elim = pow(T[lim], 4.0);
			elip = pow(T[lip], 4.0);
			eljm = pow(T[ljm], 4.0);
			eljp = pow(T[ljp], 4.0);
			el = pow(T[l], 4.0);

			grad1 = fabs((elip - elim)/(Rmed[ip] - Rmed[im]));
			grad2 = fabs((eljp - eljm)/(2.0*Rmed[i]*DTHETA));
			
			grad = grad1+grad2;
			rfld[l] = grad/(rho*kappa[l]*el);

			// if ( CPU_Rank == CPU_Highest ) {
			// 	if ( i == MaxMO_or_active-1 ) {
			// 		if ( j == 0 ) {
			// 			printf("ergrad = %g, ethgrad = %g\n", grad1, grad2);
			// 			printf("egrad = %g, rho*kappa = %g\n", grad, rho*kappa[l]);
			// 			printf("el = %g, Rfld = %g\n", el, rfld[l]);
			// 		}

			// 	}
			// }

			if ( rfld[l] <= 2.0 ) {
				lambda[l] = 2.0/(3.0 + sqrt(9.0 + (10.0*rfld[l]*rfld[l])));
			} else {
				lambda[l] = 10/((10*rfld[l]) + 9.0 + sqrt(81.0 + (180.0*rfld[l])));
			}
	
			D[l] = 4.0*STEFANK*lambda[l]*pow(T[l], 3.0)/(rho*kappa[l]);

			/* Information required to calculate the diffusion coefficients */
			// D[l] = ComputeD(rho, kappa[l], T[l], i, im, ip, l , lim, lip, ljm, ljp);
		}
	}

	if ( RadiationDebug ) {
		int check_neg = 0;
    int check_zero = 0;
		CheckField(Darr, check_neg, check_zero, "ComputeRadTransCoeffs");
 	}

 	for (i = One_or_active; i < MaxMO_or_active; i++) {
 		im = i - 1;
 		ip = i + 1;
 		gj = 2.0*Rmed[i]*DTHETA*DTHETA;
 		for (j = 0; j < ns; j++) {
 			l = j+i*ns;
			lim = j+im*ns;
			lip = j+ip*ns;

			if ( j == 0 ) {
				jm = ns - 1;
				jp = j + 1;
			} else if ( j == ns-1 ) {
				jm = j - 1;
				jp = 0;
			} else {
				jm = j - 1;
				jp = j + 1;
			}
			ljm = jm+i*ns;
			ljp = jp+i*ns;

			dr1 = Rinf[i] - Rmed[im];
			dr2 = Rmed[i] - Rmed[im];
			Dim = D[lim] + dr1*(D[l] - D[lim])/dr2;
			Dim = D[lim] + (Rinf[i] - Rmed[i-1])*(D[l] - D[lim])/(Rmed[i] - Rmed[i-1]);

			dr1 = Rinf[ip] - Rmed[i];
			dr2 = Rmed[ip] - Rmed[i];
			Dip = D[l] + dr1*(D[lip] - D[l])/dr2;
			Dip = D[l] + (Rinf[i+1] - Rmed[i])*(D[lip] - D[l])/(Rmed[i+1] - Rmed[i]);

			Djp = (D[l] + D[ljp]);
			Djm = (D[ljm] + D[l]);

			dr1 = Rinf[ip] - Rinf[i];
			dr2 = Rmed[i] - Rmed[im];
			gij = cvfac/(sigma[l]*Rmed[i]);
			gim = Rinf[i]*dr1*dr2;
			dr2 = Rmed[ip] - Rmed[i];
			gip = Rsup[i]*dr1*dr2;

			/* Coefficient for T_{i-1j}^{n+1} ie cell lim at next timestep level */
			// U1[l] = gij*Dim/gim;
			// U1[l] = (cvfac/sigma[l])*Rinf[i]*Dim/(Rmed[i]*Rmed[i])
			U1[l] = (cvfac/sigma[l])*(Rinf[i]*Dim)/(Rmed[i]*(Rsup[i] - Rinf[i])*(Rmed[i] - Rmed[i-1]));

			/* Coefficient for T_{i+1j}^{n+1} ie cell lip at next timestep level */
			// U2[l] = gij*Dip/gip;
			U2[l] = (cvfac/sigma[l])*(Rinf[i+1]*Dip)/(Rmed[i]*(Rsup[i] - Rinf[i])*(Rmed[i+1] - Rmed[i]));

			/* Coefficient for T_{ij-1}^{n+1} ie cell ljm at next timestep level */
			if ( NSEC > 1 ) {
				// U3[l] = gij*Djm/gj;
				U3[l] = (cvfac/sigma[l])*Djm/(2.0*Rmed[i]*Rmed[i]*DTHETA*DTHETA);
			} else {
				U3[l] = 0.0;
			}

			/* Coefficient for T_{ij+1}^{n+1} ie cell ljp at next timestep level */
			if ( NSEC > 1 ) {
				// U4[l] = gij*Djp/gj;
				U4[l] = (cvfac/sigma[l])*Djp/(2.0*Rmed[i]*Rmed[i]*DTHETA*DTHETA);
			} else {
				U4[l] = 0.0;
			}

			/* Coefficient for T_{ij}^{n+1} ie cell l at next timestep level */
			B[l] = 1.0 - (U1[l] + U2[l] + U3[l] + U4[l]);
 		}
 	}
	
	// Debug
 	if ( RadiationDebug ) {
 		int check_neg = 0;
    int check_zero = 0;
		CheckField(U1arr, check_neg, check_zero, "ComputeRadTransCoeffs");
		CheckField(U2arr, check_neg, check_zero, "ComputeRadTransCoeffs");
		CheckField(U3arr, check_neg, check_zero, "ComputeRadTransCoeffs");
		CheckField(U4arr, check_neg, check_zero, "ComputeRadTransCoeffs");
		CheckField(Barr, check_neg, check_zero, "ComputeRadTransCoeffs");
 	}
}

real ComputeD (rho, kappa, Temp, i, im, ip, l, lim, lip, ljm, ljp)
	// Input
	real rho, kappa, Temp;
	int i, im, ip, l, lim, lip, ljm, ljp;
{	
	// Declaration
	real  D;
	real *rfld, *lambda;

	// Assignment
	rfld = Rfld->Field;
	lambda = lambdafld->Field;

	// Function

	rfld[l] = ComputeR(rho, kappa, i, im, ip, l, lim, lip, ljm, ljp);
	if (rfld[l] <= 2.0) {
		lambda[l] = 2.0/(3.0 + sqrt(9.0+10.0*rfld[l]*rfld[l]));
	} else {
		lambda[l] = 10/(10*rfld[l] + 9.0 + sqrt(81.0 + 180.0*rfld[l]));
	}
	
	D = 4.0*Acode*Ccode*lambda[l]*pow(Temp, 3.0)/rho/kappa;

	// Output
	return D;
}

real ComputeR (rho, kappa, i, im, ip, l, lim, lip, ljm, ljp)
	// Input
	real rho, kappa;
	int i, im, ip, l, lim, lip, ljm, ljp;
{	
	// Declaration
	real *T;
	real grad1, grad2, numerator, denominator, rfld;
	real E, Eip, Eim, Ejp, Ejm;

	// Assignment
	T = Temperature->Field;

	// Function
	Eip = pow(T[lip], 4.0); /* removed Acode constant from calculation */
	Eim = pow(T[lim], 4.0);
	Ejp = pow(T[ljp], 4.0);
	Ejm = pow(T[ljm], 4.0);
	E	= pow(T[l]  , 4.0);

	grad1 = (Eip - Eim)/(Rmed[ip] - Rmed[im]);
	grad2 = (Ejp - Ejm)/(2.0*DTHETA*Rmed[i]);
	numerator = sqrt(pow(grad1, 2.0) + pow(grad2, 2.0));
	denominator = E;

	rfld = numerator/(denominator*rho*kappa);

	// Output
	return rfld;
}

void ToggleQirrRT (gas_density, dt)
	// Input
	PolarGrid *gas_density;
	real dt;
{	
	// Declaration
	int i, j, l, ns, nr;
	real *sigma, *QRT, cvfac, I;

	// Assignment
	sigma = gas_density->Field;
	QRT = QirrRT->Field;
	nr = gas_density->Nrad;
	ns = gas_density->Nsec;

	// Constants
	cvfac = dt/CV;

	// Function
#pragma omp parallel for private (j, l, I)
	for (i = 0; i < nr; i++) {
		for (j = 0; j < ns; j++) {
			l = j+i*ns;
			I = cvfac/sigma[l];
			if ( QRT[l] != 0.0 ) {
				if ( toggle_qirrrt )
					QRT[l] /= I;
				else
					QRT[l] *= I;
			}
		}
	}
	if ( toggle_qirrrt ) {
		toggle_qirrrt = 0;
	} else {
		toggle_qirrrt = 1;
	}

	// Debug
 	if ( RadiationDebug ) {
 		int check_neg = 1;
    int check_zero = 0;
 		CheckField(QirrRT, check_neg, check_zero, "ToggleQirrRT");
 	}
}

real ApplyInnerBC(T1, T2)
	// Input
	real T1, T2;
{
	// Declaration
	real Treturn, rgrad;

	// Function	
	if ( InnerBCExtrap ) {
		rgrad = (Rmed[1] - Rmed[0])/(Rmed[2] - Rmed[1]);
		Treturn = T1*(1.0+rgrad) - T2*rgrad;
	} else if ( InnerBCCons ) {
		Treturn = TINNER;
	} else if ( InnerBCGrad ) {
		Treturn = T1 - TGRADINNER*(Rmed[1] - Rmed[0]);
	} else {
		Treturn = 1.0E-10;
	}

	// Output
	return Treturn;
}

real ApplyOuterBC(T1, T2)
	// Input
	real T1, T2;
{
	// Declaration
	real Treturn, rgrad;

	// Function
	if ( OuterBCExtrap ) {
		rgrad = (Rmed[NRAD - 1] - Rmed[NRAD - 2])/(Rmed[NRAD - 2] - Rmed[NRAD - 1]);
		Treturn = T1*(1.0+rgrad) - T2*rgrad;
	} else if ( OuterBCCons ) {
		Treturn = TOUTER;
	} else if ( OuterBCGrad ) {
		Treturn = TGRADOUTER*(Rmed[NRAD-1] - Rmed[NRAD-2])+T1;
	} else {
		Treturn = 1.0E-10;
	}

	// Output
	return Treturn;
}

void ComputeNewEnergyField(gas_density, gas_energy)
	// Input
	PolarGrid *gas_density;
	PolarGrid *gas_energy;
{	
	// Declaration
	int i, j, l, nr, ns;
	real *dens, *energy, *temp;

	// Assignment
	nr = gas_density->Nrad;
	ns = gas_density->Nsec;
	dens = gas_density->Field;
	energy = gas_energy->Field;
	temp = Temperature->Field;

	// Function
#pragma omp parallel for private(j,l)
	for (i = 0; i < nr; i++) {
		for (j = 0; j < ns; j++) {
			l = j+i*ns;
			energy[l] = CV*temp[l]*dens[l];
		}
	}
	// Debug
 	if ( RadiationDebug ) {
 		int check_neg = 1;
    	int check_zero = 1;
 		CheckField(gas_energy, check_neg, check_zero, "ComputeNewEnergyField");
 	}
}


void InitGasTemperature (Density, Energy)
	// Input
	PolarGrid *Density;
  PolarGrid *Energy;
{	
	// Declaration
  int i, j, l, nr, ns;
  real *dens, *energy, *T;
  real a,b,c;

  // Assignment
  nr = Energy->Nrad;
  ns = Energy->Nsec;
  dens = Density->Field;
  energy = Energy->Field;
  T = Temperature->Field;

  // Constants
  a = 1000.0;
  b = (GlobalRmed[nr-1] - GlobalRmed[0])/2.0;
  c = b/6.0;

  // Function
#pragma omp parallel for private(j,l)
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      T[l] = a*exp(-pow(Rmed[i] - b, 2.0)/(2.0*c*c));
      energy[l] = R/MU/(ADIABATICINDEX-1.0)*T[l]*dens[l];
      EnergyMed[i] = energy[l];
    }
  }
}

void WriteRadTransInfo()
	// Input N/A
{
	// Declaration
	real tmp;
	int tmp_int;

	// Function
	MPI_Barrier (MPI_COMM_WORLD);
	NormAv /= counter;
	IterAv /= counter;
	MPI_Allreduce(&NormAv, &tmp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	NormAv = tmp/CPU_Number;
	MPI_Allreduce(&IterAv, &tmp_int, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
	IterAv = tmp_int/CPU_Number;

	// Output
	masterprint("Average Norm value = %g\n", NormAv);
	masterprint("Average No.of Iterations = %d\n", IterAv);
	NormAv = IterAv = tmp = tmp_int = 0;
	counter = 1;
}

void ResetTempSourcesSinks ()
	// Input N/A
{	// Declaration
	int i, nr, ns;
	real *Rij;

	//Assignment
	nr = TempSourcesSinks->Nrad;
	ns = TempSourcesSinks->Nsec;
	Rij = TempSourcesSinks->Field;

	//Function
#pragma omp parallel for
	for (i = 0; i < ns+nr*ns; i++) Rij[i] = 0.0;
}

real FLDConditionCFL()
	// Input N/A
{
	// Declaration
	int i, j, l, nr, ns;
	int lim, lip, ljm, ljp;
	real *U1, *U2, *U3, *U4, *B, *T, *D;
	real old_dt, new_dt = 1.0E30, factor, dTdt;

	// Assignment
	nr = Temperature->Nrad;
	ns = Temperature->Nsec;
	T  = Temperature->Field;
	U1 = U1arr->Field;
	U2 = U2arr->Field;
	U3 = U3arr->Field;
	U4 = U4arr->Field;
	B  = Barr->Field;
	D  = Darr->Field;

	// Constants
	factor = 0.1;

	// Function
	for (i = One_or_active; i < MaxMO_or_active; i++) {
		for (j = 0; j < ns; j++) {
			l = j + i*ns;
			// lim = l-ns;
			// lip = l+ns;
			// if ( j == 0 ) {
			// 	ljm = ns*(i+1)-1;
			// } else {
			// 	ljm = l-1;
			// }

			// if ( j == ns-1 ) {
			// 	ljp = i*ns;
			// } else {
			// 	ljp = l+1;
			// }

			// dTdt = U1[l]*T[lim] + U2[l]*T[lip] + U3[l]*T[ljm] + U4[l]*T[ljp] + (B[l]-1.0)*T[l];
			// dTdt = fabs(dTdt);
			// old_dt = T[l]/dTdt;
			// old_dt /= factor;
			old_dt = factor*DiffRsup[i]*DiffRsup[i]/D[l];

			if ( old_dt < new_dt ) {
				new_dt = old_dt;
			}
			old_dt = factor*Rmed[i]*Rmed[i]*DTHETA*DTHETA/D[l];

			if ( old_dt < new_dt ) {
				new_dt = old_dt;
			}
		}
	}

	// Output
	return new_dt;
}

void SubStep4_Explicit(gas_density, gas_energy, timestep)
	// Input
	PolarGrid *gas_density;
	PolarGrid *gas_energy;
	real timestep;
{
	// Declaration
	int i, j, l, nr, ns, nstep;
	int im, ip, jm, jp;
	int lim, lip, ljm, ljp;
	real *density, *T, *energy, *Tnew;
	real *U1, *U2, *U3, *U4;
	real U1T, U2T, U3T, U4T, CT;
	real *Qirrrt;
	real dt_FLD;
	real dt_fld, dt, dt_remainder = 0.0, c;
	int FLDTimeStepsCFL;

	// Assignment
	nr = gas_density->Nrad;
	ns = gas_density->Nsec;
	density = gas_density->Field;
	energy = gas_energy->Field;
	T = Temperature->Field;
	U1 = U1arr->Field;
	U2 = U2arr->Field;
	U3 = U3arr->Field;
	U4 = U4arr->Field;
	Tnew = TempGuess->Field;
	if ( RayTracingHeating ) {
		Qirrrt = QirrRT->Field;
	}

	// Function
	/* Use FLD coefficients to find how many explicit FLD sub-cycles are needed per hydro timestep */
	dt_fld = FLDConditionCFL();
	MPI_Allreduce(&dt_fld, &dt_FLD, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

	if ( dt_FLD >  timestep ) { 
    FLDTimeStepsCFL = 1;
    dt = timestep;
  } else {
    FLDTimeStepsCFL = (int)(timestep/dt_FLD);
    dt = dt_FLD;
    dt_remainder = timestep - (real)FLDTimeStepsCFL*dt_FLD;
  }
  for (i = 0; i < nr; i++) {
  	for (j = 0; j < ns; j++) {
  		l = j + i*ns;
  		Tnew[l] = 0.0;
  	}
  }
  // WriteDiskPolar(Temperature, 0);

  // masterprint("dt_hydro = %g, dt_fld = %g, remainder = %g: Number of FLD sub-cycles needed = %d\n", timestep, dt_FLD, dt_remainder, FLDTimeStepsCFL);
  /* Apply inner and outer temperature boundary conditions */
	if ( CPU_Rank == 0 ) {
		i = 0;
		for (j = 0; j < ns; j++) {
			l = j+i*ns;
			T[l] = TINNER;
			Tnew[l] = TINNER;
		}
	}
	if ( CPU_Rank == CPU_Highest ) {
		i = nr-1;
		for (j = 0; j < ns; j++) {
			l = j+i*ns;
			T[l] = TOUTER;
			Tnew[l] = TOUTER;
		}
	}

  for (nstep = 0; nstep < FLDTimeStepsCFL; nstep++) {
  	for (i = One_or_active; i < MaxMO_or_active; i++) {
  		im = i - 1;
  		ip = i + 1;
  		for (j = 0; j < ns; j++) {
  			l = j + i*ns;
  			jm = j - 1;
  			jp = j + 1;
  			if ( j == 0 ) {
  				jm = ns - 1;
  				jp = j + 1;
  			}
  			if ( j == ns-1 ) {
  				jp = 0;
  				jm = j - 1;
  			}
				
				lim = j + im*ns;
				lip = j + ip*ns;
				ljm = jm + i*ns;
				ljp = jp + i*ns;

				if ( i == 1 ) {
					U1T = U1[l]*TINNER;
				} else {
					U1T = U1[l]*T[lim];
				}

				if ( i == nr-1 ) {
					U2T = U2[l]*TOUTER;
				} else {
					U2T = U2[l]*T[lip];
				}

				U3T = U3[l]*T[ljm];
				U4T = U4[l]*T[ljp];

				CT = (U1[l] + U2[l] + U3[l] + U4[l])*T[l];

				// c = (U1[l] + U2[l] + U3[l] + U4[l]);

				// // Tnew[l] = T[l] + dt*(U1[l]*T[lim] + U2[l]*T[lip] + U3[l]*T[ljm] + U4[l]*T[ljp] + (B[l] - 1.0)*T[l]);
				// Tnew[l] = T[l] + dt*(U1[l]*T[lim] + U2[l]*T[lip] + U3[l]*T[ljm] + U4[l]*T[ljp] - c*T[l]);
				Tnew[l] = T[l] + dt*(U1T + U2T + U3T + U4T - CT);
				
  		}
  	}
  	/* Communicate active/overlap region boundaries and write new T values to T grid */
  	CommunicateFieldBoundaries(TempGuess);
  	

  	for (i = 0; i < nr; i++) {
  		for (j = 0; j < ns; j++) {
  			l = j+i*ns;
				T[l] = Tnew[l];
  		}
  	}
  	ComputeRKappa(gas_density);
  	ComputeRadTransCoeffs(gas_density, timestep);

  }

  /* Carry out one final sub-cycle with remainder timestep, to achieve time synchronicity */
  if (dt_remainder > 0.0) {
  	for (i = One_or_active; i < MaxMO_or_active; i++) {
  		im = i - 1;
  		ip = i + 1;
  		for (j = 0; j < ns; j++) {
  			l = j + i*ns;
  			if ( j == 0 ) {
  				jm = ns - 1;
  				jp = j + 1;
  			} else if ( j == ns-1 ) {
  				jp = 0;
  				jm = j - 1;
  			} else {
  				jm = j - 1;
  				jp = j + 1;
  			}
				
				lim = j + im*ns;
				lip = j + ip*ns;
				ljm = jm + i*ns;
				ljp = jp + i*ns;

				c = (U1[l] + U2[l] + U3[l] + U4[l]);

				Tnew[l] = T[l] + dt_remainder*(U1[l]*T[lim] + U2[l]*T[lip] + U3[l]*T[ljm] + U4[l]*T[ljp] - c*T[l]);
  		}
  	}

  	/* Communicate active/overlap region boundaries and write new T values to T grid */
  	CommunicateFieldBoundaries(TempGuess);

  	for (i = 0; i < nr; i++) {
  		for (j = 0; j < ns; j++) {
  			l = j+i*ns;
				T[l] = Tnew[l];
  		}
  	}
  }

  if ( RayTracingHeating ) {
  	for (i = 0; i < nr; i++) {
  		for (j = 0; j < ns; j++) {
  			l = j + i*ns;
  			T[l] += timestep*Qirrrt[l];
  		}
  	}
  }

  /* Convert new temperature field back to energy density */
  for (i = 0; i < nr; i++) {
  	for (j = 0; j < ns; j++) {
  		l = j+i*ns;
			energy[l] = T[l]*density[l]*CV;
  	}
  }
}