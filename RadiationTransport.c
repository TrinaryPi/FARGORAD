/* C Header
	* @filename        : RadiationTransport.c
	* @author          : Matthew Mutter (trinarypi)
	* @last_modified_by: trinarypi
	* @last_modified   : 2017/06/26 17:47
	* @description     :
*/
#include "mp.h"
#include "radiation.h"

static real NormAv=0;
static int IterAv=0, counter=1;
extern int FLDTimeStepsCFL;

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
	real qirrrt=0.0, qirr=0.0;
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
	/* Set Inner and Outer radii temperature (and temperature guess) values to constant boundary values */
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
		/* Instead of an iterative implicit SOR step, carry out FLD with a set of explicit sub-cycles */
		SubStep4_Explicit(gas_density, gas_energynew, timestep);
	} else {
		norm_tmp[0] = 0.0;
	 	norm_tmp[1] = 0.0;
	 	for (i = 0; i < nr; i++) {
	 		if ( ChebyshevOmega ) {
	 			omegaOpt[i] = 1.0;
	 		}
	 		for (j = 0; j < ns; j++) {
	 			l = j+i*ns;
	 			if ( ImplicitRadiative ) {
	 				if  ( RayTracingHeating ) {
		  			qirrrt = QirrRT->Field[l];
		  		}
		  		if ( Irradiation ) {
		  			qirr = Qirr->Field[l];
		  		}
		  	}
		  	Rij[l] = T[l] + qirrrt + cvfac*qirr/gas_density->Field[l];
	 			Tguess[l] = T[l];
	 			Tguess_old[l] = 0.0;
	 			norm_tmp[1] += fabs(Rij[l]);
	 		}
	 	}

	 	MPI_Allreduce(&norm_tmp[1], &norm_tmp[1], 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	 	tol = TOLERANCE*norm_tmp[1];
	 	norm = 2.0*tol;
	 	iteration = 0;

	 	while (( norm > tol ) && ( iteration < MAXITERATIONS )) {
	 		norm = 0.0;
	 		norm_tmp[0] = 0.0;

	 		/* Apply BC */
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
	 		if ( CPU_Rank == CPU_Highest ) {
	 			i = nr-2;
	#pragma omp parallel for private(l, lim, lip)
	 			for (j = 0; j < ns; j++) {
	 				l = j+i*ns;
	 				lim = l-ns;
	 				lip = l+ns;
	 				Tguess[lip] = ApplyOuterBC(Tguess[l], Tguess[lim]);
	 			}
	 		}
	 		
			/* Black-loop over all even numbered cells */
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
	 				}
	 			}
	 		}

	 		/* If specified, use Chebyshev acceleration to automatically tune the Over-relaxation parameter, omega, each half timestep */
	 		if ( ChebyshevOmega ) {
	 			for (i = 0; i < nr; i++) {
	 				if ( iteration == 0 ) {
	 					omegaOpt[i] = 1.0/(1.0 - 0.5*rhoJac[i]*rhoJac[i]);
	 				} else {
	 					omegaOpt[i] = 1.0/(1.0 - 0.25*rhoJac[i]*rhoJac[i]*omegaOpt[i]);
	 				}
	 			}
	 		}

			/* Overlap zones communication between neighbouring processors */
	 		CommunicateFieldBoundaries(TempGuess);

			/* Red-loop over all odd numbered cells */
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
	 				}
	 			}
	 		}

	 		if ( ChebyshevOmega ) {
	 			for ( i = 0; i < nr; i++) {
	 				omegaOpt[i] = 1.0/(1.0 - 0.25*rhoJac[i]*rhoJac[i]*omegaOpt[i]);
	 			}
	 		}

	 		MPI_Allreduce(&norm_tmp[0], &norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

			/* Overlap zones communication between neighbouring processors */
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

	 	for (i = 0; i < nr; i++) {
	 		for (j = 0; j < ns; j++) {
	 			l = j+i*ns;
	 			T[l] = Tguess[l];
	 			if (( ExplicitRayTracingHeating ) && (( i >= One_or_active ) || ( i < MaxMO_or_active ))) {
	 				T[l] += QirrRT->Field[l];
	 			}
	 		}
	 	}

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


void ComputeRadTransCoeffs(gas_density, dt)
	// Input
	PolarGrid *gas_density;
	real dt;
{	
	// Declaration
	int nr, ns, i, im, ip, jm, jp;
	int j, l, lim, lip, ljm, ljp;
	real *sigma, *kappa, *height, *T, *rfld, *lambda;
	real *D, *B, *U1, *U2, *U3, *U4;
	real elim, elip, eljm, eljp, el, rho;
	real grad1, grad2, grad;
	real Dip, Dim, Djm, Djp;
	real cvfac, c1;

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
	rfld = Rfld->Field;
	lambda = lambdafld->Field;
	
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
			lim = j+im*ns;
			lip = j+ip*ns;
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
			ljm = jm+i*ns;
			ljp = jp+i*ns;

			elim = pow(T[lim], 4.0);
			elip = pow(T[lip], 4.0);
			eljm = pow(T[ljm], 4.0);
			eljp = pow(T[ljp], 4.0);
			el = pow(T[l], 4.0);

			grad1 = fabs((elip - elim)/(Rmed[ip] - Rmed[im]));
			grad2 = fabs((eljp - eljm)/(2.0*Rmed[i]*DTHETA));
			grad = grad1+grad2;
			rho = c1*sigma[l]/height[l];
			rfld[l] = grad/(rho*kappa[l]*el);

			if ( rfld[l] <= 2.0 ) {
				lambda[l] = 2.0/(3.0 + sqrt(9.0 + (10.0*rfld[l]*rfld[l])));
			} else {
				lambda[l] = 10/((10*rfld[l]) + 9.0 + sqrt(81.0 + (180.0*rfld[l])));
			}
	
			D[l] = 4.0*STEFANK*lambda[l]*pow(T[l], 3.0)/(rho*kappa[l]);
			if ( D[l] < 0.0 ) {
				printf("D negative @ CPU_%d, i = %d, j = %d, lambda[l] = %g, rho[l] = %g\n", CPU_Rank, i, j, lambda[l], rho);
				printf("T[l]^3 = %g, kappa[l] = %g, STEFANK = %g\n", pow(T[l], 3.0), kappa[l], STEFANK);
			}
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

			Dim = D[lim] + (Rinf[i] - Rmed[i-1])*(D[l] - D[lim])/(Rmed[i] - Rmed[i-1]);
			Dip = D[l] + (Rinf[i+1] - Rmed[i])*(D[lip] - D[l])/(Rmed[i+1] - Rmed[i]);
			Djp = (D[l] + D[ljp]);
			Djm = (D[ljm] + D[l]);

			/* Coefficient for T_{i-1j}^{n+1} ie cell lim at next timestep level */
			U1[l] = (cvfac/sigma[l])*(Rinf[i]*Dim)/(Rmed[i]*(Rsup[i] - Rinf[i])*(Rmed[i] - Rmed[i-1]));

			/* Coefficient for T_{i+1j}^{n+1} ie cell lip at next timestep level */
			U2[l] = (cvfac/sigma[l])*(Rinf[i+1]*Dip)/(Rmed[i]*(Rsup[i] - Rinf[i])*(Rmed[i+1] - Rmed[i]));

			/* Coefficient for T_{ij-1}^{n+1} ie cell ljm at next timestep level */
			if ( NSEC > 1 ) {
				U3[l] = (cvfac/sigma[l])*Djm/(2.0*Rmed[i]*Rmed[i]*DTHETA*DTHETA);
			} else {
				U3[l] = 0.0;
			}

			/* Coefficient for T_{ij+1}^{n+1} ie cell ljp at next timestep level */
			if ( NSEC > 1 ) {
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

real FLDConditionCFL()
	// Input N/A
{
	// Declaration
	int i, j, l, ns;
	real *D;
	real old_dt, new_dt=1.0E30, factor;

	// Assignment
	ns = Darr->Nsec;
	D  = Darr->Field;

	// Constants
	factor = 0.1;

	// Function
	for (i = One_or_active; i < MaxMO_or_active; i++) {
		for (j = 0; j < ns; j++) {
			l = j + i*ns;
			old_dt = factor*DiffRsup[i]*DiffRsup[i]/fabs(D[l]);
			if ( old_dt < new_dt ) {
				new_dt = old_dt;
			}
			old_dt = factor*Rmed[i]*Rmed[i]*DTHETA*DTHETA/fabs(D[l]);
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
	real *density, *T, *energy, *Tnew, *Qirrrt;
	real *U1, *U2, *U3, *U4;
	real U1T, U2T, U3T, U4T, CT;
	real dt_FLD, dt_fld, dt, dt_remainder;

	// Assignment
	char str[256];
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
	dt_remainder = 0.0;

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
  for (nstep = 0; nstep < FLDTimeStepsCFL; nstep++) {
  	for (i = One_or_active; i < MaxMO_or_active; i++) {
  		im = i - 1;
  		ip = i + 1;
  		for (j = 0; j < ns; j++) {
  			lim = j + im*ns;
				lip = j + ip*ns;
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
				ljm = jm + i*ns;
				ljp = jp + i*ns;

				U1T = U1[l]*T[lim];
				U2T = U2[l]*T[lip];
				U3T = U3[l]*T[ljm];
				U4T = U4[l]*T[ljp];
				CT = (U1[l] + U2[l] + U3[l] + U4[l])*T[l];

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
  	if ( debug ) {
  		int check_neg = 1;
    	int check_zero = 1;
    	sprintf(str, "SubStep4 Explicit after sub-cycle %d", nstep);
    	CheckField(Temperature, check_neg, check_zero, str);
  	}
  	ComputeRKappa(gas_density);
  	ComputeRadTransCoeffs(gas_density, timestep);
  }
  /* Carry out one final sub-cycle with remainder timestep, to achieve time synchronicity */
  if ( dt_remainder > 0.0 ) {
  	for (i = One_or_active; i < MaxMO_or_active; i++) {
  		im = i - 1;
  		ip = i + 1;
  		for (j = 0; j < ns; j++) {
  			lim = j + im*ns;
				lip = j + ip*ns;
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
				ljm = jm + i*ns;
				ljp = jp + i*ns;

				U1T = U1[l]*T[lim];
				U2T = U2[l]*T[lip];
				U3T = U3[l]*T[ljm];
				U4T = U4[l]*T[ljp];
				CT = (U1[l] + U2[l] + U3[l] + U4[l])*T[l];

				Tnew[l] = T[l] + dt_remainder*(U1T + U2T + U3T + U4T - CT);
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
  // Debug
  if ( debug ) {
  	int check_neg = 1;
    int check_zero = 1;
    CheckField(Temperature, check_neg, check_zero, "SubStep4 Explicit");
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
