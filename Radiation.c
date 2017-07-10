/* C Header
	* @filename        : Radiation.c
	* @author          : Matthew Mutter (trinarypi)
	* @last_modified_by: trinarypi
	* @last_modified   : 2017/06/27 12:22
	* @description     :
*/
#include "mp.h"
#include "radiation.h"
#define LITTLE_NUMBER_FIX 1.0E-15


extern boolean BinaryOn;
extern boolean NoCFL, RadiationDebug;
extern real StarTaper;


void ComputeDiscHeight (bsys)
	// Input
	BinarySystem *bsys;
{
	// Declaration
	int i, j, l, s, nb, ns, nr;
	real *xstar, *ystar, *mb, xc, yc, dist, angle, smooth;
	real *height, *cs;
	real average, csiso, drx, dry, A, q[2], dist_smooth;

	// Assignment
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
			for (i = 0; i < nr; i++) {
				average = 0.0;
				for (j = 0; j < ns; j++) {
					l = j + i*ns;
					angle = (real)j/(real)ns*2.0*PI;
					xc = Rmed[i]*cos(angle);
					yc = Rmed[i]*sin(angle);
					height[l] = 0.0;
					for (s = 0; s < nb; s++) {
						/* Add a smoothing length to cell-binary distance, to avoid singularities */
						smooth = bsys->smooth[s]*(0.49*pow(q[s], 2.0/3.0))/(0.6*pow(q[s], 2.0/3.0) + log(1.0 + pow(q[s], 1.0/3.0)))*A;
						smooth *= smooth;
						smooth = 0.0;
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
			for (i = 0; i < nr; i++) {
				average = 0.0;
				for (j = 0; j < ns; j++) {
					l = j+i*ns;
					height[l] = cs[l]*sqrt(pow(Rmed[i], 3.0)/G/1.0); /* M = mass of stellar system = 1.0 */
					average += height[l];
				}
				MeanDiscHeight[i] = average/(real)ns;
			}
		}
	} else {
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
	if (( dist >= Rmed[NRAD-1] ) || ( dist <= Rmed[0] )) {
		smoothing = 0.0001*THICKNESSSMOOTHING;
	} else {
		smoothing = compute_aspectratio(xp, yp);
		smoothing = smoothing * THICKNESSSMOOTHING * pow(dist, 1.0+FLARINGINDEX);
	}

	// Debug
	if ( RadiationDebug ) {
		int check_neg = 1;
    int check_zero = 1;
		int flag = 0;
		int GlobalFlag = 0;
		flag = CheckValue(smoothing, check_neg, check_zero);
		MPI_Allreduce(&flag, &GlobalFlag, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
		if ( GlobalFlag != 0 ) {
			printf("Error: Non-normal value in variable disc height smoothing. Exiting.\n");
			MPI_Finalize();
			exit(GlobalFlag);
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
		while ( dangle*((real)j + 0.5) < angle ) {
			j++;
		}
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

	// Output
	return HoverR;
}


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

	// Assignment
	ns = Sigma->Nsec;
	nr = Sigma->Nrad;
	dens = Sigma->Field;
	H = DiscHeight->Field;
	rkappa = RKappaval->Field;
	qirr = Qirr->Field;
	tau = OpticalDepth->Field;
	taueff = OpticalDepthEff->Field;
	nb = IrrSources->nb;
	Rs = IrrSources->Rstar;
	Tstar = IrrSources->Tstar;
	central_source = IrrSources->CentralSource;
	xs = bsys->x;
	ys = bsys->y;

	// Constants
	disc_albedo = 0.5;
	c2 = 0.5;
	constant1 = 2.0*STEFANK;
	tau_min = 0;

	// Function
	for (i = 0; i < nr; i++) {
		for (j = 0; j < ns; j++) {
			l = j+i*ns;
			tau[l] = c2*rkappa[l]*dens[l];
			taueff[l] = (0.375*tau[l]) + (1.0/(4.0*tau[l] + tau_min)) + 0.866;
			qirr[l] = 0.0;
			for (s = 0; s < nb; s++) {
				Ts = StarTaper*Tstar[s];
				constant2 = constant1*(1.0 - disc_albedo)*pow(Ts*Ts*Rs[s], 2.0);
				if ( central_source[s] ) {
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
    int check_zero = 0;
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


void ResetTempSourcesSinks()
	// Input N/A
{	// Declaration
	int i, nr, ns;
	real *Rij;

	//Assignment
	nr = TempSourcesSinks->Nrad;
	ns = TempSourcesSinks->Nsec;
	Rij = TempSourcesSinks->Field;

	//Function
	for (i = 0; i < ns + nr*ns; i++) {
		Rij[i] = 0.0;
	}
}

