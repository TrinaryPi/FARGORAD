#include "mp.h"

extern real ScalingFactor, Abinary;
extern boolean BinaryOn, Sigma_Taper, Sigma_Cavity, CustomCooling;

/* Surface density */
real Sigma(r)
     real r;
{
  real gap_func = 1.0;
  real taper_func = 1.0;
  real gap;
  if (Sigma_Cavity == YES) {
    if (BinaryOn == YES) {
      gap = 3.0*Abinary;      /*This should be 2.5 the binary separation, a_b, in whatever units etc */
      gap_func = 1.0/(1.0+exp(10.0*(gap-r)/gap));
    } else {
      gap = 1.6*RMIN;
      gap_func = 1.0/(1.0+exp(10.0*(gap-r)/gap));
    }
  }

  /*if (r < CAVITYRADIUS) cavity = 1.0/CAVITYRATIO;
  gap_func = cavity*pow(r,-SIGMASLOPE);*/

  if (Sigma_Taper == YES) {
    taper_func = tanh((RMAX - r + 0.01*RMAX)/0.1/RMAX);
  }
  /* This is *not* a steady state */
  /* profile, if a cavity is defined. It first needs */
  /* to relax towards steady state, on a viscous time scale */
  
  return taper_func*gap_func*ScalingFactor*SIGMA0*pow(r,-SIGMASLOPE);
}

void FillSigma() {
  int i;
  for (i = 0; i < NRAD; i++) {
    SigmaMed[i] = Sigma(Rmed[i]);
    SigmaInf[i] = Sigma(Rinf[i]);
  }
}

void RefillSigma (Surfdens)
     PolarGrid *Surfdens;
{
  int i, j, nr, ns, l;
  real *field;
  real moy;
  nr = Surfdens->Nrad;
  ns = Surfdens->Nsec;
  field = Surfdens->Field;
  for (i = 0; i < nr; i++) {
    moy = 0.0;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      moy += field[l];
    }
    moy /= (real)ns;
    SigmaMed[i] = moy;
  }
  SigmaInf[0] = SigmaMed[0];
  for (i = 1; i < nr; i++) {
    SigmaInf[i] = (SigmaMed[i-1]*(Rmed[i]-Rinf[i])+\
		   SigmaMed[i]*(Rinf[i]-Rmed[i-1]))/\
      (Rmed[i]-Rmed[i-1]);
  }
}

/* Thermal energy */
real Energy(r)
     real r;
{
  real energy0;
  if (ADIABATICINDEX == 1.0) {
    fprintf (stderr, "The adiabatic index must differ from unity to initialize the gas internal energy. I must exit.\n");
    prs_exit (1);
  } else {
  	if (Sigma_Cavity == YES) {
      real gap_func = 1.0;
      real taper_func = 1.0;
      real gap;
      real sigma;
    	if (BinaryOn == YES) {
      		gap = 3.0*Abinary;      /*This should be 2.5 the binary separation, a_b, in whatever units etc */
      		gap_func = 1.0/(1.0+exp(10.0*(gap-r)/gap));
    	} else {
      		gap = 1.6*RMIN;
      		gap_func = 1.0/(1.0+exp(10.0*(gap-r)/gap));
    	}
      sigma = taper_func*gap_func*ScalingFactor*SIGMA0*pow(r,-SIGMASLOPE);
      energy0 = R/MU/(ADIABATICINDEX-1.0)*pow(ASPECTRATIO,2.0)*pow(r,-1.0+2.0*FLARINGINDEX)*sigma;
    	// energy0 = R/MU/(ADIABATICINDEX-1.0)*gap_func*SIGMA0*pow(ASPECTRATIO,2.0)*pow(r,-SIGMASLOPE-1.0+2.0*FLARINGINDEX);
  	} else {
  		energy0 = R/MU/(ADIABATICINDEX-1.0)*SIGMA0*pow(ASPECTRATIO,2.0)*pow(r,-SIGMASLOPE-1.0+2.0*FLARINGINDEX);
  	}
  	// energy0 = R/MU/(ADIABATICINDEX-1.0)*SIGMA0*pow(ASPECTRATIO,2.0)*pow(r,-SIGMASLOPE-1.0+2.0*FLARINGINDEX);
  }

  return energy0;
}

void FillEnergy() {
  int i;
  for (i = 0; i < NRAD; i++)
    EnergyMed[i] = Energy(Rmed[i]);
}


void RefillEnergy (energy)
     PolarGrid *energy;
{
  int i, j, nr, ns, l;
  real *field;
  real moy;
  nr = energy->Nrad;
  ns = energy->Nsec;
  field = energy->Field;
  for (i = 0; i < nr; i++) {
    moy = 0.0;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      moy += field[l];
    }
    moy /= (real)ns;
    EnergyMed[i] = moy;
  }
}

/* Cooling time */
real CoolingTime(r)
     real r;
{
  real ct0;
  real minFrac, maxFrac, delFrac = 0.1, RRange;
  if (CustomCooling == YES) {
    maxFrac = COOLINGTIME0;
    minFrac = maxFrac*delFrac;
    RRange = RMAX-RMIN;
    if (CPU_Master)
    	masterprint("Using variable, custom cooling time fraction\n");
    if (r < 2.0*RRange/3.0)
      ct0 = minFrac + (maxFrac - minFrac)*exp(-pow(r - (2.0/3.0)*RRange, 2.0)/2.0/pow(0.1*RRange, 2.0));
    else
      ct0 = maxFrac;

  } else {
    ct0 = COOLINGTIME0*pow(r,2.0+2.0*FLARINGINDEX);
  }
  return ct0;
}

void FillCoolingTime() {
  int i;
  for (i = 0; i < NRAD; i++)
    CoolingTimeMed[i] = CoolingTime(Rmed[i]);
}

/* Heating source term */
real Qplusinit(r)
     real r;
{
  real qp0, viscosity;
  viscosity = FViscosity(r);
  qp0 = 2.25*viscosity*SIGMA0*pow(r,-SIGMASLOPE-3.0);
  return qp0;
}

void FillQplus() {
  int i,j,l, ns;
  real *qplus;
  qplus = Qplus->Field;
  ns = Qplus->Nsec;
  for (i = 0; i < NRAD; i++){
    QplusMed[i] = Qplusinit(Rmed[i]);
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      qplus[l] = QplusMed[i];
    }
  }
    
}
