/* C Header
	* @filename        : SourceEuler.c
	* @author          : Frederic Masset
	* @last_modified_by: trinarypi
	* @last_modified   : 2017/06/28 12:11
	* @description     :
*/
#include "mp.h"
#define CFLSECURITY 0.5 			/* Maximum fraction of zone size */
															/* swept in one timestep         */

#define CVNR 1.41       			/* Shocks are spread over CVNR zones:       */
                          		/* von Neumann-Richtmyer viscosity constant */
															/* Beware of misprint in Stone and Norman's */
															/* paper : use C2^2 instead of C2           */
#define NBODYSECURITY 200.0  	/* Fraction of the smallest pair-wise period in the binary-planet system */


static PolarGrid *TemperInt;
static PolarGrid *VradNew, *VradInt;
static PolarGrid *VthetaNew, *VthetaInt;
static PolarGrid *EnergyNew, *EnergyInt;
static real timeCRASH;  
extern boolean Corotating;
extern boolean Adiabatic, DiscMassTaper;
extern boolean SelfGravity, SGZeroMode;
extern boolean ZMPlus;
real PhysicalTime=0.0, OmegaFrame, PhysicalTimeInitial;
int FirstGasStepFLAG=1;
static int AlreadyCrashed = 0, GasTimeStepsCFL;
static int RadCoolTimeStepsCFL;
extern boolean FastTransport, IsDisk, BinaryOn, HydroOn, LiveBodies;
Pair DiskOnPrimaryAcceleration;
Pair DiskOnBinaryAcceleration;
extern real mlinner, mlouter;
real SigmaSafetyFloor = 1.0E-10, EnergySafetyFloor;
extern boolean VarDiscHeight, RadCooling, Irradiation, RadTransport, TempInit;
extern boolean NoCFL, RadiativeOnly;
extern boolean RayTracingHeating, ExplicitRayTracingHeating, ImplicitRadiative;
static int crash_i, crash_j;
static real dt_hydro_av = 0.0, dt_energy_av = 0.0;
static int nsteps_av = 0;
extern boolean  AnalyticCooling;
extern boolean ExplicitRadTransport;
extern int FLDTimeStepsCFL;


boolean DetectCrash (array)
     PolarGrid *array;
{
  int i, j, l, nr, ns;
  real *ptr;
  boolean bool = NO;

  nr = array->Nrad;
  ns = array->Nsec;
  ptr= array->Field;

#pragma omp parallel for private(j,l) shared(bool)
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      if ( ptr[l] < 0.0 ) {
	      bool = YES;
        crash_i = i+IMIN;
        crash_j = j;
      }
    }
  }
  return bool;
}


void FillPolar1DArrays ()
{
  FILE *input, *output;
  int i,ii;
  real drrsep;
  float temporary;
  char InputName[256], OutputName[256];

  DTHETA = 2.0*PI/(real)NSEC;
  CV = R/MU/(ADIABATICINDEX-1.0);
  drrsep = (RMAX-RMIN)/(real)GLOBALNRAD;

  sprintf (InputName, "%s%s", OUTPUTDIR, "radii.dat");
  sprintf (OutputName, "%s%s", OUTPUTDIR, "used_rad.dat");
  input = fopen (InputName, "r");
  if ( input == NULL ) {
    mastererr ("Warning : no `radii.dat' file found. Using default.\n");
    if ( LogGrid == YES ) {
      for (i = 0; i <= GLOBALNRAD; i++) {
	      Radii[i] = RMIN*exp((real)i/(real)GLOBALNRAD*log(RMAX/RMIN));
      }
    } else {
      for (i = 0; i <= GLOBALNRAD; i++) {
	      Radii[i] = RMIN+drrsep*(real)(i);
      }
    }
  } else {
    mastererr ("Reading 'radii.dat' file.\n");
    for (i = 0; i <= GLOBALNRAD; i++) {
      fscanf (input, "%f", &temporary);
      Radii[i] = (real)temporary;
    }
  }
  for (i = 0; i < GLOBALNRAD; i++) {
    GlobalRmed[i] = 2.0/3.0*(Radii[i+1]*Radii[i+1]*Radii[i+1]-Radii[i]*Radii[i]*Radii[i]);
    GlobalRmed[i] = GlobalRmed[i] / (Radii[i+1]*Radii[i+1]-Radii[i]*Radii[i]);
  }
  for (i = 0; i <= NSEC; i++) {
    GlobalTheta[i] = (real)i/(real)NSEC*2.0*PI;
  }
  for (i = 0; i < NRAD; i++) {
    ii = i+IMIN;
    Rinf[i] = Radii[ii];
    Rsup[i] = Radii[ii+1];
    Rmed[i] = 2.0/3.0*(Rsup[i]*Rsup[i]*Rsup[i]-Rinf[i]*Rinf[i]*Rinf[i]);
    Rmed[i] = Rmed[i] / (Rsup[i]*Rsup[i]-Rinf[i]*Rinf[i]);
    Surf[i] = PI*(Rsup[i]*Rsup[i]-Rinf[i]*Rinf[i])/(real)NSEC;
    InvRmed[i] = 1.0/Rmed[i];
    InvSurf[i] = 1.0/Surf[i];
    DiffRsup[i] = Rsup[i]-Rinf[i];
    InvDiffRsup[i] = 1.0/DiffRsup[i];
    InvRinf[i] = 1.0/Rinf[i];
  }
  Rinf[NRAD]=Radii[NRAD+IMIN];
  for (i = 1; i < NRAD; i++) {
    InvDiffRmed[i] = 1.0/(Rmed[i]-Rmed[i-1]);
  }
  if ( CPU_Master ) {
    output = fopen (OutputName, "w");
    if ( output == NULL ) {
      mastererr ("Can't write %s.\nProgram stopped.\n", OutputName);
      prs_exit (1);
    }
    for (i = 0; i <= GLOBALNRAD; i++) {
      fprintf (output, "%.18g\n", Radii[i]);
    }
    fclose (output);
  }
  if ( input != NULL ) {
  	fclose (input);
  }
}


void InitEuler (Vr, Vt, Rho, Energy)
     PolarGrid *Vr, *Vt, *Rho, *Energy;
{
  InitTransport ();
  InitViscosity ();
  RhoStar      = CreatePolarGrid(NRAD, NSEC, "RhoStar");
  RhoInt       = CreatePolarGrid(NRAD, NSEC, "RhoInt");
  VradNew      = CreatePolarGrid(NRAD, NSEC, "VradNew");
  VradInt      = CreatePolarGrid(NRAD, NSEC, "VradInt");
  VthetaNew    = CreatePolarGrid(NRAD, NSEC, "VthetaNew");
  VthetaInt    = CreatePolarGrid(NRAD, NSEC, "VthetaInt");
  EnergyNew    = CreatePolarGrid(NRAD, NSEC, "EnergyNew");
  EnergyInt    = CreatePolarGrid(NRAD, NSEC, "EnergyInt");
  TemperInt    = CreatePolarGrid(NRAD, NSEC, "TemperInt");
  Potential    = CreatePolarGrid(NRAD, NSEC, "Potential");
  Pressure     = CreatePolarGrid(NRAD, NSEC, "Pressure");
  SoundSpeed   = CreatePolarGrid(NRAD, NSEC, "SoundSpeed");
  if ( TempInit == NO ) {
    Temperature  = CreatePolarGrid(NRAD, NSEC, "Temperature");
  }
  Qplus        = CreatePolarGrid(NRAD, NSEC, "Qplus");
  QDiv         = CreatePolarGrid(NRAD, NSEC, "Qdiv");
  InitComputeAccel ();
  /* Rho and Energy are already initialized: cf main.c */
  ComputeSoundSpeed (Rho, Energy);
  ComputePressureField (Rho, Energy);
  ComputeTemperatureField (Rho, Energy);
  if ( ExplicitRadTransport ) {
    BoundaryConditionsFLD(Temperature, Temperature);
    ComputeNewEnergyField(Rho, Energy);
    ComputeSoundSpeed(Rho, Energy);
    ComputePressureField(Rho, Energy);
  }
  InitGasVelocities (Vr, Vt);
}


real min2 (a, b)
     real a, b;
{
  if ( b < a ) {
  	return b;
  }
  return a;
}

real max2 (a, b)
     real a, b;
{
  if ( b > a ) {
  	return b;
  }
  return a;
}


void ActualiseGas (array, newarray)
     PolarGrid *array, *newarray;
{
  int i,j,l,ns,nr;
  real *old, *new;

  nr = array->Nrad;
  ns = array->Nsec;
  old= array->Field;
  new= newarray->Field;

#pragma omp parallel for private(j,l)
  for (i = 0; i <= nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+ns*i;
      old[l] = new[l];
    }
  }
}

void AlgoGas (force, Rho, Vrad, Vtheta, Energy, Label, sys, bsys, Ecc, TimeStep)
     Force *force;
     PolarGrid *Rho, *Vrad, *Vtheta, *Energy, *Label, *Ecc;
     PlanetarySystem *sys;
     BinarySystem *bsys;
     int TimeStep;
{
  real dt, dtemp=0.0;
  real OmegaNew, domega;
  int gastimestepcfl;
  boolean Crashed=NO;
  FirstGasStepFLAG=1;
  gastimestepcfl = 1;
  int timestep_counter = 0;
  real radcooltimestepcfl;
  real dt_rc;

  if ( Adiabatic ) {
    ComputeSoundSpeed (Rho, Energy);
    /* it is necesary to update computation of soundspeed if one uses
       alphaviscosity in FViscosity. It is not necesary in locally
       isothermal runs since cs is constant.  It is computed here for
       the needs of ConditionCFL. */ 
  }

  if ( IsDisk ) {
    CommunicateBoundaries (Rho, Vrad, Vtheta, Energy, Label);
    if ( SloppyCFL ) {
      if ( !RadiativeOnly ) {
        ComputeViscousTerms (Vrad, Vtheta, Rho);
      }
      if ( Adiabatic ) {
        ComputeQplus(Rho);
      }
      gastimestepcfl = ConditionCFL (Vrad, Vtheta, Energy, DT-dtemp, sys, bsys);
    }
  }
  MPI_Allreduce (&gastimestepcfl, &GasTimeStepsCFL, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  dt = DT / (real)GasTimeStepsCFL;

  while (dtemp < 0.999999999*DT) {
    MassTaper = PhysicalTime/(MASSTAPER*2.0*M_PI);
    MassTaper = (MassTaper > 1.0 ? 1.0 : pow(sin(MassTaper*M_PI/2.0), 2.0));
    // New Star Temperature Taper (6/01/2017)
    if (STARTAPER < dt) {
    	StarTaper = 1.0;
    } else {
    	StarTaper = PhysicalTime/(STARTAPER*2.0*M_PI); 
    	StarTaper = (StarTaper > 1.0 ? 1.0 : pow(sin(StarTaper*M_PI/2.0), 2.0));
    }
    CopyDensity(Rho);

    if ( IsDisk ) {
      CommunicateBoundaries (Rho, Vrad, Vtheta, Energy, Label);
      if ( !SloppyCFL ) {
	      gastimestepcfl = 1;
	      if ( !NoCFL ) {
          if ( !RadiativeOnly ) {
            ComputeViscousTerms (Vrad, Vtheta, Rho);
          }
          if ( Adiabatic ) {
            ComputeQplus(Rho);
          }
          gastimestepcfl = ConditionCFL (Vrad, Vtheta, Energy, DT-dtemp, sys, bsys);
        }
	      MPI_Allreduce (&gastimestepcfl, &GasTimeStepsCFL, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	      dt = (DT-dtemp)/(real)GasTimeStepsCFL;
      }
      /* The mass of the disc is slowly tapered. */
      /* This should occur over a timescale much */
      /* longer than any other in the system.    */
      if ( DiscMassTaper ) {
        DiscMassTaperCalc (Rho, dt);
      }
      AccreteOntoPlanets (Rho, Vrad, Vtheta, dt, sys);
      /* Accrete onto binary stars goes here - need to find a dynamical timescale to use for this */
      if (( BinaryOn ) && ( HydroOn)) {
        AccreteOntoStars (Rho,dt, bsys);
      }
    }
    dtemp += dt;

    DiskOnPrimaryAcceleration.x = 0.0;
    DiskOnPrimaryAcceleration.y = 0.0;
    DiskOnBinaryAcceleration.x = 0.0;
    DiskOnBinaryAcceleration.y = 0.0;

    if ( Corotating ) {
      GetPsysInfo (sys, MARK);
    }

    if ( IsDisk  ) {
      if ( HydroOn ) {
        /* Indirect term star's potential computed here */
        if ( BinaryOn ) {
          DiskOnBinaryAcceleration = ComputeIndirect (force, Rho, 0.0, bsys);
        } else {
          DiskOnPrimaryAcceleration = ComputeAccel (force, Rho, 0.0, 0.0, 0.0, 0.0);
        }
        /* Gravitational potential from star and planet(s) is computed
  	     and stored here*/
        FillForcesArrays (bsys, sys, Rho, Energy);
        /* Planets' velocities are updated here from gravitational
  	     interaction with disk */
        if ( LiveBodies ) {
          AdvanceSystemFromDisk (force, Rho, Energy, sys, dt);
          if ( BinaryOn ) {
            AdvanceBinaryFromDisk(force, Rho, Energy, bsys, dt);
          }
        }
      }
    }
    /* Planets' positions and velocities are updated from
     gravitational interaction with star and other planets */
    if ( LiveBodies ) {
      AdvanceSystemRK5 (bsys, sys, dt);
    }
    if ( BinaryOn ) {
      CalcAbin (bsys);
    }
    /* Below we correct vtheta, planet's position and velocities if we
     work in a frame non-centered on the star */
    if ( Corotating ) {
      OmegaNew = GetPsysInfo(sys, GET) / dt;
      domega = OmegaNew-OmegaFrame;
      if ( IsDisk ) {
	      CorrectVtheta (Vtheta, domega);
      }
      OmegaFrame = OmegaNew;
    }
    RotatePsys (sys, OmegaFrame*dt);
    /* Now we update gas */
    if ( IsDisk ) {
      ApplyBoundaryCondition (Vrad, Vtheta, Rho, Energy, bsys, dt);
      Crashed = DetectCrash (Rho);	/* test for negative density values */
      Crashed = DetectCrash (Energy);	/* test for negative energy values */
      if ( Crashed ) {
	      if (AlreadyCrashed == 0) {
	        timeCRASH=PhysicalTime;	/* if it appears to be the first crash */
	        fprintf (stdout,"\nCrash! at time %.12g\n", timeCRASH);
          fprintf (stdout,"Crash at r = %f (ring %d) and azimuthal cell %d \n", Rmed[crash_i], crash_i, crash_j);
          prs_exit(1);
	      }
	      AlreadyCrashed++;
	      masterprint ("c");
      } else {
      	if ( FLDTimeStepsCFL > 1) {
      		masterprint("%d.", FLDTimeStepsCFL);
      	} else {
	      	masterprint (".");
	      }
      }
     	fflush (stdout);
     	if ( ZMPlus ) {
	      compute_anisotropic_pressurecoeff (sys);
	    }
      ComputePressureField (Rho, Energy);

    	if ( HydroOn ) {
        SubStep1 (Vrad, Vtheta, Rho, dt);
        if ( debug ) {
        	int check_neg = 0;
        	int check_zero = 0;
        	CheckField(VradInt, check_neg, check_zero, "SubStep1");
        	CheckField(VthetaInt, check_neg, check_zero, "SubStep1");
        }
        SubStep2 (Rho, Energy, dt);
        if ( debug ) {
        	int check_neg = 0;
        	int check_zero = 0;
        	CheckField(VradNew, check_neg, check_zero,"SubStep2");
        	CheckField(VthetaNew, check_neg, check_zero,"SubStep2");
        }
        if ( Adiabatic ) {
        	if ( RadiativeOnly ) {
        		ActualiseGas(EnergyInt, Energy);
        	}
        	if ( debug ) {
        		int check_neg = 1;
        		int check_zero = 1;
        		CheckField(EnergyInt, check_neg, check_zero, "SubStep2");
        	}
        }
        ActualiseGas (Vrad, VradNew);
        ActualiseGas (Vtheta, VthetaNew);
    	} else {
        ActualiseGas (EnergyInt, Energy);
      }
      ApplyBoundaryCondition (Vrad, Vtheta, Rho, Energy, bsys, dt);
      
      if ( Adiabatic ) {
        ComputeSoundSpeed (Rho, EnergyInt);
        ComputeQterms(Rho, Vrad, Vtheta, EnergyInt, bsys);
        SubStep3 (Rho, dt);
        if ( RadCooling ) {
          radcooltimestepcfl = RadCoolConditionCFL(Energy);
          MPI_Allreduce (&radcooltimestepcfl, &dt_rc, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
          SubStep3_RadCool(Rho, EnergyNew, dt, dt_rc);
        }
        if (( RadTransport ) || ( RayTracingHeating )) {
          ComputeSoundSpeed (Rho, EnergyNew);
          ComputeDiscHeight (bsys);
          ComputeTemperatureField (Rho, EnergyNew);
          ComputeRKappa (Rho);
          if ( RayTracingHeating ) {
          	/* Function for calculating the midplane stellar heating */
            ComputeRayTracingHeating(Rho, bsys);
          }
          if (( !RadTransport) && ( ExplicitRayTracingHeating )) {
            SubStep4_Explicit_Irr(dt);
            ComputeNewEnergyField(Rho, EnergyNew);
          } else {
            SubStep4 (Rho, EnergyNew, bsys, dt);
          }
        }
	      ActualiseGas (Energy, EnergyNew);
      }
      if ( HydroOn ) {
        Transport (Rho, Vrad, Vtheta, Energy, Label, dt);
      }
      if ( debug ) {
      	int check_neg = 0;
      	int check_zero = 0;
      	CheckField(Vrad, check_neg, check_zero, "Transport");
      	CheckField(Vtheta, check_neg, check_zero, "Transport");
      	check_neg = 1;
      	check_zero = 1;
      	CheckField(Rho, check_neg, check_zero, "Transport");
      	if ( Adiabatic ) {
      		CheckField(Energy, check_neg, check_zero, "Transport");
      	}
      }
      ApplyBoundaryCondition (Vrad, Vtheta, Rho, Energy, bsys, dt);
      ComputeTemperatureField (Rho, Energy);
      mdcp = CircumPlanetaryMass (Rho, sys);
      exces_mdcp = mdcp - mdcp0;
    }
    PhysicalTime += dt;
    timestep_counter++;
    if ( (timestep_counter % 100) == 0 ) {
    	timestep_counter = 0;
    	masterprint("t = %g", PhysicalTime);
    }
  }
  masterprint ("\n");
}

void SubStep1 (Vrad, Vtheta, Rho, dt)
     PolarGrid *Vrad, *Vtheta, *Rho;
     real dt;
{
  int i, j, l, lim, ljm, ljp, nr, ns;
  boolean selfgravityupdate;
  extern boolean Evanescent;

  real *vrad, *vtheta, *rho;
  real *Pot, *Press;
  real *vradint, *vthetaint;
  real gradp, gradphi, vt2, dxtheta;
  real invdxtheta;
  real supp_torque=0.0;/* for imposed disk drift */

  nr = Vrad->Nrad;
  ns = Vrad->Nsec;
  rho = Rho->Field;
  vrad = Vrad->Field;
  vtheta = Vtheta->Field;
  vradint = VradInt->Field;
  vthetaint = VthetaInt->Field;
  Pot = Potential->Field;
  Press = Pressure->Field;
  /* In this substep we take into account  */
  /* the source part of Euler equations  */
  /* We evolve velocities with pressure gradients */
  /* gravitational forces and curvature terms */
#pragma omp parallel private(j,l,lim,ljm,ljp,dxtheta,vt2,gradp,gradphi,invdxtheta,supp_torque)
  {
#pragma omp for nowait
    for (i = 1; i < nr; i++) {
      for (j = 0; j < ns; j++) {
        l = j+i*ns;
        lim = l-ns;
        ljp = l+1;
        if ( j == ns-1 ) {
          ljp = i*ns;
        }
        gradp = (Press[l]-Press[lim])*2.0/(rho[l]+rho[lim])*InvDiffRmed[i];
        gradphi = (Pot[l]-Pot[lim])*InvDiffRmed[i];
        vt2 = vtheta[l]+vtheta[ljp]+vtheta[lim]+vtheta[ljp-ns];
        vt2 = vt2/4.0+Rinf[i]*OmegaFrame;
        vt2 = vt2*vt2;
        vradint[l] = vrad[l]+dt*(-gradp-gradphi+vt2*InvRinf[i]);
      }
    }
#pragma omp for
    for (i = 0; i < nr; i++) {
      supp_torque = IMPOSEDDISKDRIFT*0.5*pow(Rmed[i],-2.5+SIGMASLOPE);
      dxtheta = 2.0*PI/(real)ns*Rmed[i];
      invdxtheta = 1.0/dxtheta;
      for (j = 0; j < ns; j++) {
        l = j+i*ns;
        ljm = l-1;
        if ( j == 0 ) {
          ljm = i*ns+ns-1;
        }
        gradp = (Press[l]-Press[ljm])*2.0/(rho[l]+rho[ljm])*invdxtheta;
        if ( ZMPlus ) {
          gradp *= SG_aniso_coeff;
        }
        gradphi = (Pot[l]-Pot[ljm])*invdxtheta;
        vthetaint[l] = vtheta[l]-dt*(gradp+gradphi);
        vthetaint[l] += dt*supp_torque;
      }
    }
  }
  if ( SelfGravity ) {
    selfgravityupdate = YES;
    compute_selfgravity(Rho, VradInt, VthetaInt, dt, selfgravityupdate);
  }
  ComputeViscousTerms(VradInt, VthetaInt, Rho);
  UpdateVelocitiesWithViscosity(VradInt, VthetaInt, Rho, dt);
  if (!Evanescent) {
    ApplySubKeplerianBoundary(VthetaInt);
  }
}

void SubStep2 (Rho, Energy, dt)
     PolarGrid *Rho, *Energy;
     real dt;
{
  int i, j, l, lim, lip, ljm, ljp, nr, ns;
  real *vrad, *vtheta, *rho, *energy;
  real *vradnew, *vthetanew, *qt, *qr, *energyint;
  real dxtheta, invdxtheta;
  real dv;

  nr = Rho->Nrad;
  ns = Rho->Nsec;
  rho = Rho->Field;
  vrad = VradInt->Field;
  vtheta = VthetaInt->Field;
  qr = RhoInt->Field;
  vradnew = VradNew->Field;
  vthetanew = VthetaNew->Field;
  qt = TemperInt->Field;
  energy = Energy->Field;
  energyint = EnergyInt->Field;

#pragma omp parallel for private(j,dxtheta,l,lim,lip,ljm,ljp,dv)
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      lip = l+ns;
      ljp = l+1;
      if ( j == ns-1 ) {
        ljp = i*ns;
      }
      dv = vrad[lip]-vrad[l];
      if ( dv < 0.0 ) {
        qr[l] = CVNR*CVNR*rho[l]*dv*dv;
      } else {
        qr[l] = 0.0; 
      }
      dv = vtheta[ljp]-vtheta[l];
      if ( dv < 0.0 ) {
        qt[l] = CVNR*CVNR*rho[l]*dv*dv;
      } else {
	      qt[l] = 0.0;
      }
    }
  }
#pragma omp parallel private(l,lim,lip,ljm,ljp,j,dxtheta,invdxtheta)
  {
#pragma omp for nowait
    for (i = 1; i < nr; i++) {
      for (j = 0; j < ns; j++) {
	      l = j+i*ns;
	      lim = l-ns;
	      vradnew[l] = vrad[l]-dt*2.0/(rho[l]+rho[lim])*(qr[l]-qr[lim])*InvDiffRmed[i];
      }
    }
#pragma omp for
    for (i = 0; i < nr; i++) {
      dxtheta = 2.0*PI/(real)ns*Rmed[i];
      invdxtheta = 1.0/dxtheta;
      for (j = 0; j < ns; j++) {
	      l = j+i*ns;
	      ljm = l-1;
	      if ( j == 0 ) {
          ljm = i*ns+ns-1;
        }
	      vthetanew[l] = vtheta[l]-dt*2.0/(rho[l]+rho[ljm])*(qt[l]-qt[ljm])*invdxtheta;
      }
    }
    /* If gas disk is adiabatic, we add artificial viscosity as a source */
    /* term for advection of thermal energy polargrid */
    if (( Adiabatic ) && ( !RadiativeOnly )) {
#pragma omp for nowait
      for (i = 0; i < nr; i++) {
	      dxtheta = 2.0*PI/(real)ns*Rmed[i];
	      invdxtheta = 1.0/dxtheta;
	      for (j = 0; j < ns; j++) {
	        l = j+i*ns;
	        lip = l+ns;
	        ljp = l+1;
	        if ( j == ns-1 ) {
            ljp = i*ns;
          }
	        energyint[l] = energy[l] - \
	         	dt*qr[l]*(vrad[lip]-vrad[l])*InvDiffRsup[i] -	\
	         	dt*qt[l]*(vtheta[ljp]-vtheta[l])*invdxtheta;
	      }
      }
    }
  }
}

void SubStep3 (Rho, dt)
  // Input
  PolarGrid *Rho;
  real dt;
{ 
  // Declaration
  extern boolean Cooling;
  int i, j, l, nr, ns;
  real *dens, *pres, *energy, *energynew;
  real *divergence, *qplus, *qirr, *qdiv;
  real energypred, pressurepred, num, den;
  real qirr_temp = 0.0, div_temp = 0.0;

  // Assignment
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  dens = Rho->Field;
  pres = Pressure->Field;
  energy = EnergyInt->Field;
  energynew = EnergyNew->Field;
  divergence = DivergenceVelocity->Field;
  qplus = Qplus->Field;
  qdiv = QDiv->Field;
  if ( Irradiation ) {
    qirr = Qirr->Field;
  }

  // Function
  /* In this substep we take into account  */
  /* the source part of energy equation  */
  /* We evolve internal energy with  */
  /* compression/dilatation and heating terms */
  /* Now we can update energy with source terms from i=0 */
#pragma omp parallel for private(j, l, div_temp, qirr_temp, num, den, energypred, pressurepred)
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      if ( !RadiativeOnly ) {
        div_temp = divergence[l];
      }
      if ( !ImplicitRadiative ) {
        if ( Irradiation ) {
          qirr_temp = qirr[l];
        }
      }
      if ( Cooling ) {
        num = EnergyMed[i]*dt*dens[l]/SigmaMed[i] + CoolingTimeMed[i]*energy[l] + dt*CoolingTimeMed[i]*(qirr_temp + qplus[l]-QplusMed[i]*dens[l]/SigmaMed[i]);
        den = dt + CoolingTimeMed[i] + (ADIABATICINDEX-1.0)*dt*CoolingTimeMed[i]*div_temp;
        energynew[l] = num/den;
      } else {
        energypred = energy[l] + dt*(qirr_temp + qplus[l] - pres[l]*div_temp);
        pressurepred = (ADIABATICINDEX - 1.0)*energypred;
        qdiv[l] = 0.5*(pres[l]+pressurepred)*div_temp;
        energynew[l] = energy[l] + dt*(qirr_temp + qplus[l] - qdiv[l]);
      }
    }
  }

  // Debug
  if ( debug ) {
    int check_neg = 1;
    int check_zero = 1;
    CheckField(EnergyNew, check_neg, check_zero, "SubStep3");
  }
}


int ConditionCFL (Vrad, Vtheta, energy, deltaT, sys, bsys)
     PolarGrid *Vrad, *Vtheta, *energy;
     PlanetarySystem *sys;
	   BinarySystem *bsys;
     real deltaT;
{
  static real Vresidual[MAX1D], Vmoy[MAX1D];
  int i, j, l, ns, nr, lip, ljp;
  real invdt1, invdt2, invdt3, invdt4, cs, newdt, dt;
  int ideb, jdeb;
  real itdbg1, itdbg2, itdbg3, itdbg4, mdtdbg; /* debugging variables */
  real *vt, *vr, dxrad, dxtheta, dvr, dvt, viscr, visct;
  int nstar, nplan, ij;
  real *xs, *ys, *xp, *yp, *ms, *mp, *Ppair, Pmin, dist, dt_nbod;
  real *soundspeed;
  
  soundspeed = SoundSpeed->Field;
  ns = Vtheta->Nsec;
  nr = Vtheta->Nrad;
  vt = Vtheta->Field;
  vr = Vrad->Field;
  xs = bsys->x;
  ys = bsys->y;
  xp = sys->x;
  yp = sys->y;
  ms = bsys->mass;
  mp = sys->mass;
  nstar = bsys->nb;
  nplan = sys->nb;

  newdt = 1.0e30;
  ij = 0;
  Pmin = 1.0e30;
  Ppair = (real *)malloc (sizeof(real)*(nstar*nplan));

  for (i = 0; i < nr; i++) {
    dxrad = Rsup[i]-Rinf[i];
    dxtheta=Rmed[i]*2.0*PI/(real)ns;
    Vmoy[i] = 0.0;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      Vmoy[i] += vt[l];
    }
    Vmoy[i] /= (real)ns;
  }
  for (i = One_or_active; i < Max_or_active; i++) {
    dxrad = Rsup[i]-Rinf[i];
    dxtheta=Rmed[i]*2.0*PI/(real)ns;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      if ( FastTransport ) {
	      Vresidual[j] = vt[l]-Vmoy[i];  /* FARGO algorithm */
      } else {
	      Vresidual[j] = vt[l];	       /* Standard algorithm */
      }
    }
    Vresidual[ns]=Vresidual[0];
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      lip = l+ns;
      ljp = l+1;
      if ( j == ns-1 ) {
      	ljp = i*ns;
      }
      cs = soundspeed[l];
      invdt1 = cs/(min2(dxrad,dxtheta));
      invdt2 = fabs(vr[l])/dxrad;
      invdt3 = fabs(Vresidual[j])/dxtheta;
      dvr = vr[lip]-vr[l];
      dvt = vt[ljp]-vt[l];
      if ( dvr >= 0.0 ) {
        dvr = 1.0e-10;
      } else {
        dvr = -dvr;
      }
      if ( dvt >= 0.0 ) {
        dvt = 1.0e-10;
      } else {
        dvt = -dvt;
      }
      invdt4 = max2(dvr/dxrad,dvt/dxtheta);
      invdt4*= 4.0*CVNR*CVNR;
      if ( HydroOn ) {
        dt = CFLSECURITY/sqrt(invdt1*invdt1+invdt2*invdt2+invdt3*invdt3+invdt4*invdt4);
      } else {
        dt = CFLSECURITY/sqrt(invdt1*invdt1);
      }
      if ( dt < newdt ) {
	      newdt = dt;
	      if ( debug == YES ) {
	        ideb = i;
	        jdeb = j;
	        itdbg1 = 1.0/invdt1; itdbg2=1.0/invdt2; itdbg3=1.0/invdt3; itdbg4=1.0/invdt4;
	        mdtdbg = newdt;
	        viscr = dxrad/dvr/4.0/CVNR/CVNR;     
	        visct = dxtheta/dvt/4.0/CVNR/CVNR;
	      }
      }  
    }
  }

  if (( BinaryOn  ) && ( LiveBodies )) {
    for (i = 0; i < nstar; i++) {
      for (j = 0; j < nplan; j++) {
        dist = pow((xs[i] - xp[j]), 2.0)+pow((ys[i] - yp[j]), 2.0);
        Ppair[ij] = 2.0*PI*sqrt(dist*sqrt(dist)/G/(ms[i] + mp[j]));
        if ( Ppair[ij] < Pmin ) {
          Pmin = Ppair[ij];
        }
		    ij++;
      }
    }
    dt_nbod = Pmin/NBODYSECURITY;
	  if ( dt_nbod < newdt ) {
      newdt = dt_nbod;
    }
  }
  for (i = Zero_or_active; i < MaxMO_or_active; i++) {
    dt = 2.0*PI*CFLSECURITY/(real)NSEC/fabs(Vmoy[i]*InvRmed[i]-Vmoy[i+1]*InvRmed[i+1]);
    if ( dt < newdt ) {
    	newdt = dt;
    }
  }
  // if ( Adiabatic ) {
  //   real dt_e;
  //   dt_e = EnergyConditionCFL(energy);
  //   if ( dt_e < newdt ) {
  //     newdt = dt_e;
  //   }
  // }
  if ( deltaT < newdt ) {
    newdt = deltaT;
  }
  // if (debug == YES) {
  //   printf ("Timestep control information for CPU %d: \n", CPU_Rank);
  //   printf ("Most restrictive cell at i=%d and j=%d\n", ideb, jdeb);
  //   printf ("located at radius Rmed         : %g\n", Rmed[ideb]);
  //   printf ("Sound speed limit              : %g\n", itdbg1);
  //   printf ("Radial motion limit            : %g\n", itdbg2);
  //   printf ("Residual circular motion limit : %g\n", itdbg3);
  //   printf ("Viscosity limit                : %g\n", itdbg4);
  //   printf ("   Arise from r with limit     : %g\n", viscr);
  //   printf ("   and from theta with limit   : %g\n", visct);
  //   printf ("N-body pair limit              : %g\n", dt_nbod);
  //   printf ("Energy Equation Evolution limit: %g\n", dt_e);
  //   printf ("Limit time step for this cell  : %g\n", mdtdbg);
  //   printf ("Limit time step adopted        : %g\n", newdt);
  //   if (newdt < mdtdbg) {
  //     printf ("Discrepancy arise either from shear.\n");
  //     printf ("or from the imposed DT interval.\n");
  //   }
  // }
  return (int)(ceil(deltaT/newdt));
}


void ComputeSoundSpeed (Rho, Energy)
     PolarGrid *Rho;
     PolarGrid *Energy;
{
  int i, j, l, nr, ns;
  real *dens, *energ, *cs;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  dens = Rho->Field;
  energ = Energy->Field;
  cs = SoundSpeed->Field;

  for ( i = 0; i < nr; i++ ) {
    for ( j = 0; j < ns; j++ ) {
      l = i*ns + j;
      if ( !Adiabatic ) {
        cs[l] = AspectRatio(Rmed[i])*sqrt(G*1.0/Rmed[i])*pow(Rmed[i], FLARINGINDEX);
      } else {
        cs[l] = sqrt( ADIABATICINDEX*(ADIABATICINDEX-1.0)*energ[l]/dens[l] );
      }
    }
  }

  // Debug
  if ( debug ) {
  	int check_neg = 1;
    int check_zero = 1;
    CheckField(SoundSpeed, check_neg, check_zero, "ComputeSoundSpeed");
  }
}

void ComputePressureField (Rho, Energy)
     PolarGrid *Rho;
     PolarGrid *Energy;
{
	int i, j, l, nr, ns;
	real *dens, *pres, *energ, *cs;

	nr = Rho->Nrad;
	ns = Rho->Nsec;
	dens = Rho->Field;
	pres = Pressure->Field;
	energ = Energy->Field;
	cs = SoundSpeed->Field;

	for ( i = 0; i < nr; i++ ) {
		for ( j = 0; j < ns; j++ ) {
			l = i*ns + j;
			if ( !Adiabatic ) {
				pres[l] = dens[l]*cs[l]*cs[l]; 
	/* since SoundSpeed is not updated
	from initialization, cs remains
	axisymmetric */
			} else {
				pres[l] = (ADIABATICINDEX - 1.0)*energ[l];
			}
		}
	}
	// Debug
	if ( debug ) {
		int check_neg = 1;
  	int check_zero = 1;
  	CheckField(Pressure, check_neg, check_zero, "ComputePressureField");
	}
}


void ComputeTemperatureField (Rho, energy)
     PolarGrid *Rho, *energy;
{
  int i, j, l, nr, ns;
  real *dens, *pres, *energ, *temp;

  nr = Rho->Nrad;
  ns = Rho->Nsec;
  dens = Rho->Field;
  pres = Pressure->Field;
  energ = energy->Field;
  temp = Temperature->Field;

  for ( i = 0; i < nr; i++ ) {
    for ( j = 0; j < ns; j++ ) {
      l = i*ns + j;
      if ( !Adiabatic ) {
	      temp[l] = MU/R* pres[l]/dens[l];
      } else {
	      temp[l] = MU/R*(ADIABATICINDEX-1.0)*energ[l]/dens[l];
      }
    }
  }
  // Debug
  if ( debug ) {
  	int check_neg = 1;
    int check_zero = 1;
    CheckField(Temperature, check_neg, check_zero, "ComputeTemperatureField");
  }
}


real CircumPlanetaryMass (Rho, sys)
     PolarGrid *Rho;
     PlanetarySystem *sys;
{
  int i, j, l, ns;
  real xpl, ypl;
  real dist, mdcplocal, mdcptotal;
  real *dens, *abs, *ord;
  ns = Rho->Nsec;
  dens = Rho->Field;
  abs = CellAbscissa->Field;
  ord = CellOrdinate->Field;
  xpl = sys->x[0];
  ypl = sys->y[0];

  mdcplocal = 0.0;
  mdcptotal = 0.0;

  if ( FakeSequential && ( CPU_Rank > 0 )) {
    MPI_Recv (&mdcplocal, 1, MPI_DOUBLE, CPU_Rank-1, 0, MPI_COMM_WORLD, &stat);
  }
  for (i = Zero_or_active; i < Max_or_active; i++) {
    for (j = 0; j < ns; j++) {
      l = i*ns + j;
      dist = sqrt ((abs[l] - xpl)*(abs[l] - xpl) + \
		    (ord[l] - ypl)*(ord[l] - ypl));
      if ( dist < HillRadius ) {
        mdcplocal += Surf[i]*dens[l];
      }
    }
  }
  if ( FakeSequential ) {
     if ( CPU_Rank < CPU_Number-1 ) {
       MPI_Send (&mdcplocal, 1, MPI_DOUBLE, CPU_Rank+1, 0, MPI_COMM_WORLD);
     }
  } else {
    MPI_Allreduce (&mdcplocal, &mdcptotal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
  if ( FakeSequential ) {
    MPI_Bcast (&mdcplocal, 1, MPI_DOUBLE, CPU_Number-1, MPI_COMM_WORLD);
    mdcptotal = mdcplocal;
  }
  return mdcptotal;
}


void SolveDiscOrbits (TimeStep,Rho, Vrad, Vtheta, e_cell, per_cell)
  int TimeStep;
  PolarGrid *Rho;
  PolarGrid *Vrad;
  PolarGrid *Vtheta;
  PolarGrid *e_cell;
  PolarGrid *per_cell;
{
  int i, j, l, ns, nr;
  real *dens, *vr, *vth, *ec, *pc;
  real ed, edlocal, pd, pdlocal, Md, Mdlocal, etemp;
  real xc, yc, vxc, vyc, ang, Axc, Ayc, d, mu;
  FILE *output;
  char name[256];

  ns = Rho->Nsec;
  nr = Rho->Nrad;
  dens = Rho->Field;
  vr = Vrad->Field;
  vth = Vtheta->Field;
  ec = e_cell->Field;
  pc = per_cell->Field;
  
#pragma omp parallel for
  for (i = 0; i < (nr+1)*ns; i++) {
    ec[i] = 0.0;
    pc[i] = 0.0;
  }

  edlocal = ed = pdlocal = pd = Mdlocal = Md = 0.0;
  
  /* -- Calculates mass of whole gas disk between Rmin and Rmax -- */

  if ( FakeSequential && ( CPU_Rank > 0 )) {
    MPI_Recv (&Mdlocal, 1, MPI_DOUBLE, CPU_Rank-1, 0, MPI_COMM_WORLD, &stat);
  }

#pragma omp parallel for private (j,l)
  for (i = Zero_or_active; i < Max_or_active; i++ ) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      Mdlocal += Surf[i]*dens[l];
    }
  }
  if ( FakeSequential ) {
     if ( CPU_Rank < CPU_Number-1 ) {
       MPI_Send (&Mdlocal, 1, MPI_DOUBLE, CPU_Rank+1, 0, MPI_COMM_WORLD);
     }
  } else {
    MPI_Allreduce (&Mdlocal, &Md, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
  if ( FakeSequential ) {
    MPI_Bcast (&Mdlocal, 1, MPI_DOUBLE, CPU_Number-1, MPI_COMM_WORLD);
    Md = Mdlocal;
  }

  /* -- Calculates individual cell's orbital elements -- */

  if ( FakeSequential && ( CPU_Rank > 0 )) {
    MPI_Recv (&edlocal, 1, MPI_DOUBLE, CPU_Rank-1, 0, MPI_COMM_WORLD, &stat);
    MPI_Recv (&pdlocal, 1, MPI_DOUBLE, CPU_Rank-1, 0, MPI_COMM_WORLD, &stat);
  }

#pragma omp parallel for private(j,l, ang, xc, yc, vxc, vyc, d, mu, Axc, Ayc, etemp)
  for (i = Zero_or_active; i < Max_or_active; i++ ) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      ang = (real)j/(real)ns*2.0*PI;
      xc = Rmed[i]*cos(ang);
      yc = Rmed[i]*sin(ang);
      vxc = vr[l]*cos(ang)-vth[l]*sin(ang);
      vyc = vr[l]*sin(ang)+vth[l]*cos(ang);
      d = sqrt(xc*xc+yc*yc);

      mu = G*(Surf[i]*dens[l] + 1.0); /* Gravitational parameter */
      
      Axc = xc*vyc*vyc/mu-yc*vxc*vyc/mu-xc/d;
      Ayc = yc*vxc*vxc/mu-xc*vxc*vyc/mu-yc/d;
      etemp = Axc*Axc+Ayc*Ayc;
      if ( etemp < 0.0 ) {
      	etemp = 1.0E-20;
      }
      ec[l] = sqrt(etemp);
      if ( ec[l] < 1.0E-6 ) {
      	ec[l] = 0.0;
      }
      edlocal += dens[l]*ec[l]*Surf[i];

      if ( ec[l] != 0.0 ) {
        pc[l]=atan2(Ayc,Axc);
      } else {
        pc[l]=atan2(yc,xc);
      }
      if ( fabs(pc[l]) < 1.0E-6 ) {
      	pc[l] = 0.0;
      }
      pdlocal += dens[l]*pc[l]*Surf[i];
    } 
  }
  if ( FakeSequential ) {
     if ( CPU_Rank < CPU_Number-1 ) {
       MPI_Send (&edlocal, 1, MPI_DOUBLE, CPU_Rank+1, 0, MPI_COMM_WORLD);
       MPI_Send (&pdlocal, 1, MPI_DOUBLE, CPU_Rank+1, 0, MPI_COMM_WORLD);
     }
  } else {
    MPI_Allreduce (&edlocal, &ed, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce (&pdlocal, &pd, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
  if ( FakeSequential ) {
    MPI_Bcast (&edlocal, 1, MPI_DOUBLE, CPU_Number-1, MPI_COMM_WORLD);
    MPI_Bcast (&pdlocal, 1, MPI_DOUBLE, CPU_Number-1, MPI_COMM_WORLD);
    ed = edlocal;
    pd = pdlocal;
  }

  /* -- Average over mass of disk -- */

  ed = ed/Md;
  pd = pd/Md;

  /* -- Print out results -- */

  if ( !CPU_Master ) {
  	return;
  }
  fflush (stdout);
  sprintf (name, "%sdiscorbit.dat", OUTPUTDIR);
  output = fopen (name, "a");
  if ( output == NULL ) {
    fprintf (stderr, "Can't write 'discorbit.dat' file. Aborting.\n");
    prs_exit (1);
  }
  fprintf (output, "%d\t%#.18g\t%#.18g\t%#.18g\t%#.18g\n", TimeStep, PhysicalTime, Md, ed, pd);
  fclose (output);
  fflush (stdout);
  sprintf (name, "%sdiscmass.dat", OUTPUTDIR);
  output = fopen (name, "a");
  if ( output == NULL ) {
    fprintf (stderr, "Can't write 'discmass.dat' file. Aborting.\n");
    prs_exit (1);
  }
  fprintf (output, "%d\t%#.18g\t%#.18g\t%#.18g\t%#.18g\n", TimeStep, PhysicalTime, Md, mlinner*Surf[1], mlouter*Surf[GLOBALNRAD-1]);
  fclose (output);
  fflush (stdout);
}   


void DiscMassTaperCalc (gas_density, timestep)
  PolarGrid *gas_density;
  real timestep;
{
  int nr, ns, i, j, l;
  real *dens;

  dens = gas_density->Field;
  nr = gas_density->Nrad;
  ns = gas_density->Nsec;

  #pragma omp parallel for private(j,l)
  for (i = 0; i < nr; i++) {
    SigmaDMT[i] *= exp(-timestep/DMTTAU);
    if (Adiabatic)
      EnergyDMT[i] *= exp(-timestep/DMTTAU);
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      dens[l] *= exp(-timestep/DMTTAU);
    }
  }
}


void ApplyFloor (gas_density, gas_energy)
  // Input
  PolarGrid *gas_density;
  PolarGrid *gas_energy;
{
  // Declaration
  int nr, ns, i, j, l;
  real *dens, *energy, Tfloor, EnergySafetyFloor;

  // Assignment
  nr = gas_density->Nrad;
  ns = gas_density->Nsec;
  dens = gas_density->Field;
  energy = gas_energy->Field;

  // Constants
  Tfloor = 1.2E-5;

  // Function
#pragma omp parallel for private(j,l, EnergySafetyFloor)
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      if (dens[l] < SigmaSafetyFloor) {
        dens[l] = SigmaSafetyFloor;
      }
      EnergySafetyFloor = dens[l]*Tfloor*R/(MU*(ADIABATICINDEX-1.0));
      if (energy[l] < EnergySafetyFloor) {
        energy[l] = EnergySafetyFloor;
      }
    }
  }
}

void ComputeQplus (Rho)
  // Input
  PolarGrid *Rho;
{
  // Declaration
  int i, j, l, nr, ns;
  int lip, li2p;
  real r, rip, ri2p, qpip, qpi2p;
  real *dens;
  real *divergence, *Trr, *Trp, *Tpp, *qplus;
  real viscosity;

  // Assignment
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  dens = Rho->Field;
  divergence = DivergenceVelocity->Field;
  qplus = Qplus->Field;
  Trr = TAURR->Field;
  Trp = TAURP->Field;
  Tpp = TAUPP->Field;

  // Function
  if ( !RadiativeOnly ) {
#pragma omp for private (j, l, viscosity)
      /* We calculate the heating source term Qplus from i=1 */
      for (i = 1; i < nr; i++) { /* Trp defined from i=1 */
        viscosity = FViscosity (Rmed[i]); /* Global_Bufarray holds cs */
        for (j = 0; j < ns; j++) {
          l = j+i*ns;
          if ( viscosity != 0.0 ) {
            qplus[l] = 0.5/viscosity/dens[l]*(Trr[l]*Trr[l] + \
             Trp[l]*Trp[l] +  \
             Tpp[l]*Tpp[l] );
            qplus[l] += (2.0/9.0)*viscosity*dens[l]*divergence[l]*divergence[l];
          } else
            qplus[l] = 0.0;
        }
      }
      /* We calculate the heating source term Qplus for i=0 */
      i = 0;
      r    = Rmed[i];
      rip  = Rmed[i+1];
      ri2p = Rmed[i+2];
      for (j = 0; j < ns; j++) {
        l = j+i*ns;
        lip = l+ns;
        li2p = lip+ns;
        qpip = qplus[lip];   /* qplus(i=1,j) */
        qpi2p = qplus[li2p]; /* qplus(i=2,j) */
        if ( viscosity != 0.0 ) {
          /* power-law extrapolation */
          qplus[l] = qpip*exp( log(qpip/qpi2p) * log(r/rip) / log(rip/ri2p) );
        } else {
          qplus[l] = 0.0;
        }
      }
    } else {
#pragma omp for  
      for (i = 0; i < nr; i++) {
        for (j = 0; j < ns; j++) { 
          l = j+i*ns;
          qplus[l] = 0.0;
        }
      }
    }

    // Debug
    if ( debug ) {
      int check_neg = 1;
      int check_zero = 0;
      CheckField(Qplus, check_neg, check_zero, "SubStep3");
    }
}

real EnergyConditionCFL(Energy)
  // Input
  PolarGrid *Energy;
{
  // Declaration
  int i, j, l, nr, ns;
  real *energy, *qplus, *qirr, *pressure, *divergence;
  real dedt_abs, dedt, new_dt = 1.0E30, old_dt, factor;

  // Assignment
  nr = Energy->Nrad;
  ns = Energy->Nsec;
  energy = Energy->Field;
  if ( Irradiation ) {
    qirr = Qirr->Field;
  }
  if ( HydroOn ) {
    divergence = DivergenceVelocity->Field;
    pressure = Pressure->Field;
  }
  if ( !RadiativeOnly ) {
    qplus = Qplus->Field;
  }

  // Constants
  factor = 1.0;

  // Function
#pragma omp parallel for private(j, l, dedt, dt, newdt)
  for (i = One_or_active; i < MaxMO_or_active; i++) {
    for (j = 0; j < ns; j++) {
      l = j + i*ns;
      dedt_abs = 0.0;
      dedt = 0.0;
      if ( Irradiation ) {
        dedt_abs += fabs(qirr[l]);
        dedt += qirr[l];
      }
      if ( HydroOn) {
        dedt_abs += fabs(pressure[l]*divergence[l]);
        dedt -= pressure[l]*divergence[l];
      }
      if ( !RadiativeOnly) {
        dedt_abs += fabs(qplus[l]);
        dedt += qplus[l];
      }
      if ( dedt_abs == 0.0 ) {
      	old_dt = 1.0E30;
      } else {
      	old_dt = energy[l]/dedt_abs;
      	old_dt /= factor;
      }
      if ( old_dt < new_dt ) {
        new_dt = old_dt;
      }
    }
  }
  
  // Output
  return new_dt;
}

real RadCoolConditionCFL(Energy)
  // Input
  PolarGrid *Energy;
{
  // Declaration
  int i, j, l, nr, ns;
  real *energy, *qminus;
  real dedt_abs, new_dt=1.0E30, old_dt, factor;

  // Assignment
  nr = Energy->Nrad;
  ns = Energy->Nsec;
  energy = Energy->Field;
  qminus = Qminus->Field;

  // Constants
  factor = 1.0;

  // Function
#pragma omp parallel for private(j, l, dedt, dt, newdt)
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j + i*ns;
      dedt_abs = fabs(qminus[l]);
      old_dt = energy[l]/dedt_abs;
      old_dt /= factor;
      if ( old_dt < new_dt ) {
        new_dt = old_dt;
      }
    }
  }
  
  // Output
  return new_dt;
}


void ComputeQterms (Rho, Vrad, Vtheta, Energy, bsys)
  // Input
  PolarGrid *Rho;
  PolarGrid *Vrad;
  PolarGrid *Vtheta;
  PolarGrid *Energy;
  BinarySystem *bsys;
{
  // Function
  if (( Irradiation ) || ( RadTransport ) || ( RadCooling )) {
    ComputeDiscHeight (bsys);
    ComputeTemperatureField (Rho, EnergyInt);
    ComputeRKappa (Rho);
  }
  if ( !RadiativeOnly ) {
    ComputeViscousTerms (Vrad, Vtheta, Rho);
  }
  if ( RadTransport ) {
    ResetTempSourcesSinks();
  }
  if ( Irradiation ) {
    ComputeQirr(Rho, bsys);
  }
  if ( RadCooling ) {
    ComputeQminus(Rho);
  }
  ComputeQplus(Rho);
}


void SubStep3_RadCool(Rho, Energy, dt_hydro, dt_energy)
  // Input
  PolarGrid *Rho;
  PolarGrid *Energy;
  real dt_hydro;
  real dt_energy;
{
  // Declaration
  int i, j, l, nr, ns, nstep;
  real *density, *energy, *temperature, *qminus, *taueff;
  real dt_remainder=0.0, dt;
  real *densitycv;
  real constant1;
  real term1, term2;

  // Assignment
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  density = Rho->Field;
  energy = Energy->Field;
  temperature = Temperature->Field;
  qminus = Qminus->Field;
  taueff = OpticalDepthEff->Field;
  densitycv = (real *)malloc(nr*ns*sizeof(real));

  // Constants
  constant1 = 2.0*STEFANK;

  // Function
  /* Convert energy density to temperature */
#pragma omp parallel for private(j, l)
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      densitycv[l] = density[l]*CV;
      temperature[l] = energy[l]/densitycv[l];
    }
  }
  /* Evolve temperature field with analytical cooling term */
  if ( AnalyticCooling ) {
    for (i = 0; i < nr; i++) {
      for (j = 0; j < ns; j++) {
        l = j + i*ns;
        term1 = pow(temperature[l], -3.0);
        term2 = 3.0*constant1*dt_hydro/(densitycv[l]*taueff[l]);
        temperature[l] = pow(term1 + term2, -1.0/3.0);
      }
    }
  } else {
    /* Find number of cooling timesteps per hydro timestep */
    if ( dt_energy >  dt_hydro ) { 
      RadCoolTimeStepsCFL = 1;
      dt = dt_hydro;
    } else {
      RadCoolTimeStepsCFL = (int)(dt_hydro/dt_energy);
      dt = dt_energy;
      dt_remainder = dt_hydro - (real)(RadCoolTimeStepsCFL*dt_energy);
    }
    dt_hydro_av += dt_hydro;
    dt_energy_av += dt_energy;
    nsteps_av += RadCoolTimeStepsCFL;
    /* Carry out RadCoolTimeStepsCFL sub-cycles, with timestep dt */
    for (nstep = 0; nstep < RadCoolTimeStepsCFL; nstep++) {
  #pragma omp parallel for private(j, l)
      /* Evolve temperature field with local radiative cooling */
      for (i = 0; i < nr; i++) {
        for (j = 0; j < ns; j++) {
          l = j+i*ns;
          temperature[l] -= dt*qminus[l]/densitycv[l];
        }
      }
      /* Update Rosseland Mean Opacity with new temperature */
      ComputeRKappa(Rho);
      /* Calculate new local radiative cooling */
      ComputeQminus(Rho);
    }
    /* Update temperature to time level n+1 with last dt_remainder timestep */
    if ( dt_remainder > 0.0 ) {
  #pragma omp parallel for private(j, l)
      for (i = 0; i < nr; i++) {
        for (j = 0; j < ns; j++) {
          l = j+i*ns;
          temperature[l] -= dt_remainder*qminus[l]/densitycv[l];
        }
      }
    }
  }
  /* Convert temperature back to energy density */
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      energy[l] = temperature[l]*densitycv[l];
    }
  }
  free(densitycv);

  // Debug
  if ( debug ) {
    int check_zero = 1;
    int check_neg = 1;
    CheckField(Energy, check_neg, check_zero, "SubStep3_RadCool");
  }
}