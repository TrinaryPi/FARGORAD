/* C Header
	* @filename        : Pframeforce.c
	* @author          : Frederic Masset
	* @last_modified_by: trinarypi
	* @last_modified   : 2017/06/26 16:57
	* @description     :
*/
#include "mp.h"

extern boolean AllowAccretion, Corotating, Indirect_Term, Cooling, BinaryOn, Restart;
extern Pair DiskOnPrimaryAcceleration;
extern Pair DiskOnBinaryAcceleration;
extern real Abinary;
static Pair IndirectTerm;
static real q0[MAX1D], q1[MAX1D], PlanetMasses[MAX1D];
static real vt_int[MAX1D], vt_cent[MAX1D];
extern boolean VarDiscHeight, RadCooling, RadiativeOnly;
extern boolean StellarSmoothing;

void ComputeIndirectTerm () 
{
  if ( BinaryOn ) {
    IndirectTerm.x = -DiskOnBinaryAcceleration.x;
    IndirectTerm.y = -DiskOnBinaryAcceleration.y;
  } else {
    IndirectTerm.x = -DiskOnPrimaryAcceleration.x;
    IndirectTerm.y = -DiskOnPrimaryAcceleration.y;
  }
  if ( Indirect_Term == NO ) {
    IndirectTerm.x = 0.0;
    IndirectTerm.y = 0.0;
  }
}

/* Below : work in non-rotating frame */
/* centered on the Primary or ...Binary centre of mass*/
void FillForcesArrays (bsys, sys, Rho, Energy)
     BinarySystem *bsys;
     PlanetarySystem *sys;
     PolarGrid *Rho, *Energy;  
{
  int i, j, l, nr, ns, k, s, NbPlan, NbStar;
  real x, y, ang, d, dsmooth;
  real xs, ys, ms, M, rs;
  real xp, yp, RRoche, smooth, mp;
  real PlanD, *Pot, pot, smoothing;
  real InvPlanD3, InvD;
  real drx, dry, A;
  real Rstar[2], q[2], ssmooth[2];
  
  Pot = Potential->Field;
  nr = Potential->Nrad;
  ns = Potential->Nsec;
  NbPlan = sys->nb;
  NbStar = bsys->nb;
  M = bsys->mass[0] + bsys->mass[1];

  // Constants
  q[0] = bsys->mass[0]/bsys->mass[1];
  q[1] = 1.0/q[1];
  
  /* Indirect term star on gas here */

  ComputeIndirectTerm();

#pragma omp parallel for
  for (i = 0; i < (nr+1)*ns; i++) {
    Pot[i] = 0.0;
  }

  /* -- Gravitational potential from planet on gas -- */

  for (k = 0; k < NbPlan; k++) { 
    xp = sys->x[k];
    yp = sys->y[k];
    mp = sys->mass[k]*MassTaper;
    PlanD = sqrt(xp*xp + yp*yp);
    InvPlanD3 = 1.0/PlanD/PlanD/PlanD;
    RRoche = PlanD*pow((1.0/3.0*mp), 1.0/3.0);
    if ( RocheSmoothing ) {
      smoothing = RRoche*ROCHESMOOTHING;
    } else {
      if ( VarDiscHeight ) {
        smoothing = compute_varheight_smoothing(xp,yp);
      } else {
        smoothing = compute_smoothing (PlanD);
      }
    }
    smooth = smoothing*smoothing;
#pragma omp parallel for private(InvD, j, l, ang, x, y, d, dsmooth, pot)
    for (i = 0; i < nr; i++) {
      InvD = 1.0/Rmed[i];
      for (j = 0; j < ns; j++) {
	      l = j+i*ns;
	      ang = (real)j/(real)ns*2.0*PI;
	      x = Rmed[i]*cos(ang);
	      y = Rmed[i]*sin(ang);
	      d = (x-xp)*(x-xp)+(y-yp)*(y-yp);
	      dsmooth = sqrt(d + smooth);
        /* --  Direct term from planet -- */
	      pot = (-1.0)*G*mp/dsmooth; 
        Pot[l] += pot;
        pot = 0.0;
        if ( BinaryOn ) {
          for (s = 0; s < NbStar; s++) {
	          if (( Indirect_Term ) && ( bsys->FeelOthers[s] )) {
              xs = bsys->x[s];
              ys = bsys->y[s];
              ms = bsys->mass[s];
              d = (xs-xp)*(xs-xp) + (ys-yp)*(ys-yp);
              d = sqrt(d);
              InvPlanD3 = 1.0/d/d/d;
              /* Indirect term from planets on stars */
              pot -= G*ms*mp*((xs-xp)*x + (ys-yp)*y)*InvPlanD3;  
            }
          }
        } else {
          /* Indirect term from planet  */
          pot += G*mp*InvPlanD3*(x*xp + y*yp); 
        }
        Pot[l] += pot;
      }
    }
  }
  
  /* -- Indirect Term From Stars -- */
#pragma omp parallel for private(InvD, j, l, ang, x, y)
  for (i = 0; i < nr; i++) {
    InvD = 1.0/Rmed[i];
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      ang = (real)j/(real)ns*2.0*PI;
      x = Rmed[i]*cos(ang);
      y = Rmed[i]*sin(ang);
      Pot[l] -= IndirectTerm.x*x + IndirectTerm.y*y;
    }
  }
  
  /* -- Direct term from star(s) on gas -- */
  if ( BinaryOn ) {
    drx = bsys->x[0] - bsys->x[1];
    dry = bsys->y[0] - bsys->y[1];
    A = sqrt(drx*drx + dry*dry);
    for (s = 0; s < NbStar; s++) {
      xs = bsys->x[s];
      ys = bsys->y[s];
      rs = sqrt(xs*xs + ys*ys);
      ms = bsys->mass[s];
      if ( StellarSmoothing == NO ) {
        smoothing = 0.0;
        ssmooth[s] = 0.0;
        Rstar[s] = (0.49*pow(q[s], 2.0/3.0))/(0.6*pow(q[s], 2.0/3.0) + log(1.0 + pow(q[s], 1.0/3.0)))*A;
      } else if ( RocheSmoothing ) {
        smoothing = Abinary*pow(ms/(3.0*M), 1.0/3.0)*bsys->smooth[s];
        ssmooth[s] = smoothing*smoothing;
      } else {
        smoothing = (bsys->smooth[s]) * AspectRatio(rs) * pow(rs, 1.0+FLARINGINDEX);
        ssmooth[s] = smoothing*smoothing;
      }

#pragma omp parallel for private(InvD, j, l, ang, x, y, pot, d, dsmooth)
      for (i = 0; i < nr; i++) {
        InvD = 1.0/Rmed[i];
        for (j = 0; j < ns; j++) {
          l = j+i*ns;
          ang = (real)j/(real)ns*2.0*PI;
          x = Rmed[i]*cos(ang);
          y = Rmed[i]*sin(ang);
          d = (x-xs)*(x-xs) + (y-ys)*(y-ys);
          if ( StellarSmoothing == NO ) {
            if (sqrt(d) < Rstar[s]) {
              pot = (-1.0)*G*bsys->mass[s]*(3.0*pow(Rstar[s], 2.0) - d)/(2.0*pow(Rstar[s], 3.0));
            } else {
              pot = (-1.0)*G*bsys->mass[s]/sqrt(d);
            }
          } else {
            dsmooth = sqrt(d + ssmooth[s]);
            pot = (-1.0)*G*ms/dsmooth;
          }
          Pot[l] += pot;
        }
      }
    }
  } else {
#pragma omp parallel for private(InvD, j, l, ang, x, y, pot)
    for (i = 0; i < nr; i++) {
      InvD = 1.0/Rmed[i];
      for (j = 0; j < ns; j++) {
        l = j+i*ns;
        ang = (real)j/(real)ns*2.0*PI;
        x = Rmed[i]*cos(ang);
        y = Rmed[i]*sin(ang);
        pot = (-1.0)*G*1.0*InvD;
        Pot[l] += pot;
      }
    }
  }

  // Debug
  if ( debug ) {
  	int check_neg = 0;
    int check_zero = 0;
    CheckField(Potential, check_neg, check_zero, "");
  }
}

void AdvanceSystemFromDisk (force, Rho, Energy, sys, dt)
     Force *force;
     PlanetarySystem *sys;
     PolarGrid *Rho, *Energy;
     real dt;		       
{
  int NbPlanets, k;
  Pair gamma;
  real x, y, r, m, smoothing;
  NbPlanets = sys->nb;
  
  for (k = 0; k < NbPlanets; k++) {
    if ( sys->FeelDisk[k] ) {
      m = sys->mass[k];
      x = sys->x[k];
      y = sys->y[k];
      r = sqrt(x*x + y*y);
      if ( RocheSmoothing ) {
	      smoothing = r*pow(m/3.,1./3.)*ROCHESMOOTHING;
      } else {
        if (VarDiscHeight) {
          smoothing = compute_varheight_smoothing(x,y);
        } else {
	        smoothing = compute_smoothing (r);
        }
      }
      gamma = ComputeAccel (force, Rho, x, y, smoothing, m);
      sys->vx[k] += dt * gamma.x;
      sys->vy[k] += dt * gamma.y;
      sys->vx[k] += dt * IndirectTerm.x;
      sys->vy[k] += dt * IndirectTerm.y;
    }
  }
}

void AdvanceSystemRK5 (bsys, sys, dt)
     BinarySystem *bsys;
     PlanetarySystem *sys;
     real dt;
{
  int i, n, nplan, nstar;
  boolean *feelothers, *binfeelothers;

  nplan = sys->nb;
  nstar = bsys->nb;
  
  n = nplan + nstar;
  
  for (i = 0; i < nstar; i++) {
    q0[i] = bsys->x[i];
    q0[i+n] = bsys->y[i];
    q0[i+2*n] = bsys->vx[i];
    q0[i+3*n] = bsys->vy[i];
    PlanetMasses[i] = bsys->mass[i];
  }
  for (i = 0; i < nplan; i++) {
    q0[i+nstar] = sys->x[i];
    q0[i+nstar+n] = sys->y[i];
    q0[i+nstar+2*n] = sys->vx[i];
    q0[i+nstar+3*n] = sys->vy[i];
    PlanetMasses[i+nstar] = sys->mass[i];
  }
  binfeelothers = bsys->FeelOthers;
  feelothers = sys->FeelOthers;
  RungeKunta (q0, dt, PlanetMasses, q1, nplan, nstar, feelothers, binfeelothers);
  for (i = 0; i < bsys->nb; i++) {
      bsys->x[i] = q1[i];
      bsys->y[i] = q1[i+n];
      bsys->vx[i] = q1[i+2*n];
      bsys->vy[i] = q1[i+3*n];
  }
  for (i = 0; i < sys->nb; i++) {
      sys->x[i] = q1[i+nstar];
      sys->y[i] = q1[i+n+nstar];
      sys->vx[i] = q1[i+2*n+nstar];
      sys->vy[i] = q1[i+3*n+nstar];
  }
  /*for (i = 1-(PhysicalTime >= RELEASEDATE); i < sys->nb; i++) {
    if ( !ForcedCircular ) {
      sys->x[i] = q1[i];
      sys->y[i] = q1[i+n];
      sys->vx[i] = q1[i+2*n];
      sys->vy[i] = q1[i+3*n];
    }
    else {
      x = sys->x[i];
      y = sys->y[i];
      theta = atan2(y,x);
      vx = sys->vx[i];
      vy = sys->vy[i];
      r = sqrt(x*x + y*y);
      omega = (-y*vx + x*vy)/r/r;
      dtheta = omega*dt;
      sys->x[i] = r*cos(theta+dtheta);
      sys->y[i] = r*sin(theta+dtheta);
      sys->vx[i]= vx*cos(dtheta+theta) - vy*sin(dtheta+theta);
      sys->vy[i]= vx*sin(dtheta+theta) + vy*cos(dtheta+theta);
    }
  }
  if (PhysicalTime < RELEASEDATE) {
    x = sys->x[0];
    y = sys->y[0];
    r = sqrt(x*x+y*y);
    theta = atan2(y,x);
    rdot = (RELEASERADIUS-r)/(RELEASEDATE-PhysicalTime);
    omega = sqrt((1.+sys->mass[0])/r/r/r);
    new_r = r + rdot*dt;
    denom = r-new_r;
    if (denom != 0.0) {
      dtheta = 2.*dt*r*omega/denom*(sqrt(r/new_r)-1.);
    } else {
      dtheta = omega*dt;
    }
    vx = rdot;
    vy = new_r*sqrt((1.+sys->mass[0])/new_r/new_r/new_r);
    sys->x[0] = new_r*cos(dtheta+theta);
    sys->y[0] = new_r*sin(dtheta+theta);
    sys->vx[0]= vx*cos(dtheta+theta) - vy*sin(dtheta+theta); 
    sys->vy[0]= vx*sin(dtheta+theta) + vy*cos(dtheta+theta); 
  }*/
}

void SolveOrbits (sys)
     PlanetarySystem *sys;
{
  int i, n;
  real x, y, vx, vy;

  n = sys->nb;

  for (i = 0; i < n; i++) {
    x = sys->x[i];
    y = sys->y[i];
    vx = sys->vx[i];
    vy = sys->vy[i];
    FindOrbitalElements(x, y, vx, vy, 1.0+sys->mass[i], i);
  }
} 

real ConstructSequence (u, v, n)
     real *u, *v;
     int n;
{
  int i;
  real lapl=0.0;

  for (i = 1; i < n; i++) {
    u[i] = 2.0*v[i]-u[i-1];
  }
  for (i = 1; i < n-1; i++) {
    lapl += fabs(u[i+1]+u[i-1]-2.0*u[i]);
  }

  return lapl;
}

void InitGasDensity (Rho)
     PolarGrid *Rho;
{
  int i, j, l, nr, ns;
  real *dens;

  dens = Rho->Field;
  nr = Rho->Nrad;
  ns = Rho->Nsec;

  FillSigma();
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      dens[l] = SigmaMed[i];
    }
  }
}

void InitGasEnergy (Energy)
     PolarGrid *Energy;
{
  int i, j, l, nr, ns;
  real *energy;

  energy = Energy->Field;
  nr = Energy->Nrad;
  ns = Energy->Nsec;

  FillEnergy();
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      energy[l] = EnergyMed[i];
    }
  }
}

void InitGasVelocities (Vr, Vt)
     PolarGrid *Vr, *Vt;
{
  extern boolean SGZeroMode;
  extern boolean SelfGravity;
  int i, j, l, nr, ns;
  real *vr, *vt, *pres, *cs;
  real  r, omega, ri;
  real viscosity, t1, t2, r1, r2;

  vr  = Vr->Field;
  vt  = Vt->Field;
  nr  = Vt->Nrad;
  ns  = Vt->Nsec;
  cs = SoundSpeed->Field;
  pres = Pressure->Field;  /* Pressure is already initialized: cf initeuler in SourceEuler.c ... */

  /* --------- */
  // Initialization of azimutal velocity with exact centrifugal balance
  /* --------- */
  if ( CentrifugalBalance ) {
    /* vt_int \equiv rOmega² = grad(P)/sigma +  \partial(phi)/\partial(r)  -  acc_sg_radial */
    mpi_make1Dprofile (pres, GLOBAL_bufarray);
    /* global axisymmetric pressure field, known by all cpus*/
    for (i = 1; i < GLOBALNRAD; i++) {
      vt_int[i] = ( GLOBAL_bufarray[i] - GLOBAL_bufarray[i-1] ) /	\
        (.5*(Sigma(GlobalRmed[i])+Sigma(GlobalRmed[i-1])))/(GlobalRmed[i]-GlobalRmed[i-1]) + \
        G*(1.0/GlobalRmed[i-1]-1.0/GlobalRmed[i])/(GlobalRmed[i]-GlobalRmed[i-1]);
    }
    /* Case of a disk with self-gravity */
    if ( SelfGravity ) { // Better test with CL rigid!
      if ( !SGZeroMode ) {
        mpi_make1Dprofile(SG_Accr, GLOBAL_AxiSGAccr);
      } else {
        GLOBAL_AxiSGAccr = SG_Accr;
      }
      for (i = 1; i < GLOBALNRAD; i++) {
        vt_int[i] -= ( (Radii[i] - GlobalRmed[i-1])*GLOBAL_AxiSGAccr[i] + \
          (GlobalRmed[i] - Radii[i])*GLOBAL_AxiSGAccr[i-1] ) / (GlobalRmed[i]-GlobalRmed[i-1]);
      }
    }
    for (i = 1; i < GLOBALNRAD; i++) {
      vt_int[i] = sqrt(vt_int[i]*Radii[i])-Radii[i]*OmegaFrame;
    }
    
    t1 = vt_cent[0] = vt_int[1]+.75*(vt_int[1]-vt_int[2]);
    r1 = ConstructSequence(vt_cent, vt_int, GLOBALNRAD);
    vt_cent[0] += .25*(vt_int[1]-vt_int[2]);
    t2 = vt_cent[0];
    r2 = ConstructSequence(vt_cent, vt_int, GLOBALNRAD);
    t1 = t1-r1/(r2-r1)*(t2-t1);
    vt_cent[0] = t1;
    ConstructSequence (vt_cent, vt_int, GLOBALNRAD);
    vt_cent[GLOBALNRAD] = vt_cent[GLOBALNRAD-1];
  }
  /* --------- */
  // Initialization with self-gravity, without exact centrifugal balance
  if ( SelfGravity && !CentrifugalBalance ) {
    init_azimutalvelocity_withSG (Vt);
  }
  /* --------- */
  if ( ViscosityAlpha ) {
    mpi_make1Dprofile (cs, GLOBAL_bufarray);
  }
  /* We calculate here the cooling time radial profile (Theo.c) */
  if ( Cooling ) {
    FillCoolingTime();
  }
    /* To fill qplus, one requires to calculate viscosity, hence cs if
       one uses alpha-viscosity */
  if (( Cooling ) || ( RadCooling )) {
    FillQplus();
  }
  
  for (i = 0; i <= nr; i++) {
    if ( i == nr ) {
      r = Rmed[nr-1];
      ri= Rinf[nr-1];
    } else {
      r = Rmed[i];
      ri= Rinf[i];
    }
    viscosity = FViscosity (r);
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      /* --------- */
      if ( !SelfGravity ) {
	      omega = sqrt(G*1.0/r/r/r);
	      vt[l] = omega*r*sqrt(1.0-pow(MeanDiscHeight[i]/r,2.0)*			\
			   pow(r,2.0*FLARINGINDEX)*			\
			   (1.+SIGMASLOPE-2.0*FLARINGINDEX) );
      }
      /* --------- */
      vt[l] -= OmegaFrame*r;
      if ( CentrifugalBalance ) {
	      vt[l] = vt_cent[i+IMIN];
      }
      if ( i == nr ) {
	      vr[l] = 0.0;
      } else {
	      vr[l] = IMPOSEDDISKDRIFT*SIGMA0/SigmaInf[i]/ri;
	      if ( ViscosityAlpha ) {
	        vr[l] -= 3.0*viscosity/r*(-SIGMASLOPE+2.0*FLARINGINDEX+1.0);
	      } else {
	        vr[l] -= 3.0*viscosity/r*(-SIGMASLOPE+.5);
	      }
      }
    }
  }
  for (j = 0; j < ns; j++) {
    vr[j] = vr[j+ns*nr] = 0.0;
  }
  // Commented out for generic simualtions, to match Crida2006 (i think) runs uncomment
  // if ( RadiativeOnly ) {
  //   for (i = 0; i < nr; i++) {
  //     for (j = 0; j < ns; j++) {
  //       l = j+i*ns;
  //       vr[l] = 0;
  //       vt[l] = pow(Rmed[i], -0.5);
  //     }
  //   }
  // }
}


