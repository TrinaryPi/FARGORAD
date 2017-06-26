/* C Header
	* @filename        : binary.c
	* @author          : Matthew Mutter (trinarypi)
	* @last_modified_by: trinarypi
	* @last_modified   : 2017/06/26 16:58
	* @description     :
*/
#include "mp.h"

static real Xbinary, Ybinary, VXbinary, VYbinary, Mbinary;
extern boolean OpenInner, NonReflecting, OuterSourceMass, Evanescent, ClosedInner, MassFlow;
extern boolean SelfGravity, SGZeroMode, Adiabatic, Indirect_Term;
extern Pair DiskOnPrimaryAcceleration;
extern Pair DiskOnBinaryAcceleration;
extern Pair bin_acc;
static Pair IndirectTerm;
extern int dimfxy;
extern real Abinary;

/* This file includes functions for initialising the binary system (if switched on) */


/* Allocates the memory and variables required for the binary system */
BinarySystem *AllocBinarySystem (nb)
int nb;
{
  real *mass, *x, *y, *vx, *vy, *smooth, *acc;
  boolean *feeldisk, *binfeelothers;
  int i;
  BinarySystem *bsys;
  bsys  = (BinarySystem *)malloc (sizeof(BinarySystem));
  if (bsys == NULL) {
    fprintf (stderr, "Not enough memory.\n");
    prs_exit (1);
  }
  x    = (real *)malloc (sizeof(real)*(nb+1));
  y    = (real *)malloc (sizeof(real)*(nb+1));
  vy   = (real *)malloc (sizeof(real)*(nb+1));
  vx   = (real *)malloc (sizeof(real)*(nb+1));
  mass = (real *)malloc (sizeof(real)*(nb+1));
  smooth = (real *)malloc (sizeof(real)*(nb+1));
  acc = (real *)malloc (sizeof(real)*(nb+1));
  if ((x == NULL) || (y == NULL) || (vx == NULL) || (vy == NULL)  || (mass == NULL) || (acc == NULL)) {
    fprintf (stderr, "Not enough memory.\n");
    prs_exit (1);
  }
  feeldisk = (boolean *)malloc (sizeof(real)*(nb+1));
  binfeelothers = (boolean *)malloc (sizeof(real)*(nb+1));
  if ((feeldisk == NULL) || (binfeelothers == NULL)) {
    fprintf (stderr, "Not enough memory.\n");
    prs_exit (1);
  }
  bsys->x = x;
  bsys->y = y;
  bsys->vx = vx;
  bsys->vy = vy;
  bsys->mass = mass;
  bsys->smooth = smooth;
  bsys->acc = acc;
  bsys->FeelDisk = feeldisk;
  bsys->FeelOthers = binfeelothers;
  for (i = 0; i < nb; i++) {
    x[i] = y[i] = vx[i] = vy[i] = mass[i] = smooth[i]= acc[i] = 0.0;
    feeldisk[i] = binfeelothers[i] = YES;
  }
  return bsys;
}


/*Initialises the variables of the binary system stars, with the input parameters, ready for simulation */
BinarySystem *InitBinarySystem (filename)
char *filename;
{
  FILE *input;
  char s[512], nm[512], test1[512], test2[512], *s1;
  BinarySystem *bsys;
  int i = 0, nb;
  float mass, abin, ecc, per, rsmooth, accret;
  boolean feeldis, binfeelothers;
  nb = FindNumberOfPlanets (filename);
  if (CPU_Master)
    printf ("%d star(s) found.\n", nb);
  bsys = AllocBinarySystem (nb);
  input = fopen (filename, "r");
  bsys->nb = nb;
  while (fgets(s, 510, input) != NULL) {
    sscanf(s, "%s ", nm);
    if (isalpha(s[0])) {
      s1 = s + strlen(nm);
      sscanf(s1 + strspn(s1, "\t :=>_"), "%f %f %f %f %f %f %s %s", &abin, &mass, &ecc, &per, &rsmooth, &accret, test1, test2);
      Abinary = abin;
      bsys->mass[i] = (real)mass;
      bsys->smooth[i] = (real)rsmooth;
      feeldis = binfeelothers = YES;
      if (tolower(*test1) == 'n') feeldis = NO;
      /*
      if ( (feeldis == YES) && (ForcedCircular == YES) ) {
	masterprint ("Careful: there is a contradiction between FeelDisk = Yes in your planet configuration file, and ForcedCircular = Yes in your parameter file. I decided to put FeelDisk = No. Please check this is what you really want to mean.");
	feeldis = NO;
      }
      */
      if (tolower(*test2) == 'n') binfeelothers = NO;
      bsys->x[i] = (real)(1.0-mass)*abin*cos(per)*(1.0-ecc);
      bsys->y[i] = 0.0;
      bsys->vy[i] = (real)pow(abin,-3.0/2.0)*cos(per)*abin*(1.0-mass)*(1.0+ecc)/sqrt(1.0-ecc*ecc);
      bsys->vx[i] = 0.0;
      bsys->acc[i] = accret;
      bsys->FeelDisk[i] = feeldis;
      bsys->FeelOthers[i] = binfeelothers;
      i++;
    }
  }
  return bsys;
}


/* Clears the binary system parameters and arrays */
void FreeBinary (bsys)
     BinarySystem *bsys;
{
  free (bsys->x);
  free (bsys->vx);
  free (bsys->y);
  free (bsys->vy);
  free (bsys->mass);
  free (bsys->acc);
  free (bsys->FeelOthers);
  free (bsys->FeelDisk);
  free (bsys);
}


boolean CountStars (sys)
	BinarySystem *sys;
{
	int nb;
	boolean bool=YES;
	nb = sys->nb;
	if (nb < 2) {
		printf ("Tried to find a binary system, two stars not found\n");
		printf ("Reverting to a single star system...\n");
		bool = NO;
	}
	return bool;
}


/* Prints the stellar properties at start of runtime */
void ListStars (bsys)
     BinarySystem *bsys;
{
  int nb;
  int i;
  nb = bsys->nb;
  if (!CPU_Master) return;
  for (i = 0; i < nb; i++) {
    printf ("Star number %d\n", i);
    printf ("---------------\n");
    printf ("x = %.10f\ty = %.10f\n", bsys->x[i], bsys->y[i]);
    printf ("vx = %.10f\tvy = %.10f\n", bsys->vx[i], bsys->vy[i]);
    if (bsys->acc[i] == 0)
      printf ("Non-accreting. \n");
    else
      printf ("Accretion Time = %.10f\n", 1.0/(bsys->acc[i]));
    if (bsys->FeelDisk[i] == YES) {
      printf ("Feels the disk potential\n");
    } else {
      printf ("Doesn't feel the disk potential\n");
    }
    if (bsys->FeelOthers[i] == YES) {
      printf ("Feels the other planets potential\n");
    } else {
      printf ("Doesn't feel the other planets potential\n");
    }
    if (Indirect_Term == YES) {
      printf ("Removing Indirect Term from system\n");
    } else {
      printf ("No Indirect Term calculated\n");
    }
    printf ("\n");
  }
}


void WriteBinaryFile (TimeStep, n)
     int TimeStep;
     int n;
{
  FILE *output;
  char name[256];
  if (!CPU_Master) return;
  printf ("Updating 'binary%d.dat'...", n);
  fflush (stdout);
  sprintf (name, "%sbinary%d.dat", OUTPUTDIR, n);
  output = fopen (name, "a");
  if (output == NULL) {
    fprintf (stderr, "Can't write 'binary%d.dat' file. Aborting.\n", n);
    prs_exit (1);
  }
  fprintf (output, "%d\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\n", TimeStep, Xbinary, Ybinary, VXbinary, VYbinary, Mbinary);
  fclose (output);
  printf ("done\n");
  fflush (stdout);
}

void WriteBinaryAccretionFile (TimeStep)
     int TimeStep;
{
  FILE *output;
  char OutputName[256];
  if (!CPU_Master) return;
  fflush(stdout);
  sprintf (OutputName, "%saccretion.dat", OUTPUTDIR);
  output = fopen (OutputName, "a");
  if (output == NULL) {
    fprintf (stderr, "Can't write 'accretion.dat' file. Just thought you should know...");
  }
  fprintf (output, "%d\t%#.12g\t%#.12g\n", TimeStep, bin_acc.x, bin_acc.y );
  fclose (output);
  fflush (stdout);
}


void WriteBinarySystemFile (bsys, t)
     BinarySystem *bsys;
     int t;
{
  int i, n;
  n = bsys->nb;
  for (i = 0; i < n; i++) {
    Xbinary = bsys->x[i];
    Ybinary = bsys->y[i];
    VXbinary = bsys->vx[i];
    VYbinary = bsys->vy[i];
    Mbinary = bsys->mass[i];
    WriteBinaryFile (t, i);
  }
}


real GetFromBinaryFile (TimeStep, column, n)
     int TimeStep, column, n;
{
  FILE *input;
  char name[256];
  char testline[256];
  int time;
  char *pt;
  double value;
  sprintf (name, "%sbinary%d.dat", OUTPUTDIR, n);
  input = fopen (name, "r");
  if (input == NULL) {
    mastererr ("Can't read 'binary%d.dat' file. Aborting restart.\n",n);
    prs_exit (1);
  }
  if (column < 2) {
    mastererr ("Invalid column number in 'binary%d.dat'. Aborting restart.\n",n);
    prs_exit (1);
  }
  do {
    pt = fgets (testline, 255, input);
    sscanf (testline, "%d", &time);
  } while ((time != TimeStep) && (pt != NULL));
  if (pt == NULL) {
    mastererr ("Can't read entry %d in 'binary%d.dat' file. Aborting restart.\n", TimeStep,n);
    prs_exit (1);
  }
  fclose (input);
  pt = testline;
  while (column > 1) {
    pt += strspn(pt, "eE0123456789-.");
    pt += strspn(pt, "\t :=>_");
    column--;
  }
  sscanf (pt, "%lf", &value);
  return (real)value;
}

real GetFromBinaryAccretionFile (TimeStep, column)
     int TimeStep, column;
{
  FILE *input;
  char InputName[256];
  char testline[256];
  int time;
  char *pt;
  double value;
  sprintf (InputName, "%saccretion.dat", OUTPUTDIR);
  input = fopen (InputName, "r");
  if (input == NULL) {
    mastererr ("Can't read 'accretion.dat' file. Aborting restart.\n");
    prs_exit (1);
  }
  if (column < 2) {
    mastererr ("Invalid column number in 'accretion.dat'. Aborting restart.\n");
    prs_exit (1);
  }
  do {
    pt = fgets (testline, 255, input);
    sscanf (testline, "%d", &time);
  } while ((time != TimeStep) && (pt != NULL));
  if (pt == NULL) {
    mastererr ("Can't read entry %d in 'accretion.dat' file. Aborting restart.\n", TimeStep);
    prs_exit (1);
  }
  fclose (input);
  pt = testline;
  while (column > 1) {
    pt += strspn(pt, "eE0123456789-.");
    pt += strspn(pt, "\t :=>_");
    column--;
  }
  sscanf (pt, "%lf", &value);
  return (real)value;
}


void RestartBinarySystem (timestep, bsys)
     BinarySystem *bsys;
     int timestep;
{
  int k;
  for (k = 0; k < bsys->nb; k++) {
    bsys->x[k] = GetFromBinaryFile (timestep, 2, k);
    bsys->y[k] = GetFromBinaryFile (timestep, 3, k);
    bsys->vx[k] = GetFromBinaryFile (timestep, 4, k);
    bsys->vy[k] = GetFromBinaryFile (timestep, 5, k);
    bsys->mass[k] = GetFromBinaryFile (timestep, 6, k);
  }
}

Pair ComputeIndirect (force, Rho, mass, bsys)
     Force *force;
     PolarGrid *Rho;
     real mass;
     BinarySystem *bsys; 
{
  int k;
  real ms, M, xs, ys, rs, smoothing;
  Pair acceleration;
    M = bsys->mass[0] + bsys->mass[1];
    acceleration.x = 0.0;
    acceleration.y = 0.0;
    for (k = 0; k < bsys->nb; k++ ) {
      xs = bsys->x[k];
      ys = bsys->y[k];
      rs = sqrt(xs*xs + ys*ys);
      ms = bsys->mass[k];
      if (RocheSmoothing)
        smoothing = Abinary*pow(ms/(3.0*M), 1.0/3.0)*(bsys->smooth[k]);
      else
        smoothing = (bsys->smooth[k]) * AspectRatio(rs) * pow(rs, 1.0+FLARINGINDEX);
      if (bsys->FeelDisk[k] == YES) {
  		  ComputeForce (force, Rho, xs, ys, smoothing, mass, dimfxy);
  		  if (ExcludeHill) {
    		  acceleration.x += ms*(force->fx_ex_inner) + ms*(force->fx_ex_outer);
    		  acceleration.y += ms*(force->fy_ex_inner) + ms*(force->fy_ex_outer);
        
  		  } else {
    		  acceleration.x += ms*(force->fx_inner) + ms*(force->fx_outer);
    		  acceleration.y += ms*(force->fy_inner) + ms*(force->fy_outer);
    	  }
	    }
 	  }	
  return acceleration;
}

void AdvanceBinaryFromDisk (force, Rho, Energy, bsys, dt)
     Force *force;
     BinarySystem *bsys;
     PolarGrid *Rho, *Energy;
     real dt;		       
{
  int nstar, k;
  Pair gamma;
  real xs, ys, rs, ms, M, smoothing;
  nstar = bsys->nb;

  M = bsys->mass[0] + bsys->mass[1];
  for (k = 0; k < nstar; k++) {
    if (bsys->FeelDisk[k] == YES) {
      ms = bsys->mass[k];
      xs = bsys->x[k];
      ys = bsys->y[k];
      rs = sqrt(xs*xs + ys*ys);
      if (RocheSmoothing)
        smoothing = Abinary*pow(ms/(3.0*M), 1.0/3.0)*bsys->smooth[k];
      else
        smoothing = (bsys->smooth[k]) * AspectRatio(rs) * pow(rs, 1.0+FLARINGINDEX);
      gamma = ComputeAccel (force, Rho, xs, ys, smoothing, ms);
      bsys->vx[k] += dt * gamma.x;
      bsys->vy[k] += dt * gamma.y;
      bsys->vx[k] += dt * IndirectTerm.x;
      bsys->vy[k] += dt * IndirectTerm.y;
    }
  }
}

void SolveBinOrbits (bsys)
     BinarySystem *bsys;
{ 
  real x0, y0, vx0, vy0, x1, y1, vx1, vy1, m0, m1;
  real vx, vy, x, y, v2, r;
  real a, e, e_temp, V, mu, omega;
  real ex, ey;

  x0  = bsys->x[0];
  y0  = bsys->y[0];
  vx0 = bsys->vx[0];
  vy0 = bsys->vy[0];
  x1  = bsys->x[1];
  y1  = bsys->y[1];
  vx1 = bsys->vx[1];
  vy1 = bsys->vy[1];

  x  = x0-x1;
  y  = y0-y1;
  r  = sqrt(x*x + y*y);
  vx = vx0-vx1;
  vy = vy0-vy1;
  v2 = vx*vx + vy*vy;

  m0 = bsys->mass[0];
  m1 = bsys->mass[1];
  mu = G*(m0 + m1);

  a = 1.0/((2.0/r)-(v2/mu));
  Abinary = a;

  ex = vy*vy*x/mu - vy*vx*y/mu - x/r;
  ey = vx*vx*y/mu - vx*vy*x/mu - y/r;
  e_temp = ex*ex + ey*ey;
  if (e_temp < 0.0) e_temp = 0.0;
  e = sqrt(e_temp);
  if (e < 1.0E-6) e = 0.0;
  
  if (e != 0.0) {
    V = acos((a*(1.0-e*e)/r-1.0)/e);
  } else {
    V = 0.0;
  }

  if (e != 0.0) {
    omega=atan2(ey,ex);
  } else {
    omega=atan2(y,x);
  }

  if (abs(omega) < 1.0E-6) omega = 0.0;
  
  FILE *output;
  char name[256];
  if (CPU_Rank != CPU_Highest) return;
  sprintf (name, "%sbinorbit.dat", OUTPUTDIR);
  output = fopen (name, "a");
  if (output == NULL) {
    message ("Can't open 'binorbit.dat'. Exited.\n");
    prs_exit (1);
  }
  fprintf (output, "%.16g\t%.16g\t%.16g\t%.16g\t%.16g\n", PhysicalTime, e, a, V, omega);
  fclose (output);
}

void CalcAbin (bsys) 
    BinarySystem *bsys;
{
  real x0, y0, vx0, vy0, x1, y1, vx1, vy1, m0, m1, mu;
  real vx, vy, x, y, v2, r;

  x0 = bsys->x[0];
  y0 = bsys->y[0];
  vx0 = bsys->vx[0];
  vy0 = bsys->vy[0];
  x1 = bsys->x[1];
  y1 = bsys->y[1];
  vx1 = bsys->vx[1];
  vy1 = bsys->vy[1];
  m0 = bsys->mass[0];
  m1 = bsys->mass[1];

  x = x0-x1;
  y = y0-y1;
  vx = vx0-vx1;
  vy = vy0-vy1;
  r = sqrt(x*x + y*y);
  v2 = vx*vx + vy*vy;
  mu = G*(m0 + m1);
  Abinary = 1.0/((2.0/r)-(v2/mu));
}

void AccreteOntoStars (Rho, dt, bsys)
     real dt;
     PolarGrid *Rho;
     BinarySystem *bsys;
{
  real RRoche, Rstar, distance, dx, dy, deltaM, angle, temp, binsep;
  int i_min,i_max, j_min, j_max, i, j, l, jf, ns, nr, k;
  real Xstar, Ystar, Mstar, q;
  real facc, facc1, facc2, frac1, frac2; /* We adopt the same notations as W. Kley */
  real *dens, *abs, *ord, binacc[2];
  real xc, yc;
  real dMstar;
  nr     = Rho->Nrad;
  ns     = Rho->Nsec;
  dens   = Rho->Field;
  abs    = CellAbscissa->Field;
  ord    = CellOrdinate->Field;
  for (k=0; k < bsys->nb; k++) {
    if (bsys->acc[k] > 1e-10) {
      dMstar = binacc[k] = 0.0;
      /* Hereafter : initialization of W. Kley's parameters */
      facc = dt*(bsys->acc[k]);
      facc1 = 1.0/3.0*facc;
      facc2 = 2.0/3.0*facc;
      frac1 = 0.75;
      frac2 = 0.45;
      /* W. Kley's parameters initialization finished */
      Xstar = bsys->x[k];
      Ystar = bsys->y[k];
      binsep = sqrt(((bsys->x[0])-(bsys->x[1]))*((bsys->x[0])-(bsys->x[1]))
        +((bsys->y[0])-(bsys->y[1]))*((bsys->y[0])-(bsys->y[1])));
      Mstar = bsys->mass[k];
      Rstar = sqrt(Xstar*Xstar+Ystar*Ystar);
      /* Think about what value to ACTUALLY have here! Roche Lobe or Hill Sphere */
      /* This is more applicable to the hill sphere
      RRoche = pow((1.0/3.0*Mplanet),(1.0/3.0))*Rplanet; */
      q = Mstar/(1-Mstar);
      RRoche = binsep*(0.38 + 0.2*log(q));
      /* Central mass is 1.0 */
      i_min=0;
      i_max=nr-1;
      while ((Rsup[i_min] < Rstar-RRoche) && (i_min < nr)) i_min++;
      while ((Rinf[i_max] > Rstar+RRoche) && (i_max > 0)) i_max--;
      angle = atan2 (Ystar, Xstar);
      j_min =(int)((real)ns/2.0/PI*(angle - 2.0*RRoche/Rstar));
      j_max =(int)((real)ns/2.0/PI*(angle + 2.0*RRoche/Rstar));

#pragma omp parallel for private(j,jf,l,xc,yc,dx,dy,distance,deltaM) shared(dMstar)
      for (i = i_min; i <= i_max; i++) {
        for (j = j_min; j <= j_max; j++) {
          jf = j;
          while (jf <  0)  jf += ns;
          while (jf >= ns) jf -= ns;
          l   = jf+i*ns;
          xc = abs[l];
          yc = ord[l];
          dx = Xstar-xc;
          dy = Ystar-yc;
          distance = sqrt(dx*dx+dy*dy);
            if (distance < frac1*RRoche) {
            deltaM = facc1*dens[l]*Surf[i];
            if (i < Zero_or_active) deltaM = 0.0;
            if (i >= Max_or_active) deltaM = 0.0;
            dens[l] *= (1.0 - facc1);
#pragma omp atomic
            dMstar     += deltaM;
          }
          if (distance < frac2*RRoche) {
            deltaM = facc2*dens[l]*Surf[i];
            if (i < Zero_or_active) deltaM = 0.0;
            if (i >= Max_or_active) deltaM = 0.0;
            dens[l] *= (1.0 - facc2);
#pragma omp atomic
            dMstar     += deltaM;
          }
        }
      }
      MPI_Allreduce (&dMstar, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      dMstar = temp;
      
    } else {
      dMstar = 0.0;

    }
    binacc[k] = dMstar;
  }
  bin_acc.x += binacc[0];
  bin_acc.y += binacc[1];
}
