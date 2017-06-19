#include "mp.h"
#define RKACCURACY 1e-1

static real k1[MAX1D], k2[MAX1D], k3[MAX1D], k4[MAX1D], k5[MAX1D], k6[MAX1D];
extern boolean Indirect_Term;

void DerivMotionRK5 (q_init, masses, deriv, np, ns, dt, feelothers, binfeelothers)
     real *q_init, *deriv, *masses, dt;
     boolean *feelothers, *binfeelothers;
     int np, ns;
{
  real *x,*y,*vx,*vy, dist, acomx, acomy;
  real *derivx, *derivy, *derivvx, *derivvy;
  int i, j, n, ip, jp;
  n = np + ns;
  x = q_init;
  y = x+n;
  vx = y+n;
  vy = vx+n;
  derivx = deriv;
  derivy = derivx+n;
  derivvx = derivy+n;
  derivvy = derivvx+n;
  //printf ("Calling DerivMotionRK5\n");
  acomx = 0.0;
  acomy = 0.0;

  /* Evolution of Stars */

  if (BinaryOn == YES) {
    for (i = 0; i < ns; i++) {
      for (j = 0; j < np; j++) {
        if ((Indirect_Term == YES) && (binfeelothers[i] == YES)) {
          jp = j+ns;
          dist = (x[i]-x[jp])*(x[i]-x[jp])+(y[i]-y[jp])*(y[i]-y[jp]);
          dist = sqrt(dist);
          acomx -= G*masses[i]*masses[jp]*(x[jp]-x[i])/dist/dist/dist;    // Indirect term from planets on the binary
          acomy -= G*masses[i]*masses[jp]*(y[jp]-y[i])/dist/dist/dist;
        }

      }

    }
    for (i = 0; i < ns; i++) {
      derivx[i] = vx[i];
      derivy[i] = vy[i];
      derivvx[i] = acomx;
      derivvy[i] = acomy;
      for (j = 0; j < ns; j++) {
        if (j != i) {
          dist = (x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j]);
          dist = sqrt(dist);
          derivvx[i] -= G*masses[j]*(x[i]-x[j])/dist/dist/dist;           // Direct term from other star
          derivvy[i] -= G*masses[j]*(y[i]-y[j])/dist/dist/dist;
        }

      }
      for (j = 0; j < np; j++){
        if (binfeelothers[i] == YES){
          jp = j+ns;
          dist = (x[i]-x[jp])*(x[i]-x[jp])+(y[i]-y[jp])*(y[i]-y[jp]);
          dist = sqrt(dist);
          derivvx[i] -= G*masses[jp]*(x[i]-x[jp])/dist/dist/dist;         // Direct term from planets
          derivvy[i] -= G*masses[jp]*(y[i]-y[jp])/dist/dist/dist;
        }

      }

    }

    /* Evolution of the planets */

    for (i = 0; i < np; i++){
      ip = i+ns;
      derivx[ip] = vx[ip];
      derivy[ip] = vy[ip];
      derivvx[ip] = acomx;
      derivvy[ip] = acomy;
      for (j = 0; j < np; j++) {
        jp = j+ns;
        if ((jp != ip) && (feelothers[i] == YES)) {
          dist = (x[ip]-x[jp])*(x[ip]-x[jp])+(y[ip]-y[jp])*(y[ip]-y[jp]);
          dist = sqrt(dist);
          derivvx[ip] -= G*masses[jp]*(x[ip]-x[jp])/dist/dist/dist;       //Direct term from other planets
          derivvy[ip] -= G*masses[jp]*(y[ip]-y[jp])/dist/dist/dist;
        }
      }
      for (j = 0; j < ns; j++) {
        dist = (x[ip]-x[j])*(x[ip]-x[j])+(y[ip]-y[j])*(y[ip]-y[j]);
        dist = sqrt(dist);
        derivvx[ip] -= G*masses[j]*(x[ip]-x[j])/dist/dist/dist;           //Direct term from stars
        derivvy[ip] -= G*masses[j]*(y[ip]-y[j])/dist/dist/dist;
      }

    }

  } else {
    for (i = 0; i < np; i++) {
      derivx[i] = vx[i];
      derivy[i] = vy[i];
      dist = (x[i]*x[i])+(y[i]*y[i]);
      dist = sqrt(dist);
      derivvx[i] = -G*1.0*x[i]/dist/dist/dist;
      derivvy[i] = -G*1.0*y[i]/dist/dist/dist;                            //Direct term from star
      for (j = 0; j < np; j++) {
        if (Indirect_Term == YES) {
          dist = (x[j]*x[j])+(y[j]*y[j]);
          dist = sqrt(dist);
	        derivvx[i] -= G*1.0*masses[j]*x[j]/dist/dist/dist;              //Indirect term from planet on primary star
          derivvy[i] -= G*1.0*masses[j]*y[j]/dist/dist/dist;
        }
        if ((j != i) && (feelothers[i] == YES)) {
	       dist = (x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j]);
	       dist = sqrt(dist);
	       derivvx[i] -= G*masses[j]*(x[i]-x[j])/dist/dist/dist;            //Direct term from other planets
	       derivvy[i] -= G*masses[j]*(y[i]-y[j])/dist/dist/dist;
        }

      }

    }

  }
  for (i = 0; i < 4*n; i++)
      deriv[i] *= dt;
}

void TranslatePlanetRK5 (qold, c1, c2, c3, c4, c5, qnew, n)
     real *qold, *qnew;
     real c1, c2, c3, c4, c5;
     int n;
{
  int i;
  //printf ("Calling TranslatePlanetRK5\n");
  for (i = 0; i < 4*n; i++)
    qnew[i] = qold[i]+c1*k1[i]+c2*k2[i]+c3*k3[i]+c4*k4[i]+c5*k5[i];
}

void RungeKunta (q0, dt, masses, q1, np, ns, feelothers, binfeelothers)
     real *q0, *q1;
     real dt, *masses;
     boolean *feelothers, *binfeelothers;
     int np, ns;
{
  int i, n;
  n = ns + np;
  real timestep;
  timestep = dt;
  //printf ("Calling RungeKunta\n");
  for (i = 0; i < n*4; i++) {
    k1[i] = k2[i] = k3[i] = k4[i] = k5[i] = k6[i];
  }
  DerivMotionRK5 (q0, masses, k1, np, ns, timestep, feelothers, binfeelothers);
  TranslatePlanetRK5 (q0, 0.2, 0.0, 0.0, 0.0, 0.0, q1, n);
  DerivMotionRK5 (q1, masses, k2, np, ns, timestep, feelothers, binfeelothers);
  TranslatePlanetRK5 (q0, 0.075, 0.225, 0.0, 0.0, 0.0, q1, n);
  DerivMotionRK5 (q1, masses, k3, np, ns, timestep, feelothers, binfeelothers);
  TranslatePlanetRK5 (q0, 0.3, -0.9, 1.2, 0.0, 0.0, q1, n);
  DerivMotionRK5 (q1, masses, k4, np, ns, timestep, feelothers, binfeelothers);
  TranslatePlanetRK5 (q0, -11.0/54.0, 2.5, -70.0/27.0, 35.0/27.0, 0.0, q1, n);
  DerivMotionRK5 (q1, masses, k5, np, ns, timestep, feelothers, binfeelothers);
  TranslatePlanetRK5 (q0, 1631.0/55296.0, 175.0/512.0, 575.0/13824.0, 44275.0/110592.0, 253.0/4096.0, q1, n);
  DerivMotionRK5 (q1, masses, k6, np, ns, timestep, feelothers, binfeelothers);
  for (i = 0; i < 4*n; i++) {
    q1[i]=q0[i]+37.0/378.0*k1[i]+250.0/621.0*k3[i]+125.0/594.0*k4[i]+512.0/1771.0*k6[i];
  }
}
 
