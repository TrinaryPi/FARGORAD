/* C Header
	* @filename        : types.h
	* @author          : Matthew Mutter (trinarypi)
	* @last_modified_by: trinarypi
	* @last_modified   : 2017/06/26 16:59
	* @description     :
*/
#include <sys/times.h>

typedef int     boolean;
typedef double	real;

struct torque_set {
  real          InnerDisk;
  real          OuterDisk;
  real		ExcludeHill;
  real          Total;
};

typedef struct torque_set TorqueSet;

struct triplet {
  real            x;
  real            y;
  real		z;
};

typedef struct triplet Triplet;

struct pair {
  real            x;
  real            y;
};

typedef struct pair Pair;

struct force {
  real fx_inner;
  real fy_inner;
  real fx_ex_inner;
  real fy_ex_inner;
  real fx_outer;
  real fy_outer;
  real fx_ex_outer;
  real fy_ex_outer;
  real *GlobalForce;
};

typedef struct force Force;

struct polargrid {
  int             Nrad;
  int             Nsec;
  real           *Field;
  char           *Name;
};

typedef struct polargrid PolarGrid;

#define		YES	1
#define		NO	0
#define		REAL	1
#define		INT	0
#define		STRING  2
#define 	SINE	0
#define		COSINE	1
#define		ABSCISSA	0
#define		ORDINATE	1
#define		HEIGHT		2
#define		INF		0
#define 	SUP		1
#define         GET             0
#define         MARK            1
#define         FREQUENCY       2
#define         COM_DENSITY     0
#define         COM_VRAD        1
#define         COM_VTHETA      2

#define		MAX1D	16384

struct param {
  char name[80];
  int  type;
  char *variable;
  int read;
  int necessary;
};

typedef struct param Param;

struct timeprocess {
  char name[80];
  clock_t clicks;
};

typedef struct timeprocess TimeProcess;

struct planetary_system {
  int nb;			/* Number of planets */
  real *mass;			/* their masses */
  real *x;			/* their coordinates */
  real *y;
  real *vx;			/* their velocities */
  real *vy;
  real *acc;			/* Their accretion times^-1 */
  char **name;			/* their names */
  boolean *FeelDisk;		/* will "migrate" ? */
  boolean *FeelOthers;		/* will feel other planets ? */
};

typedef struct planetary_system PlanetarySystem;

/* Structure Needed for binary star system */
struct binary_system {
  int nb;     /* Number of stars */
  real *mass;     /* their masses */
  real *x;      /* their coordinates */
  real *y;
  real *vx;     /* their velocities */
  real *vy;
  real *smooth;
  char **name;      /* their names */
  boolean *FeelDisk;    /* will "migrate" ? */
  boolean *FeelOthers;    /* will feel other planets ? */
  real *acc;
};

typedef struct binary_system BinarySystem;

struct rad_kappa {
  int Ni;
  real *Kappa0;
  real *Alpha;
  real *Beta;
  real *Tmax;
  real *Constant;
  real *Exponent;
};

typedef struct rad_kappa OpacTable;

struct stellar_irr {
  int nb;
  real *Rstar;
  real *Tstar;
  boolean *CentralSource;
};

typedef struct stellar_irr IrradiationSource;

// struct intersections {
//   real *x;
//   real *y;
//   int n;
// };

// typedef struct intersections ints;


// struct ray_struct {
//   real x;
//   real y;
//   real dx;
//   real dy;
//   real dr;
//   real r;
//   real th;
//   real diff;
//   real length;
//   real dtau;
//   real tau;
//   int env;
//   ints *Ints;
// };

struct ray_struct {
  real x;
  real y;
  real dx;
  real dy;
  real dr;
  real r;
  real th;
  real diff;
  real length;
  real separation;
  real dtau;
  real tau;
  real x_int[2];
  real y_int[2];
  int env;
  int n_int;
};

typedef struct ray_struct RayStruct;
 
