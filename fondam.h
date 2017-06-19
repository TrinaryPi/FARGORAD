#define		G	1.0
#define	      PI4	12.56637061435917295376
#define	       PI	3.14159265358979323844
#define   CPUOVERLAP    5	/* Zeus-like overlap kernel. 2:transport; 2: source, 1:viscous stress  */
#define      MU         1.0  /* Mean molecular weight */
#define      R          1.0  /* Universal Gas Constant in code units */
//  	STEFANK		15621.414985793097


// Unit conversion factors for code units to cgs units //

#define 	LCGS  		1.4959E13
#define 	MCGS  		1.9891E33*MSTARSOLAR
#define 	GCGS  		6.6738E-8
#define 	TCGS  		pow(GCGS*MCGS/LCGS/LCGS/LCGS, -0.5)
#define 	VCGS  		LCGS/TCGS
#define 	ECGS  		MCGS*pow(LCGS/TCGS, 2.0)
#define 	E2CGS 		ECGS/LCGS/LCGS
#define 	R_MUCGS 	3.522881356E7
#define 	TEMPCGS 	E2CGS/MCGS*LCGS*LCGS/R_MUCGS
#define  	STEFANKCGS	5.67037321E-5
#define		STEFANK		STEFANKCGS/(pow(R_MUCGS, 4.0)*pow(GCGS, -2.5)* pow(MCGS, -1.5)* pow(LCGS, -0.5))
#define		CCGS		2.998E10
#define		Ccode 			CCGS*TCGS/LCGS
#define		ACGS		(4.0*STEFANKCGS)/CCGS
#define		Acode			ACGS/(ECGS*pow(LCGS, -3.0)*pow(TEMPCGS, -4.0))

// STEFANKCGS = STEFANK/(pow(R_MUCGS,-4)*pow(GCGS,2.5)* pow(MCGS,1.5)* pow(LCGS,0.5))
