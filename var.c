#define __LOCAL
#include "mp.h"
#undef __LOCAL

void
InitVariables()
{
  var("DT", &DT, REAL, YES, "1.");
  var("SIGMA0", &SIGMA0, REAL, YES, "173.");
  var("NINTERM", &NINTERM, INT, YES, "10.");
  var("NTOT", &NTOT, INT, YES, "1501.");
  var("OUTPUTDIR", OUTPUTDIR, STRING, YES, "~masset");
  var("INNERBOUNDARY", OPENINNERBOUNDARY, STRING, NO, "WALL");
  var("OUTERBOUNDARY", OPENOUTERBOUNDARY, STRING, NO, "WALL");
  var("LABELADVECTION", ADVLABEL, STRING, NO, "NO");
  var("TRANSPORT", TRANSPORT, STRING, NO, "FAST");
  var("PLANETCONFIG", PLANETCONFIG, STRING, NO, "Systems/SolarSystem.cfg");
  var("STARCONFIG", STARCONFIG, STRING, NO, "in/binary0.cfg");
  var("MASSTAPER", &MASSTAPER, REAL, NO, "0.0000001");
  var("MMAX", &MMAX, REAL, NO, "0.001");
  var("CONFIGMP", CONFIGMP, STRING, NO, "NO");
  var("CONFIGPOS", CONFIGPOS, STRING, NO, "NO");
  var("SIGMATAPER", SIGMATAPER, STRING, NO, "NO");
  var("SIGMACAVITY", SIGMACAVITY, STRING, NO, "NO");
  var("RADIALSPACING", GRIDSPACING, STRING, NO, "ARITHMETIC");
  var("NRAD", &NRAD, INT, YES, "64.0");
  var("NSEC", &NSEC, INT, YES, "64.0");
  var("RMIN", &RMIN, REAL, YES, "1.0");
  var("RMAX", &RMAX, REAL, YES, "1.0");
  var("THICKNESSSMOOTHING", &THICKNESSSMOOTHING, REAL, NO, "0.0");
  var("ROCHESMOOTHING", &ROCHESMOOTHING, REAL, NO, "0.0");
  var("ASPECTRATIO", &ASPECTRATIO, REAL, YES, "0.05");
  var("VISCOSITY", &VISCOSITY, REAL, NO, "0.0");
  var("ALPHAVISCOSITY", &ALPHAVISCOSITY, REAL, NO, "0.0");  
  var("SIGMASLOPE", &SIGMASLOPE, REAL, YES, "0.0");  
  var("RELEASERADIUS", &RELEASERADIUS, REAL, NO, "0.0");  
  var("RELEASEDATE", &RELEASEDATE, REAL, NO, "0.0");  
  var("OMEGAFRAME", &OMEGAFRAME, REAL, NO, "0.0");
  var("DISK", DISK, STRING, NO, "YES");
  //new - to turn all HD effects and functions off
  var("HYDRO", HYDRO, STRING, NO, "YES");
  //new - set purely Keplerian, but fully-Radiative disc - no transport etc.
  var("RADIATIVEONLY", RADIATIVEONLY, STRING, NO, "NO");
  var("FRAME", FRAME, STRING, NO, "FIXED");
  var("OUTERSOURCEMASS", OUTERSOURCEMASS, STRING, NO, "NO");
  var("WRITEDENSITY", WRITEDENSITY, STRING, NO, "YES");
  var("WRITEVELOCITY", WRITEVELOCITY, STRING, NO, "YES");
  var("WRITEENERGY", WRITEENERGY, STRING, NO, "NO");
  var("WRITETEMPERATURE", WRITETEMPERATURE, STRING, NO, "NO");
  var("WRITEDIVV", WRITEDIVV, STRING, NO, "NO");
  var("WRITEQPLUS", WRITEQPLUS, STRING, NO, "NO");
  var("WRITEECCENTRCITY", WRITEECCENTRICITY, STRING, NO, "NO");
  var("WRITEPERICENTRE", WRITEPERICENTRE, STRING, NO, "NO");
  var("INDIRECTTERM", INDIRECTTERM, STRING, NO, "YES");
  var("EXCLUDEHILL", EXCLUDEHILL, STRING, NO, "NO");
  var("IMPOSEDDISKDRIFT", &IMPOSEDDISKDRIFT, REAL, NO, "0.0");
  var("FLARINGINDEX", &FLARINGINDEX, REAL, NO, "0.0");
  var("ECCENTRICITY", &ECCENTRICITY, REAL, NO, "0.0");
  var("CAVITYRADIUS", &CAVITYRADIUS, REAL, NO, "0.0");
  var("CAVITYRATIO", &CAVITYRATIO, REAL, NO, "1.0");
  var("CAVITYWIDTH", &CAVITYWIDTH, REAL, NO, "1.0");
  var("TRANSITIONRADIUS", &TRANSITIONRADIUS, REAL, NO, "0.0");
  var("TRANSITIONRATIO", &TRANSITIONRATIO, REAL, NO, "1.0");
  var("TRANSITIONWIDTH", &TRANSITIONWIDTH, REAL, NO, "1.0");
  var("LAMBDADOUBLING", &LAMBDADOUBLING, REAL, NO, "0.0");
  var("SELFGRAVITY", SELFGRAVITY, STRING, NO, "NO");
  var("CICPLANET", CICPLANET, STRING, NO, "NO");
  var("FORCEDCIRCULAR", FORCEDCIRCULAR, STRING, NO, "NO");
  var("ZMPLUS", ZMPLUS, STRING, NO, "NO");
  var("ADIABATIC", ADIABATIC, STRING, NO, "NO");
  var("ADIABATICINDEX", &ADIABATICINDEX, REAL, YES, "1.4");
  var("COOLING", COOLING, STRING, NO, "NO");
  var("COOLINGTIME0", &COOLINGTIME0, REAL, NO, "6.28");
  // new - to change the constant fraction cooling time in the disc to a stepped function
  var("CUSTOMCOOLING", CUSTOMCOOLING, STRING, NO, "NO");
  // new - these variables set the tapering of the disc mass with time
  var("DISCMASSTAPER", DISCMASSTAPER, STRING, NO, "NO"); 
  var("DMTTAU", &DMTTAU, REAL, NO, "0.0");
  var("DMTOUTPUTSTART", &DMTOUTPUTSTART, REAL, NO, "0.0");
  // Options required in input file for the new modules
  var("RADCOOLING", RADCOOLING, STRING, NO, "NO");
  var("IRRADIATION", IRRADIATION, STRING, NO, "NO");
  var("RADTRANSPORT", RADTRANSPORT, STRING, NO, "NO");
  var("RAYTRACINGHEATING", RAYTRACINGHEATING, STRING, NO, "NO");
  var("MSTARSOLAR", &MSTARSOLAR, REAL, NO, "1.0");
  var("WRITEDISCHEIGHT", WRITEDISCHEIGHT, STRING, NO, "NO");
  var("WRITEQMINUS", WRITEQMINUS, STRING, NO, "NO");
  var("WRITEKAPPA", WRITEKAPPA, STRING, NO, "NO");
  var("WRITECOEFFS", WRITECOEFFS, STRING, NO, "NO");
  var("WRITERFLD", WRITERFLD, STRING, NO, "NO");
  var("WRITEQIRR", WRITEQIRR, STRING, NO, "NO");
  var("WRITERTQIRR", WRITERTQIRR, STRING, NO, "NO");
  var("WRITETAUS", WRITETAUS, STRING, NO, "NO");
  var("NOCFL", NOCFL, STRING, NO, "NO");
  var("WRITERESIDUAL", WRITERESIDUAL, STRING, NO, "NO");
  var("IRRADIATIONCONFIG", IRRADIATIONCONFIG, STRING, NO, "in/irradiationsources.cfg");
  var("OPTIMALOMEGA", OPTIMALOMEGA, STRING, NO, "NO");
  var("SOROMEGA", &SOROMEGA, REAL, NO, "1.0");
  var("TOLERANCE", &TOLERANCE, REAL, NO, "1E-8");
  var("MAXITERATIONS", &MAXITERATIONS, INT, NO, "10000.0");
  var("INNERTEMPBC", INNERTEMPBC, STRING, NO, "Constant");
  var("OUTERTEMPBC", OUTERTEMPBC, STRING, NO, "Constant");
  var("TINNER", &TINNER, REAL, NO, "100.0");
  var("TOUTER", &TOUTER, REAL, NO, "100.0");
  var("TGRADINNER", &TGRADINNER, REAL, NO, "0.0");
  var("TGRADOUTER", &TGRADOUTER, REAL, NO, "0.0");
  var("REVERSETEMPERATUREINIT", REVERSETEMPERATUREINIT, STRING, NO,"NO");
  var("NORMRESIDUAL", NORMRESIDUAL, STRING, NO, "NO");
  var("NORMRELATIVE", NORMRELATIVE, STRING, NO, "NO");
  var("RADIATIONDEBUG", RADIATIONDEBUG, STRING, NO, "NO");
  var("OPACITYTABLE", OPACITYTABLE, STRING, NO, "Bell");
  var("OPACITYSMOOTHING", OPACITYSMOOTHING, STRING, NO, "NO");
  // new - Star Temperature Taper (6/01/2017)
  var("STARTAPER", &STARTAPER, REAL, NO, "0.0000001");
  // new - Constant factor between Rosseland opacity and visible opacity (i.e. Kley2014 etc) (6/01/2017)
  var("BITSCHSKAPPA", BITSCHSKAPPA, STRING, NO, "NO");
  // new - Explicit Raytracing heating step (outside of SOR step) (6/01/2017)
  var("RAYTRACINGHEATING", RAYTRACINGHEATING, STRING, NO, "NO"); // can be IMPLICIT, EXPLICIT, NO;
  // new - Variables for two source ray-tracing (04/04/2017)
  var("EMPTYCAVITY", EMPTYCAVITY, STRING, NO, "YES");
  var("TAUCEILING", &TAUCEILING, REAL, NO, "100000");
  var("RTPRECISION", &RTPRECISION, REAL, NO, "2.0");
  var("HUBENYSKAPPA", HUBENYSKAPPA, STRING, NO, "NO");

  var("STELLARSMOOTHING", STELLARSMOOTHING, STRING, NO, "YES");
  var("OUTFLOWBETA", &OUTFLOWBETA, INT, NO, "5");
  var("NONKEPBOUNDARIES", NONKEPBOUNDARIES, STRING, NO, "NO");
  var("IMPLICITRAD", IMPLICITRAD, STRING, NO, "NO");
  var("LIVEBODIES", LIVEBODIES, STRING, NO, "YES");
  var("WRITEOPTICALDEPTHS", WRITEOPTICALDEPTHS, STRING, NO, "NO");
  var("ANALYTICCOOLING", ANALYTICCOOLING, STRING, NO, "NO");
}
