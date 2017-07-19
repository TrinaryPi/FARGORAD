/* C Header
	* @filename        : Interpret.c
	* @author          : Frederic Masset
	* @last_modified_by: trinarypi
	* @last_modified   : 2017/06/26 14:56
	* @description     :
*/
#include "mp.h"
#define MAXVARIABLES 500

extern int      begin_i;
extern boolean  OpenInner, RadiationDebug, debug;
static Param    VariableSet[MAXVARIABLES];
static int      VariableIndex = 0;
static int	FirstStep = YES;
static clock_t  First, Preceeding, Current, FirstUser, CurrentUser, PreceedingUser;
static long	Ticks;
boolean         FastTransport = YES, GuidingCenter = NO, Indirect_Term = YES;
boolean         IsDisk = YES, HydroOn = YES,  Corotating = NO;
boolean         NonReflecting = NO, Evanescent = NO;
boolean         ClosedInner = NO, ViscousInner = NO, EvanescentInner = NO;
boolean         OuterSourceMass = NO, EvanescentOuter = NO, OpenOuter = NO;
boolean         Write_Density = YES, Write_Velocity = YES, Write_Energy = NO, Sigma_Taper = NO, Sigma_Cavity;
boolean         Write_Eccentricity = NO, Write_Pericentre = NO;
boolean         Write_Temperature = NO, Write_DivV = NO, Write_Qplus = NO;
boolean         SelfGravity = NO, SGZeroMode = NO, ZMPlus = NO;
boolean         Adiabatic = NO, Cooling = NO;
boolean         CICPlanet = NO, ForcedCircular = NO;
boolean         ConfigMp = NO, ConfigPos = NO;
boolean         DiscMassTaper = NO;
boolean         CustomCooling = NO;
boolean         RadCooling = NO, Irradiation = NO, RadTransport = NO, VarDiscHeight = NO, ImplicitRayTracingHeating = NO, ExplicitRayTracingHeating = NO, RayTracingHeating = NO;
boolean         Write_DiscHeight = NO, Write_Qminus = NO, Write_Qirr = NO, Write_QirrRT = NO, Write_Kappa = NO, Write_Coeffs = NO, Write_Residual = NO, Write_Rfld = NO, Write_Taus;
boolean         InnerBCCons = NO, InnerBCGrad = NO, InnerBCExtrap = NO, OuterBCCons = NO, OuterBCGrad = NO, OuterBCExtrap = NO;
boolean         TempInit = NO;
boolean         Residual_Max = NO, Residual_Diff = NO, Relative_Max = NO, Relative_Diff = NO;
boolean         OptimalOmega = NO, ChebyshevOmega = NO;
boolean         NoCFL = NO;
boolean         LinPap1985_Opacity = NO;
boolean         BellLin1994_Opacity = NO;
boolean         PowerLaw_Opacity = NO;
boolean         RadiativeOnly = NO;
// new - Constant factor between Rosseland opacity and visible opacity (i.e. Bitsch2013 etc) (6/01/2017)
boolean         BitschSKappa = NO;
boolean         HubenySKappa = NO;
// new - For ray-tracing, is there material outside of the grid domain?
boolean         EmptyCavity = YES;
boolean         StellarSmoothing = YES;
boolean         NonKepBoundaries = NO;
boolean         ImplicitRadiative = NO;
boolean         LiveBodies = YES;
boolean         OpacitySmoothing = NO;
boolean         Write_OpticalDepths = NO;
boolean         ExplicitRadTransport = NO;
boolean         AnalyticCooling = NO;
boolean         Constant_Opacity = NO;
boolean					Relative_Source = NO;


void var(name, ptr, type, necessary, deflt)
     char           *name;
     char           *ptr;
     int             type;
     int             necessary;
     char           *deflt;
{
  real valuer;
  int valuei;
  float temp;

  sscanf (deflt, "%f", &temp);
  valuer = (real) (temp);
  valuei = (int) valuer;
  strcpy(VariableSet[VariableIndex].name, name);
  VariableSet[VariableIndex].variable = ptr;
  VariableSet[VariableIndex].type = type;
  VariableSet[VariableIndex].necessary = necessary;
  VariableSet[VariableIndex].read = NO;
  if (necessary == NO) {
    if (type == INT) {
      *((int *) ptr) = valuei;
    } else if (type == REAL) {
      *((real *) ptr) = valuer;
    } else if (type == STRING) {
      strcpy (ptr, deflt);
    }
  }
  VariableIndex++;
}

void ReadVariables(filename)
     char *filename;
{
  char            nm[300], s[350],stringval[290];
  char           *s1;
  float           temp;
  real            valuer;
  int             i, found, valuei, success, type;
  int            *ptri;
  real           *ptrr;
  FILE           *input;
  
  InitVariables();
  input = fopen(filename, "r");
  if ( input == NULL ) {
    mastererr ("Unable to read '%s'. Program stopped.\n",filename);
    prs_exit(1);
  }
  mastererr ("Reading parameters file '%s'.\n", filename);
  while ( fgets(s, 349, input) != NULL ) {
    success = sscanf(s, "%s ", nm);
    if (( nm[0] != '#' ) && ( success == 1 )) {  /* # begins a comment line */
      s1 = s + strlen(nm);
      sscanf(s1 + strspn(s1, "\t :=>_"), "%f", &temp);
      sscanf(s1 + strspn(s1, "\t :=>_"), "%289s ", stringval);
      valuer = (real) temp;
      valuei = (int) temp;
      for (i = 0; i < strlen(nm); i++) {
        nm[i] = (char) toupper(nm[i]);
      }
      found = NO;
      for (i = 0; i < VariableIndex; i++) {
        if ( strcmp(nm, VariableSet[i].name) == 0 ) {
          if ( VariableSet[i].read == YES ) {
            mastererr("Warning : %s defined more than once.\n", nm);
          }
          found = YES;
          VariableSet[i].read = YES;
          ptri = (int *) (VariableSet[i].variable);
          ptrr = (real *) (VariableSet[i].variable);
          if ( VariableSet[i].type == INT ) {
            *ptri = valuei;
          } else if ( VariableSet[i].type == REAL ) {
            *ptrr = valuer;
          } else if ( VariableSet[i].type == STRING ) {
            strcpy (VariableSet[i].variable, stringval);
          }
        }
      }
      if ( found == NO ) {
        mastererr("Warning : variable %s defined but non-existent in code.\n", nm);
      }
    }
  }
  found = NO;
  for (i = 0; i < VariableIndex; i++) {
    if (( VariableSet[i].read == NO ) && ( VariableSet[i].necessary == YES )) {
      if ( found == NO ) {
        mastererr("Fatal error : undefined mandatory variable(s):\n");
        found = YES;
      }
      mastererr("%s\n", VariableSet[i].name);
    }
    if (found == YES) {
      prs_exit(1);
    }
  }
  found = NO;
  for (i = 0; i < VariableIndex; i++) {
    if ( VariableSet[i].read == NO ) {
      if ( found == NO ) {
        mastererr("Secondary variables omitted :\n");
        found = YES;
      }
      if ( (type = VariableSet[i].type) == REAL ) {
        mastererr("%s ;\t Default Value : %.5g\n", VariableSet[i].name, *((real *) VariableSet[i].variable));
      }
      if ( type == INT ) {
        mastererr("%s ;\t Default Value : %d\n", VariableSet[i].name, *((int *) VariableSet[i].variable));
      }
      if ( type == STRING ) {
        mastererr("%s ;\t Default Value : %s\n", VariableSet[i].name, VariableSet[i].variable);
      }
    }
  }

  /* Inner and Outer Hydro Boundary Conditions */
  if (( *OPENINNERBOUNDARY == 'O' ) || ( *OPENINNERBOUNDARY == 'o' )) {
    OpenInner = YES;
  }
  if (( *OPENINNERBOUNDARY == 'C' ) || ( *OPENINNERBOUNDARY == 'c' )) {
    ClosedInner = YES;
  }
  if (( *OPENINNERBOUNDARY == 'N' ) || ( *OPENINNERBOUNDARY == 'n' )) {
    NonReflecting = YES;
  }
  if (( *OPENINNERBOUNDARY == 'E' ) || ( *OPENINNERBOUNDARY == 'e' )) {
    EvanescentInner = YES;
  }
  if (( *OPENINNERBOUNDARY == 'V' ) || ( *OPENINNERBOUNDARY == 'v' )) {
    ViscousInner = YES;
  }
  if (( *OPENOUTERBOUNDARY == 'O' ) || ( *OPENOUTERBOUNDARY == 'o' )) {
    OpenOuter = YES;
  }
  if (( *OPENOUTERBOUNDARY == 'W' ) || ( *OPENOUTERBOUNDARY == 'w' )) {
    EvanescentOuter = YES;
  }
  if (( EvanescentInner ) && (EvanescentOuter )) {
    Evanescent = YES;
    EvanescentInner = NO;
    EvanescentOuter = NO;
  }
  if (( *NONKEPBOUNDARIES == 'Y' ) || ( *NONKEPBOUNDARIES == 'y' )) {
    NonKepBoundaries = YES;
  }
  if (( *OUTERSOURCEMASS == 'y' ) || ( *OUTERSOURCEMASS == 'Y' )) {
    OuterSourceMass = YES;
  }

  /* Hydro Numerical Method Parameters */
  if (( *DISK == 'N' ) || ( *DISK == 'n' )) {
    IsDisk = NO;
  }
  if (( *HYDRO == 'N' ) || ( *HYDRO == 'n' )) {
    HydroOn = NO;
  }
  if (( *LIVEBODIES == 'N' ) || ( *LIVEBODIES == 'n' )) {
    LiveBodies = NO;
  }
  if (( *FRAME == 'C' ) || ( *FRAME == 'c' )) {
    Corotating = YES;
  }
  if (( *FRAME == 'G' ) || ( *FRAME == 'g' )) {
    Corotating = YES;
    GuidingCenter = YES;
  }
  if (( *NOCFL == 'Y' ) || ( *NOCFL == 'y' )) {
    NoCFL = YES;
  }
  if (( *INDIRECTTERM == 'N' ) || ( *INDIRECTTERM == 'n' )) {
    Indirect_Term = NO;
  }
  if (( *CONFIGMP == 'Y' ) || ( *CONFIGMP == 'y' )) {
    ConfigMp = YES;
  }
  if (( *CONFIGPOS == 'Y' ) || ( *CONFIGPOS == 'y' )) {
    ConfigPos = YES;
  }
  if (( *EXCLUDEHILL == 'Y' ) || ( *EXCLUDEHILL == 'y' )) {
    ExcludeHill = YES;
  }
  if (( *CICPLANET == 'Y' ) || ( *CICPLANET == 'y' )) {
    CICPlanet = YES;
  }
  if (( *FORCEDCIRCULAR == 'Y' ) || ( *FORCEDCIRCULAR == 'y' )) {
    ForcedCircular = YES;
  }
  if (( ALPHAVISCOSITY != 0.0 ) && ( VISCOSITY != 0.0 )) {
    mastererr ("You cannot use at the same time\n");
    mastererr ("VISCOSITY and ALPHAVISCOSITY.\n");
    mastererr ("Edit the parameter file so as to remove\n");
    mastererr ("one of these variables and run again.\n");
    prs_exit (1);
  }
  if ( ALPHAVISCOSITY != 0.0 ) {
    ViscosityAlpha = YES;
    masterprint ("Viscosity is of alpha type\n");
  }
  if (( *STELLARSMOOTHING == 'N' ) || ( *STELLARSMOOTHING == 'n' )) {
    StellarSmoothing = NO;
  }
  if (( THICKNESSSMOOTHING != 0.0 ) && ( ROCHESMOOTHING != 0.0 )) {
    mastererr ("You cannot use at the same time\n");
    mastererr ("`ThicknessSmoothing' and `RocheSmoothing'.\n");
    mastererr ("Edit the parameter file so as to remove\n");
    mastererr ("one of these variables and run again.\n");
    prs_exit (1);
  }
  if (( THICKNESSSMOOTHING <= 0.0 ) && ( ROCHESMOOTHING <= 0.0 )) {
    mastererr ("A non-vanishing potential smoothing length is required.\n");
    mastererr ("Please use either of the following variables:\n");
    mastererr ("`ThicknessSmoothing' *or* `RocheSmoothing'.\n");
    mastererr ("before launching the run again.\n");
    prs_exit (1);
  }
  if ( ROCHESMOOTHING != 0.0 ) {
    RocheSmoothing = YES;
    masterprint ("Planet potential smoothing scales with their Hill sphere.\n");
  }
  if ((( ROCHESMOOTHING != 0.0 ) || ( THICKNESSSMOOTHING != 0.0 )) && ( StellarSmoothing == NO )) {
    masterprint("The option has been chosen to treat the stellar potential as non-point source. \n");
    masterprint("Thickness smoothing or Roche smoothing will not be used. \n");
    masterprint("Please check if this is what you require. \n");
  }
  if (( *DISCMASSTAPER == 'Y' ) || ( *DISCMASSTAPER == 'y')) {
    DiscMassTaper = YES;
  }
  if ((( DMTTAU != 0.0 ) || ( DMTOUTPUTSTART != 0.0 )) && ( DiscMassTaper == NO )) {
    masterprint("DiscMassTaper = NO, so DMTTSTART and/or DMTTAU do not need a non-zero value. \n");
    masterprint("Please check the input parameter file to see if this is the set-up you want");
  }

  /* Equation of State Parameters */
  if (( *ADIABATIC == 'Y' ) || ( *ADIABATIC == 'y' )) {
    Adiabatic = YES;
    Write_Temperature = YES;
  }
  if (( Adiabatic ) && ( ADIABATICINDEX == 1 )) {
    masterprint ("You cannot have Adiabatic = YES and AdiabatcIndex = 1. I decided to put Adiabatic = No, to simulate a locally isothermal equation of state. Please check that it is what you really wanted to do!\n");
    Adiabatic = NO;
  }
  if (( *COOLING == 'Y' ) || ( *COOLING == 'y' )) {
    Cooling = YES;
  }
  if (( *CUSTOMCOOLING == 'Y' ) || ( *CUSTOMCOOLING == 'y' )) {
    Cooling = YES;
    CustomCooling = YES;
  }
  if ((*REVERSETEMPERATUREINIT == 'Y') || (*REVERSETEMPERATUREINIT == 'y')) {
    TempInit = YES;
    if (Adiabatic == YES) {
      masterprint("This option is used for testing\n");
      masterprint("Please check you want to override the default Temperature and Energy initialisation.\n");
    }
  }

  /* Self-Gravity Module Parameters */
  if (( *SELFGRAVITY == 'Y' ) || ( *SELFGRAVITY == 'y' )) {
    SelfGravity = YES;
  }
  if (( *SELFGRAVITY == 'Z' ) || ( *SELFGRAVITY == 'z' )) {
    SelfGravity = YES;
    SGZeroMode = YES;
  }
  if (( *ZMPLUS == 'Y' ) || ( *ZMPLUS == 'y' )) ZMPlus = YES;
  if (( ZMPlus ) && ( !SGZeroMode )) {
    masterprint ("This is not very meaningfull to involve the anisotropic pressure model (ZMPlus=Yes) without taking into account the axisymmetric component of the disk self-gravity. I decided to put ZMPlus = No. Please check again!");
    ZMPlus = NO;
  }

  /* Disk Grid/Set-up Parameters */
  if (( *GRIDSPACING == 'L' ) || ( *GRIDSPACING == 'l' )) {
    LogGrid = YES;
  }
  if (( *SIGMATAPER == 'Y' ) || ( *SIGMATAPER == 'y' )) {
    Sigma_Taper = YES;
  }
  if (( *SIGMACAVITY == 'Y' ) || ( *SIGMACAVITY == 'y' )) {
    Sigma_Cavity = YES;
  }

  /* General Radiation Module Parameters */
  if (( *RADIATIVEONLY == 'Y' ) || ( *RADIATIVEONLY == 'y' )) {
    RadiativeOnly = YES;
  }
  if (( *RADCOOLING == 'Y' ) || ( *RADCOOLING == 'y' )) {
    RadCooling = YES;
  }
  if (( RadCooling ) && (( Cooling ) || ( CustomCooling ))) {
    masterprint("There are too many types/sources of cooling present. I am turning radiative cooling off. Please recheck input parameter files! \n");
    RadCooling = NO;
  }
  if (( *ANALYTICCOOLING == 'Y' ) || ( *ANALYTICCOOLING == 'y' )) {
    AnalyticCooling = YES;
  }
  if (( *IRRADIATION == 'Y' ) || ( *IRRADIATION == 'y' )) {
    Irradiation = YES;
  }
  if (( *RAYTRACINGHEATING == 'I' ) || ( *RAYTRACINGHEATING == 'i' )) {
    ImplicitRayTracingHeating = YES;
    RayTracingHeating = YES;
  }
  if (( *RAYTRACINGHEATING == 'E' ) || ( *RAYTRACINGHEATING == 'e' )) {
    ExplicitRayTracingHeating = YES;
    RayTracingHeating = YES;
  }
  if (( *EMPTYCAVITY == 'N' ) || ( *EMPTYCAVITY == 'N' )) {
    EmptyCavity = NO;
  }
  
  /* Opacity Parameters */
  if (( *OPACITYTABLE == 'B' ) || ( *OPACITYTABLE == 'b' )) {
    BellLin1994_Opacity = YES;
  }
  if (( *OPACITYTABLE == 'L' ) || ( *OPACITYTABLE == 'l' )) {
    LinPap1985_Opacity = YES;
  }
  if (( *OPACITYTABLE == 'P' ) || ( *OPACITYTABLE == 'p' )) {
    PowerLaw_Opacity = YES;
  }
  if (( *OPACITYTABLE == 'C' ) || ( *OPACITYTABLE == 'c' )) {
    Constant_Opacity = YES;
  }
  if (( *OPACITYSMOOTHING == 'Y' ) || ( *OPACITYSMOOTHING == 'y' )) {
    OpacitySmoothing = YES;
  }
  if (( *BITSCHSKAPPA == 'Y' ) || ( *BITSCHSKAPPA == 'y' )) {
    BitschSKappa = YES;
  }
  if (( *HUBENYSKAPPA == 'Y' ) || ( *HUBENYSKAPPA == 'y' )) {
    HubenySKappa = YES;
  }
  

  /* Radiative Diffusion (FLD/SOR) Set-up */
  if (( *RADTRANSPORT == 'Y' ) || ( *RADTRANSPORT == 'y' )) {
    RadTransport = YES;
  }
  if (( *RADTRANSPORT == 'E' ) || ( *RADTRANSPORT == 'e' )) {
    RadTransport = YES;
    ExplicitRadTransport = YES;
  }
  if (( *RADTRANSPORT == 'I' ) || ( *RADTRANSPORT == 'i' )) {
    RadTransport = YES;
    ExplicitRadTransport = NO;
  }
  if (( *INNERTEMPBC == 'C' ) || ( *INNERTEMPBC == 'c' )) {
    InnerBCCons = YES;
  }
  if (( *INNERTEMPBC == 'G' ) || ( *INNERTEMPBC == 'g' )) {
    InnerBCGrad = YES;
  }
  if (( *INNERTEMPBC == 'E' ) || ( *INNERTEMPBC == 'e' )) {
    InnerBCExtrap = YES;
  }
  if (( *OUTERTEMPBC == 'C' ) || ( *OUTERTEMPBC == 'c' )) {
    OuterBCCons = YES;
  }
  if (( *OUTERTEMPBC == 'G' ) || ( *OUTERTEMPBC == 'g' )) {
    OuterBCGrad = YES;
  }
  if (( *OUTERTEMPBC == 'E' ) || ( *OUTERTEMPBC == 'e' )) {
    OuterBCExtrap = YES;
  }
  if (( RadTransport ) && (( InnerBCCons == NO ) && ( InnerBCGrad == NO ) && ( InnerBCExtrap == NO ))) {
    masterprint("You have not set any value for a boundary condition at the inner disc edge\n");
    masterprint("I will continue, but check this is what you want\n");
  }
  if (( RadTransport ) && (( OuterBCCons == NO ) && ( OuterBCGrad == NO ) && ( OuterBCExtrap == NO ))) {
    masterprint("You have not set any value for a boundary condition at the outer disc edge\n");
    masterprint("I will continue, but check this is what you want\n");
  }
  if (( *OPTIMALOMEGA == 'Y' ) || ( *OPTIMALOMEGA == 'y' )) {
    OptimalOmega = YES;
  }
  if (( *OPTIMALOMEGA == 'C' ) || ( *OPTIMALOMEGA == 'c' )) {
    ChebyshevOmega = YES;
  }
  if (( *NORMRESIDUAL == 'M' ) || ( *NORMRESIDUAL == 'm' )) {
    Residual_Max = YES;
  }
  if (( *NORMRELATIVE == 'M' ) || ( *NORMRELATIVE == 'm' )) {
    Relative_Max = YES;
  }
  if (( *NORMRESIDUAL == 'D' ) || ( *NORMRESIDUAL == 'd' )) {
    Residual_Diff = YES;
  }
  if (( *NORMRELATIVE == 'D' ) || ( *NORMRELATIVE == 'd' )) {
    Relative_Diff = YES;
  }
  if (( Residual_Max == NO ) && ( Residual_Diff == NO ) && ( Relative_Max == NO ) && ( Relative_Diff == NO )) {
    Relative_Source = YES;
  }
  if (( RadTransport ) && ((( *IMPLICITRAD == 'Y' ) || ( *IMPLICITRAD == 'Y' )))) {
    ImplicitRadiative = YES;
    if ( ExplicitRayTracingHeating ) {
      ExplicitRayTracingHeating = NO;
      ImplicitRayTracingHeating = YES;
    }
  }
  if ( ExplicitRadTransport ) {
    if ( ImplicitRadiative ) {
      ImplicitRadiative = NO;
    }
    if ( ImplicitRayTracingHeating ) {
      ImplicitRayTracingHeating = NO;
      ExplicitRayTracingHeating = YES;
    }
  }
  if (( RadCooling ) || ( RadTransport ) || ( Irradiation ) || ( RayTracingHeating )) {
    VarDiscHeight = YES;
  }

  /* Consistency Checks and Debugging */
  if (( Adiabatic == NO ) && (( RadCooling ) || ( RadTransport ) || ( Irradiation ) || ( RayTracingHeating ))) {
    masterprint("Radiative cooling, transport, ray-trace heating or irradiation make no sense with an isothermal EoS. Please recheck input parameter files! \n");
    RadCooling = NO;
    RadTransport = NO;
    Irradiation = NO;
    RayTracingHeating = NO;
    VarDiscHeight = NO;
  }
  if (( debug ) && (( RadCooling ) || ( RadTransport ) || ( Irradiation ) || ( RayTracingHeating ))) {
  	RadiationDebug = YES;
  }

  /* Output Parameters */
  if (( *WRITEVELOCITY == 'N' ) || ( *WRITEVELOCITY == 'n' )) {
    Write_Velocity = NO;
  }
  if (( *WRITEDENSITY == 'N' ) || ( *WRITEDENSITY == 'n' )) {
    Write_Density = NO;
  }
  if (( *WRITEENERGY == 'Y' ) || ( *WRITEENERGY == 'y' )) {
    Write_Energy = YES;
  }
  if (( *WRITETEMPERATURE == 'Y' ) || ( *WRITETEMPERATURE == 'y' )) {
    Write_Temperature = YES;
  }
  if (( *WRITEDIVV == 'Y' ) || ( *WRITEDIVV == 'y' )) {
    Write_DivV = YES;
  }
  if (( *WRITEQPLUS == 'Y' ) || ( *WRITEQPLUS == 'y' )) {
    Write_Qplus = YES;
  }
  if (( *WRITEECCENTRICITY == 'Y' ) || ( *WRITEECCENTRICITY == 'y' )) {
    Write_Eccentricity = YES;
  }
  if (( *WRITEPERICENTRE == 'Y' ) || ( *WRITEPERICENTRE == 'y' )) {
    Write_Pericentre = YES;
  }
  if (( VarDiscHeight ) && (( *WRITEDISCHEIGHT == 'Y' ) || ( *WRITEDISCHEIGHT == 'y' ))) {
    Write_DiscHeight = YES;
  }
  if ((( *WRITEQMINUS == 'Y' ) || ( *WRITEQMINUS == 'y' )) && ( RadCooling )) {
    Write_Qminus = YES;
  }
  if ((( *WRITEQIRR == 'Y' ) || ( *WRITEQIRR == 'y' )) && ( Irradiation )) {
    Write_Qirr = YES;
  }
  if ((( *WRITERTQIRR == 'Y' ) || ( *WRITERTQIRR == 'y' )) && ( RayTracingHeating )) {
    Write_QirrRT = YES;
  }
  if (( *WRITEKAPPA == 'Y' ) || ( *WRITEKAPPA == 'y' )) {
    Write_Kappa = YES;
  }
  if ((( *WRITECOEFFS == 'Y' ) || ( *WRITECOEFFS == 'y' )) && ( RadTransport )) {
    Write_Coeffs = YES;
  }
  if ((( *WRITERESIDUAL == 'Y' ) || ( *WRITERESIDUAL == 'y' )) && ( RadTransport )) {
    Write_Residual = YES;
  }
  if ((( *WRITERFLD == 'Y' ) || ( *WRITERFLD == 'y' )) && ( RadTransport )) {
    Write_Rfld = YES;
  }
  if (( *WRITETAUS == 'Y' ) || ( *WRITETAUS == 'y' )) {
    Write_Taus = YES;
  }
  if ((( *WRITEOPTICALDEPTHS == 'Y' ) || ( *WRITEOPTICALDEPTHS == 'y' )) && (( Irradiation ) || ( RadCooling ))) {
    Write_OpticalDepths = YES;
  }
  if (( *WRITEENERGY == 'N' ) || ( *WRITEENERGY == 'n' )) {
    Write_Energy = NO;
  }
  if ( OverridesOutputdir ) {
    sprintf (OUTPUTDIR, "%s", NewOutputdir);
  }
  /* Add a trailing slash to OUTPUTDIR if needed */
  if ( *(OUTPUTDIR+strlen(OUTPUTDIR)-1) != '/' ) {
    strcat (OUTPUTDIR, "/");
  }
}

void PrintUsage (execname)
     char *execname;
{
  mastererr("Usage : %s [-abcdeimnptvz] [-(0-9)] [-s number] [-f scaling] parameters file\n", execname);
  mastererr("\n-a : Monitor mass and angular momentum at each timestep\n");
  mastererr("-b : Adjust azimuthal velocity to impose strict centrifugal balance at t=0\n");
  mastererr("-c : Sloppy CFL condition (checked at each DT, not at each timestep)\n");
  mastererr("-d : Print some debugging information on 'stdout' at each timestep\n");
  mastererr("-e : Activate EU test problem torque file output\n");
  mastererr("-f : Scale density array by 'scaling'. Useful to increase/decrease\n");
  mastererr("     disk surface density after a restart, for instance.            \n");
  mastererr("-i : tabulate Sigma profile as given by restart files\n");
  mastererr("-m : Merge output files from different CPUs\n");
  mastererr("-n : Disable simulation. The program just reads parameters file\n");
  mastererr("-o : Overrides output directory of input file.\n");
  mastererr("-p : Give profiling information at each time step\n");
  mastererr("-s : Restart simulation, taking #'number' files as initial conditions\n");
  mastererr("-t : Monitor CPU time usage at each time step\n");
  mastererr("-v : Verbose mode. Tells everything about parameters file\n");
  mastererr("-z : fake sequential built when evaluating sums on HD meshes\n");
  mastererr("-(0-9) : only write initial (or restart) HD meshes,\n");
  mastererr("     proceed to the next nth output and exit\n");
  mastererr("     This option must stand alone on one switch (-va -4 is legal, -v4a is not)\n");
  prs_exit (1);
}

real TellNbOrbits (time)
     real time;
{
  return time/2.0/PI*sqrt(G*1.0/1.0/1.0/1.0);
}

real TellNbOutputs (time)
     real time;
{
  return (time/DT/NINTERM);
}

void TellEverything ()
{
  real temp, nbfileoutput;
  if (!CPU_Master) return;
  printf ("\nDisc properties:\n");
  printf ("----------------\n");
  printf ("Inner Radius          : %g\n", RMIN);
  printf ("Outer Radius          : %g\n", RMAX);
  printf ("Aspect Ratio          : %g\n", ASPECTRATIO);
  printf ("VKep at inner edge    : %.3g\n", sqrt(G*1.0*(1.-0.0)/RMIN));
  printf ("VKep at outer edge    : %.3g\n", sqrt(G*1.0/RMAX));
  temp=2.0*PI*SIGMA0/(2.0-SIGMASLOPE)*(pow(RMAX,2.0-SIGMASLOPE) - pow(RMIN,2.0-SIGMASLOPE));	/* correct this and what follows... */
  printf ("Initial Disk Mass             : %g\n", temp);
  temp=2.0*PI*SIGMA0/(2.0-SIGMASLOPE)*(1.0 - pow(RMIN,2.0-SIGMASLOPE));
  printf ("Initial Mass inner to r=1.0  : %g \n", temp);
  temp=2.0*PI*SIGMA0/(2.0-SIGMASLOPE)*(pow(RMAX,2.0-SIGMASLOPE) - 1.0);
  printf ("Initial Mass outer to r=1.0  : %g \n", temp);
  printf ("Travelling time for acoustic density waves :\n");
  temp = 2.0/3.0/ASPECTRATIO*(pow(RMAX,1.5)-pow(RMIN,1.5)); // Any mentions of ASPECTRATIO from this point are probably wrong if there is an adiabatic EoS or any radiative stuff running!
  printf (" * From Rmin to Rmax  : %.2g = %.2f orbits ~ %.1f outputs\n", temp, TellNbOrbits(temp), TellNbOutputs(temp));
  temp = 2.0/3.0/ASPECTRATIO*(pow(RMAX,1.5)-pow(1.0,1.5));
  printf (" * From r=1.0 to Rmax: %.2g = %.2f orbits ~ %.1f outputs\n", temp, TellNbOrbits(temp), TellNbOutputs(temp));
  temp = 2.0/3.0/ASPECTRATIO*(pow(1.0,1.5)-pow(RMIN,1.5));
  printf (" * From r=1.0 to Rmin: %.2g = %.2f orbits ~ %.1f outputs\n", temp, TellNbOrbits(temp), TellNbOutputs(temp));
  temp = 2.0*PI*sqrt(RMIN*RMIN*RMIN/G/1.0);
  printf ("Orbital time at Rmin  : %.3g ~ %.2f outputs\n", temp, TellNbOutputs(temp));
  temp = 2.0*PI*sqrt(RMAX*RMAX*RMAX/G/1.0);
  printf ("Orbital time at Rmax  : %.3g ~ %.2f outputs\n", temp, TellNbOutputs(temp));
  printf ("Sound speed :\n");
  printf (" * At unit radius     : %.3g\n", ASPECTRATIO*sqrt(G*1.0));
  printf (" * At outer edge      : %.3g\n", ASPECTRATIO*sqrt(G*1.0/RMAX));
  printf (" * At inner edge      : %.3g\n", ASPECTRATIO*sqrt(G*1.0/RMIN));
  printf ("\nGrid properties:\n");
  printf ("----------------\n");
  printf ("Number of rings       : %d\n", NRAD);
  printf ("Number of sectors     : %d\n", NSEC);
  printf ("Total cells           : %d\n", NRAD*NSEC);
  printf ("\nOutputs properties:\n");
  printf ("-------------------\n");
  printf ("Time increment between outputs : %.3f = %.3f orbits\n", NINTERM*DT, TellNbOrbits(NINTERM*DT));
  printf ("At each output #i, the following files are written:\n");
  printf ("gasdens[i].dat : %d bytes\n",(int)(GLOBALNRAD*NSEC*sizeof(real)));
  printf ("gasvrad[i].dat : %d bytes\n",(int)(GLOBALNRAD*NSEC*sizeof(real)));
  printf ("gasvtheta[i].dat : %d bytes\n",(int)(GLOBALNRAD*NSEC*sizeof(real)));
  if ( Adiabatic ) {
    printf ("gasTemperature[i].dat : %d bytes\n",(int)(GLOBALNRAD*NSEC*sizeof(real)));
  }
  if ( AdvecteLabel ) {
    printf ("gaslabel[i].dat : %d bytes\n",(int)(GLOBALNRAD*NSEC*sizeof(real)));
  }
  printf ("There will be in total %d outputs\n", NTOT/NINTERM);
  printf ("(which correspond to an elapsed time = %.3f or to %.2f orbits)\n", NTOT*DT, TellNbOrbits(NTOT*DT));
  nbfileoutput = 3.0;
  if ( Adiabatic ) {
    nbfileoutput += 1.0;
  }
  if ( AdvecteLabel ) {
    nbfileoutput += 1.0;
  }
  temp =nbfileoutput*GLOBALNRAD*NSEC*sizeof(real);
  temp *= (real)NTOT/(real)NINTERM;
  temp /= 1024.0*1024.0;
  printf ("So the code will produce ~%.2f Mbytes of data\n", temp);
  printf ("Check (eg by issuing a 'df' command) that you have enough disk space,\n");
  printf ("otherwise you will get a system full and the code will stop.\n");
  fflush (stdout);
}

void GiveTimeInfo (number)
     int number;
{
  struct tms buffer;
  real total, last, mean, totalu;

  Current = times (&buffer);
  CurrentUser = buffer.tms_utime;
  if ( FirstStep ) {
    First = Current;
    FirstUser = CurrentUser;
    fprintf (stderr, "Time counters initialized\n");
    FirstStep = NO;
    Ticks = sysconf (_SC_CLK_TCK);
  } else {
    total = (real)(Current - First)/Ticks;
    totalu= (real)(CurrentUser-FirstUser)/Ticks;
    last  = (real)(CurrentUser - PreceedingUser)/Ticks;
    number -= begin_i/NINTERM;
    mean  = totalu / number;
    fprintf (stderr, "Total Real Time elapsed    : %.3f s\n", total);
    fprintf (stderr, "Total CPU Time of process  : %.3f s (%.1f %%)\n", totalu, 100.*totalu/total);
    fprintf (stderr, "CPU Time since last time step : %.3f s\n", last);
    fprintf (stderr, "Mean CPU Time between time steps : %.3f s\n", mean);
    fprintf (stderr, "CPU Load on last time step : %.1f %% \n", (real)(CurrentUser-PreceedingUser)/(real)(Current-Preceeding)*100.);
    
  }	
  PreceedingUser = CurrentUser;
  Preceeding = Current;
}

void InitSpecificTime (profiling, process_name, title)
     boolean profiling;
     TimeProcess *process_name;
     char *title;
{
  struct tms buffer;

  if (profiling == NO) {
    return;
  }
  Ticks = sysconf (_SC_CLK_TCK);
  times (&buffer);
  process_name->clicks = buffer.tms_utime;
  strcpy (process_name->name, title);
}

void GiveSpecificTime (profiling, process_name)
     boolean profiling;
     TimeProcess process_name;
{
  struct tms buffer;
  long ticks;
  real t;

  if (profiling == NO) {
    return;
  }
  Ticks = sysconf (_SC_CLK_TCK);
  times (&buffer);
  ticks = buffer.tms_utime - process_name.clicks;
  t = (real)ticks / (real)Ticks;
  fprintf (stderr, "Time spent in %s : %.3f s\n", process_name.name, t);
}

