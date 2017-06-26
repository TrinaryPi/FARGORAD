/* C Header
	* @filename        : main.c
	* @author          : Matthew Mutter (trinarypi)
	* @last_modified_by: trinarypi
	* @last_modified   : 2017/06/26 16:58
	* @description     :
*/
#include "mp.h"

boolean         TimeToWrite, Restart=NO, OpenInner=NO, BinaryOn=NO, DiscElem=NO, RadiationDebug=NO, PreInitialisation;
int             TimeStep = 0, begin_i = 0, NbRestart = 0, verbose = NO;
int             dimfxy = 11;
static int      InnerOutputCounter=0, StillWriteOneOutput;
extern real     LostMass;
extern boolean  Corotating;
extern boolean  SelfGravity, SGZeroMode, Adiabatic, DiscMassTaper;
real            ScalingFactor = 1.0;
real            Abinary;
Pair            bin_acc, acc_rate;
extern boolean  RadCooling, Irradiation, RadTransport, RayTracingHeating;
extern boolean  TempInit;

int
main(argc, argv)
     int argc;
     char *argv[];
{
  PolarGrid   *gas_density;
  PolarGrid   *gas_v_rad;
  PolarGrid   *gas_v_theta;
  PolarGrid   *gas_energy; 
  PolarGrid   *gas_label;
  PolarGrid   *gas_e_cell;
  PolarGrid   *gas_per_cell;

  int          i;
  real         foostep = 0.;
  boolean      disable = NO, TimeInfo = NO, Profiling = NO;
  boolean      Stockholm = NO, updatevelocities = NO;
  TimeProcess  t_Hydro;
  char         ParameterFile[256];
  PlanetarySystem *sys;
  BinarySystem *bsys;
  Force *force;
  int FirstRestartOutput = 0;
  PreInitialisation = YES;
  mlinner = 0.0;
  mlouter = 0.0;
  
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &CPU_Rank);
  MPI_Comm_size (MPI_COMM_WORLD, &CPU_Number);
  CPU_Master = (CPU_Rank == 0 ? 1 : 0);
  setfpe ();  /* Control behavior for floating point
		 exceptions trapping (default is not to do anything) */
  if (argc == 1) {
    PrintUsage (argv[0]);
  }
  strcpy (ParameterFile, "");
  for (i = 1; i < argc; i++) {
    if (*(argv[i]) == '-') {
      if (strspn (argv[i], "-secndovtpfamzibkqr0123456789") != strlen (argv[i])){
	      PrintUsage (argv[0]);
      }
      if ( strchr (argv[i], 'n') ) {
        disable = YES;
      }
      if ( strchr (argv[i], 'e') ) {
        Stockholm = YES;
      }
      if ( strchr (argv[i], 'v') ) {
        verbose = YES;
      }
      if ( strchr (argv[i], 't') ) {
        TimeInfo = YES;
      }
      if ( strchr (argv[i], 'c') ) {
        SloppyCFL = YES;
      }
      if ( strchr (argv[i], 'p') ) {
        Profiling = YES;
      }
      if ( strchr (argv[i], 'd') ) {
        debug = YES;
      }
      if ( strchr (argv[i], 'b') ) {
        CentrifugalBalance = YES;
      }
      if ( strchr (argv[i], 'q') ) {
        BinaryOn = YES;
      }
      if ( strchr (argv[i], 'k') ) {
        DiscElem = YES;
      }
      if ( strchr (argv[i], 'm') ) {
        Merge = YES;
      }
      if ( strchr (argv[i], 'a') ) {
        MonitorIntegral = YES;
      }
      if ( strchr (argv[i], 'z') ) {
        FakeSequential = YES;
      }
      if ( strchr (argv[i], 'r') ) {
        RadiationDebug = YES;
      }
      if ( strchr (argv[i], 'i') ) {
	      StoreSigma = YES; 
	      if ( Adiabatic ) {
	        StoreEnergy = YES;
        }
      }
      if ( strchr (argv[i], '0') ) {
	      OnlyInit = YES;
      }
      if (( argv[i][1] >= '1' ) && ( argv[i][1] <= '9' )) {
	      GotoNextOutput = YES;
	      StillWriteOneOutput = (int)(argv[i][1]-'0');
      }
      if ( strchr (argv[i], 's') ) {
	      Restart = YES;
	      FirstRestartOutput = 1;
	      i++;
	      NbRestart = atoi(argv[i]);
	      if ( (NbRestart < 0) ) {
	        masterprint ("Incorrect restart number\n");
	        PrintUsage (argv[0]);
	      }
      }
      if ( strchr (argv[i], 'o') ) {
	      OverridesOutputdir = YES;
	      i++;
	      sprintf (NewOutputdir, "%s", argv[i]);
      } else {
	      if ( strchr (argv[i], 'f') ) {
          i++;
          ScalingFactor = atof(argv[i]);
          masterprint ("Scaling factor = %g\n", ScalingFactor);
	        if ( ScalingFactor <= 0) {
	          masterprint ("Incorrect scaling factor\n");
	          PrintUsage (argv[0]);
	        }
	      }
      }
    } else {
      strcpy (ParameterFile, argv[i]);
    }
  }
  
  if (( StoreSigma || StoreEnergy ) && !Restart ) {
    mastererr ("You cannot use tabulated surface density\n");
    mastererr ("or surface internal energy in a non-restart run.\n");
    mastererr ("Aborted\n");
    prs_exit (0);
  }
  if ( ParameterFile[0] == 0 ) {
    PrintUsage (argv[0]);
  }
  ReadVariables (ParameterFile);

  if ( Restart == NO ) {
    EmptyTargetFolder();
  }

  SplitDomain ();
  if ( verbose ) {
    TellEverything ();
  }
  if ( disable ) {
    prs_exit (0);
  }
  DumpSources (argc, argv);
  masterprint ("Allocating arrays...");
  fflush (stdout);
  gas_density        = CreatePolarGrid(NRAD, NSEC, "dens");
  gas_v_rad          = CreatePolarGrid(NRAD, NSEC, "vrad");
  gas_v_theta        = CreatePolarGrid(NRAD, NSEC, "vtheta");
  gas_energy         = CreatePolarGrid(NRAD, NSEC, "energy");
  gas_label          = CreatePolarGrid(NRAD, NSEC, "label");
  gas_e_cell         = CreatePolarGrid(NRAD, NSEC, "eccen");
  gas_per_cell       = CreatePolarGrid(NRAD, NSEC, "perihelion");
  masterprint ("done.\n");
  masterprint ("Allocating 1D arrays...");
  FillPolar1DArrays ();
  masterprint ("done.\n");
  masterprint ("Allocating force arrays...");
  force = AllocateForce (dimfxy);
  masterprint ("done.\n");
  
  /* Here planets are initialized feeling star potential but they do
     not feel disk potential */
  masterprint ("Initialising planet and star data structures and systems...");
  sys = InitPlanetarySystem (PLANETCONFIG);
  bsys = InitBinarySystem (STARCONFIG);
  masterprint ("done.\n");

  if ( BinaryOn ) {
    BinaryOn = CountStars (bsys);
  } else {
    bsys->nb = 0;
  }
  
  /* Gas density initialization */
  InitGasDensity (gas_density);
  /* If energy equation is taken into account, we initialize the gas
     thermal energy */
  if ( Adiabatic ) {
    if ( TempInit ) {
      masterprint("Initialising the temperature profile before the energy profile\n");
      Temperature  = CreatePolarGrid(NRAD, NSEC, "Temperature");
      InitGasTemperature (gas_density, gas_energy);
    } else {
      InitGasEnergy (gas_energy);
    }
  }

  if ( SelfGravity ) {
    /* If SelfGravity = YES or Z, planets are initialized feeling disk
       potential. Only the surface density is required to calculate
       the radial self-gravity acceleration. The disk radial and
       azimuthal velocities are not updated */
    compute_selfgravity (gas_density, gas_v_rad, gas_v_theta, foostep, updatevelocities);
    init_planetarysys_withSG (sys);
  }
  
  OmegaFrame = OMEGAFRAME;
  if ( Corotating ) {
    OmegaFrame = GetPsysInfo (sys, FREQUENCY);
  }

  if ( debug ) {
  	PrintBooleanUsage();
  }
  ListEOSParams();

  /* Only gas velocities remain to be initialized */
  Initialization (gas_density, gas_v_rad, gas_v_theta, gas_energy, gas_label, sys);

  if ( DiscMassTaper ) {
    for (i = 0; i < NRAD; i++) {
      SigmaDMT[i] = SigmaMed[i];
      EnergyDMT[i] = EnergyMed[i]; 
    }
  }
  
  if ( BinaryOn ) {
    ListStars(bsys);
  }
  ListPlanets(sys);

  WriteSimVariableFile();

  /* Initial gas_density is used to compute the circumplanetary mass
     with initial density field */
  mdcp0 = CircumPlanetaryMass(gas_density, sys);
  bin_acc.x = 0.0;
  bin_acc.y = 0.0;
  acc_rate.x = 0.0;
  acc_rate.y = 0.0;

  if ( DiscElem ) {
    SolveDiscOrbits(TimeStep, gas_density, gas_v_rad, gas_v_theta, gas_e_cell, gas_per_cell);
  }
  SolveOrbits(sys);
  if ( BinaryOn ) {
    SolveBinOrbits(bsys);
  }
  
  if ( Restart ) {
    begin_i = NbRestart * NINTERM;
    RestartPlanetarySystem (NbRestart, sys);
    if ( BinaryOn ) {
      RestartBinarySystem (NbRestart, bsys);
      bin_acc.x = GetFromBinaryAccretionFile(NbRestart, 2);
      bin_acc.y = GetFromBinaryAccretionFile(NbRestart, 3);
    }
    LostMass = GetfromPlanetFile (NbRestart, 7, 0); /* 0 refers to planet #0 */
    PhysicalTime  = GetfromPlanetFile (NbRestart, 8, 0);
    OmegaFrame  = GetfromPlanetFile (NbRestart, 9, 0);
  } else {
    /* We initialize 'planet[i].dat' file */
    EmptyPlanetSystemFile (sys);
  }
  if ( MonitorIntegral ) {
    CheckMomentumConservation (gas_density, gas_v_theta, sys);
  }
  PhysicalTimeInitial = PhysicalTime;
  MultiplyPolarGridbyConstant (gas_density, ScalingFactor);

  /* Initialisation of the fields required for Radiation Modules */
  if (( Irradiation ) || ( RadTransport ) || ( RadCooling ) || ( RayTracingHeating )) {
    StarTaper = (STARTAPER > 1.0 ? 0 : 1.0);
    InitialiseRadiationModule (gas_density, bsys);
  }

  PreInitialisation = NO;

  for (i = begin_i; i <= NTOT; i++) {
    InnerOutputCounter++;
    if ( InnerOutputCounter == 1 ) {
      InnerOutputCounter = 0;
      WriteBigPlanetSystemFile (sys, TimeStep);
      UpdateLog (force, sys, gas_density, gas_energy, TimeStep, PhysicalTime, dimfxy);
      if ( Stockholm ) {
	      UpdateLogStockholm (sys, gas_density, gas_energy, TimeStep, PhysicalTime);
      }
    }
    if ( NINTERM * (TimeStep = (i / NINTERM) ) == i ) {
      /* Outputs are done here */
      TimeToWrite = YES;
      SendOutput (TimeStep, gas_density, gas_v_rad, gas_v_theta, gas_energy, gas_label, gas_e_cell);
      if ( FirstRestartOutput == 0 ) {
      	WritePlanetSystemFile (sys, TimeStep);
      	if ( BinaryOn ) {
        	WriteBinarySystemFile (bsys, TimeStep);
        	WriteBinaryAccretionFile (TimeStep);
      	}
      }
      if (( OnlyInit ) || (( GotoNextOutput ) && ( !StillWriteOneOutput ))) {
	      MPI_Finalize();
	      return 0;
      }
      FirstRestartOutput = 0;
      StillWriteOneOutput--;
      if ( TimeInfo ) {
        /* Time monitoring is done here */
        GiveTimeInfo (TimeStep);
      }
    } else {
      TimeToWrite = NO;
    }
    
    /* Algorithm loop begins here */
    /***********************/
    /* Hydrodynamical Part */
    /***********************/
    InitSpecificTime (Profiling, &t_Hydro, "Eulerian Hydro algorithms");
    AlgoGas (force, gas_density, gas_v_rad, gas_v_theta, gas_energy, gas_label, sys, bsys, gas_e_cell, TimeStep);
    GiveSpecificTime (Profiling, t_Hydro);
    if ( DiscElem ) {
      SolveDiscOrbits (TimeStep, gas_density, gas_v_rad, gas_v_theta, gas_e_cell, gas_per_cell);
    }
    if ( BinaryOn ) {
      SolveBinOrbits (bsys);
    }
    SolveOrbits (sys);
    
    if ( MonitorIntegral ) {
      CheckMomentumConservation (gas_density, gas_v_theta, sys);
      masterprint ("Gas Momentum   : %.18g\n", GasMomentum (gas_density, gas_v_theta));
      masterprint ("Gas total Mass : %.18g\n", GasTotalMass (gas_density));
      masterprint ("Gas total Energy : %.18g\n", GasTotalEnergy (gas_density, gas_v_rad, gas_v_theta, gas_energy));
    }
  }
  FreePlanetary (sys);
  FreeBinary (bsys);
  FreeForce (force);
  if (( SelfGravity ) && (!SGZeroMode )) {
    rfftwnd_mpi_destroy_plan(SGP_fftplan_forward);
    rfftwnd_mpi_destroy_plan(SGP_fftplan_backward);
  }
  MPI_Finalize ();
  return 0;
}
