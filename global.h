/* C Header
	* @filename        : global.h
	* @author          : Frederic Masset
	* @last_modified_by: trinarypi
	* @last_modified   : 2017/06/26 16:58
	* @description     :
*/
int CPU_Rank;
int CPU_Number;
boolean CPU_Master;
int CPU_Next, CPU_Prev, CPU_Highest;
/* ------------------------------------- */
/* Variables specific to fftw mesh split */
/* ------------------------------------- */
int CPU_Friend, CPU_NoFriend;
real *dens_friend;
real *SGP_buffft_Accr_friend, *SGP_buffft_Acct_friend;
real *ffttohydro_transfer, *ffttohydro_transfer_friend;
int local_Nx, local_i_start, local_i_start_friend, total_local_size_friend, local_Nx_friend;
int local_Ny_after_transpose, local_j_start_after_transpose;
int total_local_size, ifront, Zero_or_active_friend;
int transfer_size, transfer_size_friend;
int hydro_totalsize, active_hydro_totalsize, active_hydro_totalsize_friend;  
/* ------------------------------------- */
int IMIN;
int IMAX;
int Zero_or_active;
int Max_or_active;
int One_or_active;
int MaxMO_or_active;		/* MO: Minus One */
int GLOBALNRAD;
real Rinf[MAX1D], Rsup[MAX1D], Rmed[MAX1D], Surf[MAX1D];
real InvRmed[MAX1D], InvSurf[MAX1D], InvDiffRmed[MAX1D];
real InvDiffRsup[MAX1D], InvRinf[MAX1D], Radii[MAX1D], GlobalRmed[MAX1D];
real SigmaMed[MAX1D], SigmaInf[MAX1D], MassTaper;
real EnergyMed[MAX1D], QplusMed[MAX1D], CoolingTimeMed[MAX1D];
real OmegaFrame, PhysicalTime, PhysicalTimeInitial;
int TimeStep;
real HillRadius, mdcp, mdcp0, exces_mdcp;
real GLOBAL_bufarray[MAX1D];
boolean Merge, AdvecteLabel, FakeSequential, MonitorIntegral, debug, OnlyInit;
boolean	GotoNextOutput, StoreSigma, StoreEnergy, ViscosityAlpha, RocheSmoothing;
boolean CentrifugalBalance, ExcludeHill, SloppyCFL, BinaryOn;
MPI_Status stat;
PolarGrid *CellAbscissa, *CellOrdinate;
PolarGrid *RhoStar, *RhoInt, *Potential, *Pressure, *SoundSpeed, *Temperature;
PolarGrid *DivergenceVelocity, *TAURR, *TAUPP, *TAURP, *Qplus;
boolean LogGrid;
boolean OverridesOutputdir;
char NewOutputdir[1024];
real SigmaDMT[MAX1D], EnergyDMT[MAX1D], mlinner, mlouter;
/* --------------------------------------------------- */
/* Variables required by new radiation/cooling modules */
/* --------------------------------------------------- */
real DTHETA, CV;
PolarGrid *Qminus, *DiscHeight, *RKappaval, *Qirr, *QirrRT;
PolarGrid *Darr, *Barr, *U1arr, *U2arr, *U3arr, *U4arr, *Residual, *Rfld, *Tau_cell, *Tau_grid;
OpacTable *RKappa;
IrradiationSource *IrrSources;
real QminusMed[MAX1D], MeanDiscHeight[MAX1D], omegaOpt[MAX1D], rhoJac[MAX1D];
// New Star Temperature Taper (6/01/2017)
real StarTaper;
boolean RadiationDebug, PreInitialisation;
// new - Global variables for two source ray-traced heating - (04/04/2017)
int global_real_size;
RayStruct *ray;
real *Global_sigma, *Global_rkappa, *Global_height, *Global_temperature, *Global_qrt;
real DiffRsup[MAX1D], GlobalTheta[MAX1D];

PolarGrid *QDiv, *TempSourcesSinks, *OpticalDepth, *OpticalDepthEff, *lambdafld;
PolarGrid *TempGuess, *TempGuess_old, *Density;
int FLDTimeStepsCFL;