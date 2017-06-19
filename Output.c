#include "mp.h"

static real     Xplanet, Yplanet, VXplanet, VYplanet, MplanetVirtual;
extern real     LostMass, OmegaFrame;
extern boolean  Write_Density, Write_Velocity, Write_Energy, IsDisk, Write_Eccentricity, Write_Pericentre;
extern boolean  Write_Temperature, Write_DivV, Write_Qplus, DiscElem;
extern boolean  ConfigMp, ConfigPos, BinaryOn, RadTransport;
extern boolean  Write_DiscHeight, Write_Qminus, Write_Kappa, Write_Qirr, Write_Residual, Write_Coeffs, Write_Rfld, Write_QirrRT, Write_Taus;
extern boolean  Write_OpticalDepths;

void EmptyPlanetSystemFile (sys)
     PlanetarySystem *sys;
{
  FILE *output;
  char name[256];
  int i, n;
  n = sys->nb;
  if (!CPU_Master) return;
  for (i = 0; i < n; i++) {
    sprintf (name, "%splanet%d.dat", OUTPUTDIR, i);
    output = fopen (name, "w");
    if (output == NULL) {
      fprintf (stderr, "Can't write %s file. Aborting.\n", name);
      prs_exit (1);
    }
    fclose (output);
  }
}

void WritePlanetFile (TimeStep, n)
     int TimeStep;
     int n;
{
  FILE *output;
  char name[256];
  if (!CPU_Master) return;
  printf ("Updating 'planet%d.dat'...", n);
  fflush (stdout);
  sprintf (name, "%splanet%d.dat", OUTPUTDIR, n);
  output = fopen (name, "a");
  if (output == NULL) {
    fprintf (stderr, "Can't write 'planet%d.dat' file. Aborting.\n", n);
    prs_exit (1);
  }
  fprintf (output, "%d\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\n", TimeStep, Xplanet, Yplanet, VXplanet, VYplanet, MplanetVirtual, LostMass, PhysicalTime, OmegaFrame, mdcp, exces_mdcp);
  fclose (output);
  printf ("done\n");
  fflush (stdout);
}

void WritePlanetSystemFile (sys, t)
     PlanetarySystem *sys;
     int t;
{
  int i, n;
  n = sys->nb;
  for (i = 0; i < n; i++) {
    Xplanet = sys->x[i];
    Yplanet = sys->y[i];
    VXplanet = sys->vx[i];
    VYplanet = sys->vy[i];
    MplanetVirtual = sys->mass[i];
    WritePlanetFile (t, i);
  }
}
   

void WriteBigPlanetFile (TimeStep, n)
     int TimeStep;
     int n;
{
  FILE *output;
  char name[256];
  if (!CPU_Master) return;
  fflush (stdout);
  sprintf (name, "%sbigplanet%d.dat", OUTPUTDIR, n);
  output = fopen (name, "a");
  if (output == NULL) {
    fprintf (stderr, "Can't write 'bigplanet.dat' file. Aborting.\n");
    prs_exit (1);
  }
  fprintf (output, "%d\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\n", TimeStep, Xplanet, Yplanet, VXplanet, VYplanet, MplanetVirtual, LostMass, PhysicalTime, OmegaFrame, mdcp, exces_mdcp);
  fclose (output);
  fflush(stdout);
}

void WriteBigPlanetSystemFile (sys, t)
     PlanetarySystem *sys;
     int t;
{
  int i, n;
  n = sys->nb;
  for (i = 0; i < n; i++) {
    Xplanet = sys->x[i];
    Yplanet = sys->y[i];
    VXplanet = sys->vx[i];
    VYplanet = sys->vy[i];
    MplanetVirtual = sys->mass[i];
    WriteBigPlanetFile (t, i);
  }
}

real GetfromPlanetFile (TimeStep, column, n)
     int TimeStep, column, n;
{
  FILE *input;
  char name[256];
  char testline[256];
  int time;
  char *pt;
  double value;
  sprintf (name, "%splanet%d.dat", OUTPUTDIR, n);
  input = fopen (name, "r");
  if (input == NULL) {
    mastererr ("Can't read 'planet%d.dat' file. Aborting restart.\n",n);
    prs_exit (1);
  }
  if (column < 2) {
    mastererr ("Invalid column number in 'planet%d.dat'. Aborting restart.\n",n);
    prs_exit (1);
  }
  do {
    pt = fgets (testline, 255, input);
    sscanf (testline, "%d", &time);
  } while ((time != TimeStep) && (pt != NULL));
  if (pt == NULL) {
    mastererr ("Can't read entry %d in 'planet%d.dat' file. Aborting restart.\n", TimeStep,n);
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

void RestartPlanetarySystem (timestep, sys)
     PlanetarySystem *sys;
     int timestep;
{
  int k;
  for (k = 0; k < sys->nb; k++) {
    if (ConfigPos == NO) {
      sys->x[k] = GetfromPlanetFile (timestep, 2, k);
      sys->y[k] = GetfromPlanetFile (timestep, 3, k);
      sys->vx[k] = GetfromPlanetFile (timestep, 4, k);
      sys->vy[k] = GetfromPlanetFile (timestep, 5, k);
    }
    if (ConfigMp == NO)
      sys->mass[k] = GetfromPlanetFile (timestep, 6, k);
  }
}

void WriteDiskPolar(array, number)
     PolarGrid 	*array;
     int 	number;
{
  int           Nr, Ns;
  FILE          *dump;
  char 		name[256];
  real 		*ptr;
  ptr = array->Field;
  if (CPU_Master)
    sprintf (name, "%s%s%d.dat", OUTPUTDIR, array->Name, number);
  else
    sprintf (name, "%s%s%d.dat.%05d", OUTPUTDIR, array->Name, number, CPU_Rank);
  Nr = array->Nrad;
  Ns = array->Nsec;
  dump = fopen(name, "w");
  if (dump == NULL) {
    fprintf(stderr, "Unable to open '%s'.\n", name);
    prs_exit(1);
  }
  masterprint ("Writing '%s%d.dat'...", array->Name, number);
  fflush (stdout);
  MPI_Barrier (MPI_COMM_WORLD);
/* We strip the first CPUOVERLAP rings if the current CPU is not the 
   innermost one */
  if (CPU_Rank > 0) {
    ptr += CPUOVERLAP*Ns;
    Nr -=CPUOVERLAP ;
  }
/* We strip the last CPUOVERLAP rings if the current CPU is not the outermost
   one, equal to CPU_Highest in all cases */
  if (CPU_Rank != CPU_Highest) {
    Nr -=CPUOVERLAP;
  }
  fwrite (ptr, sizeof(real), Nr*Ns,dump);
  fclose(dump);
  fprintf(stdout, "%d/", CPU_Rank);  
  fflush(stdout);
  MPI_Barrier (MPI_COMM_WORLD);
  masterprint("done\n");
}

void WriteDim () {	  
  char filename[256];
  FILE 	*dim;
  if (!CPU_Master) return;
  sprintf (filename, "%sdims.dat", OUTPUTDIR);
  if ((dim = fopen (filename, "w")) == NULL) {
    fprintf (stderr, "Unable to open %s. Program stopped\n", filename);
    prs_exit (1);
  }
  fprintf (dim,"%d\t%d\t\t%d\t%d\t%f\t%d\t%d\t%d\n",0,0,0,0,RMAX, NTOT/NINTERM, GLOBALNRAD, NSEC);
  fclose (dim);
}

void SendOutput (index, dens, gasvr, gasvt, gasenerg, label, e_cell)
     int          index;
     PolarGrid   *dens, *gasvr, *gasvt, *label, *gasenerg, *e_cell;
{
  if (CPU_Master)
    printf ("\n*** OUTPUT %d ***\n", index);
  if (IsDisk == YES) {
      if (Write_Density == YES) WriteDiskPolar (dens, index);
      if (Write_Velocity == YES) {
	       WriteDiskPolar (gasvr, index);
	       WriteDiskPolar (gasvt, index);
      }
      if (Write_Eccentricity == YES) WriteDiskPolar (e_cell, index);
      if (Write_Energy == YES) WriteDiskPolar (gasenerg, index);
      if (Write_Temperature == YES) WriteDiskPolar (Temperature, index);
      if (Write_DivV == YES) WriteDiskPolar (DivergenceVelocity, index);
      if (Write_Qplus == YES)  WriteDiskPolar (Qplus, index);
      if (AdvecteLabel == YES) WriteDiskPolar (label, index);
      if (Write_Qminus == YES) WriteDiskPolar (Qminus, index);
      if (Write_Qirr == YES) WriteDiskPolar (Qirr, index);
      if (Write_QirrRT == YES) WriteDiskPolar (QirrRT, index);
      if (Write_DiscHeight == YES) WriteDiskPolar (DiscHeight, index);
      if (Write_Kappa == YES) WriteDiskPolar (RKappaval, index);
      if (Write_Coeffs == YES) {
        WriteDiskPolar (Darr, index);
        WriteDiskPolar (Barr, index);
        WriteDiskPolar (U1arr, index);
        WriteDiskPolar (U2arr, index);
        WriteDiskPolar (U3arr, index);
        WriteDiskPolar (U4arr, index);
      }
      if (Write_Taus == YES) {
        WriteDiskPolar(Tau_cell, index);
        WriteDiskPolar(Tau_grid, index);
      }
      if (Write_Residual == YES) WriteDiskPolar (Residual, index);
      if (Write_Rfld == YES) {
        WriteDiskPolar(Rfld, index);
        WriteDiskPolar(lambdafld, index);
      }
      if (Write_OpticalDepths == YES) {
        WriteDiskPolar(OpticalDepth, index);
        WriteDiskPolar(OpticalDepthEff, index);
      }

      MPI_Barrier (MPI_COMM_WORLD);
      if (Merge && (CPU_Number > 1)) merge (index);
      if (RadTransport == YES)
        WriteRadTransInfo();
  }
}

void WriteSimVariableFile ()

{
  int lg_grid = 0, bin_on = 0;
  FILE *output;
  char name[256];
  if (!CPU_Master) return;
  if (BinaryOn == YES) bin_on = 1;
  if (LogGrid == YES) lg_grid = 1;
  fflush (stdout);
  sprintf (name, "%ssimvars.dat", OUTPUTDIR);
  output = fopen (name, "w");
  if (output == NULL) {
    fprintf (stderr, "Can't write 'simvars.dat' file. Aborting.\n");
    prs_exit (1);
  }
  fprintf (output, "%f\t%f\t%d\t%d\t%d\t%d\n", RMIN, RMAX, GLOBALNRAD, NSEC, lg_grid, bin_on);
  fclose (output);
  fflush(stdout); 
}
		

	

 
