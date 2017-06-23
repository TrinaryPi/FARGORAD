#include "mp.h"

void CheckField(testField, checkNegative, checkZero, note_string)
	// Input
	PolarGrid *testField;
	int checkNegative;
	int	checkZero;
	char *note_string;
{
	// Declaration
	int		nr, ns, i, j, l;
	int		flagNonFinite, flagNegative, flagZero, flag1=0, flag2=0, flag3=0;
	real	*fieldvals;

	// Assignment
	nr = testField->Nrad;
	ns = testField->Nsec;
	fieldvals = testField->Field;

	// Function
	for (i = 0; i < nr; i++) {
		for (j = 0; j < ns; j++) {
			l = j+i*ns;
			flag1 = (isfinite(fieldvals[l]) != 1 ? 1 : flag1);
   		if ( checkNegative == 1 ) {
   			flag2 = (fieldvals[l] < 0.0 ? 1 : flag2);
			}
			if ( checkZero == 1 ) {
   			flag3 = (fieldvals[l] == 0.0 ? 1 : flag3);
			}
		}
	}
	// Output
	MPI_Allreduce (&flag1, &flagNonFinite, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	if ( flagNonFinite != 0 ) {
    masterprint("Error: Non-finite value in %s (%s). Exiting.\n", testField->Name, note_string);
    if ( RadiationDebug ) {
			DumpRadiationFields(testField);
		}
   	MPI_Finalize();
		exit(flagNonFinite);
	}
  if ( checkNegative == 1 ) {
		MPI_Allreduce (&flag2, &flagNegative, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
		if ( flagNegative != 0 ) {
			masterprint("Error: Negative value in %s (%s). Exiting.\n", testField->Name, note_string);
			if ( RadiationDebug ) {
				DumpRadiationFields(testField);
			}
			MPI_Finalize();
			exit(flagNegative);
		}
	}
	if ( checkZero == 1 ) {
		MPI_Allreduce (&flag3, &flagZero, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
		if ( flagZero != 0 ) {
			masterprint("Error: Zero value in %s (%s). Exiting.\n", testField->Name, note_string);
			if ( RadiationDebug ) {
				DumpRadiationFields(testField);
			}
			MPI_Finalize();
			exit(flagZero);
		}
	}
}

int CheckValue(Value, checkNegative, checkZero)
	// Input
	real Value;
	int checkNegative;
	int	checkZero;
{
	// Declaration
	int flag=0, tmp=0;

	// Function (and Output)
	tmp = isfinite(Value);
	if ( tmp != 1 ) {
		flag = 1;
	}
	if ( checkNegative == 1 ) {
		if (Value < 0.0) {
			flag = 1;
		}
	}
	if (checkZero == 1) {
		if (Value == 0.0) {
			flag = 1;
		}
	}
	return flag;
}

void PrintBooleanUsage()
	// Input N/A
{
	// Declaration
	extern boolean RadCooling, Irradiation, RadTransport, RayTracingHeating, ExplicitRayTracingHeating, Cooling, CustomCooling;
	extern boolean VarDiscHeight, BinaryOn, Adiabatic, DiscMassTaper, TempInit;
	extern boolean SelfGravity, SGZeroMode, ZMPlus;
	extern boolean FastTransport, IsDisk, HydroOn, NoCFL, RadiativeOnly;
	extern boolean RadiationDebug, ViscosityAlpha;

	// Function (and Output)
	if ( CPU_Rank == 0 ) {
		printf("Adiabatic Equation of State     : %s\n", Adiabatic ? "YES" : "NO");
		printf("  Cooling                       : %s\n", Cooling ? "YES" : "NO");
		printf("  Custom Cooling                : %s\n", CustomCooling ? "YES" : "NO");
		printf("  Radiative Cooling             : %s\n", RadCooling ? "YES" : "NO");
		printf("  Irradiation (surface)         : %s\n", Irradiation ? "YES" : "NO");
		printf("  Radiation Transport (fld)     : %s\n", RadTransport ? "YES" : "NO");
		printf("  Ray Tracing Heating           : %s\n", RayTracingHeating ? "YES" : "NO");
		printf("  Ray Tracing Heating (explicit): %s\n", ExplicitRayTracingHeating ? "YES" : "NO");
		printf("  Radiation Debug               : %s\n", RadiationDebug ? "YES" : "NO");
		printf("Variable Disc Height            : %s\n", VarDiscHeight ? "YES" : "NO");
		printf("Temperature Initialisation      : %s\n", TempInit ? "YES" : "NO");
		printf("Binary On                       : %s\n", BinaryOn ? "YES" : "NO");
		printf("Simulate Disc                   : %s\n", IsDisk ? "YES" : "NO");
		printf("HD processes on                 : %s\n", HydroOn ? "YES" : "NO");
		printf("Radiative processes only        : %s\n", RadiativeOnly ? "YES" : "NO");
		printf("No CFL criterion for dt         : %s\n", NoCFL ? "YES" : "NO");
		printf("Fast Transport Algorithm        : %s\n", FastTransport ? "YES" : "NO");
		printf("Alpha Viscosity                 : %s\n", ViscosityAlpha ? "YES" : "NO");
		printf("Disc mass taper                 : %s\n", DiscMassTaper ? "YES" : "NO");
		printf("Self Gravity                    : %s\n", SelfGravity ? "YES" : "NO");
		printf("  Self Gravity - Zero Mode      : %s\n", SGZeroMode ? "YES" : "NO");
		printf("  Self Gravity ZMPlus           : %s\n", ZMPlus ? "YES" : "NO");
	}
}

void DumpRadiationFields(Field)
	PolarGrid *Field;
{
	// Declaration
	extern boolean RadCooling, Irradiation, RayTracingHeating, RadTransport, VarDiscHeight, PreInitialisation;

	// Function (and Output)
	if ( PreInitialisation == NO ) {
		if ( RadCooling ) {
		WriteDiskPolar(Qminus, 9999);
		}
		if ( Irradiation ) {
			WriteDiskPolar(Qirr, 9999);
		}
		if ( RayTracingHeating ) {
			WriteDiskPolar(QirrRT, 9999);
		}
		if ( RadTransport ) {
			WriteDiskPolar(Rfld, 9999);
			WriteDiskPolar(lambdafld, 9999);
			WriteDiskPolar(Darr, 9999);
			WriteDiskPolar(Barr, 9999);
			WriteDiskPolar(U1arr, 9999);
			WriteDiskPolar(U2arr, 9999);
			WriteDiskPolar(U3arr, 9999);
			WriteDiskPolar(U4arr, 9999);
			WriteDiskPolar(TempGuess, 9999);
		}
		if (( Irradiation ) || ( RadCooling )) {
			WriteDiskPolar(OpticalDepth, 9999);
			WriteDiskPolar(OpticalDepthEff, 9999);
		}
		if ( VarDiscHeight ) {
			WriteDiskPolar(DiscHeight, 9999);
		}
		WriteDiskPolar(RKappaval, 9999);
	}
	WriteDiskPolar(SigmaGlobal, 9999);
	WriteDiskPolar(Temperature, 9999);
	WriteDiskPolar(Field, 9999);
	WriteDiskPolar(Qplus, 9999);
	WriteDiskPolar(QDiv, 9999);
}
