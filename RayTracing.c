
#include "mp.h"

extern boolean BinaryOn, RayTracingHeating, ExplicitRayTracingHeating, ImplicitRayTracingHeating;
extern real StarTaper;
extern boolean	BitschSKappa, HubenySKappa, EmptyCavity, RadiationDebug;

static real *SendInnerBoundary;
static real *SendOuterBoundary;
static real *RecvInnerBoundary;
static real *RecvOuterBoundary;
static int allocated_radcomm = 0, allocated_globalfields = 0;

void InitRayTracing ()
	// Input N/A
{
	// Declaration
	int nstar;

	// Assignment
	nstar = IrrSources->nb;

	// Function
	QirrRT 		= CreatePolarGrid(NRAD, NSEC, "QirrRT");
	global_real_size = GLOBALNRAD*NSEC*sizeof(real);
	
	masterprint("Initialising parameters for RT irradiative heating...");
	if ( nstar > 1 ) {
		masterprint("two sources found...\n");
		if ( CPU_Rank == 0) {
			ray = AllocRayStructure();
			InitRayStructure(0.0, 0.0);
			PrintRayInfo();
		}
	} else {
		masterprint("one source found...");
		Tau_cell  = CreatePolarGrid(NRAD, NSEC, "Tau_cell");
		Tau_grid  = CreatePolarGrid(NRAD, NSEC, "Tau_grid");
	}
	masterprint("Done.\n\n");

}

void InitRayStructure (xstart, ystart)
	//Input
	real xstart;
	real ystart;
{
	// Declaration
	real x[2]={0.0}, y[2]={0.0};
	// real x0, y0;

	// x0 = ray->x;
	// y0 = ray->y;
	// x0 = xstart;
	// y0 = ystart;

	// Function
	ray->x = xstart;
	ray->y = ystart;
	ray->dx = 0.0;
	ray->dy = 0.0;
	ray->dr = 0.0;
	ray->r = 0.0;
	ray->th = 0.0;
	ray->diff = 0.0;
	ray->length = 0.0;
	ray->separation = 0.0;
	ray->tau = 0.0;
	ray->dtau = 0.0;
	ray->env = 0;
	ray->n_int = 0;

	memcpy(ray->x_int, x, sizeof x);
	memcpy(ray->y_int, y, sizeof y);
}


RayStruct *AllocRayStructure ()
	// Input N/A
{	
	// Declaration
	RayStruct *array;

	// Function
	array  = (RayStruct *)malloc (sizeof(RayStruct));

	if ( array == NULL ) {
		fprintf (stderr, "Not enough memory.\n");
    	prs_exit (1);
	}

	// Output
	return array;
}


void ComputeRayTracingHeating (gas_density, bsys)
	// Input
	PolarGrid *gas_density;
	BinarySystem *bsys;
{	
	// Declaration
	int nstar;

	//Assignment
	nstar = IrrSources->nb;

	//Function
	if ( nstar == 1 ) {
		ComputeSingleSourceRT(gas_density);
	} else {
		ComputeBinarySourceRT(gas_density, bsys);
	}
}

void ComputeBinarySourceRT (gas_density, bsys)
	// Input
	PolarGrid *gas_density;
	BinarySystem *bsys;
{
	// Declaration
	FILE	*dump;
  char	name[256];
	int ns, nstar, ncp, n;
	int i, im, ip, ib;
	int j, jm, jp;
	int l, lijm, lijp, lip, ljp;
	int s, s2;
	int blocked;
	int *continue_in_i = NULL;
	int first_active, last_active, active_size, ffas; /* ff = four-field */
	int *sizes = NULL, *ff_sizes = NULL, *displs = NULL, *ff_displs = NULL;
	int shift_address, start_address, ff_start_address, sizes_bytes;
	int sgn_diff_x, sgn_diff_y;
	real *sigma, *rkappa, *height, *temperature, *QRT;
	real Cavity_sigma, Cavity_rkappa, Cavity_height, Cavity_temperature;
	real Ray_sigma, Ray_rkappa, Ray_height, Ray_temperature;
	real *Send_FieldBuffer, *Recv_FieldBuffer;
	real *Tstar, *Rstar;
	real rho;
	real *xb, *yb, rb;
	real c1, dth, starcons;
	real angle, xi, yi;
	real diff_x, abs_diff_x;
	real diff_y, abs_diff_y;
	real ratio_xy, denom, dx, dy, dray;
	real r2r, rr1, r2r1, th2th, thth1;
	real skappa, term;
	real t_elapsed;
	time_t t_start, t_end;


	// Assignment
	ns 			= gas_density->Nsec;
	sigma 		= gas_density->Field;
	rkappa 		= RKappaval->Field;
	temperature = Temperature->Field;
	height 		= DiscHeight->Field;
	QRT 		= QirrRT->Field;
	nstar		= IrrSources->nb;
	Rstar 		= IrrSources->Rstar;
	Tstar 		= IrrSources->Tstar;
	xb			= bsys->x;
	yb			= bsys->y;

	/* Allocate memory for global Sigma, Kappa_R, Height, Temperature (four-fields)
	and Q_irr^RT fields */
	masterprint("AllocateGlobalFIelds...");
	if (( CPU_Rank == 0 ) && ( allocated_globalfields == 0 )) {
		AllocateGlobalFields();
	}
	masterprint("done.\n");

	// Constants
	c1 = 0.5;
	dth = 2.0*PI/(real)ns;
	first_active = (IMIN + Zero_or_active)*ns;
	last_active = (IMIN + Max_or_active)*ns - 1;
	active_size = (last_active - first_active + 1);
	ffas = 4*active_size;

	// Function

	/* Gather array and four-field array sizes onto root */
	if ( CPU_Rank == 0 ) {
		continue_in_i = (int *)malloc(sizeof(int)* ns);
		sizes 		= (int *)malloc(sizeof(int) * CPU_Number);
		ff_sizes 	= (int *)malloc(sizeof(int) * CPU_Number);
		displs 		= (int *)malloc(sizeof(int) * CPU_Number);
		ff_displs 	= (int *)malloc(sizeof(int) * CPU_Number);
	}
	masterprint("Gather field sizes to master...");
	MPI_Gather(&active_size, 1, MPI_INT, sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Gather(&ffas, 1, MPI_INT, ff_sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);
	masterprint("done.\n");

	/* Allocate memory for four-field send buffer,
	and copy active zone regions of the fields to this buffer  */
	masterprint("Create Send_FieldBuffer and copy local fields into it...");
	Send_FieldBuffer = (real *)malloc((ffas)*sizeof(real));
	memcpy(&Send_FieldBuffer[0*active_size], &sigma[Zero_or_active*ns], active_size*sizeof(real));
	memcpy(&Send_FieldBuffer[1*active_size], &rkappa[Zero_or_active*ns], active_size*sizeof(real));
	memcpy(&Send_FieldBuffer[2*active_size], &height[Zero_or_active*ns], active_size*sizeof(real));
	memcpy(&Send_FieldBuffer[3*active_size], &temperature[Zero_or_active*ns], active_size*sizeof(real));
	masterprint("done.\n");

	/* On Master CPU create receive buffer,
	for cell and grid tau buffers respectively.
	Then prepare array of start displacement values (stride length)
	to send correct array portions back to CPUs */
	masterprint("Create Recv_FieldBuffer and create displacement and stride lengths...");
	if ( CPU_Rank == 0 ) {
		Recv_FieldBuffer = (real *)malloc((4*global_real_size));
		/* Create array of stride lengths */
		ff_displs[0] = 0;
		displs[0] = 0;
		for (i = 1; i < CPU_Number; i++) {
			displs[i] = displs[i-1] + sizes[i-1];
			ff_displs[i] = ff_displs[i-1] + ff_sizes[i-1];
		}
	}
	masterprint("done.\n");

	/* Gather four-field values to buffer
	on Master CPU. Split fields from buffer into
	separate Global fields on Master CPU */
	masterprint("Gather all local fields to global field buffer on master...");
	MPI_Gatherv(Send_FieldBuffer, ffas, MPI_DOUBLE, Recv_FieldBuffer, ff_sizes, ff_displs, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	if ( CPU_Rank == 0 ) {
		for (ncp = 0; ncp < CPU_Number; ncp++) {
			shift_address = ff_sizes[ncp]/4;
			start_address = displs[ncp];
			sizes_bytes = sizes[ncp]*sizeof(real);

			ff_start_address = ff_displs[ncp];
			memcpy(&Global_sigma[start_address], &Recv_FieldBuffer[ff_start_address], sizes_bytes);

			ff_start_address += shift_address;
			memcpy(&Global_rkappa[start_address], &Recv_FieldBuffer[ff_start_address], sizes_bytes);

			ff_start_address += shift_address;
			memcpy(&Global_height[start_address], &Recv_FieldBuffer[ff_start_address], sizes_bytes);

			ff_start_address += shift_address;
			memcpy(&Global_temperature[start_address], &Recv_FieldBuffer[ff_start_address], sizes_bytes);
		}
		free(Recv_FieldBuffer);
	}
	free(Send_FieldBuffer);
	masterprint("Done.\n");

	masterprint("Starting ray-traced heating calculation.\n\n");
	if ( CPU_Rank == 0 ) {
		masterprint("Initialising heating field to 0...");
		for (i = 0; i < GLOBALNRAD; i++) {
			for (j = 0; j < GLOBALNRAD; j++) {
				l=j+i*ns;
				Global_qrt[l] = 0.0;
			}
		}
		masterprint("done.\n");
		// memset(Global_qrt, 0.0, global_real_size);
		masterprint("Setting field values of cavity material...");
		if ( EmptyCavity == NO ) {
			Cavity_sigma = Cavity_rkappa = Cavity_height = Cavity_temperature = 0.0;
			for (j = 0; j < ns; j ++) {
				Cavity_sigma += Global_sigma[j];
				Cavity_rkappa += Global_rkappa[j];
				Cavity_height += Global_height[j];
				Cavity_temperature += Global_temperature[j];
			}
			Cavity_sigma /= ns;
			Cavity_rkappa /= ns;
			Cavity_height /= ns;
			Cavity_temperature /= ns;
		} else {
			Cavity_sigma = 0.0;
			Cavity_rkappa = 0.0;
			Cavity_height = 0.0;
			Cavity_temperature = 0.0;
		}
		masterprint("done.\n");

		for (s = 0; s < nstar; s++) {
			masterprint("Setting up stellar parameters for stat %d...", s+1);
			Tstar[s] =  StarTaper*(IrrSources->Tstar[s]);
			starcons = STEFANK*pow(Tstar[s]*Tstar[s]*Rstar[s], 2.0);
			rb = hypot(xb[s], yb[s]);

			for (j = 0; j < ns; j++) {
				continue_in_i[j] = 1;
			}
			masterprint("done.\n");
			// memset(blocked, 0, sizeof(int)*NSEC*GLOBALNRAD);
			masterprint("Calculating star %d ray-traced heating...", s+1);
			for (i = 0; i < GLOBALNRAD; i++) {
				for (j = 0; j < ns; j++) {
					l = j + i*ns;
					// angle = 2.0*PI*((real)j+0.5)/(real)ns;
					angle = GlobalTheta[j];
					xi = GlobalRmed[i]*cos(angle);
					yi = GlobalRmed[i]*sin(angle);
					diff_x = xi - xb[s];
					diff_y = yi - yb[s];
					
					abs_diff_x = fabs(diff_x);
					abs_diff_y = fabs(diff_y);
					sgn_diff_x = sign2(diff_x);
					sgn_diff_y = sign2(diff_y);
					ratio_xy = abs_diff_x/abs_diff_y;
					denom = sqrt(ratio_xy*ratio_xy + 1.0);

					InitRayStructure(xb[s], yb[s]);

					ray->r = hypot(ray->x, ray->y);
					ray->th = atan3(ray->y, ray->x);
					ray->diff = hypot(diff_x, diff_y);
					ray->separation = hypot(diff_x, diff_y);

					ib = 0;
					n = 0;
					while ( rb > GlobalRmed[ib] ) {
				    	ib = ib + 1;
					}

					if ( s == 0 ) {
						s2 = 1;
					} else {
						s2 = 0;
					}

					if ( continue_in_i[j] == 0 ) {
						ray->diff = 0;
						ray->env = -1;
					} else {
						blocked = CircleLineProjection(xi, yi, xb[s2], yb[s2], Rstar[s2]);
						/* If ray intersects other star (ray->n_int == 1), ray is blocked
						   and heating in target cell is 0 from star in question */
						if ( blocked != 0 ) {
							ray->diff = 0;
							ray->env = -1;
						}
					}

					/* If not, advance ray through disc, updating
					   the optical depth along its length */
					while (( ray->diff > 0 ) || ( ray->env > 0 )) {

						/* find local minimum grid spacing */
						if ( ib == GLOBALNRAD-1 ) {
							dray = GlobalRmed[GLOBALNRAD-1] - GlobalRmed[GLOBALNRAD-2];
						} else {
							dray = RTPRECISION*(GlobalRmed[ib+1] - GlobalRmed[ib]);
						}
					
						/* use this factor to calculate x and y increments */
						if (abs_diff_x == 0) {
							dx = 0.0;
							dy = dray*sgn_diff_y;
						} else if (abs_diff_y == 0) {
							dx = dray*sgn_diff_x;
							dy = 0.0;
						} else {
							dx = dray*ratio_xy/denom;
							dy = dx/ratio_xy;
							dx = sgn_diff_x*dx;
							dy = sgn_diff_y*dy;
						}

						/* use dx and dy and advance ray in direction */
						IterateRay(dx, dy, dray, xi, yi);
						n++;
						if (( ray->length > ray->separation ) && ( ray->env != 0 )) {
							masterprint("Ray length greater than source->target cell separation. Exiting\n");
							prs_exit();
						}

						/* update tau from previous values of tau and delta tau */
						ray->tau += ray->dtau;

						if (( ray->tau > TAUCEILING ) && (GlobalRmed[i] >= rb)) {
							ray->diff = 0;
							ray->env = -1;
							continue_in_i[j] = 0;
						}

						if ( ray->env == 0 ) {
							/* ray advanced to target cell
							   relevant field values are 
							   the on-grid values */
							if ( ray->length <= 1.5*Rstar[s] ) {
								/* ray is within radius of the star */
								Ray_sigma = 0.0;
								Ray_rkappa = 0.0;
								Ray_height = 0.0;
								Ray_temperature = 0.0;
							} else {
								Ray_sigma = Global_sigma[l];
								Ray_rkappa = Global_rkappa[l];
								Ray_height = Global_temperature[l];
								Ray_temperature = Global_height[l];
							}
						} else if ( ray->env == 2 ) {
							/* ray has crossed the cavity to the grid edge */
							Ray_sigma = Cavity_sigma;
							Ray_rkappa = Cavity_rkappa;
							Ray_height = Cavity_height;
							Ray_temperature = Cavity_temperature;
						} else if (( ray->length <= 1.5*Rstar[s] ) || ( ray->env < 0 )) {
							/*ray is within radius of the star */
							Ray_sigma = 0.0;
							Ray_rkappa = 0.0;
							Ray_height = 0.0;
							Ray_temperature = 0.0;
						} else {
							/* find correct r indices, depending on ray environment */
							if ( ray->env == 3 ) {
								im = 0;
								ip = 0;
								r2r = 1;
								rr1 = 1;
							} else {
								ip = 1;
								while ( ray->r > GlobalRmed[ip] ) {
									ip++;
								}
								im = ip - 1;
								r2r1 = GlobalRmed[ip] - GlobalRmed[im];
								r2r = (GlobalRmed[ip] - ray->r)/r2r1;
								rr1 = (ray->r - GlobalRmed[im])/r2r1;
							}
							ib = ip;

							jm = (int)((ray->th)/dth);
							jp = jm + 1;
							th2th = (GlobalTheta[jp] - ray->th)/dth;
							thth1 = (ray->th - GlobalTheta[jm])/dth;
							if ( jp == NSEC ) {
								jp = 0;
							}

							lijm = jm+im*ns;
							lip = jm+ip*ns;
							ljp = jp+im*ns;
							lijp = jp+ip*ns;

							Ray_sigma = th2th*(r2r*Global_sigma[lijm] + rr1*Global_sigma[lip]) + thth1*(r2r*Global_sigma[ljp] + rr1*Global_sigma[lijp]);
							Ray_rkappa = th2th*(r2r*Global_rkappa[lijm] + rr1*Global_rkappa[lip]) + thth1*(r2r*Global_rkappa[ljp] + rr1*Global_rkappa[lijp]);
							Ray_height = th2th*(r2r*Global_height[lijm] + rr1*Global_height[lip]) + thth1*(r2r*Global_height[ljp] + rr1*Global_height[lijp]);
							Ray_temperature = th2th*(r2r*Global_temperature[lijm] + rr1*Global_temperature[lip]) + thth1*(r2r*Global_temperature[ljp] + rr1*Global_temperature[lijp]);

							if (( Ray_sigma <= 0 ) || ( Ray_height <= 0 ) || ( Ray_temperature <= 0 ) || ( Ray_rkappa <= 0 )) {
								if (( i != 0 ) && (EmptyCavity != 0)){
									masterprint("Negative/zero interpolated values @ i = %d, j = %d, s = %d\n", i, j, s);
									prs_exit(1);
								}
							}
						}

						if (( Ray_sigma == 0 ) || ( Ray_height == 0 ) || ( Ray_temperature == 0 ) || ( Ray_rkappa == 0 )) {
							ray->dtau = 0.0;
						} else {
							rho = c1*Ray_sigma/Ray_height;
							skappa = ComputeSkappa(Ray_sigma, Ray_rkappa, Ray_temperature, Tstar[s]);
							ray->dtau = rho*skappa*ray->dr;
						}

						if (ray->env > 1) {
							ray->env = ray->env - 1;
						}
					}

					term = ComputeQRT(starcons, ray->length, ray->tau, ray->dtau, ray->dr, ray->env);
          			Global_qrt[l] += term;
          			Global_qrt[l] /= sigma[l]*CV;
				}
			}
			masterprint("done.\n");
		}
		time(&t_end);
		t_elapsed = difftime(t_end, t_start);
		// printf("Two-source ray-traced irradiative heating calculation took %f seconds\n", t_elapsed);

		Send_FieldBuffer = (real *)malloc(global_real_size);
		memcpy(Send_FieldBuffer, Global_qrt, global_real_size);
	}
	masterprint("Finished calculating ray-traced heating.\n");

	/* Send grid tau values to Slaves, copy buffer to 
		 grid tau array */
	masterprint("Create Recv_FieldBuffer for heating field and scatter global field on master to local cpu fields...\n");
	Recv_FieldBuffer = (real *)malloc((active_size)*sizeof(real));
	MPI_Scatterv(Send_FieldBuffer, sizes, displs, MPI_DOUBLE, Recv_FieldBuffer, active_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	memcpy(&QRT[Zero_or_active*ns], Recv_FieldBuffer, active_size*sizeof(real));
	if ( CPU_Rank == 0 ) {
		free(Send_FieldBuffer);
	}
	masterprint("Done.\n");
	masterprint("Free allocated memory...\n");
	free(Recv_FieldBuffer);
	free(continue_in_i);
	free(sizes);
	free(displs);
	free(ff_displs);
	free(ff_sizes);
	masterprint("done.\n");
	masterprint("Copy field overlap zones to neighbouring processors...");
	CommunicateFieldBoundaries(QirrRT);
	masterprint("done.\n");

	// Debug
	if ( RadiationDebug ) {
		int check_neg = 1;
    	int check_zero = 0;
		CheckField(QirrRT, check_neg, check_zero, "ComputeBinarySourceRT");
	}
}

int CircleLineProjection(xi, yi, xc, yc, Rc)
	// Input
	real xi;
	real yi;
	real xc;
	real yc;
	real Rc;
{
	// Declaration
	real seg_vx, seg_vy, seg_v;
	real pt_vx, pt_vy;
	real proj_vx, proj_vy, proj_v;
	real closest_x, closest_y;
	real dist_vx, dist_vy, dist_v;
	real x, y;

	// Assignment
	x = ray->x;
	y = ray->y;

	// Function
	seg_vx = xi - x;
	seg_vy = yi - y;
	seg_v = hypot(seg_vx, seg_vy);

	pt_vx = xc - x;
	pt_vy = yc - y;

	proj_v = (pt_vx*seg_vx/seg_v) + (pt_vy*seg_vy/seg_v);
	proj_vx = proj_v*seg_vx/seg_v;
	proj_vy = proj_v*seg_vy/seg_v;

	if ( proj_v <= 0 ) {
		closest_x = x;
		closest_y = y;
	} else if ( proj_v > seg_v ) {
		closest_x = xi;
		closest_y = yi;
	} else {
		closest_x = x + proj_vx;
		closest_y = y + proj_vy;
	}

	dist_vx = xc - closest_x;
	dist_vy = yc - closest_y;
	dist_v = hypot(dist_vx, dist_vy);

	if ( dist_v <= Rc ) {
		// Output
		return 1; /* Line Segment between start of ray and target cell intersects radius of other star */
	} else {
		// Output
		return 0; /* Line Segment between start of ray and target cell doesn't intersect radius of other star */
	}
}

void IterateRay(dx, dy, dray, xi, yi)
	// Input
	real dx;
	real dy;
	real dray;
	real xi;
	real yi;
{
	// Declaration
	int index;
	real xtemp, ytemp, rtemp;

	if ( ray->diff <= 1.5*dray ) {
		ray->dx = xi - ray->x;
		ray->dy = yi - ray->y;
		ray->dr = hypot(ray->dx, ray->dy);
		ray->diff = 0.0;
		ray->length += ray->dr;
		ray->x = xi;
		ray->y = yi;
		ray->r = hypot(xi, yi);
		ray->th = atan3(yi, xi);
		ray->env = 0.0;
		return;
	}

	if (( ray->r < RMIN ) || ( ray->env == 2 )) {
		index = 1;
		if ( ray->env != 2 ) {
			CircleLineIntersection(xi, yi, RMIN);
			index = 0;
		}
		ray->dx = ray->x_int[index] - ray->x;
		ray->dy = ray->y_int[index] - ray->y;
		ray->dr = hypot(ray->dx, ray->dy);
		ray->length += ray->dr;
		ray->x = ray->x_int[index];
		ray->y = ray->y_int[index];
		ray->diff = hypot((xi - ray->x), (yi - ray->y));
		ray->r = RMIN;
		ray->th = atan3(ray->y, ray->x);
		ray->env = 2;
		return;
	}

	xtemp = ray->x + dx;
	ytemp = ray->y + dy;
	rtemp = hypot(xtemp, ytemp);
	
	if ( rtemp < RMIN ) {
		CircleLineIntersection(xi, yi, RMIN);
		ray->dx = ray->x_int[0] - ray->x;
		ray->dy = ray->y_int[0] - ray->y;
		ray->dr = hypot(ray->dx, ray->dy);
		ray->length += ray->dr;
		ray->x = ray->x_int[0];
		ray->y = ray->y_int[0];
		ray->diff = hypot((xi - ray->x), (yi - ray->y));
		ray->r = RMIN;
		ray->th = atan3(ray->y, ray->x);
		ray->env = 3;

		return;
	}

	ray->x = xtemp;
	ray->y = ytemp;
	ray->dx = dx;
	ray->dy = dy;
	ray->dr = hypot(ray->dx, ray->dy);
	ray->r = rtemp;
	ray->th = atan3(ray->y, ray->x);
	ray->diff = hypot((xi - ray->x), (yi - ray->y));
	ray->length += ray->dr;
	ray->env = 1;

	return;
}

void CircleLineIntersection(xi, yi, Rc)
	// Input
	real xi;
	real yi;
	real Rc;
{
	// Declaration
	real dx, dy, dr2, D, D2;
	real discriminant;
	real abs_dy;
	real xint[2], yint[2];
	int sgn_dx, sgn_dy;

	// Function
	dx = xi - ray->x;
	dy = yi - ray->y;
	dr2 = dx*dx + dy*dy;
	D = ray->x*yi - xi*ray->y;
	D2 = D*D;
	discriminant = Rc*Rc*dr2 - D2;

	memset(ray->x_int, 0.0, 2*sizeof(real));
	memset(ray->y_int, 0.0, 2*sizeof(real));
	ray->n_int = 0;
	
	if (discriminant > 0) {
		sgn_dx = sign2(dx);
		sgn_dy = sign2(dy);
		abs_dy = fabs(dy);

		discriminant = sqrt(discriminant);

		xint[0] = (D*dy + sgn_dy*dx*discriminant)/dr2;
		yint[0] = (-D*dx + abs_dy*discriminant)/dr2;
		xint[1] = (D*dy - sgn_dy*dx*discriminant)/dr2;
		yint[1] = (-D*dx - abs_dy*discriminant)/dr2;

		if (ray->r < Rc) {
			if ( sgn_dx == sign2(xint[0]-ray->x) ) {
				ray->x_int[0] = xint[0];
			} else {
				ray->x_int[0] = xint[1];
			}

			if ( sgn_dy == sign2(yint[0]-ray->y) ) {
				ray->y_int[0] = yint[0];
			} else {
				ray->y_int[0] = yint[1];
			}
			ray->n_int = 1;
		} else {
			if ( dx == 0 ) {
				// ray->x_int = xint;
				memcpy(ray->x_int, xint, 2*sizeof(real));
			} else {
				if ( fabs(xint[0] - ray->x) <= fabs(xint[1] - ray->x) ) {
					// ray->x_int = xint;
					memcpy(ray->x_int, xint, 2*sizeof(real));
				} else {
					ray->x_int[1] = xint[0];
					ray->x_int[0] = xint[1];
				}
			}
			if (dy == 0) {
				// ray->y_int = yint;
				memcpy(ray->y_int, yint, 2*sizeof(real));
			} else {
				if ( fabs(yint[0] - ray->y) <= fabs(yint[1] - ray->y) ) {
					// ray->y_int = yint;
					memcpy(ray->y_int, yint, 2*sizeof(real));
				} else {
					ray->y_int[1] = yint[0];
					ray->y_int[0] = yint[1];
				}
			}
			ray->n_int = 2;
		}
	}

	// Output
	return;
}

void ComputeSingleSourceRT (gas_density)
	// Input
	PolarGrid *gas_density;
{
	// Declaration
	int nr, ns, i, j, l, lim, first_active, last_active, active_size;
	int *sizes = NULL, *displs = NULL;
	real *cellTau, *gridTau, *QRT, *H, *Sigma, *Rkappa, *T;
	real rho, Skappa;
	real Tstar, Rs, Ts, Fs, dr;
	real *Recv_cellTauBuffer = NULL, *Recv_gridTauBuffer;
	real *Send_gridTauBuffer = NULL, *Send_cellTauBuffer;
	real c1;

	// Assignment
	nr 			= gas_density->Nrad;
	ns 			= gas_density->Nsec;
	Sigma 		= gas_density->Field;
	Rkappa 		= RKappaval->Field;
	T 			= Temperature->Field;
	H 			= DiscHeight->Field;
	cellTau 	= Tau_cell->Field;
	gridTau 	= Tau_grid->Field;
	QRT 		= QirrRT->Field;
	Rs 			= IrrSources->Rstar[0];
	Tstar 		= IrrSources->Tstar[0];

	// Constants
	Ts = Tstar*StarTaper;
	Fs = STEFANK*pow(Ts*Ts*Rs, 2.0);
	c1 = 0.5;
	first_active = (IMIN + Zero_or_active)*ns;
	last_active = (IMIN + Max_or_active)*ns - 1;
	active_size = (last_active - first_active + 1);

	// Function
	/* Calculation of tau in each individual cell */
 #pragma omp parallel for private (dr, j, l, Skappa, tau, taueff, rho)
	for (i = 0; i < nr; i++) {
		dr = Rsup[i] - Rinf[i];
		for (j = 0; j < ns; j++) {
			l = j+i*ns;
			/* Choose whether to have a constant factor between \kappa_R and \kappa_
				 , say 10 or 1/10 (inline with Bitsch 2014), or use the prescription for \kappa_
				 used in Dobbs-Dixon 2011, where it depends on a relationship between T and /kappa_R */
			Skappa = ComputeSkappa(Sigma[l], Rkappa[l], T[l], Tstar);
			/* Commented out
			if ( NoCFL )
				H = Rmed[i]*ASPECTRATIO;
			else
				H = DiscHeight->Field[l]; */
			rho = c1*Sigma[l]/H[l];
			cellTau[l] = Skappa*rho*dr;
		}
	}

	/* Gather array sizes onto root */
	if ( CPU_Rank == 0 ) {
		sizes  = (int *)malloc(sizeof(int) * CPU_Number);
		displs = (int *)malloc(sizeof(int) * CPU_Number);
	}
	MPI_Gather(&active_size, 1, MPI_INT, sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);

	/* Trim local arrays ready for sending and 
	   prime receive buffer with cell tau values. */
	Send_cellTauBuffer = (real *)malloc((active_size)*sizeof(real));
	Recv_gridTauBuffer = (real *)malloc((active_size)*sizeof(real));
	memcpy(Send_cellTauBuffer, &cellTau[Zero_or_active*ns], active_size*sizeof(real));

	/* On Master CPU create receive and send buffers,
		 for cell and grid tau buffers respectively.
		 Then prepare array of start displacement values (stride length)
		 to send correct array portions back to CPUs */
	if ( CPU_Rank == 0 ) {
		Recv_cellTauBuffer = (real *)malloc(global_real_size);
		Send_gridTauBuffer = (real *)malloc(global_real_size);
		// Create array of stride lengths
		displs[0] = 0;
		for (i = 1; i < CPU_Number; i++)
			displs[i] = displs[i-1] + sizes[i-1];
	}

	/* Gather cell tau values to Master, and then
	   free the memory allocated to send buffers */
	
	MPI_Gatherv(Send_cellTauBuffer, active_size, MPI_DOUBLE, Recv_cellTauBuffer, sizes, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD );

  /* Calculate grid tau values of entire disc on Master CPU */

	if ( CPU_Rank == 0 ) {
		i = 0;
		for (j = 0; j < ns; j++) {
			l = j+i*ns;
			Send_gridTauBuffer[j] = 0.0;
		}

		for (i = 1; i < GLOBALNRAD; i++) {
			for (j = 0; j < ns; j++) {
				l = j+i*ns;
				lim = l-ns;
				Send_gridTauBuffer[l] = Send_gridTauBuffer[lim] + Recv_cellTauBuffer[lim];
			}
		}
	}

	/* Send grid tau values to Slaves, copy buffer to 
		 grid tau array */
	MPI_Scatterv(Send_gridTauBuffer, sizes, displs, MPI_DOUBLE, Recv_gridTauBuffer, active_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	memcpy(&gridTau[Zero_or_active*ns], Recv_gridTauBuffer, active_size*sizeof(real));

	CommunicateFieldBoundaries(Tau_grid);

	/* Calculate value of Q_{irr}^{RT} in each individual cell */
 #pragma omp parallel for private (j, l)
	for (i = 0; i < nr; i++) {
		dr = Rsup[i] - Rinf[i];
		for (j = 0; j < ns; j++) {
			l = j+i*ns;
			QRT[l] = ComputeQRT(Fs, Rmed[i], gridTau[l], cellTau[l], dr, 0);
		}
	}

	/* Free Send and Receive buffers
		 and stride length arrays */
	free(Send_cellTauBuffer);
	free(Recv_cellTauBuffer);
	free(Recv_gridTauBuffer);
	free(Send_gridTauBuffer);
	free(sizes);
	free(displs);

	// Debug
 	if ( RadiationDebug ) {
 		int check_neg = 1;
    	int check_zero = 0;
 		CheckField(Tau_cell, check_neg, check_zero, "ComputeSingleSourceRT");
 		CheckField(Tau_grid, check_neg, check_zero, "ComputeSingleSourceRT");
 		CheckField(QirrRT, check_neg, check_zero, "ComputeSingleSourceRT");
 	}
}

void AllocateGlobalFields ()
	// Input N/A
{
	// Function
	Global_sigma = (real *)malloc (global_real_size);
	Global_rkappa = (real *)malloc (global_real_size);
	Global_height = (real *)malloc (global_real_size);
	Global_temperature = (real *)malloc (global_real_size);
	Global_qrt = (real *)malloc (global_real_size);

	if (( Global_sigma       == NULL )||\
	  	( Global_rkappa      == NULL )||\
	  	( Global_height      == NULL )||\
	  	( Global_temperature == NULL )||\
	  	( Global_qrt         == NULL )) {
		fprintf (stderr, "CPU %d didn't have enough memory to allocate global fields.\n", CPU_Rank);
		prs_exit(0);
	}
	allocated_globalfields = 1;
}

void AllocateRadComm ()
	// Input N/A
{
	// Declaration
	real size_com;

	// Constant
	size_com = NSEC * CPUOVERLAP;

	// Function
	SendInnerBoundary = malloc (size_com * sizeof(real));
	SendOuterBoundary = malloc (size_com * sizeof(real));
	RecvInnerBoundary = malloc (size_com * sizeof(real));
	RecvOuterBoundary = malloc (size_com * sizeof(real));

	if (( SendInnerBoundary == NULL )||\
	  	( SendOuterBoundary == NULL )||\
	  	( RecvInnerBoundary == NULL )||\
	  	( RecvOuterBoundary == NULL )) {
		fprintf (stderr, "CPU %d didn't have enough memory to allocate communicators.\n", CPU_Rank);
		prs_exit(0);
	}
	allocated_radcomm = 1;
}

void CommunicateFieldBoundaries (field)
	// Input
  PolarGrid *field;
{	
	// Declaration
	MPI_Request req1, req2, req3, req4;
	int l, oo, o, nr;
	real size_com;

	// Assignment
	nr = field->Nrad;

	//Constants
	size_com = NSEC * CPUOVERLAP;
	l = CPUOVERLAP*NSEC;
	oo = (nr-CPUOVERLAP)*NSEC;
	o = (nr-2*CPUOVERLAP)*NSEC;

	// Function
	if (allocated_radcomm == NO)
		AllocateRadComm ();

	memcpy (SendInnerBoundary, field->Field+l, l*sizeof(real));
	memcpy (SendOuterBoundary, field->Field+o, l*sizeof(real));
	if ( CPU_Rank%2 == 0 ) {
	if ( CPU_Rank > 0 ) {
	  MPI_Isend (SendInnerBoundary, size_com, MPI_DOUBLE, CPU_Prev, 0, MPI_COMM_WORLD, &req1);
	  MPI_Irecv (RecvInnerBoundary, size_com, MPI_DOUBLE, CPU_Prev, 0, MPI_COMM_WORLD, &req2);
	}
	if ( CPU_Rank != CPU_Highest ) {
	  MPI_Isend (SendOuterBoundary, size_com, MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD, &req3);
	  MPI_Irecv (RecvOuterBoundary, size_com, MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD, &req4);
	}
	} else {
	if ( CPU_Rank != CPU_Highest ) {
	  MPI_Irecv (RecvOuterBoundary, size_com, MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD, &req3);
	  MPI_Isend (SendOuterBoundary, size_com, MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD, &req4);
	}
	if ( CPU_Rank > 0 ) {
	  MPI_Irecv (RecvInnerBoundary, size_com, MPI_DOUBLE, CPU_Prev, 0, MPI_COMM_WORLD, &req1);
	  MPI_Isend (SendInnerBoundary, size_com, MPI_DOUBLE, CPU_Prev, 0, MPI_COMM_WORLD, &req2);
	}
	}
	if ( CPU_Rank > 0 ) {
	MPI_Wait (&req1, &stat);
	MPI_Wait (&req2, &stat);
	memcpy (field->Field, RecvInnerBoundary, l*sizeof(real));
	}
	if ( CPU_Rank != CPU_Highest ) {
	MPI_Wait (&req3, &stat);
	MPI_Wait (&req4, &stat);
	memcpy (field->Field+oo, RecvOuterBoundary, l*sizeof(real));
	}
}

real atan3(y, x)
	// Input
	real y;
	real x;
{	
	// Declaration
	real angle, absx, absy;
	if ((x == 0.0) && (y == 0.0)) {
		angle = 0.0;
		return angle;
	}

	// Function
	if (x < 0.0) {
		absx = -x;
	} else {
		absx = x;
	}

	if (y < 0.0){
		absy = -y;
	} else {
		absy = y;
	}

	if (absy - absx == absy) {
		if (y < 0.0){
			angle = PI;
		} else {
			angle = 0.0;
		}
		return angle;
	}

	if (absx - absy == absx){
		angle = 0.0;
	} else {
		angle = atan(y/x);
	}

	if (x < 0.0) {
		angle = angle + PI;
	} else {
		if (y < 0.0) {
			angle = angle + 2.0*PI;
		}
	}

	// Output
	return angle;
}

int sign2(val)
	// Input
	real val;
{
	// Declaration
	int sign;

	// Function
	if (val >= 0) {
		sign = 1;
	} else {
		sign = -1;
	}

	// Output
	return sign;
}

real ComputeSkappa(S, Rk, T, Tstar)
	// Input
	real S;
	real Rk;
	real T;
	real Tstar;
{
	// Declaration
	real tau, taueff, Teff, W, RtoS_k, c1;
	real Sk;

	// Constants
	W = 2.2E-3;
	RtoS_k = 10.0;
	c1 = 0.5;

	// Function
	if ( BitschSKappa ) {
		Sk = RtoS_k*Rk;
	} else if ( HubenySKappa ) {
		tau = c1*Rk*S;
		taueff = (0.375*tau) + (0.25/tau) + 0.866;
		tau += (1.0/sqrt(3.0));
		Teff = pow(T, 4.0)/taueff;
		Sk = (T*Rk)/((3.0*Teff*tau/4.0) + (W*pow(Tstar, 4.0)));
	} else {
		Sk = (Rk/W)*pow(T/Tstar, 4.0);
	}

	// Output
	return Sk;
}

real ComputeQRT(starcons, radius, tau, dtau, dr, flag)
	// Input
	real starcons;
	real radius;
	real tau;
	real dtau;
	real dr;
	int flag;

{
	// Declaration
	real qrt;

	// Function
	if (( flag < 0 ) || ( tau > TAUCEILING )) {
		qrt = 0.0;
	} else if (dtau < 1) {
		qrt = starcons*pow(radius, -2.0)*exp((-1.0)*tau)*dtau/dr;
	} else {
		qrt = starcons*pow(radius, -2.0)*exp((-1.0)*tau)*((1.0 - exp((-1.0)*dtau))/dr);
	}

	// Output
	return qrt;
}

void PrintRayInfo ()
	// Input N/A
{
	masterprint("############################################################################\n");
	masterprint("# ray.x = %f, ray.y = %f, ray.r = %f, ray.th = %f.\n", ray->x, ray->y, ray->r, ray->th);
	masterprint("# ray.dx = %f, ray.dy = %f, ray.dr = %f.\n", ray->dx, ray->dy, ray->dr);
	masterprint("# ray.diff = %f, ray.length = %f, ray.sep = %f, ray.env = %d.\n", ray->diff, ray->length, ray->separation, ray->env);
	masterprint("# ray.x_int[0] = %f, ray.x_int[1] = %f.\n", ray->x_int[0], ray->x_int[1]);
	masterprint("# ray.y_int[0] = %f, ray.y_int[1] = %f.\n", ray->y_int[0], ray->y_int[1]);
	masterprint("# ray.n_int = %d\n", ray->n_int);
	masterprint("############################################################################\n");
}

