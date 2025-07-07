#include <iostream>
#include <fstream>
#include <cmath>
#include "main.h"
#include "models.h"
#include "deform.h"
#include "output.h"
#include "solvers.h"



/********** GLOBAL VARIABLES *********************/

Domain		D;
Body		B1, B2;

MM_density			RHO;
MM_ls_viscosity		LSV;
MM_gn_viscosity		GNV;
MM_th_conductivity	THC;
MM_ht_capacity		HCP;

/*************************************************/

void herzian_parameters(Body B1, Body B2, int kval, double W, double& a, double& b, double& ph)
{
/* calculate the herzian contact parameters */

	double RX, RY, R, Estar, k, eps, f;
	double r1x, r1y, r2x, r2y, YM1, YM2, PR1, PR2;

	r1x = B1.Rx;
	r1y = B1.Ry;
	YM1 = B1.YM;
	PR1 = B1.PR;

	r2x = B2.Rx;
	r2y = B2.Ry;
	YM2 = B2.YM;
	PR2 = B2.PR;

	if (kval == 1) { /* point contact*/

		RX = pow(1 / r1x + 1 / r2x, -1);
		RY = pow(1 / r1y + 1 / r2y, -1);

		R = pow(1 / RX + 1 / RY, -1);

		Estar = pow(0.5 * ((1 - PR1 * PR1) / YM1 + (1 - PR2 * PR2) / YM2), -1);

		k = 1.0339 * pow(RY / RX, 0.636);
		eps = 1.0003 + 0.5968 / (RY / RX);
		f = 1.5277 + 0.6023 * log(RY / RX);

		a = pow((6 * k * k * eps * W * R) / (pi * Estar), 1 / 3.0);

		b = pow((6 * eps * W * R) / (pi * k * Estar), 1 / 3.0);

		ph = (3 * W) / (2 * pi * a * b);

		/* delta = f * (4.5 / (eps * R) * pow(W / (pi * k * Estar), 1 / 3.0)); */
	}
	else { /* line contact */

		/* line contact formulas here */

	}

} /* proc */

void elbody_initialize(int N, double Rx, double Ry, double YM, double PR, double K, double Cp, 
	double vel[2], double Ts, double IT1, string IT2, Body& B)
{
/* elastic body initialization */
	
	B.Num = N;
	B.Rx = Rx;
	B.Ry = Ry;
	B.YM = YM;
	B.PR = PR;
	B.YMe = YM / (1 - PR * PR);
	B.K = K;
	B.Cp = Cp;
	B.Ts = Ts;
	B.vel[0] = vel[0];
	B.vel[1] = vel[1];
	B.IT1 = IT1;
	B.IT2 = IT2;

}/* proc */

void domain_initialize(int dkey, Body B1, Body B2, int rs, double a, double b, double kx, double ky, double W, double ph, 
	double rhoR, double muR, double KR, double CpR, double T0, int pkey, string IT1, Domain& D)
{
	/* domain initialization */

	int i, j, l, Nx, Ny, Nz, Nz0, Nz1, Nz_curr;
	double dx, dy, z1_dim, z2_dim, h, a1, b1, x_curr, y_curr, p_curr;
	double xmin, xmax, ymin, ymax;
	double r1x, r1y, r2x, r2y, h0, q;


	dx = pow(1 / 2.0, rs);	/* dimensionless meshsize for symmetry part of area */
	dy = dx;				/* we use only equidistant points in 2D-grid! */

	/* dimensionless boundaries */
	a1 = kx;
	b1 = ky;

	/* number of points along x- and y-axis in XY-area */
	Nx = static_cast<int>(floor(2 * a1 / dx) + 1);
	Ny = static_cast<int>(floor(2 * b1 / dy) + 1);

	
	/* number of poitns along z axis outside refine subarea*/
	Nz0 = 5;					/* for 3/8-integtation rule (Nz-1) must be a multiple of 3 !! */				

	/* boundaries for subarea with refine z-discretization (dimensionless form) */
	xmin = -1;
	xmax = 1;
	ymin = -1;
	ymax = 1;
	
	/* number of poitns along z axis inside refine subarea */
	Nz1 = 10;					/* for 3/8-integtation rule (Nz-1) must be a multiple of 3 !! */				

	/* get bodies geometry */
	r1x = B1.Rx;
	r2x = B2.Rx;
	r1y = B1.Ry;
	r2y = B2.Ry;
	h0  = B2.IT1;


	/* write general information about domain */
	D.dkey = dkey;
	D.Nx = Nx;
	D.Ny = Ny;
	D.Nz0 = Nz0; /* for 3/8-integtation rule (Nz-1) must be a multiple of 3 !! */
	
	D.RF.x1  = xmin;
	D.RF.x2  = xmax;
	D.RF.y1  = ymin;
	D.RF.y2	 = ymax;
	D.RF.val = Nz1;

	D.dx = dx;
	D.dy = dy;
	D.a = a;
	D.b = b;
	D.kx = kx;
	D.ky = ky;
	D.W  = W;
	D.ph = ph;
	D.h0 = B2.IT1;
	D.rhoR = rhoR;
	D.muR = muR;
	D.KR = KR;
	D.CpR = CpR;
	D.T0 = T0;
	D.IT1 = IT1;

	if (dkey == 2) { /* 2D-domain (rectangular area) initialization */

		D.Area = new Point * [Nx];

		for (i = 0; i < Nx; ++i) {

			D.Area[i] = new Point[Ny];
			for (j = 0; j < Ny; ++j) {

				x_curr = -(a * kx) + i * (dx * a);
				y_curr = -(b * ky) + j * (dy * b);

				/* z-coords for surfaces of bodies in dimension form */
				z1_dim = -1 / (2 * r1x) * x_curr * x_curr - 1 / (2 * r1y) * y_curr * y_curr;
				z2_dim = 1 / (2 * r2x) * x_curr * x_curr + 1 / (2 * r2y) * y_curr * y_curr + h0;

				h = abs(z1_dim - z2_dim); /* film thickness, m */

				D.Area[i][j].z1 = z1_dim / h;
				D.Area[i][j].z2 = z2_dim / h;
				D.Area[i][j].h = h;

				/* set number of points along z axis in refined subarea */
				if ((x_curr >= xmin) && (x_curr <= xmax) && (y_curr >= ymin) && (y_curr <= ymax))
					D.Area[i][j].Nz = Nz1;
				else D.Area[i][j].Nz = Nz0;

				/* set initial pressure: uniform (q=0); herzian (q=0.5) */
				if (pkey == 0) q = 0; else q = 0.5;

				if (( (x_curr * x_curr) / (a * a) + (y_curr * y_curr) / (b * b) - 1) >= 0) D.Area[i][j].p = 0;
				else D.Area[i][j].p = pow((1 - (x_curr * x_curr) / (a * a) - (y_curr * y_curr) / (b * b)), q);
			
				/* set initial temperature, density, viscosity (uniform, dimensionless)  */
				/* create zero vectors for further storage of shear stress, velocity comp. */
				Nz_curr = D.Area[i][j].Nz;
				for (l = 0; l < Nz_curr; ++l) {
					D.Area[i][j].T.push_back(1.0);
					D.Area[i][j].Told.push_back(1.0);
					D.Area[i][j].rho.push_back(1.0);
					D.Area[i][j].mu.push_back(1.0);
					D.Area[i][j].th.push_back(1.0);
					D.Area[i][j].hc.push_back(1.0);
					D.Area[i][j].nu.push_back(0.0);
					D.Area[i][j].tau.push_back(0.0);
					D.Area[i][j].u.push_back(0.0);
					D.Area[i][j].v.push_back(0.0);
				}
			} /* j loop */
		} /* i loop */
	} /* if */
	else { /* 1-D domain (line segment) initialization */


	} /* else */
} /* proc */

void domain_finalize(Domain& D)
{
/* freeing memory from dynamic array data at the end of the program */

	int Nx;

	Nx = D.Nx;

	for (int i = 0; i < Nx; ++i) 
		delete[] D.Area[i];
			
	delete[] D.Area;

} /* proc */


int main() {

	double a, b, ph;



	/**************** INPUT DATA *********************/

	int		rs = 1;

	double	kx = 1.0;
	double	ky = 1.0;

	double	r1x = 1e-2;	/* m */
	double	r1y = 1e-2;	/* m */
	double	r2x = 1e-2;	/* m */
	double	r2y = 1e-2;	/* m */

	double	YM1 = 2e11;	/* Pa */
	double	PR1 = 0.3;
	double  K1 = 45.0;	/* W / (m K) */
	double  CP1 = 480.0;	/* J / (kg K) */
	double  Ts1 = 295.0;	/* K */

	double	vel1[2];	/* vector of surface velocity for body 1 */
	
	vel1[0] = 1.0;		/* m/s */
	vel1[1] = 0.0;		/* m/s */

	double	YM2 = 2e11;	/* Pa */
	double	PR2 = 0.3;
	double  K2 = 45.0;	/* W / (m K) */
	double  CP2 = 480.0;	/* J / (kg K) */
	double  Ts2 = 295.0;	/* K */

	double	vel2[2];	/* vector of surface velocity for body 1 */

	vel2[0] = 0.0;		/* m/s */
	vel2[1] = 0.0;		/* m/s */

	double	W = 50;		/* normal contact force, N */

	double	h0 = 1e-5;	/* any small, but non-zero value, m */

	int		pkey = 1;	/* elliptical profile for initial pressure */

	double	rhoR = 852.9; /* kg / m^3 */
	double	muR = 0.41055;	/* Pa s */
	
	double	KR = 0.12;	/* W / (m K) */
	double	CpR = 1670.0;	/* J / (kg K) */
	double	T0 = 295.0;	/* K */

	/*****************************************************/


	/**************** INITIALIZATION *********************/

	int dtype = 2; /* point contact, 2D-domain */


	/* elastic bodies initialization*/
	elbody_initialize(1, r1x, r1y, YM1, PR1, K1, CP1, vel1, Ts1, 0, "This is bottom body", B1);
	elbody_initialize(2, r2x, r2y, YM2, PR2, K2, CP2, vel2, Ts2, h0, "This is top body", B2);
		
	/* get point contact parameters (herzian model) */
	if (dtype == 2)
		herzian_parameters(B1, B2, 1, W, a, b, ph);
	else {/* line contact */ }
	
	/* initialize domain for TEHL (point contact) calculation */
	if (dtype == 2)
		domain_initialize(2, B1, B2, rs, a, b, kx, ky, W, ph, rhoR, muR, KR, CpR, T0, pkey,
		"This is 2D-domain for TEHL (point contact) calculations with herzian initial pressure field ", D);
	else {/* line contact */ }

	/*****************************************************/

	/******* MODELS FOR LUBRICANT PROPERTIES *************/

	/* set density model */
	RHO.mkey = 1;			/* Tait (80W-90 Gear oil, Bair, 2019. p.167) */
	RHO.arg1 = 852.9;		/* rhoR  */
	RHO.arg2 = 293;			/* TR	 */
	RHO.arg3 = 4.873e-3;	/* betaK */
	RHO.arg4 = 7.187e-4;	/* aV	 */
	RHO.arg5 = 10.65;		/* K0	 */
	RHO.arg6 = 7.232e9;		/* K00	 */

	/* set low-shear viscosity model */
	LSV.mkey = 1;			/* Doolitle (80W-90 Gear oil, Bair, 2019. p.167) */
	LSV.arg1 = 0.08474;		/* muR */
	LSV.arg2 = 323;			/* TR  */
	LSV.arg3 = 4.336;		/* B   */
	LSV.arg4 = 0.6988;		/* R0  */
	LSV.arg5 = -7.679e-4;	/* eps */
	LSV.arg6 = 0.0;			/* -   */
	LSV.arg7 = 0.0;			/* -   */
	LSV.arg8 = 0.0;			/* -   */

	/* set thermal conductivity model */
	THC.mkey = 1;			/* Power-law (T9 oil, Habchi, 2010. p.1845) */
	THC.arg1 = 0.053;		/* B */
	THC.arg2 = 0.026;		/* C */
	THC.arg3 = 7.6;			/* s */
	THC.arg3 = 295;			/* TR */

	/* set heat capacity model */
	HCP.mkey = 1;			/* Linear (T9 oil, Habchi, 2010. p.1845) */
	HCP.arg1 = 1.17e6;		/* C0 */
	HCP.arg2 = 0.39e6;		/* m  */
	HCP.arg3 = 0.0;			/* -  */
	HCP.arg4 = 0.0;			/* -  */
	HCP.arg5 = 0.0;			/* -  */
	HCP.arg6 = 0.0;			/* -  */
	HCP.arg7 = 295;			/* TR */

	/* set rheology model (shear stress - viscosity)  */
	GNV.mkey = 1;			/* Ree-Eyring ( ) */
	GNV.arg1 = 3e6;			/* tau_0 */
	GNV.arg2 = 0.0;			/* -  */
	GNV.arg3 = 0.0;			/* -  */
	GNV.arg4 = 0.0;			/* -  */
	GNV.arg5 = 0.0;			/* -  */
	GNV.arg6 = 0.0;			/* -  */

	/*****************************************************/
	
	/* calculation of lubricant physical properties in domain */
	
	
	calc_physical_fields(D, RHO, LSV, GNV,THC, HCP, B1.vel, B2.vel);

	/* calculation of the bodies elastic deformation in domain */
	calc_deform(D, B1, B2);




	/******************* OUTPUT RESULTS ******************/

	output_data(D, B1, B2, 0, 0, '-', 0, 0, "GENERAL_INF");
	output_data(D, B1, B2, 0, 1, '-', 0, 0, "B1-SURF_DL");
	output_data(D, B1, B2, 0, 1, '-', 0, 1, "B1-SURF_DM");
	output_data(D, B1, B2, 0, 2, '-', 0, 0, "B2-SURF_DL");
	output_data(D, B1, B2, 0, 2, '-', 0, 1, "B2-SURF_DM");

	output_data(D, B1, B2, 0, 3, '-', 0, 0, "Z-RESOL");

	output_data(D, B1, B2, 0, 4, '-', 0  , 1, "THK_DM");
	output_data(D, B1, B2, 0, 4, 'X', 0.5, 1, "THK-YZ_DM");
	output_data(D, B1, B2, 0, 4, 'Y', 0.5, 1, "THK_XZ-DM");

	output_data(D, B1, B2, 1, 1, '-', 0,   0, "DF1_DL");
	output_data(D, B1, B2, 1, 1, 'X', 0.5, 0, "DF1-YZ_DL");
	output_data(D, B1, B2, 1, 1, 'Y', 0.5, 0, "DF1_XZ-DL");
	output_data(D, B1, B2, 1, 1, '-', 0,   1, "DF1_DM");
	output_data(D, B1, B2, 1, 1, 'X', 0.5, 1, "DF1-YZ_DM");
	output_data(D, B1, B2, 1, 1, 'Y', 0.5, 1, "DF1_XZ-DM");

	output_data(D, B1, B2, 1, 2, '-', 0,   0, "DF2_DL");
	output_data(D, B1, B2, 1, 2, 'X', 0.5, 0, "DF2-YZ_DL");
	output_data(D, B1, B2, 1, 2, 'Y', 0.5, 0, "DF2_XZ-DL");
	output_data(D, B1, B2, 1, 2, '-', 0,   1, "DF2_DM");
	output_data(D, B1, B2, 1, 2, 'X', 0.5, 1, "DF2-YZ_DM");
	output_data(D, B1, B2, 1, 2, 'Y', 0.5, 1, "DF2_XZ-DM");

	output_data(D, B1, B2, 1, 3, '-', 0,   0, "PRS_DL");
	output_data(D, B1, B2, 1, 3, 'X', 0.5, 0, "PRS-YZ_DL");
	output_data(D, B1, B2, 1, 3, 'Y', 0.5, 0, "PRS-XZ_DL");
	output_data(D, B1, B2, 1, 3, '-', 0,   1, "PRS_DM");
	output_data(D, B1, B2, 1, 3, 'X', 0.5, 1, "PRS-YZ_DM");
	output_data(D, B1, B2, 1, 3, 'Y', 0.5, 1, "PRS-XZ_DM");

	output_data(D, B1, B2, 1, 4, '-', 0,   0, "RHO_DL");
	output_data(D, B1, B2, 1, 4, '-', 0,   1, "RHO_DM");

	output_data(D, B1, B2, 1, 5, '-', 0, 0, "LSV_DL");
	output_data(D, B1, B2, 1, 5, '-', 0, 1, "LSV_DM");

	output_data(D, B1, B2, 1, 6, '-', 0, 0, "TAU_DL");
	output_data(D, B1, B2, 1, 6, '-', 0, 1, "TAU_DM");

	output_data(D, B1, B2, 1, 7, '-', 0, 0, "GNV_DL");
	output_data(D, B1, B2, 1, 7, '-', 0, 1, "GNV_DM");

	/* free CPU memory */
	domain_finalize(D);

} /* end of program */









