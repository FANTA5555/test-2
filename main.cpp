#include <iostream>
#include <fstream>
#include <cmath>
#include <format>
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

void herzian_parameters(Body B1, Body B2, int kval, double W, double len, double& a, double& b, double& ph, double& delta)
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
	
	RX = pow(1 / r1x + 1 / r2x, -1);
	RY = pow(1 / r1y + 1 / r2y, -1);
	R = pow(1 / RX + 1 / RY, -1);

	/* equvivalent (for reduced case) YM */
	Estar = pow(0.5 * ((1 - PR1 * PR1) / YM1 + (1 - PR2 * PR2) / YM2), -1);
	
	if (kval == 1) { /* point contact */

		k	= 1.0339 * pow(RY / RX, 0.636);
		eps = 1.0003 + 0.5968 / (RY / RX);
		f   = 1.5277 + 0.6023 * log(RY / RX);

		// Harris, 1991

		a = pow((6.0 * k * k * eps * W * R) / (pi * Estar), 1 / 3.0);

		b = pow((6.0 * eps * W * R) / (pi * k * Estar), 1 / 3.0);

		ph = (3 * W) / (2 * pi * a * b);

		delta = f * pow(4.5 / (eps * R) * pow(W / (pi * k * Estar),2), (1 / 3.0));
	}
	else { /* line contact */

		a = pow(8 * W * RX / (len * pi * Estar), 0.5);
		
		ph = (2 * W) / (len * pi * a); 

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

void domain_initialize(int dkey, Body B1, Body B2, int rsx, int rsy, double a, double b, double kx, double ky, double W, double ph, double delta,
	double rhoR, double muR, double KR, double CpR, double T0, int pkey, string IT1, Domain& D)
{
	/* domain initialization */

	int i, j, l, Nx, Ny, Nz, Nz0, Nz1, Nz_curr, nnodes;
	double dx, dy, h, a1, b1, x_curr, y_curr, p_curr;
	double xmin, xmax, ymin, ymax;
	double r1x, r1y, r2x, r2y, q, RX, RY;
	double Estar, alpha, um, vm, k, x, y, z, par;
	
	/* dimensionless meshsize for symmetry part of area */
	dx = pow(1 / 2.0, rsx);		
	dy = pow(1 / 2.0, rsy);	

	/* dimensionless area boundaries (scale factors) */
	a1 = kx;
	b1 = ky;

	/* number of points along x- and y-axis in XY-area */
	Nx = static_cast<int>(floor(2 * a1 / dx) + 1);
	Ny = static_cast<int>(floor(2 * b1 / dy) + 1);

	/* number of poitns along z axis outside refine subarea */
	Nz0 = 5; // for 3/8-integtation rule (Nz-1) must be a multiple of 3 !!				

	/* boundaries for subarea with refine z-discretization (dimensionless form) */
	xmin = -0.25;
	xmax =  0.25;
	ymin = -0.25;
	ymax =  0.25;
	
	/* number of poitns along z axis inside refine subarea */
	Nz1 = 5;				

	/* get bodies geometry */
	r1x = B1.Rx;
	r2x = B2.Rx;
	r1y = B1.Ry;
	r2y = B2.Ry;
	
	/* write general information about domain */
	D.dkey = dkey;
	D.Nx = Nx;
	D.Ny = Ny;
	D.Nz0 = Nz0; 
	
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
	D.W  = W; // for line contact W = W/l here
	D.ph = ph;
	D.delta = delta;
	D.rhoR = rhoR;
	D.muR = muR;
	D.etaR = muR; // etaR = muR
	D.KR = KR;
	D.CpR = CpR;
	D.T0 = T0;
	D.IT1 = IT1;

	// parameters for reduced geometry
	RX = pow(1 / B1.Rx + 1 / B2.Rx, -1);
	RY = pow(1 / B1.Ry + 1 / B2.Ry, -1);
	Estar = pow(0.5 * ((1 - B1.PR * B1.PR) / B1.YM + (1 - B2.PR * B2.PR) / B2.YM), -1);
	
	D.rRX = RX;
	D.rRY = RY;
	D.rYm = Estar;

	um = (B1.vel[0] + B2.vel[0]) / 2;
	vm = (B1.vel[1] + B2.vel[1]) / 2;
	
	// piezocoefficient of viscosity
	alpha = 2e-8; // we use here a constant value, but in future we have to use it as variable
	
	std::fstream extfile_p;
	std::fstream extfile_h;

	if (dkey == 2) { /* 2D-domain (rectangular area) initialization */

		D.Area = new Point * [Nx];

		for (i = 0; i < Nx; ++i) {

			D.Area[i] = new Point[Ny];
			for (j = 0; j < Ny; ++j) {

				x_curr = -(a * kx) + i * (dx * a); /* m */
				y_curr = -(b * ky) + j * (dy * b); /* m */

				/* z-coords for surfaces of bodies in dimension form */
				
				/* original geometry, m */
				D.Area[i][j].z1_n = -1 / (2 * r1x) * x_curr * x_curr - 1 / (2 * r1y) * y_curr * y_curr;
				D.Area[i][j].z1_n =  1 / (2 * r2x) * x_curr * x_curr + 1 / (2 * r2y) * y_curr * y_curr + D.h0;

				/* reduced geometry, m */
				D.Area[i][j].z1 = 0.0; // in future, we can throw this out for memory saving
				D.Area[i][j].z2 = 1 / (2 * RX) * x_curr * x_curr + 1 / (2 * RY) * y_curr * y_curr;

				/* initial film thickness (before surface deformation), m */
				D.Area[i][j].h = D.Area[i][j].z2 + D.h0;
				
				/* set number of points along z axis in refined subarea */
				if ((x_curr >= xmin) && (x_curr <= xmax) && (y_curr >= ymin) && (y_curr <= ymax))
					D.Area[i][j].Nz = Nz1;
				else D.Area[i][j].Nz = Nz0;

				/* set initial dimensionless pressure */
				switch (pkey) {
					case 0: { // zero initialization
					D.Area[i][j].p = 0.0;
					break;
				}
					
					case 3: { // initalize from file
						if ((i == 0) && (j == 0)) {

							extfile_p.open("PC-PRS_DL_r" + std::to_string(D.rsx - 1) + "_init.DAT", std::ios::in);
							extfile_p >> x >> y >> par >> z; // read the first string
							nnodes = static_cast<int>(par);
							
							extfile_h.open("PC-THK_DL_r" + std::to_string(D.rsx - 1) + "_init.DAT", std::ios::in);
							extfile_h >> x >> y >> par >> z; // read the first string

							if ((!extfile_p) || (!extfile_h)) {
								std::cout << "INITIALIZATION FILES NOT FOUND, TRY ANOTHER INITIALIZATION METHOD" << "\n";
								std::exit(1);
							}
							
							D.h0 = z;

							if ( (par != D.Nx * D.Ny) || (x != kx) || (y != ky) ){
								std::cout << "INITIALIZATION ERROR, THE GRID NODES OR VALUES OF THE AREA SCALE FACTORS ARE INCORRECT" << "\n";
								std::exit(2);
							}
						}
									 
								extfile_p >> x >> y >> par;
								D.Area[i][j].p = par;
								extfile_h>> x >> y >> par;
								D.Area[i][j].h = par;
							
						break;
					} // case 3
						
					default: { // uniform or elliptic initialization
						q = 0;
						if (pkey == 2) q = 0.5;
					
						if (( (x_curr * x_curr) / (a * a) + (y_curr * y_curr) / (b * b) - 1) >= 0) D.Area[i][j].p = 0;
						else D.Area[i][j].p = pow((1 - (x_curr * x_curr) / (a * a) - (y_curr * y_curr) / (b * b)), q);
					} // default case
				} //switch
			
				//init_from_file(D, 1);


				/* set initial temperature, density, viscosity (uniform, dimensionless)  */
				/* create zero vectors for further storage of shear stress, velocity comp. */
				Nz_curr = D.Area[i][j].Nz;
				for (l = 0; l < Nz_curr; ++l) {
					D.Area[i][j].T.push_back(1.0);
					D.Area[i][j].Told.push_back(1.0);
					D.Area[i][j].rho.push_back(rhoR);
					D.Area[i][j].mu.push_back(muR);
					D.Area[i][j].th.push_back(KR);
					D.Area[i][j].hc.push_back(CpR);
					D.Area[i][j].eta.push_back(muR);
					D.Area[i][j].tau.push_back(0.0);
					D.Area[i][j].u.push_back(0.0);
					D.Area[i][j].v.push_back(0.0);
				}
					
				/* set initial deformation (they are equal to zero) */
				D.Area[i][j].d1 = 0.0;
				D.Area[i][j].d2 = 0.0;

				/* initialize the storage for pressure from previous iteration */
				D.Area[i][j].pold = D.Area[i][j].p;

				D.Area[i][j].dp = 0.0; // initialize pressure additions
			
			} /* j loop */
		
		} /* i loop */

		/* slide-to-roll ratios in domain */
		D.SRR[0] = 2 * (B2.vel[0] - B1.vel[0]) / (B2.vel[0] + B1.vel[0]);
		if ((B1.vel[1] == 0) && (B2.vel[1] == 0)) 
			D.SRR[1] = 0.0;
		else 
			D.SRR[1] = 2 * (B2.vel[1] - B1.vel[1]) / (B2.vel[1] + B1.vel[1]);

		/* Dowson-Higginson dimensionless parameters */
		D.W_HD = D.W / (2 * Estar * RX * RX);
		D.G_HD = 2 * alpha * Estar;
		D.U_HD = D.muR * pow(um * um + vm * vm, 0.5) / (2 * Estar * RX);

		k = 1.0339 * pow(RY / RX, 0.636); //ellipticity ratio

		/* The central and minimum film thickness (m), Hamrock, 2004 */
		D.Hc = 2.69 * pow(D.W_HD, -0.067) * pow(D.U_HD, 0.67) * pow(D.G_HD, 0.53) *
		(1 - 0.61 * exp(-0.73 * k)) * pow((RX / a),2);

		D.Hm = 3.63 * pow(D.W_HD, -0.073) * pow(D.U_HD, 0.68) * pow(D.G_HD, 0.49) *
		(1 - exp(-0.68 * k)) * pow((RX / a), 2);
		
		if (extfile_p) extfile_p.close();
		if (extfile_h) extfile_h.close();

	} /* if */
	else { /* 1-D domain (line segment) initialization */
		
		D.Area = new Point * [Nx];

		/* geometry, pressure and physical properties initialization */
		for (i = 0; i < Nx; ++i) {

			D.Area[i] = new Point[1];
			x_curr = -(a * kx) + i * (dx * a); /* m */

			/* original geometry, m */
			D.Area[i][0].z1_n = -1 / (2 * r1x) * x_curr * x_curr;
			D.Area[i][0].z2_n =  1 / (2 * r2x) * x_curr * x_curr + D.h0;

			/* reduced geometry, m */
			D.Area[i][0].z1 = 0.0; // may be we can throw this out from domain for memory saving
			D.Area[i][0].z2 = 1 / (2 * RX) * x_curr * x_curr;

			D.Area[i][0].h = D.Area[i][0].z2 + D.h0; // thickness before surface deformation

			/* set number of points along z axis in refined subarea */
			if ((x_curr >= xmin) && (x_curr <= xmax)) D.Area[i][0].Nz = Nz1;
			else D.Area[i][0].Nz = Nz0;
			
			/* set initial pressure */
			switch (pkey) {
					
				case 0: { // zero initialization
					D.Area[i][0].p = 0.0;
				break;
				}
				
				case 3: { // initalize from file
					if (i == 0) {

						extfile_p.open("LC-PRS_DL_r" + std::to_string(D.rsx - 1) + "_init.DAT", std::ios::in);
						extfile_p >> x >> par >> z; // read the first string
						nnodes = static_cast<int>(par);

						extfile_h.open("LC-THK_DL_r" + std::to_string(D.rsx - 1) + "_init.DAT", std::ios::in);
						extfile_h >> x >> par >> z; // read the first string

						if ((!extfile_p) || (!extfile_h)) {
							std::cout << "INITIALIZATION FILES NOT FOUND, TRY ANOTHER INITIALIZATION METHOD" << "\n";
							std::exit(1);
						}

						D.h0 = z;

						if ((par != D.Nx) || (x != kx)) {
							std::cout << "INITIALIZATION ERROR, THE GRID NODES OR VALUE OF THE AREA SCALE FACTOR IS INCORRECT" << "\n";
							std::exit(2);
						}
					}

					extfile_p >> x >> par;
					D.Area[i][0].p = par;
					extfile_h >> x >> par;
					D.Area[i][0].h = par;

					break;
				} // case 3


					default: { // uniform or elliptic initialization
					q = 0;
					if (pkey == 2) q = 0.5;

					if (((x_curr * x_curr) / (a * a) - 1) >= 0) D.Area[i][0].p = 0;
					else D.Area[i][0].p = pow((1 - (x_curr * x_curr) / (a * a)), q);
									
					}
				}

				/* set initial temperature, density, viscosity (uniform, dimensionless)  */
				/* create zero vectors for further storage of shear stress, velocity comp. */
				Nz_curr = D.Area[i][0].Nz;
				for (l = 0; l < Nz_curr; ++l) {
					D.Area[i][0].T.push_back(1.0);
					D.Area[i][0].Told.push_back(1.0);
					D.Area[i][0].rho.push_back(rhoR);
					D.Area[i][0].mu.push_back(muR);
					D.Area[i][0].th.push_back(KR);
					D.Area[i][0].hc.push_back(CpR);
					D.Area[i][0].eta.push_back(muR);
					D.Area[i][0].tau.push_back(0.0);
					D.Area[i][0].u.push_back(0.0);
					D.Area[i][0].v.push_back(0.0);
				}
					
				/*set initial deformation (they are equal to zero) */
				D.Area[i][0].d1 = 0.0;
				D.Area[i][0].d2 = 0.0;

				/* initialize the storage for pressure from previous iteration */
				D.Area[i][0].pold = D.Area[i][0].p;
	
		} /* i loop */
		
		/* slide-to-roll ratios in domain */
		D.SRR[0] = 2 * (B2.vel[0] - B1.vel[0]) / (B2.vel[0] + B1.vel[0]);
		
		/* Dowson-Higginson dimensionless parameters */
		D.W_HD = D.W / (2 * Estar * RX);
		D.G_HD = 2 * alpha * Estar;
		D.U_HD = D.muR * um / (2 * Estar * RX);

		/* The central and minimum film thickness (m), Hamrock, 2004 */
		k = 1.0339 * pow(RY / RX, 0.636);
		D.Hc = pow((RX / a), 2) * (2.69 * pow(D.U_HD, 0.67) * pow(D.G_HD, 0.53) * pow(D.W, -0.067) *
			   (1 - 0.61 * exp(-0.73 * k)));
		D.Hm = pow((RX / a), 2) * (3.63 * pow(D.U_HD, 0.68) * pow(D.G_HD, 0.49) * pow(D.W, -0.073) *
			   (1 - 1.00 * exp(-0.68 * k)));

		if (extfile_p) extfile_p.close();
		if (extfile_h) extfile_h.close();
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


void prsmodif(Domain& D)
{
	for (int i = 0; i < D.Nx; ++i)
		for (int j = 0; j < D.Ny; ++j) {
			
			if (abs(D.Area[i][j].p) < 1e-6) {
				D.Area[i][j].p = D.Area[i][j].pold * 1;
				D.Area[i][j].pold = D.Area[i][j].p;
			}
		}
}


void init_from_file(Domain& D, int kval) {
// kval: 1 - pressure, 2 - thickness, 3 - temperature
	
	int i, j, nnodes;
	double x, y, par;

	string fname;

	
	std::fstream extfile;
		
	switch (kval) {

	case 1: // pressure 
	{
		fname = "LC-bPRS_DL.DAT";
		extfile.open(fname, std::ios::in);
		break;
	} // case 1

	case 2: // thickness 
	{
		fname = "LC-bTHK_DL.DAT";
		extfile.open(fname, std::ios::in);

	} // case 2
	break;
	} // switch
	
	extfile >> par; // read the number of grid nodes
	nnodes = static_cast<int>(par);
	
	if (par != D.Nx * D.Ny) 
		std::cout << "INITIALIZATION ERROR, WRONG NUMBER OF GRID NODES" << "\n";
	else {

		i = 0; j = 0;

		while (extfile >> x >> y >> par)
		{
			std::cout << x << "\t" << y << "\t" << par << "\n";

			switch (kval) {

			case 1: // pressure 
			{
				D.Area[i][j].p = par;
				break;
			} // case 1

			case 2: // thickness 
			{
				D.Area[i][j].h = par;
				break;
			} // case 2
			
			break;
			} // switch

			i = i + 1; 
			j = j + 1;
		} // while
	} // else
} //proc



int main() {

	double a, b, ph, delta;
	double rsW_curr, rsT_curr;
	int lc;
	double rel_d;
	double pint_curr, err_curr, vel1i, vel2i, dh, dh_new;
	string dtprfx;

	
	/*****************************************************/
	/****************** INPUT DATA ***********************/
	/*****************************************************/
	
	int		rsx = 3;	// resolution in x-direction  (LC - 3..10; PC - 3..5)
	int		rsy = 3;	/* resolution in y-direction */

	double	kx = 2.0;	// x-scale factor for LC and PC (for PC must be kx=ky)
	double	ky = 2.0;	// y-scale factor  for LC this parameter is not used

	double	r1x = 10e-3;	/* m */
	double	r1y = 10e-3;	/* m */

	double	r2x = 10e-3;		/* m */
	double	r2y = 10e-3;		/* m */

	double	YM1 = 2e11;		/* Pa */
	double	PR1 = 0.3;
	double  K1 = 45.0;		/* W / (m K) */
	double  CP1 = 480.0;	/* J / (kg K) */
	double  Ts1 = 295.0;	/* K */

	double	vel1[2];	/* vector of surface velocity for body 1 */

	vel1[0] = 2.12;		/* m/s */
	vel1[1] = 0.0;		/* m/s */

	double	YM2 = 2e11;	/* Pa */
	double	PR2 = 0.3;
	double  K2 = 45.0;	/* W / (m K) */
	double  CP2 = 480.0;	/* J / (kg K) */
	double  Ts2 = 295.0;	/* K */

	double	vel2[2];	/* vector of surface velocity for body 1 */

	vel2[0] = 2.12;		/* m/s */
	vel2[1] = 0.0;		/* m/s */

	double	F =	10;		/* normal contact force, N or N/m */
	
	double len = 1;	/* length for line contact problem */

	int		pkey = 3;	/* initial pressure profile: 0-zero; 1-uniform; 2-elliptical; 3-from file */

	/* reference values */
	
	double	rhoR = 829.0;	/* kg / m^3 */
		
	double	muR = 0.0847;	/* Pa s */
	//double	muR = 0.028;
	
	double	KR = 0.12;		/* W / (m K) */
	double	CpR = 1670.0;	/* J / (kg K) */
	double	T0 = 295.0;		/* K */

	/* folder for output files, must be created */
	std::string pth = "E:\\_work\\EHL-prg\\TEHL\\TEHL-v2\\_OUTPUT\\";
	D.pth = pth;

	/****************** END OF SECTION *******************/



	/*****************************************************/
	/**************** INITIALIZATION *********************/
	/*****************************************************/
	
	int dtype = 2; /* 1 - line contact; 2 - point contact, 2D-domain */

	D.rsx = rsx;
	D.rsy = rsy;

	/* ELASTIC BODIES initialization */
	elbody_initialize(1, r1x, r1y, YM1, PR1, K1, CP1, vel1, Ts1, 0.0, "This is bottom body", B1);
	elbody_initialize(2, r2x, r2y, YM2, PR2, K2, CP2, vel2, Ts2, 0.0, "This is top body", B2);

	/* CONTACT parameters (herzian model) */
	{
		if (dtype == 2) /* point contact */
			herzian_parameters(B1, B2, 1, F, 0.0, a, b, ph, delta);
		else
			/* line contact */
			herzian_parameters(B1, B2, 0, F, len, a, b, ph, delta);
	}
	
	/* PHYSICAL MODELS initialization */
	{
		/* set DENSITY model */

			// mkey=0 - constant 
			// mkey=1 - Dowson & Higginson
			// mkey=2 - Tait 
			// mkey=3 - Murnaghan 
			// mkey=4 - Dowson & Higginson var 1

		RHO.mkey = 2;

		switch (RHO.mkey) {

		case 0: {

			RHO.IT1 = "CONSTANT VALUE";

			RHO.arg1 = 873.0;		/* rho0  */
			RHO.narg1 = "RHO_0\t";

			break;
		}

		case 1: {

			RHO.IT1 = "DOWSON & HIGGINSON";

			RHO.arg1 = 873.0;		/* rhoO */
			RHO.arg2 = 0.59e-9;		/* C1   */
			RHO.arg3 = 1.7e-9;		/* C2   */
			RHO.arg4 = 0.0;			/* -    */
			RHO.arg5 = 0.0;			/* -    */
			RHO.arg6 = 0.0;			/* -    */

			RHO.narg1 = "RHO_0\t";
			RHO.narg2 = "C1\t\t";
			RHO.narg3 = "C2\t\t";
			
			break;
		}

		case 2: {

			RHO.IT1 = "TAIT EQUATION";

			// Mineral Oil 23D, Encyclopedia of Tribology,2013. p.3528-3529
			RHO.arg1 = 829;				/* rhoR  */
			RHO.arg2 = 273.15;			/* TR	 */
			RHO.arg3 = 5.567e-3;		/* betaK */
			RHO.arg4 = 8.25e-4;			/* aV	 */
			RHO.arg5 = 10.341;			/* K0	 */
			RHO.arg6 = 7.334e9;			/* K00	 */

			RHO.narg1 = "RHO_R\t";
			RHO.narg2 = "T_R\t\t";
			RHO.narg3 = "BETA_K\t";
			RHO.narg4 = "aV\t\t";
			RHO.narg5 = "K_0\t\t";
			RHO.narg6 = "K_00\t";

			break;
		}

		case 3: {

			RHO.IT1 = "MURNAGHAN EQUATION";

			// Shell T9 oil, Habchi, 2018. p.423
			RHO.arg1 = 872.0;			/* rhoR  */
			RHO.arg2 = 298.0;			/* TR	 */
			RHO.arg3 = 6.09e-3;			/* betaK */
			RHO.arg4 = 7.734e-4;		/* aV	 */
			RHO.arg5 = 10.545;			/* K0	 */
			RHO.arg6 = 9.234e9;			/* K00	 */

			RHO.narg1 = "RHO_R\t";
			RHO.narg2 = "T_R\t\t";
			RHO.narg3 = "BETA_K\t";
			RHO.narg4 = "aV\t\t";
			RHO.narg5 = "K_0\t\t";
			RHO.narg6 = "K_00\t";

			break;
		}

		case 4: {

			RHO.IT1 = "DOWSON & HIGGINSON VAR 1";

			RHO.arg1 = 873.0;		/* rhoO */
			RHO.arg2 = 5.9e8;		/* C   */
			RHO.arg3 = 1.34;		/* pcoef   */
			RHO.arg4 = 0.0;			/* -    */
			RHO.arg5 = 0.0;			/* -    */
			RHO.arg6 = 0.0;			/* -    */

			RHO.narg1 = "RHO_0\t";
			RHO.narg2 = "C\t\t";
			RHO.narg3 = "P_COEF\t";

			break;
		}

		} /* switch */

		// set LOW-SHEAR VISCOSITY model

			// mkey=0 - constant 
			// mkey=1 - Barus 
			// mkey=2 - Roeland 
			// mkey=3 - Johari & Whalley 
			// mkey=4 - Doolitle 
			// mkey=5 - Yasutomi (improved)

		LSV.mkey = 2;

		switch (LSV.mkey) {

		case 0: {
			LSV.IT1 = "CONSTANT VALUE";

			LSV.arg1 = 0.08474;		/* mu0   */
			LSV.narg1 = "MU_0\t";

			break;
		}

		case 1: {
			LSV.IT1 = "BARUS EQUATION";

			LSV.arg1 = 0.08474;		/* muR   */
			//LSV.arg1 = 0.028;		/* muR   */
			LSV.arg2 = 2.1e-8;		/* alpha */
			LSV.arg3 = 0.0;			/* -	 */
			LSV.arg4 = 0.0;			/* -     */
			LSV.arg5 = 0.0;			/* -     */
			LSV.arg6 = 0.0;			/* -     */
			LSV.arg7 = 0.0;			/* -     */
			LSV.arg8 = 0.0;			/* -     */

			LSV.narg1 = "MU_R\t";
			LSV.narg2 = "ALPHA\t";

			break;
		}

		case 2: {
			LSV.IT1 = "ROELAND EQUATION";
			// Mineral Oil 23D, Encyclopedia of Tribology,2013. p.3535
			LSV.arg1 = 0.0847;		/* muR   */
			//LSV.arg1 = 0.028;		/* muR   */
			LSV.arg2 = 0.0631e-3;	/* muP	 */
			LSV.arg3 = -196e6;		/* Pp	 */
			LSV.arg4 = 273.0;		/* TR    */
			LSV.arg5 = 138.0;		/* Tinf  */
			LSV.arg6 = 0.79;		/* Z     */
			LSV.arg7 = 1.09;		/* S     */
			LSV.arg8 = 0.0;			/* -     */

			LSV.narg1 = "MU_R\t";
			LSV.narg2 = "MU_P\t";
			LSV.narg3 = "P_P\t\t";
			LSV.narg4 = "TR\t\t";
			LSV.narg5 = "TINF\t";
			LSV.narg6 = "Z\t\t";
			LSV.narg7 = "S\t\t";

			break;
		}

		case 3: {

			LSV.IT1 = "JOHARI & WHALLEY EQUATION";
			// Mineral Oil 23D, Encyclopedia of Tribology, 2013. p.3535
			LSV.arg1 = 2.78e-14;	/* mu_inf */
			LSV.arg2 = 2.73e9;		/* p_inf  */
			LSV.arg3 = 27.5;		/* CF     */
			LSV.arg4 = 0.0;			/* -	  */
			LSV.arg5 = 0.0;			/* -	  */
			LSV.arg6 = 0.0;			/* -	  */
			LSV.arg7 = 0.0;			/* -	  */
			LSV.arg8 = 311.0;		/* T_ref  */
			// there is no T_ref in equation, but it define a set of model constant
			// this model has the following conditon: p < p_inf

			LSV.narg1 = "MU_INF\t";
			LSV.narg2 = "P_INF\t";
			LSV.narg3 = "CF\t\t";
			LSV.narg8 = "T_REF\t";

			break;
		}

		case 4: {

			LSV.IT1 = "DOOLITLE EQUATION";
			// 80W-90 Gear oil, Bair, 2019. p.167 
			LSV.arg1 = 34.8e-3;		/* muR */
			LSV.arg2 = 323;			/* TR  */
			LSV.arg3 = 4.336;		/* B   */
			LSV.arg4 = 0.6988;		/* R0  */
			LSV.arg5 = -7.679e-4;	/* eps */
			LSV.arg6 = 0.0;			/* -   */
			LSV.arg7 = 0.0;			/* -   */
			LSV.arg8 = 0.0;			/* -   */

			LSV.narg1 = "MU_R\t";
			LSV.narg2 = "TR\t\t";
			LSV.narg3 = "B\t\t";
			LSV.narg4 = "R0\t\t";
			LSV.narg5 = "EPS_C\t";

			break;
		}

		case 5: {

			LSV.IT1 = "YASUTOMI (IMPROVED) EQUATION";
			// 80W-90 Gear oil, Bair, 2019
			LSV.arg1 = 1e12;		/* muG, Pa s */
			LSV.arg2 = 201.45;		/* TG0, K    */
			LSV.arg3 = 483.75;		/* A1, K     */
			LSV.arg4 = 0.455e9;		/* A2, Pa    */
			LSV.arg5 = 7.606e9;		/* b1, Pa    */
			LSV.arg6 = -0.4193;		/* b2        */
			LSV.arg7 = 16.372;		/* C1, K     */
			LSV.arg8 = 303.14;		/* C2, K     */

			LSV.narg1 = "MU_G\t";
			LSV.narg2 = "TG0\t\t";
			LSV.narg3 = "A1\t\t";
			LSV.narg4 = "A2\t\t";
			LSV.narg5 = "b1\t\t";
			LSV.narg6 = "b2\t\t";
			LSV.narg7 = "C1\t\t";
			LSV.narg8 = "C2\t\t";

			break;
		}

		} /* switch */

		/* set THERMAL CONDUCTIVITY model */

			// mkey = 0 - constant
			// mkey = 1 - pressure depended
			// mkey = 2 - power-law 

		THC.mkey = 1;

		switch (THC.mkey) {

		case 0: {
			THC.IT1 = "CONSTANT VALUE";

			THC.arg1 = 0.135;		/* K0 */
			THC.narg1 = "K0\t\t";

			break;
		}

		case 1: {

			THC.IT1 = "PRESSURE DEPENDED";

			// Paraphinic mineral oil, Larsson, 2000)
			THC.arg1 = 0.137;		/* lambda0 */
			THC.arg2 = 1.72;		/* c1  */
			THC.arg3 = 0.54;		/* c2  */
			THC.arg4 = 0.0;			/* -  */
			THC.arg5 = 0.0;			/* -  */

			THC.narg1 = "LAMBDA_0";
			THC.narg2 = "C1\t\t";
			THC.narg3 = "C2\t\t";

			break;
		}

		case 2: {

			HCP.IT1 = "POWER-LAW EQUATION";

			// Shell T9 oil, Habchi, 2018. p.427 
			HCP.arg1 = 0.053;		/* B */
			HCP.arg2 = 0.026;		/* C */
			HCP.arg3 = 7.66;		/* s */
			HCP.arg4 = 0.0;			/* - */
			HCP.arg5 = 0.0;			/* - */
			HCP.arg6 = 0.0;			/* - */
			HCP.arg7 = 0.0;			/* - */
			HCP.arg8 = 0.0;			/* - */

			HCP.narg1 = "B\t\t";
			HCP.narg2 = "C\t\t";
			HCP.narg3 = "s\t\t";

			break;
		}

		} /* switch */

		/* set HEAT CAPACITY model */

			// mkey=0 - constant 
			// mkey=1 - linear
			// mkey=2 - power-law 

		HCP.mkey = 1;

		switch (HCP.mkey) {

		case 0: {
			HCP.IT1 = "CONSTANT VALUE";

			HCP.arg1 = 1.71e-6;		/* rho*Cp0 */
			HCP.narg1 = "RHO*Cp_0\t";

			break;
		}

		case 1: {

			HCP.IT1 = "LINEAR EQUATION";

			// Shell T9 oil, Habchi, 2018. p.428)
			HCP.arg1 = 1.17e6;		/* C0 */
			HCP.arg2 = 0.39e6;		/* m  */
			HCP.arg3 = 0.0;			/* -  */
			HCP.arg4 = 0.0;			/* -  */
			HCP.arg5 = 0.0;			/* -  */
			HCP.arg6 = 0.0;			/* -  */
			HCP.arg7 = 298.0;		/* TR */
			
			HCP.narg1 = "C0\t\t";
			HCP.narg2 = "m\t\t";
			HCP.narg7 = "TR\t\t";

			break;
		}

		case 2: {

			HCP.IT1 = "PRESSURE DEPENDED";

			// Paraphinic mineral oil, Larsson, 2000 
			HCP.arg1 = 0.47;		/* k1 */
			HCP.arg2 = 0.81;		/* k2 */
			HCP.arg3 = 9.3e-4;		/* beta0 */
			HCP.arg4 = 1.4;			/* b1 */
			HCP.arg5 = -0.51;		/* b2 */
			HCP.arg6 = 295;			/* T0 */
			HCP.arg7 = 1.71e-6;		/* rhoCp0  */
			
			HCP.narg1 = "k1\t\t";
			HCP.narg2 = "k2\t\t";
			HCP.narg3 = "beta0\t";
			HCP.narg4 = "b1\\tt";
			HCP.narg5 = "b2\t\t";
			HCP.narg6 = "T0\t\t";
			HCP.narg7 = "RHO*CP0\t";
			break;
		}

		} /* switch */


		/* set RHEOLOGY MODEL (shear stress - viscosity)  */

			// mkey=0 - Newton behavior 
			// mkey=1 - Ree-Eyring
			// mkey=2 - Carreau
			// mkey=3 - Carreau-Yasuda (double modified)

		GNV.mkey = 0;

		switch (GNV.mkey) {

		case 0: {
			GNV.IT1 = "NEWTON BEHAVIOR";

			GNV.arg1 = 0.0;		/* - */
			GNV.arg2 = 0.0;		/* -  */
			GNV.arg3 = 0.0;		/* -  */
			GNV.arg4 = 0.0;		/* -  */
			GNV.arg5 = 0.0;		/* -  */
			GNV.arg6 = 0.0;		/* -  */
			GNV.arg7 = 0.0;		/* -  */

			break;
		}

		case 1: {

			GNV.IT1 = "REE-EYRING EQUATION";

			GNV.arg1 = 10e6;		/* tau0 */
			GNV.arg2 = 0.0;		/* -  */
			GNV.arg3 = 0.0;		/* -  */
			GNV.arg4 = 0.0;		/* -  */
			GNV.arg5 = 0.0;		/* -  */
			GNV.arg6 = 0.0;		/* -  */
			GNV.arg7 = 0.0;		/* -  */

			GNV.narg1 = "TAU_0\t";

			break;
		}

		case 2: {

			GNV.IT1 = "CARREAU EQUATION";

			GNV.arg1 = 0.0;		/* - */
			GNV.arg2 = 0.0;		/* -  */
			GNV.arg3 = 0.0;		/* -  */
			GNV.arg4 = 0.0;		/* -  */
			GNV.arg5 = 0.0;		/* -  */
			GNV.arg6 = 0.0;		/* -  */
			GNV.arg7 = 0.0;		/* -  */

			GNV.narg1 = "  \t";

			break;
		}

		case 3: {

			GNV.IT1 = "CARREAU-YASUDA (DBL. MOD.) EQUATION";

			GNV.arg1 = 0.0;		/* - */
			GNV.arg2 = 0.0;		/* -  */
			GNV.arg3 = 0.0;		/* -  */
			GNV.arg4 = 0.0;		/* -  */
			GNV.arg5 = 0.0;		/* -  */
			GNV.arg6 = 0.0;		/* -  */
			GNV.arg7 = 0.0;		/* -  */

			GNV.narg1 = "  \t";

			break;
		}

		} /* switch */


	/* write models to domian */
		D.RHO = RHO;
		D.LSV = LSV;
		D.GNV = GNV;
		D.THC = THC;
		D.HCP = HCP;
	}

	/* INITIALIZE DOMAIN */
	
	if (dtype == 2) { //point contact
		domain_initialize(2, B1, B2, rsx, rsy, a, b, kx, ky, F, ph, delta, rhoR, muR, KR, CpR, T0, pkey,
			"This is 2D-domain for TEHL (point contact) calculations with herzian initial pressure field", D);
		
		/* write the initial pressure profile */
		//output_data(D, B1, B2, 1, 3, '-', 0.0, 0, pth + "PC_PRS0_DL");
		//output_data(D, B1, B2, 1, 3, 'Y', 0.0, 0, pth + "PC_PRS0-XZ_DL");

		/* write the initial profile of reduced geometry for TEHL model */
		//output_data(D, B1, B2, 0, 2, '-', 0.0, 1, pth + "PC_B2-SURF_DM");
		//output_data(D, B1, B2, 0, 2, 'Y', 0.0, 1, pth + "PC_B2-SURF-XZ_DM");

		/* write initial deformations and film thickness */
		//calc_deform(D, B1, B2);
		//output_data(D, B1, B2, 0, 4, '-', 1, 0, pth + "PC_THK_DL");
		//output_data(D, B1, B2, 0, 4, 'Y', 1, 0, pth + "PC_THK-XZ_DL");
		//output_data(D, B1, B2, 0, 4, '-', 1, 1, pth + "PC_THK_DM");
		//output_data(D, B1, B2, 0, 4, 'Y', 1, 1, pth + "PC_THK-XZ_DM");

		/* write the real profile of reduced geometry for TEHL model */
		//output_data(D, B1, B2, 1, 2, '-', 0.0, 1, pth + "PC_B2-REAL_DM");
		//output_data(D, B1, B2, 1, 2, 'Y', 0.0, 1, pth + "PC_B2-REAL-XZ_DM");
	}
	else { // line contact

		len = 1.0; //length for line contact problem 
		domain_initialize(1, B1, B2, rsx, rsy, a, b, kx, ky, F / len, ph, delta, rhoR, muR, KR, CpR, T0, pkey,
			"This is 1D-domain for TEHL (line contact) calculations with herzian initial pressure field", D);
	
		/* write the initial pressure profile */
		//output_data(D, B1, B2, 1, 3, '-', 0.0, 0, pth + "LC_PRS0_DL");
		
		/* write the initial profile of reduced geometry for TEHL model */
		//output_data(D, B1, B2, 0, 2, '-', 0.0, 1, pth + "LC_B2-SURF_DM");

		/* write initial deformations and film thickness */
		//calc_deform(D, B1, B2);
		//output_data(D, B1, B2, 0, 4, '-', 1, 0, pth + "LC_THK_DL");
		//output_data(D, B1, B2, 0, 4, '-', 1, 1, pth + "LC_THK_DM");

		/* write the real profile of reduced geometry for TEHL model */
		//output_data(D, B1, B2, 1, 2, '-', 0.0, 1, pth + "LC_B2-REAL_DM");
	}
	
	/* write general information about current TEHL calculation */
	//output_data(D, B1, B2, 0, 0, '-', 1, 0, pth + "GENERAL-INF");

	/****************** END OF SECTION *******************/

	
	


	/*****************************************************/
	/*********** TEHL CALCULATION MAIN CYCLE *************/
	/*****************************************************/
		
	switch (D.dkey) {
	
	case 1: { /* LINE CONTACT */
	
		/* pressure and film thickness calculation */
		calc_pressure_field_1D(D, B1, B2);
				
		/* temperature calculation on body surfaces  */
		//calc_sftemp_1D(D, B1, B2);

		/* calculation of the new temperature field in domain */
		//calc_temperature_field_1D(D, B1, B2, 0);

		/* tangential forces and friction coefficient calculation */
		calc_fric_1D(D, B1, B2);

		dumpvar1D(D, 1, 1, 0); // 2X-interpolate and write to file the pressure field (for initialization)
		dumpvar1D(D, 2, 1, 0); // 2X-interpolate and write to file the thickness field (for initialization)
		dumpvar1D(D, 1, 1, 1); // 2X-interpolate and write to file the pressure field (for plot the surface)
		dumpvar1D(D, 2, 1, 1); // 2X-interpolate and write to file the thickness field (for plot the surface)

		/* calculation the heat flux on body surfaces */
		//calc_sfhflux_1D(D, B1, B2);


		break;
	}// case 1

	case 2: { /* POINT CONTACT */
		//calc_pressure_field_2D(D, B1, B2);
		//output_data(D, B1, B2, 1, 3, '-', 1, 0, pth + "PC-0-" + "PRS_DL");
		
		calc_pressure_field_2D(D, B1, B2);
		calc_fric_2D(D, B1, B2);
		
		dumpvar2D(D, 1, 1, 0); // 2X-interpolate and write to file the pressure field (for initialization)
		dumpvar2D(D, 2, 1, 0); // 2X-interpolate and write to file the thickness field (for initialization)
		dumpvar2D(D, 1, 1, 1); // 2X-interpolate and write to file the pressure field (for plot the surface)
		dumpvar2D(D, 2, 1, 1); // 2X-interpolate and write to file the thickness field (for plot the surface)

		

		break;
	}// case 2

	} // main switch

	/****************** END OF SECTION *******************/



	/*****************************************************/
	/******************* OUTPUT RESULTS ******************/
	/*****************************************************/
	
	/* set output pattern */

	bool isGEOMTR = true;	/* domain geometry (B1, B2, Z, DF1, DF2) */
	bool isDEFORM = true;	/* deformed body surfaces, reduced geometry (B1real, B2real) */
	
	bool isPROPS1 = true;	/* density, viscosity (low-shear) */
	bool isPROPS2 = true;	/* shear stress, viscosity (generalized) */
	bool isPROPS3 = false;	/* thermal conductivity, heat capacity */
	
	bool isPRESSR = true;	/* pressure field */
	bool isTEMPER = false;	/* temperature field */
	bool isVELCTY = true;	/* velocity field */

	bool isFRICTN = true;	/* friction in domain (distributed and integrated results) */
	bool isHFLUXS = false;	/* heating in domain (distributed and integrated results) */

	
	double spx_all = 0.0; /* DL-position for section in YZ plane */
	double spy_all = 0.0; /* DL-position for section in XZ plane */

	double spx = 0.0; 
	double spy = 0.0;

	bool isall = true; /* if true, then we will use spx_all, spy_all for all items */


	/* write general information about current TEHL calculation */
	output_data(D, B1, B2, 0, 0, '-', 1, 0, pth + "GENERAL-INF");


	if (D.dkey == 2) dtprfx = "PC_"; else dtprfx = "LC_";

	if (isGEOMTR) {
		if (isall = true) {
			spx = spx_all;
			spy = spy_all;
		}

		else {
			spx = 0.25; /* DL-position for section in YZ plane */
			spy = 0.25; /* DL-position for section in XZ plane */
		}

		output_data(D, B1, B2, 0, 1, '-',   1, 1, pth + dtprfx + "B1-SURF_DM");
		if (D.dkey == 2) {
			output_data(D, B1, B2, 0, 1, 'X', spx, 1, pth + dtprfx + "B1-SURF-YZ_DM");
			output_data(D, B1, B2, 0, 1, 'Y', spy, 1, pth + dtprfx + "B1-SURF-XZ_DM");
		}
		
		output_data(D, B1, B2, 0, 2, '-', 1  , 1, pth + dtprfx + "B2-SURF_DM");
		if (D.dkey == 2) {
			output_data(D, B1, B2, 0, 2, 'X', spx, 1, pth + dtprfx + "B2-SURF-YZ_DM");
			output_data(D, B1, B2, 0, 2, 'Y', spy, 1, pth + dtprfx + "B2-SURF-XZ_DM");
		}
		
		output_data(D, B1, B2, 0, 3, '-', 1  , 0, pth + dtprfx + "B2-REDU_DL");
		if (D.dkey == 2) {
			output_data(D, B1, B2, 0, 3, 'X', spx, 0, pth + dtprfx + "B2-REDU-YZ_DL");
			output_data(D, B1, B2, 0, 3, 'Y', spy, 0, pth + dtprfx + "B2-REDU-XZ_DL");
		}
		
		output_data(D, B1, B2, 0, 3, '-', 1  , 1, pth + dtprfx + "B2-REDU_DM");
		if (D.dkey == 2) {
			output_data(D, B1, B2, 0, 3, 'X', spx, 1, pth + dtprfx + "B2-REDU-YZ_DM");
			output_data(D, B1, B2, 0, 3, 'Y', spy, 1, pth + dtprfx + "B2-REDU-XZ_DM");
		}

		output_data(D, B1, B2, 0, 4, '-', 1, 0, pth + dtprfx + "Z-RESOL");
	}

	if (isDEFORM) {
		
		if (isall = true) {
			spx = spx_all;
			spy = spy_all;
		}

		else {
			spx = 0.25; /* DL-position for section in YZ plane */
			spy = 0.25; /* DL-position for section in XZ plane */
		}
				
		output_data(D, B1, B2, 1, 2, '-', 1  , 0, pth + dtprfx + "B2-REAL_DL");
		if (D.dkey == 2) {
			output_data(D, B1, B2, 1, 2, 'X', spx, 0, pth + dtprfx + "B2-REAL-YZ_DL");
			output_data(D, B1, B2, 1, 2, 'Y', spy, 0, pth + dtprfx + "B2-REAL-XZ_DL");
		}
		
		output_data(D, B1, B2, 1, 2, '-', 1, 1, pth + dtprfx + "B2-REAL_DM");
		if (D.dkey == 2) {
			output_data(D, B1, B2, 1, 2, 'X', spx, 1, pth + dtprfx + "B2-REAL-YZ_DM");
			output_data(D, B1, B2, 1, 2, 'Y', spy, 1, pth + dtprfx + "B2-REAL-XZ_DM");
		}
	
		output_data(D, B1, B2, 0, 5, '-', 1, 0, pth + dtprfx + "THK_DL");
		if (D.dkey == 2) {
			output_data(D, B1, B2, 0, 5, 'X', spx, 0, pth + dtprfx + "THK-YZ_DL");
			output_data(D, B1, B2, 0, 5, 'Y', spy, 0, pth + dtprfx + "THK-XZ_DL");
		}

		output_data(D, B1, B2, 0, 5, '-', 1, 1, pth + dtprfx + "THK_DM");
		if (D.dkey == 2) {
			output_data(D, B1, B2, 0, 5, 'X', spx, 1, pth + dtprfx + "THK-YZ_DM");
			output_data(D, B1, B2, 0, 5, 'Y', spy, 1, pth + dtprfx + "THK-XZ_DM");
		}

		output_data(D, B1, B2, 0, 7, '-', 1, 0, pth + dtprfx + "DF2_DL");
		if (D.dkey == 2) {
			output_data(D, B1, B2, 0, 7, 'X', spx, 0, pth + dtprfx + "DF2-YZ_DL");
			output_data(D, B1, B2, 0, 7, 'Y', spy, 0, pth + dtprfx + "DF2-XZ_DL");
		}

		output_data(D, B1, B2, 0, 7, '-', 1, 1, pth + dtprfx + "DF2_DM");
		if (D.dkey == 2) {
			output_data(D, B1, B2, 0, 7, 'X', spx, 1, pth + dtprfx + "DF2-YZ_DM");
			output_data(D, B1, B2, 0, 7, 'Y', spy, 1, pth + dtprfx + "DF2-XZ_DM");
		}
	}

	if (isPROPS1) {

		if (isall = true) {
			spx = spx_all;
			spy = spy_all;
		}

		else {
			spx = 0.25; /* DL-position for section in YZ plane */
			spy = 0.25; /* DL-position for section in XZ plane */
		}

		output_data(D, B1, B2, 1, 4, '-', 1, 0, pth + dtprfx + "RHO_DL");
		if (D.dkey == 2) {
			output_data(D, B1, B2, 1, 4, 'X', spx, 0, pth + dtprfx + "RHO-YZ_DL");
			output_data(D, B1, B2, 1, 4, 'Y', spy, 0, pth + dtprfx + "RHO-XZ_DL");
		}
		
		output_data(D, B1, B2, 1, 4, '-', 1, 1, pth + dtprfx + "RHO_DM");
		if (D.dkey == 2) {
			output_data(D, B1, B2, 1, 4, 'X', spx, 1, pth + dtprfx + "RHO-YZ_DM");
			output_data(D, B1, B2, 1, 4, 'Y', spy, 1, pth + dtprfx + "RHO-XZ_DM");
		}
		
		output_data(D, B1, B2, 1, 5, '-', 1, 0, pth + dtprfx + "LSV_DL");
		if (D.dkey == 2) {
			output_data(D, B1, B2, 1, 5, 'X', spx, 0, pth + dtprfx + "LSV-YZ_DL");
			output_data(D, B1, B2, 1, 5, 'Y', spy, 0, pth + dtprfx + "LSV-XZ_DL");
		}
		
		output_data(D, B1, B2, 1, 5, '-', 1, 1, pth + dtprfx + "LSV_DM");
		if (D.dkey == 2) {
			output_data(D, B1, B2, 1, 5, 'X', spx, 1, pth + dtprfx + "LSV-YZ_DM");
			output_data(D, B1, B2, 1, 5, 'Y', spy, 1, pth + dtprfx + "LSV-XZ_DM");
		}
	}

	if (isPROPS2) {

		if (isall = true) {
			spx = spx_all;
			spy = spy_all;
		}

		else {
			spx = 0.25; /* DL-position for section in YZ plane */
			spy = 0.25; /* DL-position for section in XZ plane */
		}

		output_data(D, B1, B2, 1, 6, '-', 1, 0, pth + dtprfx + "TAU_DL");
		if (D.dkey == 2) {
			output_data(D, B1, B2, 1, 6, 'X', spx, 0, pth + dtprfx + "TAU-YZ_DL");
			output_data(D, B1, B2, 1, 6, 'Y', spy, 0, pth + dtprfx + "TAU-XZ_DL");
		}

		output_data(D, B1, B2, 1, 6, '-', 1, 1, pth + dtprfx + "TAU_DM");
		if (D.dkey == 2) {
			output_data(D, B1, B2, 1, 6, 'X', spx, 1, pth + dtprfx + "TAU-YZ_DM");
			output_data(D, B1, B2, 1, 6, 'Y', spy, 1, pth + dtprfx + "TAU-XZ_DM");
		}

		output_data(D, B1, B2, 1, 7, '-', 1, 0, pth + dtprfx + "GNV_DL");
		if (D.dkey == 2) {
			output_data(D, B1, B2, 1, 7, 'X', spx, 0, pth + dtprfx + "GNV-YZ_DL");
			output_data(D, B1, B2, 1, 7, 'Y', spy, 0, pth + dtprfx + "GNV-XZ_DL");
		}

		output_data(D, B1, B2, 1, 7, '-', 1, 1, pth + dtprfx + "GNV_DM");
		if (D.dkey == 2) {
			output_data(D, B1, B2, 1, 7, 'X', spx, 1, pth + dtprfx + "GNV-YZ_DM");
			output_data(D, B1, B2, 1, 7, 'Y', spy, 1, pth + dtprfx + "GNV-XZ_DM");
		}
	}

	if (isPROPS3) {

		if (isall = true) {
			spx = spx_all;
			spy = spy_all;
		}

		else {
			spx = 0.25; /* DL-position for section in YZ plane */
			spy = 0.25; /* DL-position for section in XZ plane */
		}
			

	}

	if (isPRESSR) {

		if (isall = true) {
			spx = spx_all;
			spy = spy_all;
		}

		else {
			spx = 0.25; /* DL-position for section in YZ plane */
			spy = 0.25; /* DL-position for section in XZ plane */
		}

		output_data(D, B1, B2, 1, 3, '-', 1, 0, pth + dtprfx + "PRS_DL");
		if (D.dkey == 2) {
			output_data(D, B1, B2, 1, 3, 'X', spx, 0, pth + dtprfx + "PRS-YZ_DL");
			output_data(D, B1, B2, 1, 3, 'Y', spy, 0, pth + dtprfx + "PRS-XZ_DL");
		}
		
		output_data(D, B1, B2, 1, 3, '-', 1, 1, pth + dtprfx + "PRS_DM");
		if (D.dkey == 2) {
			output_data(D, B1, B2, 1, 3, 'X', spx, 1, pth + dtprfx + "PRS-YZ_DM");
			output_data(D, B1, B2, 1, 3, 'Y', spy, 1, pth + dtprfx + "PRS-XZ_DM");
		}
	}
	
	if (isTEMPER) {

		if (isall = true) {
			spx = spx_all;
			spy = spy_all;
		}

		else {
			spx = 0.25; /* DL-position for section in YZ plane */
			spy = 0.25; /* DL-position for section in XZ plane */
		}


	}

	if (isVELCTY) {

		if (isall = true) {
			spx = spx_all;
			spy = spy_all;
		}

		else {
			spx = 0.25; /* DL-position for section in YZ plane */
			spy = 0.25; /* DL-position for section in XZ plane */
		}

		output_data(D, B1, B2, 1, 8, '-', 1, 1, pth + dtprfx + "VEL-U_DM");
		if (D.dkey == 2) {
			output_data(D, B1, B2, 1, 8, 'X', spx, 1, pth + dtprfx + "VEL-U-YZ_DM");
			output_data(D, B1, B2, 1, 8, 'Y', spy, 1, pth + dtprfx + "VEL-U-XZ_DM");

			output_data(D, B1, B2, 1, 9, '-', 1, 1  , pth + dtprfx + "VEL-V_DM");
			output_data(D, B1, B2, 1, 9, 'X', spx, 1, pth + dtprfx + "VEL-V-YZ_DM");
			output_data(D, B1, B2, 1, 9, 'Y', spy, 1, pth + dtprfx + "VEL-V-XZ_DM");

		}




	}
	
	if (isFRICTN) {

		if (isall = true) {
			spx = spx_all;
			spy = spy_all;
		}

		else {
			spx = 0.25; /* DL-position for section in YZ plane */
			spy = 0.25; /* DL-position for section in XZ plane */
		}

		output_data(D, B1, B2, 2, 1, '-', 1, 0, pth + dtprfx + "FRICTION");
	}

	if (isHFLUXS) {

		if (isall = true) {
			spx = spx_all;
			spy = spy_all;
		}

		else {
			spx = 0.25; /* DL-position for section in YZ plane */
			spy = 0.25; /* DL-position for section in XZ plane */
		}


	}


	/* output graphics */
	if (D.dkey == 1) {
		//output_graph(D, 1);
	}
	else 
		//output_graph(D, 2);

	/* free CPU memory */
	domain_finalize(D);

} /* end of program */









