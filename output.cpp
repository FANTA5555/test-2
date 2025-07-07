#include <iostream>
#include <format>
#include "main.h"
//#include <matplot/matplot.h>
#include <set>
#include <thread>
#include <vector>

/*

void output_graph(Domain& D, int dtype) {

	int i, j, N;
	double par1, par2, par3;
	string fnamep, fnameh;

	using namespace matplot;

	std::fstream extfile_p;
	std::fstream extfile_h;

	switch (dtype) {
	case 1: {
		N = D.Nx;

		std::vector<double> X;
		std::vector<double> Y;
		std::vector<double> Zp;
		std::vector<double> Zh;

		std::fstream extfile_p;
		std::fstream extfile_h;

		fnamep = D.pth + "LC_PRS_DL.TXT";
		fnameh = D.pth + "LC_THK_DL.TXT";

		extfile_p.open(fnamep, std::ios::in);
		extfile_h.open(fnameh, std::ios::in);

		for (int i = 0; i < N; ++i) {

			extfile_p >> par1 >> par2;

			X.push_back(par1);
			Zp.push_back(par2);

			extfile_h >> par1 >> par2;
			Zh.push_back(par2);

		}

		extfile_p.close();
		extfile_h.close();

		auto p1 = plot(X, Zp, "-o");
		p1->line_width(4);
		xlabel("X");
		ylabel("PRESSURE");

		hold(on);

		auto p2 = plot(X, Zh, "-o");
		p2->use_y2(true);
		p2->line_width(4);
		y2label("FILM THICKNESS");

		show();

		break;
	} // case 1

	case 2: {
		N = D.Nx;
		std::vector<std::vector<double>> X(N);
		std::vector<std::vector<double>> Y(N);
		std::vector<std::vector<double>> Zp(N);
		std::vector<std::vector<double>> Zh(N);

		std::fstream extfile_p;
		std::fstream extfile_h;

		fnamep = D.pth + "PC_PRS_DL.TXT";
		fnameh = D.pth + "PC_THK_DL.TXT";

		extfile_p.open(fnamep, std::ios::in);
		extfile_h.open(fnameh, std::ios::in);

		for (i = 0; i < N; ++i)
			for (j = 0; j < N; ++j) {
				
				extfile_p >> par1 >> par2 >> par3;

				X[i].push_back(par1);
				Y[i].push_back(par2);
				Zp[i].push_back(par3);

				extfile_h >> par1 >> par2 >> par3;
				Zh[i].push_back(par3);
			}

		extfile_p.close();
		extfile_h.close();


			//subplot(2, 1, 0);
			surf(X, Y, Zp)->face_alpha(1.0).edge_color("none").lighting(true).primary(0.8f);
			//view(115, 335);
			grid(off);
			box(off);
			xlabel("X");
			ylabel("Y");
			zlabel("PRESSURE");
			show();
						
			subplot(2, 1, 1);
			surf(X, Y, Zh)->face_alpha(1.0).edge_color("none").lighting(true).primary(0.8f);
			grid(off);
			box(off);
			view(320, 60);
			xlabel("X");
			ylabel("Y");
			zlabel("FILM THICKNESS");
			
			show();

		break;
	} // case 2
	} // switch
}// proc
*/


void output_data(Domain D, Body B1, Body B2, int enttype, int entnum, char fixvar, double item1, int kval, string fname)
{
	/********** OUTPUT DOMAIN DATA TO EXTERNAL FILE *********/
	/* D		- domain									*/
	/* enttype 	- number of information type:				*/
	/*				  0 - geometry;	1 - mechanic; 			*/
	/*				  2 - thermal;	3 - results;			*/
	/*				  4 - service information.				*/
	/* entnum 	- number of entity to write in file			*/
	/* fixvar	- name of fixed variable (X or Y)			*/
	/* item1	- value for fixed variable (X or Y)			*/
	/* kval		- key for dimension/dimensionless form		*/
	/* fname	- filename without extension (*.TXT)		*/
	/* Values for <entnum> and other information can be		*/
	/* found in reference materials to this code			*/
	/********************************************************/

	int		i, j, l, Nx, Ny, Nz, Nz_curr,c;
	double	z1_curr, z2_curr, dx, dy, a, b, a1, b1, YMe1, YMe2, ph, RX, RY, Eeqv;
	double	c1, c2, c3, c4, h, d1, d2, C;
	double	rhoR, muR, tauR, etaR, alpha_HD, mu0_HD, um, vm;
	double	tolX, tolY;
	double W_HD, G_HD, U_HD, M_M, L_M, g_H, g_V, g_E, k, THK;
	string	pstfx;
	int jend;


	Nx = D.Nx;
	Ny = D.Ny;
	dx = D.dx;
	dy = D.dy;
	a = D.a;
	b = D.b;
	a1 = 1 * D.kx;
	b1 = 1 * D.ky;

	RX = pow(1 / B1.Rx + 1 / B2.Rx, -1);
	RY = pow(1 / B1.Ry + 1 / B2.Ry, -1);
	Eeqv = 1 / ((1 - B1.PR * B1.PR) / B1.YM + (1 - B2.PR * B2.PR) / B2.YM);

	YMe1 = B1.YMe;
	YMe2 = B2.YMe;

	ph = D.ph; /* Pa */
	rhoR = D.rhoR;
	muR = D.muR;
	etaR = D.etaR;

	tolX = dx / 2; /* tolerance for comparing item1 and fixed X (YZ plane results) */
	tolY = dy / 2; /* tolerance for comparing item1 and fixed Y (XZ plane results) */

	std::fstream extfile;
	std::string	 stritem1 = std::to_string(item1);
	stritem1.erase(5); /* 5 symbols for "pos" remains */

	if (item1 != 1) 
		pstfx = "_pos=" + stritem1; 
	else pstfx = "";

	extfile.open(fname + pstfx + ".TXT", std::ios::out);


	if ((enttype == 0) && (entnum == 0)) // general model information 

		{
			extfile << "========================================================" << "\n";
			extfile << "======== TEHL MODEL. GENERAL OUTPUT INFORMATION ========" << "\n";
			extfile << "========================================================" << "\n\n";

			extfile << "	TYPE OF CALCULATED EHL CONTACT:	";

			if (D.dkey == 1) {
				extfile << "\tLINE CONTACT";
				extfile << "\n\n";
				extfile << "********************************************************" << "\n";
				extfile << "************ HERZIAN CONTACT PARAMETERS ****************" << "\n";
				extfile << "********************************************************" << "\n\n";
				extfile << "	CONTACT LENGTH:\n";
				extfile << "		a\t\t\t=\t" << std::format("{:.5e}", D.a) << "\t\tm" << "\n\n";
				extfile << "	NORMAL FORCE (PER UNIT LENGTH):\n";
				extfile << "		W\t\t\t=\t" << std::format("{:.2f}", D.W) << "\t\t\tN/m" << "\n\n";
				extfile << "	MAX PRESSURE:\n";
				extfile << "		PH\t\t\t=\t" << std::format("{:.5e}", D.ph) << "\t\tPa" << "\n\n";
			}

			else {
				extfile << "\tPOINT CONTACT";
				extfile << "\n\n";
				extfile << "********************************************************" << "\n";
				extfile << "************ HERZIAN CONTACT PARAMETERS ****************" << "\n";
				extfile << "********************************************************" << "\n\n";
				extfile << "	SEMI-AXES OF CONTACT ELLIPSE:\n";
				extfile << "		A\t\t\t=\t" << std::format("{:.5e}", D.a) << "\t\tm" << "\n";
				extfile << "		B\t\t\t=\t" << std::format("{:.5e}", D.b) << "\t\tm" << "\n";
				extfile << "		RATIO\t\t=\t" << std::format("{:.5f}", D.a / D.b) << "\n\n";
				extfile << "	NORMAL CONTACT FORCE:\n";
				extfile << "		FN\t\t\t=\t" << std::format("{:.2f}", D.W) << "\t\t\tN" << "\n\n";
				extfile << "	MAX CONTACT PRESSURE:\n";
				extfile << "		PH\t\t\t=\t" << std::format("{:.5e}", D.ph) << "\t\tPa" << "\n\n";
			}

			extfile << "********************************************************" << "\n";
			extfile << "************ COORDINATE SYSTEM DESCRIPTION: ************" << "\n";
			extfile << "********************************************************" << "\n\n";
			extfile << "	THE CS ORIGIN IS PLACED ON BODY 1 (BOTTOM)" << "\n\n";
			extfile << "	X - AXIS: HORIZONTAL LINE" << "\n";
			extfile << "	Z - AXIS: VERTICAL LINE" << "\n";

			if (D.dkey == 2) extfile << "	Y - AXIS: AWAY FROM US" << "\n\n";
			else extfile << "\n";

			extfile << "********************************************************" << "\n";
			extfile << "******************* DOMAIN GEOMETRY ********************" << "\n";
			extfile << "********************************************************" << "\n\n";
			extfile << "	AREA BOUNDARIES:	" << "\n";
			extfile << "		X = [" << -a1 << " , " << a1 << "] | ";
			extfile << "[" << std::format("{:.5e}", -a1 * a) << " , " << std::format("{:.5e}", a1 * a) << "], m" << "\n";

			if (D.dkey == 2) {
				extfile << "		Y = [" << -b1 << " , " << b1 << "] | ";
				extfile << "[" << std::format("{:.5e}", -b1 * b) << " , " << std::format("{:.5e}", b1 * b) << "], m" << "\n\n";
			}
			else extfile << "\n";

			extfile << "	SCALE FACTORS TO ENLARGE COMPUTATIONAL AREA:	" << "\n";
			extfile << "		KX = " << std::format("{:.2f}", D.kx) << "\n";

			if (D.dkey == 2) extfile << "		KY = " << std::format("{:.2f}", D.ky) << "\n\n";
			else extfile << "\n";

			extfile << "	THE NUMBER OF GRID POINTS IN CONTACT AREA:" << "\n";
			extfile << "		NX = " << D.Nx << "\n";

			if (D.dkey == 2) extfile << "		NY = " << D.Ny << "\n\n";
			else extfile << "\n";

			extfile << "	THE DISTANCE BETWEEN GRID POINTS:" << "\n";
			extfile << "		DX = " << std::format("{:.5f}", D.dx) << " | ";
			extfile << std::format("{:.5e}", D.dx * a) << " m" << "\n";

			if (D.dkey == 2) {
				extfile << "		DY = " << std::format("{:.5f}", D.dy) << " | ";
				extfile << std::format("{:.5e}", D.dy * a) << " m" << "\n\n";
			}
			else extfile << "\n";

			extfile << "	SUBAREA WITH REFINED Z-RESOLUTION:	" << "\n";
			extfile << "		X = [" << D.RF.x1 << " , " << D.RF.x2 << "] | ";
			extfile << "[" << std::format("{:.5e}", D.RF.x1 * a) << " , " << std::format("{:.5e}", D.RF.x2 * a) << "], m" << "\n";

			if (D.dkey == 2) {
				extfile << "		Y = [" << D.RF.y1 << " , " << D.RF.y2 << "] | ";
				extfile << "[" << std::format("{:.5e}", D.RF.y1 * a) << " , " << std::format("{:.5e}", D.RF.y2 * a) << "], m" << "\n\n";
			}
			else extfile << "\n";

			extfile << "	Z-RESOLUTION IN DOMAIN:	" << "\n";
			extfile << "		NZ_NORMAL	=\t" << D.Nz0 << "\n";
			extfile << "		NZ_REFINED	=\t" << D.RF.val << "\n\n";

			extfile << "********************************************************" << "\n";
			extfile << "******************* BODY 1 (BOTTOM) ********************" << "\n";
			extfile << "********************************************************" << "\n\n";
			extfile << "	RADII OF CURVATURE:" << "\n";
			extfile << "		RX\t\t\t\t\t\t=\t" << std::format("{:.3f}", B1.Rx) << "\t\tm" << "\n";

			if (D.dkey == 2) {
				extfile << "		RY\t\t\t\t\t\t=\t" << std::format("{:.3f}", B1.Ry) << "\t\tm" << "\n\n";
			}
			else extfile << "\n";

			extfile << "	ELASTIC PROPERTIES:" << "\n";
			extfile << "		Young modulus\t\t\t=\t" << std::format("{:.2e}", B1.YM) << "\tPa" << "\n"
				<< "		Poisson ratio\t\t\t=\t" << std::format("{:.3f}", B1.PR) << "\n\n";
			extfile << "	THERMAL PROPERTIES:" << "\n";
			extfile << "		Thermal conductivity\t=\t" << std::format("{:.2f}", B1.K) << "\t\tW/(m K) " << "\n"
				<< "		Heat capacity\t\t\t=\t" << std::format("{:.2f}", B1.Cp) << "\t\tJ/(kg K)" << "\n\n";

			extfile << "********************************************************" << "\n";
			extfile << "********************* BODY 2 (TOP) *********************" << "\n";
			extfile << "********************************************************" << "\n\n";
			extfile << "	RADII OF CURVATURE:" << "\n";
			extfile << "		RX\t\t\t\t\t\t=\t" << std::format("{:.3f}", B2.Rx) << "\t\tm" << "\n";

			if (D.dkey == 2) {
				extfile << "		RY\t\t\t\t\t\t=\t" << std::format("{:.3f}", B1.Ry) << "\t\tm" << "\n\n";
			}
			else extfile << "\n";

			extfile << "	ELASTIC PROPERTIES:" << "\n";
			extfile << "		Young modulus\t\t\t=\t" << std::format("{:.2e}", B2.YM) << "\tPa" << "\n"
				<< "		Poisson ratio\t\t\t=\t" << std::format("{:.3f}", B2.PR) << "\n\n";
			extfile << "	THERMAL PROPERTIES:" << "\n";
			extfile << "		Thermal conductivity\t=\t" << std::format("{:.2f}", B2.K) << "\t\tW/(m K) " << "\n"
				<< "		Heat capacity\t\t\t=\t" << std::format("{:.2f}", B2.Cp) << "\t\tJ/(kg K)" << "\n\n";

			extfile << "********************************************************" << "\n";
			extfile << "***************** REDUCED GEOMETRY *********************" << "\n";
			extfile << "********************************************************" << "\n\n";
			extfile << "	RADII OF CURVATURE:" << "\n";
			extfile << "		RX_REDU\t\t\t\t\t=\t" << std::format("{:.3f}", D.rRX) << "\t\tm" << "\n";
			
			if (D.dkey == 2) 
				extfile << "		RY_REDU\t\t\t\t\t=\t" << std::format("{:.3f}", D.rRY) << "\t\tm" << "\n";
			extfile << "\n";
			extfile << "	YOUNG MODULUS:" << "\n";
			extfile << "		YM_REDU\t\t\t\t\t=\t" << std::format("{:.2e}", D.rYm) << "\tPa" << "\n";
			extfile << "\n";
			
			extfile << "********************************************************" << "\n";
			extfile << "******************** BODIES KINEMATIC ******************" << "\n";
			extfile << "********************************************************" << "\n\n";
			extfile << "	BODY 1 VELOCITY:" << "\n";
			extfile << "		U (X-COMP)\t\t\t\t\t=\t" << std::format("{:.2f}", B1.vel[0]) << "\t\tm/s" << "\n";

			if (D.dkey == 2) {
				extfile << "		V (Y-COMP)\t\t\t\t\t=\t" << std::format("{:.2f}", B1.vel[1]) << "\t\tm/s" << "\n\n";

			}
			else extfile << "\n";
			extfile << "	BODY 2 VELOCITY:" << "\n";
			extfile << "		U (X-COMP)\t\t\t\t\t=\t" << std::format("{:.2f}", B2.vel[0]) << "\t\tm/s" << "\n";

			if (D.dkey == 2) {
				extfile << "		V (Y-COMP)\t\t\t\t\t=\t" << std::format("{:.2f}", B2.vel[1]) << "\t\tm/s" << "\n\n";

			}
			else extfile << "\n";
			extfile << "	SLIDE-TO-ROLL RATIO:" << "\n";
			extfile << "		SRR_X\t\t\t\t\t\t=\t" << std::format("{:.2f}", D.SRR[0]) << "\n";
			if (D.dkey == 2)
				extfile << "		SRR_Y\t\t\t\t\t\t=\t" << std::format("{:.2f}", D.SRR[1]) << "\n";
			extfile << "\n";
						
			extfile << "********************************************************" << "\n";
			extfile << "***************** LUBRICANT PROPERTIES *****************" << "\n";
			extfile << "********************************************************" << "\n\n";

			extfile << "	DENSITY (RHO)" << "\n";
			extfile << "		MODEL NAME:\t\t\t\t" << D.RHO.IT1 << "\n";
			extfile << "		MODEL PARAMETERS:" << "\n";
			extfile << "\t\t\t\t" << D.RHO.narg1 << "\t=\t" << std::format("{:.2f}", D.RHO.arg1) << "\n";

			if (!D.RHO.narg2.empty())
				extfile << "\t\t\t\t" << D.RHO.narg2 << "\t=\t" << std::format("{:.5e}", D.RHO.arg2) << "\n";
			if (!D.RHO.narg3.empty())
				extfile << "\t\t\t\t" << D.RHO.narg3 << "\t=\t" << std::format("{:.5e}", D.RHO.arg3) << "\n";
			if (!D.RHO.narg4.empty())
				extfile << "\t\t\t\t" << D.RHO.narg4 << "\t=\t" << std::format("{:.5e}", D.RHO.arg4) << "\n";
			if (!D.RHO.narg5.empty())
				extfile << "\t\t\t\t" << D.RHO.narg5 << "\t=\t" << std::format("{:.5e}", D.RHO.arg5) << "\n";
			if (!D.RHO.narg6.empty())
				extfile << "\t\t\t\t" << D.RHO.narg6 << "\t=\t" << std::format("{:.5e}", D.RHO.arg6) << "\n";

			extfile << "\n";
			extfile << "\t\tREFERENCE VALUE:" << "\n";
			extfile << "\t\t\t\tRHO_R\t\t=\t" << std::format("{:.2f}", D.rhoR) << "\tkg/(m^3) " << "\n\n\n";

			extfile << "	LOW-SHEAR VISCOSITY (MU)" << "\n";
			extfile << "		MODEL NAME:\t\t\t\t" << D.LSV.IT1 << "\n";
			extfile << "		MODEL PARAMETERS:" << "\n";
			extfile << "\t\t\t\t" << D.LSV.narg1 << "\t=\t" << std::format("{:.5e}", D.LSV.arg1) << "\n";

			if (!D.LSV.narg2.empty())
				extfile << "\t\t\t\t" << D.LSV.narg2 << "\t=\t" << std::format("{:.5e}", D.LSV.arg2) << "\n";
			if (!D.LSV.narg3.empty())
				extfile << "\t\t\t\t" << D.LSV.narg3 << "\t=\t" << std::format("{:.5e}", D.LSV.arg3) << "\n";
			if (!D.LSV.narg4.empty())
				extfile << "\t\t\t\t" << D.LSV.narg4 << "\t=\t" << std::format("{:.5e}", D.LSV.arg4) << "\n";
			if (!D.LSV.narg5.empty())
				extfile << "\t\t\t\t" << D.LSV.narg5 << "\t=\t" << std::format("{:.5e}", D.LSV.arg5) << "\n";
			if (!D.LSV.narg6.empty())
				extfile << "\t\t\t\t" << D.LSV.narg6 << "\t=\t" << std::format("{:.5e}", D.LSV.arg6) << "\n";
			if (!D.LSV.narg7.empty())
				extfile << "\t\t\t\t" << D.LSV.narg7 << "\t=\t" << std::format("{:.5e}", D.LSV.arg7) << "\n";
			if (!D.LSV.narg8.empty())
				extfile << "\t\t\t\t" << D.LSV.narg8 << "\t=\t" << std::format("{:.5e}", D.LSV.arg8) << "\n\n";

			extfile << "\n";
			extfile << "\t\tREFERENCE VALUE:" << "\n";
			extfile << "\t\t\t\tMU_R\t\t=\t" << std::format("{:.5e}", D.muR) << "\tPa s " << "\n\n\n";

			extfile << "	GENERALIZED VISCOSITY (ETA)" << "\n";
			extfile << "		MODEL NAME:\t\t\t\t" << D.GNV.IT1 << "\n";
			extfile << "		MODEL PARAMETERS:" << "\n";
			if (!D.GNV.narg1.empty())
				extfile << "\t\t\t\t" << D.GNV.narg1 << "\t=\t" << std::format("{:.5e}", D.GNV.arg1) << "\n";

			if (!D.GNV.narg2.empty())
				extfile << "\t\t\t\t" << D.GNV.narg2 << "\t=\t" << std::format("{:.5e}", D.GNV.arg2) << "\n";
			if (!D.GNV.narg3.empty())
				extfile << "\t\t\t\t" << D.GNV.narg3 << "\t=\t" << std::format("{:.5e}", D.GNV.arg3) << "\n";
			if (!D.GNV.narg4.empty())
				extfile << "\t\t\t\t" << D.GNV.narg4 << "\t=\t" << std::format("{:.5e}", D.GNV.arg4) << "\n";
			if (!D.GNV.narg5.empty())
				extfile << "\t\t\t\t" << D.GNV.narg5 << "\t=\t" << std::format("{:.5e}", D.GNV.arg5) << "\n";
			if (!D.GNV.narg6.empty())
				extfile << "\t\t\t\t" << D.GNV.narg6 << "\t=\t" << std::format("{:.5e}", D.GNV.arg6) << "\n";

			extfile << "\n";
			extfile << "\t\tREFERENCE VALUE:" << "\n";
			extfile << "\t\t\t\tETA_R\t\t=\t" << std::format("{:.5e}", D.etaR) << "\tPa s " << "\n\n\n";

			extfile << "	THERMAL CONDUCTIVITY (K)" << "\n";
			extfile << "		MODEL NAME:\t\t\t\t" << D.THC.IT1 << "\n";
			extfile << "		MODEL PARAMETERS:" << "\n";
			extfile << "\t\t\t\t" << D.THC.narg1 << "\t=\t" << std::format("{:.5e}", D.THC.arg1) << "\n";

			if (!D.THC.narg2.empty())
				extfile << "\t\t\t\t" << D.THC.narg2 << "\t=\t" << std::format("{:.5e}", D.THC.arg2) << "\n";
			if (!D.THC.narg3.empty())
				extfile << "\t\t\t\t" << D.THC.narg3 << "\t=\t" << std::format("{:.5e}", D.THC.arg3) << "\n";
			if (!D.THC.narg4.empty())
				extfile << "\t\t\t\t" << D.THC.narg4 << "\t=\t" << std::format("{:.5e}", D.THC.arg4) << "\n";
			if (!D.THC.narg5.empty())
				extfile << "\t\t\t\t" << D.THC.narg5 << "\t=\t" << std::format("{:.5e}", D.THC.arg5) << "\n";

			extfile << "\n";
			extfile << "\t\tREFERENCE VALUE:" << "\n";
			extfile << "\t\t\t\tK_R\t\t\t=\t" << std::format("{:.5e}", D.KR) << "\tW/(m K) " << "\n\n\n";

			extfile << "	HEAT CAPACITY (RHO*Cp)" << "\n";
			extfile << "		MODEL NAME:\t\t\t\t" << D.HCP.IT1 << "\n";
			extfile << "		MODEL PARAMETERS:" << "\n";
			extfile << "\t\t\t\t" << D.HCP.narg1 << "\t=\t" << std::format("{:.5e}", D.HCP.arg1) << "\n";

			if (!D.HCP.narg2.empty())
				extfile << "\t\t\t\t" << D.HCP.narg2 << "\t=\t" << std::format("{:.5e}", D.HCP.arg2) << "\n";
			if (!D.HCP.narg3.empty())
				extfile << "\t\t\t\t" << D.HCP.narg3 << "\t=\t" << std::format("{:.5e}", D.HCP.arg3) << "\n";
			if (!D.HCP.narg4.empty())
				extfile << "\t\t\t\t" << D.HCP.narg4 << "\t=\t" << std::format("{:.5e}", D.HCP.arg4) << "\n";
			if (!D.HCP.narg5.empty())
				extfile << "\t\t\t\t" << D.HCP.narg5 << "\t=\t" << std::format("{:.5e}", D.HCP.arg5) << "\n";
			if (!D.HCP.narg6.empty())
				extfile << "\t\t\t\t" << D.HCP.narg6 << "\t=\t" << std::format("{:.5e}", D.HCP.arg6) << "\n";
			if (!D.HCP.narg7.empty())
				extfile << "\t\t\t\t" << D.HCP.narg7 << "\t=\t" << std::format("{:.5e}", D.HCP.arg7) << "\n";

			extfile << "\n";
			extfile << "\t\tREFERENCE VALUE:" << "\n";
			extfile << "\t\t\t\tRHO*Cp_R\t=\t" << std::format("{:.5e}", D.CpR) << "\tJ/(K m^3) " << "\n\n";

			/* set alpha for Moes parameters */
			switch (D.LSV.mkey) {
			case 0: /* constant */
				alpha_HD = 0.0;
				break;
			case 1: /* Barus */
				alpha_HD = D.LSV.arg2;
				break;
			case 2: /* Roeland */
				alpha_HD = 1e-8;
				break;
			case 3: /* J & W */
				alpha_HD = 1e-8;
				break;
			case 4: /* Doolitle */
				alpha_HD = 1e-8;
				break;
			case 5: /* Yasutomi (modified) */
				alpha_HD = 1e-8;
				break;
			}

			um = 0.5 * (B1.vel[0] + B2.vel[0]);
			vm = 0.5 * (B1.vel[1] + B2.vel[1]);

			c = static_cast<int>((Nx - 1) * 0.5);

			if (D.dkey == 2) { /* dimensionless groups for point contact */

				/* dimensionless thickness in center of contact area */
				THK = D.Area[c][c].h * (RX / (D.a * D.a));
				
				// according to Habchi, 2018
				W_HD = D.W / (2 * Eeqv * RX * RX);
				G_HD = 2 * alpha_HD * Eeqv;
				U_HD = D.muR * pow(um * um + vm * vm, 0.5) / (2 * Eeqv * RX);
				M_M = W_HD * pow((2 * U_HD), -0.75); 
				L_M = G_HD * pow((2 * U_HD), 0.25);

				g_H = pow(W_HD / U_HD, 2) * THK;
				g_V = pow(G_HD * pow(W_HD, 3) / (U_HD * U_HD), 2);
				g_E = pow(W_HD, (8.0 / 3.0)) * U_HD * U_HD;
				k = 1.0339 * pow(RY / RX, 0.636);

			}

			else { /* dimensionless groups for line contact */

				/* dimensionless thickness in center of contact area */
				THK = D.Area[c][0].h * (RX / (D.a * D.a));

				W_HD = D.W / (2 * Eeqv * RX);
				G_HD = 2 * alpha_HD * Eeqv;
				U_HD = D.muR * um / (2 * Eeqv * RX);
				M_M = W_HD * pow((2 * U_HD), -0.5);
				L_M = G_HD * pow((2 * U_HD), 0.25);

				g_H = (W_HD / U_HD) * THK;
				g_V = G_HD * pow(W_HD, 1.5) / pow(U_HD, 0.5);
				g_E = W_HD / pow(U_HD, 0.5);

			}

			extfile << "********************************************************" << "\n";
			extfile << "***************** DIMENSIONLESS GROUPS *****************" << "\n";
			extfile << "********************************************************" << "\n\n";
			extfile << "	HAMROCK & DOWSON:\n";
			extfile << "		W\t\t\t=\t" << std::format("{:.5e}", W_HD) << "\n";
			extfile << "		G\t\t\t=\t" << std::format("{:.5e}", G_HD) << "\n";
			extfile << "		U\t\t\t=\t" << std::format("{:.5e}", U_HD) << "\n\n";
			extfile << "	MOES:\n";
			extfile << "		M\t\t\t=\t" << std::format("{:.1f}", M_M) << "\n";
			extfile << "		L\t\t\t=\t" << std::format("{:.1f}", L_M) << "\n\n";
			extfile << "	G-PARAMETERS:\n";
			extfile << "		G_H\t\t\t=\t" << std::format("{:.5e}", g_H) << "\n";
			extfile << "		G_V\t\t\t=\t" << std::format("{:.5e}", g_V) << "\n";
			extfile << "		G_E\t\t\t=\t" << std::format("{:.5e}", g_E) << "\n\n";

			extfile << "********************************************************" << "\n";
			extfile << "******* FILM THICKNESSES BY EMPIRICAL FORMULAS  ********" << "\n";
			extfile << "********************************************************" << "\n\n";
			extfile << "	HAMROCK & DOWSON:\n";
			extfile << "		H_C\t\t\t=\t" << std::format("{:.5e}", D.Hc) <<  "(" <<
				std::format("{:.5e}", D.Hc * (a * a) / RX) << ",m)" << "\n";
			extfile << "		H_M\t\t\t=\t" << std::format("{:.5e}", D.Hm) << "(" <<
				std::format("{:.5e}", D.Hm * (a * a) / RX) << ",m)" << "\n";
			


		} /* GGI out */
		
	else  
		if ((enttype == 2) && (entnum == 1)) { // Friction results (integrated)
		
			extfile << "FRICTION ON BODY SURFACE:\n\n";
			extfile << "\tBODY 1 (BOTTOM)\n";
			extfile << "\t\tTANGENTIAL FORCE\t\t\t\t=\t" << std::format("{:.5e}", B1.Ft) << "\tN\n";
			extfile << "\t\tFRICTION COEFFICIENT\t\t\t=\t" << std::format("{:.5e}", B1.fc) << "\n";
			extfile << "\n\n";
			extfile << "\tBODY 2 (TOP)\n";
			extfile << "\t\tTANGENTIAL FORCE\t\t\t\t=\t" << std::format("{:.5e}", B2.Ft) << "\tN\n";
			extfile << "\t\tFRICTION COEFFICIENT\t\t\t=\t" << std::format("{:.5e}", B2.fc) << "\n";
		
		}
		
		else { //distributed results

		jend = 1; /* line contact is set by default */
		if (D.dkey == 2) jend = Ny;

		for (i = 0; i < Nx; ++i)
			for (j = 0; j < jend; ++j) {

				/* current point coordinates */
				/* remember that we set the position in item1 through X and Y for DL- and DM-output!! */
				c1 = -a1 + i * dx;		/* dimensionless X */
				c2 = -b1 + j * dy;		/* dimensionless Y */

				//h = (D.Area[i][j].z2 - D.Area[i][j].z1) + (D.Area[i][j].d2 + D.Area[i][j].d1) + D.h0;	/* m */

				h = D.Area[i][j].h;
;				z1_curr = (D.Area[i][j].z1 - D.Area[i][j].d1) / h;			/* start Z-coord */
				z2_curr = (D.Area[i][j].z2 + D.Area[i][j].d2 + D.h0) / h;	/* end Z-coord */

				Nz_curr = D.Area[i][j].Nz;

				if (enttype == 0) {

					switch (entnum) {

					case 1: /* B1 surface (nominal) */

						switch (kval) {

						case 0: /* DL-values */

							if (fixvar == '-')
								if (D.dkey == 2) extfile << std::format("{:.5e}", c1) << "\t" << std::format("{:.5e}", c2) << "\t" << std::format("{:.5e}", D.Area[i][j].z1_n / h) << "\n";
							     else extfile << std::format("{:.5e}", c1) << "\t" << std::format("{:.5e}", D.Area[i][j].z1_n / h) << "\n";
							else if (fixvar == 'X') {
								if (abs(c1 - item1) <= tolX)
									extfile << std::format("{:.5e}", c2) << "\t" << std::format("{:.5e}", D.Area[i][j].z1_n / h) << "\n";
							}
							else if (fixvar == 'Y') {
								if (abs(c2 - item1) <= tolY)
									extfile << std::format("{:.5e}", c1) << "\t" << std::format("{:.5e}", D.Area[i][j].z1_n / h) << "\n";
							}
							break;

						case 1: /* DM-values */

							if (fixvar == '-')
								if (D.dkey == 2) extfile << std::format("{:.5e}", c1 * a) << "\t" << std::format("{:.5e}", c2 * b) << "\t" << std::format("{:.5e}", D.Area[i][j].z1_n) << "\n";
								else extfile << std::format("{:.5e}", c1 * a) << "\t" << std::format("{:.5e}", D.Area[i][j].z1_n) << "\n";
								else if (fixvar == 'X') {
								if (abs(c1 - item1) <= tolX)
									extfile << std::format("{:.5e}", c2 * b) << "\t" << std::format("{:.5e}", D.Area[i][j].z1_n) << "\n";
							}
							else if (fixvar == 'Y') {
								if (abs(c2 - item1) <= tolY)
									extfile << std::format("{:.5e}", c1 * a) << "\t" << std::format("{:.5e}", D.Area[i][j].z1_n) << "\n";
							}
							break;
						} /* kval switch */
						break;

					case 2: /* B2 surface (nominal) */

						switch (kval) {

						case 0: /* DL-values */

							if (fixvar == '-')
								if (D.dkey == 2) extfile << std::format("{:.5e}", c1) << "\t" << std::format("{:.5e}", c2) << "\t" << std::format("{:.5e}", D.Area[i][j].z2_n / h) << "\n";
								else extfile << std::format("{:.5e}", c1) << "\t" << std::format("{:.5e}", D.Area[i][j].z2_n / h) << "\n";
								else if (fixvar == 'X') {
								if (abs(c1 - item1) <= tolX)
									extfile << std::format("{:.5e}", c2) << "\t" << std::format("{:.5e}", D.Area[i][j].z2_n / h) << "\n";
							}
							else if (fixvar == 'Y') {
								if (abs(c2 - item1) <= tolY)
									extfile << std::format("{:.5e}", c1) << "\t" << std::format("{:.5e}", D.Area[i][j].z2_n / h) << "\n";
							}
							break;

						case 1: /* DM-values */

							if (fixvar == '-')
								if (D.dkey == 2) extfile << std::format("{:.5e}", c1 * a) << "\t" << std::format("{:.5e}", c2 * b) << "\t" << std::format("{:.5e}", D.Area[i][j].z2_n) << "\n";
								else extfile << std::format("{:.5e}", c1 * a) << "\t" << std::format("{:.5e}", D.Area[i][j].z2_n) << "\n";
								else if (fixvar == 'X') {
								if (abs(c1 - item1) <= tolX)
									extfile << std::format("{:.5e}", c2 * b) << "\t" << std::format("{:.5e}", D.Area[i][j].z2_n) << "\n";
							}
							else if (fixvar == 'Y') {
								if (abs(c2 - item1) <= tolY)
									extfile << std::format("{:.5e}", c1 * a) << "\t" << std::format("{:.5e}", D.Area[i][j].z2_n) << "\n";
							}
							break;

						} /* kval switch */
						break;

					case 3: /* B2 surface (reduced) */

						switch (kval) {

						case 0: /* DL-values */

							if (fixvar == '-')
								if (D.dkey == 2) extfile << std::format("{:.5e}", c1) << "\t" << std::format("{:.5e}", c2) << "\t" << std::format("{:.5e}", D.Area[i][j].z2 / h) << "\n";
								else extfile << std::format("{:.5e}", c1) << "\t" << std::format("{:.5e}", D.Area[i][j].z2 / h) << "\n";
							else if (fixvar == 'X') {
								if (abs(c1 - item1) <= tolX)
									extfile << std::format("{:.5e}", c2) << "\t" << std::format("{:.5e}", D.Area[i][j].z2 / h) << "\n";
							}
							else if (fixvar == 'Y') {
								if (abs(c2 - item1) <= tolY)
									extfile << std::format("{:.5e}", c1) << "\t" << std::format("{:.5e}", D.Area[i][j].z2 / h) << "\n";
							}
							break;

						case 1: /* DM-values */

							if (fixvar == '-')
								if (D.dkey == 2) extfile << std::format("{:.5e}", c1 * a) << "\t" << std::format("{:.5e}", c2 * b) << "\t" << std::format("{:.5e}", D.Area[i][j].z2) << "\n";
								else extfile << std::format("{:.5e}", c1 * a) << "\t" << std::format("{:.5e}", D.Area[i][j].z2) << "\n";
							else if (fixvar == 'X') {
								if (abs(c1 - item1) <= tolX)
									extfile << std::format("{:.5e}", c2 * b) << "\t" << std::format("{:.5e}", D.Area[i][j].z2) << "\n";
							}
							else if (fixvar == 'Y') {
								if (abs(c2 - item1) <= tolY)
									extfile << std::format("{:.5e}", c1 * a) << "\t" << std::format("{:.5e}", D.Area[i][j].z2) << "\n";
							}
							break;

						} /* kval switch */
						break;

					case 4: /* Z resolutions */
					{
						if (D.dkey == 2) extfile << std::format("{:.5e}", c1) << "\t" << std::format("{:.5e}", c2) << "\t"
							<< D.Area[i][j].Nz << "\n";
						else extfile << std::format("{:.5e}", c1) << "\t" << D.Area[i][j].Nz << "\n";
						break;
					}

					case 5: /* Thickness (h0 and deformations are included) */

						switch (kval) {

						case 0: /* DL-values */

							if (fixvar == '-')
								if (D.dkey == 2) extfile << std::format("{:.5e}", c1) << "\t" << std::format("{:.5e}", c2) << "\t"
									<< std::format("{:.5e}", (RX / (a * a)) * h) << "\n";
									else extfile << std::format("{:.5e}", c1) << "\t" << std::format("{:.5e}", (RX / (a * a)) * h) << "\n";
							else if (fixvar == 'X') {
								if (abs(c1 - item1) <= tolX)
									extfile << std::format("{:.5e}", c2) << "\t"
									<< std::format("{:.5e}", (RX / (a * a)) * h) << "\n";
							}
							else if (fixvar == 'Y') {
								if (abs(c2 - item1) <= tolY)
									extfile << std::format("{:.5e}", c1) << "\t"
									<< std::format("{:.5e}", (RX / (a * a)) * h) << "\n";
							}
							break;

						case 1: /* DM-values */

							if (fixvar == '-')
								if (D.dkey == 2) extfile << std::format("{:.5e}", c1 * a) << "\t" << std::format("{:.5e}", c2 * b) << "\t"
								<< std::format("{:.5e}", h) << "\n";
								else extfile << std::format("{:.5e}", c1 * a) << "\t" << std::format("{:.5e}", h) << "\n";
							else if (fixvar == 'X') {
								if (abs(c1 - item1) <= tolX)
									extfile << std::format("{:.5e}", c2 * b) << "\t"
									<< std::format("{:.5e}", h) << "\n";
							}
							else if (fixvar == 'Y') {
								if (abs(c2 - item1) <= tolY)
									extfile << std::format("{:.5e}", c1 * a) << "\t"
									<< std::format("{:.5e}", h) << "\n";
							}
							break;
						} /* kval switch */
						break;

					case 6: /* Deformations - B1 (reduced geometry) */

						switch (kval) {

						case 0: /* DL-values */

							if (fixvar == '-') {
								if (D.dkey == 2) extfile << std::format("{:.5e}", c1) << "\t" << std::format("{:.5e}", c2) << "\t"
									<< std::format("{:.5e}", (RX / (a * a)) * D.Area[i][j].d1) << "\n";
								else extfile << std::format("{:.5e}", c1) << "\t" << std::format("{:.5e}", (RX / (a * a)) * D.Area[i][j].d1) << "\n";
							}
							else if (fixvar == 'X') {
								if (abs(c1 - item1) <= tolX)
									extfile << std::format("{:.5e}", c2) << "\t"
									<< std::format("{:.5e}", (RX / (a * a)) * D.Area[i][j].d1) << "\n";
							}
							else if (fixvar == 'Y') {
								if (abs(c2 - item1) <= tolY)
									extfile << std::format("{:.5e}", c1) << "\t"
									<< std::format("{:.5e}", (RX / (a * a)) * D.Area[i][j].d1) << "\n";
							}
							break;

						case 1: /* DM-values */

							if (fixvar == '-') {
								if (D.dkey == 2) extfile << std::format("{:.5e}", c1 * a) << "\t" << std::format("{:.5e}", c2 * b) << "\t"
									<< std::format("{:.5e}", D.Area[i][j].d1) << "\n";
								else extfile << std::format("{:.5e}", c1 * a) << "\t" << std::format("{:.5e}", D.Area[i][j].d1) << "\n";
							}
							else if (fixvar == 'X') {
								if (abs(c1 - item1) <= tolX)
									extfile << std::format("{:.5e}", c2 * b) << "\t"
									<< std::format("{:.5e}", D.Area[i][j].d1) << "\n";
							}
							else if (fixvar == 'Y') {
								if (abs(c2 - item1) <= tolY)
									extfile << std::format("{:.5e}", c1 * a) << "\t"
									<< std::format("{:.5e}", D.Area[i][j].d1) << "\n";
							}
							break;
						} /* kval switch */
						break;

					case 7: /* Deformations - B2 (reduced geometry) */

						switch (kval) {

						case 0: /* DL-values */

							if (fixvar == '-')
								if (D.dkey == 2) extfile << std::format("{:.5e}", c1) << "\t" << std::format("{:.5e}", c2) << "\t"
								<< std::format("{:.5e}", (RX / (a * a)) * D.Area[i][j].d2) << "\n";
								else extfile << std::format("{:.5e}", c1) << "\t" << std::format("{:.5e}", (RX / (a * a)) * D.Area[i][j].d2) << "\n";
							else if (fixvar == 'X') {
								if (abs(c1 - item1) <= tolX)
									extfile << std::format("{:.5e}", c2) << "\t"
									<< std::format("{:.5e}", (RX / (a * a)) * D.Area[i][j].d2) << "\n";
							}
							else if (fixvar == 'Y') {
								if (abs(c2 - item1) <= tolY)
									extfile << std::format("{:.5e}", c1) << "\t"
									<< std::format("{:.5e}", (RX / (a * a)) * D.Area[i][j].d2) << "\n";
							}
							break;

						case 1: /* DM-values */

							if (fixvar == '-')
								if (D.dkey == 2) extfile << std::format("{:.5e}", c1 * a) << "\t" << std::format("{:.5e}", c2 * b) << "\t"
								<< std::format("{:.5e}", D.Area[i][j].d2) << "\n";
								else extfile << std::format("{:.5e}", c1 * a) << "\t" << std::format("{:.5e}", D.Area[i][j].d2) << "\n";
							else if (fixvar == 'X') {
								if (abs(c1 - item1) <= tolX)
									extfile << std::format("{:.5e}", c2 * b) << "\t"
									<< std::format("{:.5e}", D.Area[i][j].d2) << "\n";
							}
							else if (fixvar == 'Y') {
								if (abs(c2 - item1) <= tolY)
									extfile << std::format("{:.5e}", c1 * a) << "\t"
									<< std::format("{:.5e}", D.Area[i][j].d2) << "\n";
							}
							break;
						} /* kval switch */
						break;
					} // entnum switch
				} // end if for (enttype)

				if (enttype == 1) {

					switch (entnum) {

					case 0: break; /* reserved */

					case 1: /* B1 real profile (deformations are included) */

						switch (kval) {

						case 0: /* DL-values */

							if (fixvar == '-')
								if (D.dkey == 2) extfile << std::format("{:.5e}", c1) << "\t" << std::format("{:.5e}", c2)
									<< "\t" << std::format("{:.5e}", (-1) * (RX / (a * a)) * D.Area[i][j].d1) << "\n";
								else extfile << std::format("{:.5e}", c1) << "\t" << std::format("{:.5e}", (-1) * (RX / (a * a)) * D.Area[i][j].d1) << "\n";
							else if (fixvar == 'X') {
								if (abs(c1 - item1) <= tolX)
									extfile << std::format("{:.5e}", c2) << "\t"
									<< std::format("{:.5e}", (-1) * (RX / (a * a)) * D.Area[i][j].d1) << "\n";
							}
							else if (fixvar == 'Y') {
								if (abs(c2 - item1) <= tolY)
									extfile << std::format("{:.5e}", c1) << "\t"
									<< std::format("{:.5e}", (-1) * (RX / (a * a)) * D.Area[i][j].d1) << "\n";
							}
							break;

						case 1: /* DM-values */

							if (fixvar == '-')
								if (D.dkey == 2) extfile << std::format("{:.5e}", c1 * a) << "\t" << std::format("{:.5e}", c2 * b) << "\t"
									<< std::format("{:.5e}", (-1) * D.Area[i][j].d1) << "\n";
								else extfile << std::format("{:.5e}", c1 * a) << "\t" << std::format("{:.5e}", (-1) * D.Area[i][j].d1) << "\n";
							else if (fixvar == 'X') {
								if (abs(c1 - item1) <= tolX)
									extfile << std::format("{:.5e}", c2 * b) << "\t"
									<< std::format("{:.5e}", (-1) * D.Area[i][j].d1) << "\n";
							}
							else if (fixvar == 'Y') {
								if (abs(c2 - item1) <= tolY)
									extfile << std::format("{:.5e}", c1 * a) << "\t"
									<< std::format("{:.5e}", (-1) * D.Area[i][j].d1) << "\n";
							}
							break;
						} /* kval switch */
						break; /* entnum case */

					case 2: /* B2 REAL PROFILE (h0 and deformations are included) */

						switch (kval) {

						case 0: /* DL-values */

							if (fixvar == '-')
								if (D.dkey == 2) extfile << std::format("{:.5e}", c1) << "\t" << std::format("{:.5e}", c2)
									<< "\t" << std::format("{:.5e}", z2_curr + (RX / (a * a)) * (D.Area[i][j].d2 + D.h0)) << "\n";
								else extfile << std::format("{:.5e}", c1) << "\t" << std::format("{:.5e}", z2_curr + (RX / (a * a)) * (D.Area[i][j].d2 + D.h0)) << "\n";
							else if (fixvar == 'X') {
								if (abs(c1 - item1) <= tolX)
									extfile << std::format("{:.5e}", c2) << "\t"
									<< std::format("{:.5e}", z2_curr + (RX / (a * a)) * (D.Area[i][j].d2 + D.h0)) << "\n";
							}
							else if (fixvar == 'Y') {
								if (abs(c2 - item1) <= tolY)
									extfile << std::format("{:.5e}", c1) << "\t"
									<< std::format("{:.5e}", z2_curr + (RX / (a * a)) * (D.Area[i][j].d2 + D.h0)) << "\n";
							}
							break;

						case 1: /* DM-values */

							if (fixvar == '-')
								if (D.dkey == 2) extfile << std::format("{:.5e}", c1 * a) << "\t" << std::format("{:.5e}", c2 * b) << "\t"
									<< std::format("{:.5e}", D.Area[i][j].z2 + D.Area[i][j].d2 + D.h0) << "\n";
								else extfile << std::format("{:.5e}", c1 * a) << "\t" << std::format("{:.5e}", D.Area[i][j].z2 + D.Area[i][j].d2 + D.h0) << "\n";
							else if (fixvar == 'X') {
								if (abs(c1 - item1) <= tolX)
									extfile << std::format("{:.5e}", c2 * b) << "\t"
									<< std::format("{:.5e}", D.Area[i][j].z2 + D.Area[i][j].d2 + D.h0) << "\n";
							}
							else if (fixvar == 'Y') {
								if (abs(c2 - item1) <= tolY)
									extfile << std::format("{:.5e}", c1 * a) << "\t"
									<< std::format("{:.5e}", D.Area[i][j].z2 + D.Area[i][j].d2 + D.h0) << "\n";
							}
							break;
						} /* kval switch */
						break; /* entnum case */

					case 3: /* PRESSURE */

						switch (kval) {

						case 0: /* DL-values */

							if (fixvar == '-') {
								if (D.dkey == 2) extfile << std::format("{:.5e}", c1) << "\t" << std::format("{:.5e}", c2) << "\t"
									<< std::format("{:.5e}", D.Area[i][j].p) << "\n";
								else extfile << std::format("{:.5e}", c1) << "\t" << std::format("{:.5e}", D.Area[i][j].p) << "\n";
							}
							else if (fixvar == 'X') {
								if (abs(c1 - item1) <= tolX)
									extfile << std::format("{:.5e}", c2) << "\t" <<
									std::format("{:.5e}", D.Area[i][j].p) << "\n";
							}
							else if (fixvar == 'Y') {
								if (abs(c2 - item1) <= tolY)
									extfile << std::format("{:.5e}", c1) << "\t"
									<< std::format("{:.5e}", D.Area[i][j].p) << "\n";
							}
							break;

						case 1: /* DM-values */

							if (fixvar == '-')
								if (D.dkey == 2) extfile << std::format("{:.5e}", c1 * a) << "\t" << std::format("{:.5e}", c2 * b) << "\t"
									<< std::format("{:.5e}", D.Area[i][j].p * ph) << "\n";
								else extfile << std::format("{:.5e}", c1 * a) << "\t" << std::format("{:.5e}", D.Area[i][j].p * ph) << "\n";
							else if (fixvar == 'X') {
								if (abs(c1 - item1) <= tolX)
									extfile << std::format("{:.5e}", c2 * b) << "\t"
									<< std::format("{:.5e}", D.Area[i][j].p * ph) << "\n";
							}
							else if (fixvar == 'Y') {
								if (abs(c2 - item1) <= tolY)
									extfile << std::format("{:.5e}", c1 * a) << "\t"
									<< std::format("{:.5e}", D.Area[i][j].p * ph) << "\n";
							}
							break;
						} /* kval switch */
						break;/* entnum case */

					case 4: /* DENSITY */

						switch (kval) {

						case 0: /* DL-values */ {

							for (l = 0; l < Nz_curr; ++l) {

								c3 = z1_curr + (1.0 * l / (Nz_curr - 1));

								if (fixvar == '-')
									if (D.dkey == 2) extfile << std::format("{:.5e}", c1) << "\t" << std::format("{:.5e}", c2) << "\t"
									<< std::format("{:.5e}", c3) << "\t" << std::format("{:.5e}", D.Area[i][j].rho[l] / rhoR) << "\n";
									else extfile << std::format("{:.5e}", c1) << "\t" << std::format("{:.5e}", c3) << "\t" 
										<< std::format("{:.5e}", D.Area[i][j].rho[l] / rhoR) << "\n";

								else if (fixvar == 'X') {
									if (abs(c1 - item1) <= tolX)
										extfile << std::format("{:.5e}", c2) << "\t" << std::format("{:.5e}", c3) << "\t"
										<< std::format("{:.5e}", D.Area[i][j].rho[l] / rhoR) << "\n";
								}
								else if (fixvar == 'Y') {
									if (abs(c2 - item1) <= tolY)
										extfile << std::format("{:.5e}", c1) << "\t" << std::format("{:.5e}", c3) << "\t"
										<< std::format("{:.5e}", D.Area[i][j].rho[l] / rhoR) << "\n";
								}
							}
							break;
						}/* kval case */

						case 1: /* DM-values */ {

							h = (D.Area[i][j].z2 - D.Area[i][j].z1) +
								(D.Area[i][j].d2 + D.Area[i][j].d1) + D.h0;	/* m */

							for (l = 0; l < Nz_curr; ++l) {

								c3 = z1_curr + (1.0 * l / (Nz_curr - 1));

								if (fixvar == '-')
									if (D.dkey == 2) extfile << std::format("{:.5e}", c1 * a) << "\t" << std::format("{:.5e}", c2 * b) << "\t"
									<< std::format("{:.5e}", c3 * h) << "\t" << std::format("{:.5e}", D.Area[i][j].rho[l]) << "\n";
									else extfile << std::format("{:.5e}", c1 * a) << "\t" << std::format("{:.5e}", c3 * h) << "\t" 
										<< std::format("{:.5e}", D.Area[i][j].rho[l]) << "\n";

								else if (fixvar == 'X') {
									if (abs(c1 - item1) <= tolX)
										extfile << std::format("{:.5e}", c2 * b) << "\t" << std::format("{:.5e}", c3 * h) << "\t"
										<< std::format("{:.5e}", D.Area[i][j].rho[l]) << "\n";
								}
								else if (fixvar == 'Y') {
									if (abs(c2 - item1) <= tolY)
										extfile << std::format("{:.5e}", c1 * a) << "\t" << std::format("{:.5e}", c3 * h) << "\t"
										<< std::format("{:.5e}", D.Area[i][j].rho[l]) << "\n";
								}
							}
							break;
						}/* kval case */

						}/* kval switch */

						break;/* entnum case */

					case 5: /* VISCOSITY (LOW-SHEAR) */

						switch (kval) {

						case 0: /* DL-values */ {

							for (l = 0; l < Nz_curr; ++l) {

								c3 = z1_curr + (1.0 * l / (Nz_curr - 1));

								if (fixvar == '-')
									if (D.dkey == 2) extfile << std::format("{:.5e}", c1) << "\t" << std::format("{:.5e}", c2) << "\t" << std::format("{:.5e}", c3) << "\t"
									<< std::format("{:.5e}", D.Area[i][j].mu[l] / muR) << "\n";
									else extfile << std::format("{:.5e}", c1) << "\t" << std::format("{:.5e}", c3) << "\t"
									<< std::format("{:.5e}", D.Area[i][j].mu[l] / muR) << "\n";
								else if (fixvar == 'X') {
									if (abs(c1 - item1) <= tolX)
										extfile << std::format("{:.5e}", c2) << "\t" << std::format("{:.5e}", c3) << "\t"
										<< std::format("{:.5e}", D.Area[i][j].mu[l] / muR) << "\n";
								}
								else if (fixvar == 'Y') {
									if (abs(c2 - item1) <= tolY)
										extfile << std::format("{:.5e}", c1) << "\t" << std::format("{:.5e}", c3) << "\t"
										<< std::format("{:.5e}", D.Area[i][j].mu[l] / muR) << "\n";
								}
							}
							break;

						}/* kval case */

						case 1: /* DM-values */ {

							h = (D.Area[i][j].z2 - D.Area[i][j].z1) +
								(D.Area[i][j].d2 + D.Area[i][j].d1) + D.h0;	/* m */

							for (l = 0; l < Nz_curr; ++l) {

								c3 = z1_curr + (1.0 * l / (Nz_curr - 1));

								if (fixvar == '-')
									if (D.dkey == 2) extfile << std::format("{:.5e}", c1 * a) << "\t" << std::format("{:.5e}", c2 * b) << "\t" << std::format("{:.5e}", c3 * h) << "\t"
									<< std::format("{:.5e}", D.Area[i][j].mu[l]) << "\n";
									else extfile << std::format("{:.5e}", c1 * a) << "\t" << std::format("{:.5e}", c3 * h) << "\t"
										<< std::format("{:.5e}", D.Area[i][j].mu[l]) << "\n";
								else if (fixvar == 'X') {
									if (abs(c1 - item1) <= tolX)
										extfile << std::format("{:.5e}", c2 * b) << "\t" << std::format("{:.5e}", c3 * h) << "\t"
										<< std::format("{:.5e}", D.Area[i][j].mu[l]) << "\n";
								}
								else if (fixvar == 'Y') {
									if (abs(c2 - item1) <= tolY)
										extfile << std::format("{:.5e}", c1 * a) << "\t" << std::format("{:.5e}", c3 * h) << "\t"
										<< std::format("{:.5e}", D.Area[i][j].mu[l]) << "\n";
								}
							}
							break;

						}/* kval case */

						} /* kval switch */
						break;/* entnum case */

					case 6: /* SHEAR STRESS */

						switch (kval) {

						case 0: /* DL-values */ {

							for (l = 0; l < Nz_curr; ++l) {

								c3 = z1_curr + (1.0 * l / (Nz_curr - 1));

								if (fixvar == '-')
									if (D.dkey == 2) extfile << std::format("{:.5e}", c1) << "\t" << std::format("{:.5e}", c2) << "\t"
									<< std::format("{:.5e}", c3) << "\t" << std::format("{:.5e}", D.Area[i][j].tau[l] / ph) << "\n";
									else extfile << std::format("{:.5e}", c1) << "\t" << std::format("{:.5e}", c3) << "\t" << std::format("{:.5e}", D.Area[i][j].tau[l] / ph) << "\n";

									else if (fixvar == 'X') {
									if (abs(c1 - item1) <= tolX)
										extfile << std::format("{:.5e}", c2) << "\t" << std::format("{:.5e}", c3) << "\t"
										<< std::format("{:.5e}", D.Area[i][j].tau[l] / ph) << "\n";
								}
								else if (fixvar == 'Y') {
									if (abs(c2 - item1) <= tolY)
										extfile << std::format("{:.5e}", c1) << "\t" << std::format("{:.5e}", c3) << "\t"
										<< std::format("{:.5e}", D.Area[i][j].tau[l] / ph) << "\n";
								}
							}
							break;

						} /* kval case*/

						case 1: /* DM-values */ {

							h = (D.Area[i][j].z2 - D.Area[i][j].z1) +
								(D.Area[i][j].d2 + D.Area[i][j].d1) + D.h0;	/* m */

							for (l = 0; l < Nz_curr; ++l) {

								c3 = z1_curr + (1.0 * l / (Nz_curr - 1));

								if (fixvar == '-')
									if (D.dkey == 2) extfile << std::format("{:.5e}", c1 * a) << "\t" << std::format("{:.5e}", c2 * b) << "\t" << std::format("{:.5e}", c3 * h) << "\t"
									<< std::format("{:.5e}", D.Area[i][j].tau[l]) << "\n";
									else extfile << std::format("{:.5e}", c1 * a) << "\t" << std::format("{:.5e}", c3 * h) << "\t"
										<< std::format("{:.5e}", D.Area[i][j].tau[l]) << "\n";
								else if (fixvar == 'X') {
									if (abs(c1 - item1) <= tolX)
										extfile << std::format("{:.5e}", c2 * b) << "\t" << std::format("{:.5e}", c3 * h) << "\t"
										<< std::format("{:.5e}", D.Area[i][j].tau[l]) << "\n";
								}
								else if (fixvar == 'Y') {
									if (abs(c2 - item1) <= tolY)
										extfile << std::format("{:.5e}", c1 * a) << "\t" << std::format("{:.5e}", c3 * h) << "\t"
										<< std::format("{:.5e}", D.Area[i][j].tau[l]) << "\n";
								}
							}
							break;

						} /* kval case */

						} /* kval switch */
						break;/* entnum case */

					case 7: /* GENERALIZED VISCOSITY */

						switch (kval) {

						case 0: /* DL-values */ {

							for (l = 0; l < Nz_curr; ++l) {

								c3 = z1_curr + (1.0 * l / (Nz_curr - 1));

								if (fixvar == '-')
									if (D.dkey == 2) extfile << std::format("{:.5e}", c1) << "\t" << std::format("{:.5e}", c2) << "\t" << std::format("{:.5e}", c3) << "\t"
									<< std::format("{:.5e}", D.Area[i][j].eta[l] / etaR) << "\n";
									else extfile << std::format("{:.5e}", c1) << "\t" << std::format("{:.5e}", c3) << "\t"
										<< std::format("{:.5e}", D.Area[i][j].eta[l] / etaR) << "\n";
								else if (fixvar == 'X') {
									if (abs(c1 - item1) <= tolX)
										extfile << std::format("{:.5e}", c2) << "\t" << std::format("{:.5e}", c3) << "\t"
										<< std::format("{:.5e}", D.Area[i][j].eta[l] / etaR) << "\n";
								}
								else if (fixvar == 'Y') {
									if (abs(c2 - item1) <= tolY)
										extfile << std::format("{:.5e}", c1) << "\t" << std::format("{:.5e}", c3) << "\t"
										<< std::format("{:.5e}", D.Area[i][j].eta[l] / etaR) << "\n";
								}
							}
							break;
						} /* kval case */

						case 1: /* DM-values */ {

							h = (D.Area[i][j].z2 - D.Area[i][j].z1) +
								(D.Area[i][j].d2 + D.Area[i][j].d1) + D.h0;	/* m */

							for (l = 0; l < Nz_curr; ++l) {

								c3 = z1_curr + (1.0 * l / (Nz_curr - 1));

								if (fixvar == '-')
									if (D.dkey == 2) extfile << std::format("{:.5e}", c1 * a) << "\t" << std::format("{:.5e}", c2 * b) << "\t" << std::format("{:.5e}", c3 * h) << "\t"
									<< std::format("{:.5e}", D.Area[i][j].eta[l]) << "\n";
									else extfile << std::format("{:.5e}", c1 * a) << "\t" << std::format("{:.5e}", c3 * h) << "\t"
										<< std::format("{:.5e}", D.Area[i][j].eta[l]) << "\n";
								else if (fixvar == 'X') {
									if (abs(c1 - item1) <= tolX)
										extfile << std::format("{:.5e}", c2 * b) << "\t" << std::format("{:.5e}", c3 * h) << "\t"
										<< std::format("{:.5e}", D.Area[i][j].eta[l]) << "\n";
								}
								else if (fixvar == 'Y') {
									if (abs(c2 - item1) <= tolY)
										extfile << std::format("{:.5e}", c1 * a) << "\t" << std::format("{:.5e}", c3 * h) << "\t"
										<< std::format("{:.5e}", D.Area[i][j].eta[l]) << "\n";
								}
							}
							break;

						} /* kval case */

						} /* kval switch */
						break;/* entnum case */

					case 8: /* VELOCITY (U-component) */

						switch (kval) {
						case 0: /* DL-values */ {

							break;
						}// kval case

						case 1: /* DM-values */ {

							h = (D.Area[i][j].z2 - D.Area[i][j].z1) +
								(D.Area[i][j].d2 + D.Area[i][j].d1) + D.h0;	/* m */

							for (l = 0; l < Nz_curr; ++l) {

								c3 = z1_curr + (1.0 * l / (Nz_curr - 1));

								if (fixvar == '-')
									if (D.dkey == 2) extfile << std::format("{:.5e}", c1 * a) << "\t" << std::format("{:.5e}", c2 * b) << "\t" << std::format("{:.5e}", c3 * h) << "\t"
										<< std::format("{:.5e}", D.Area[i][j].u[l]) << "\n";
									else extfile << std::format("{:.5e}", c1 * a) << "\t" << std::format("{:.5e}", c3 * h) << "\t"
										<< std::format("{:.5e}", D.Area[i][j].u[l]) << "\n";
								else if (fixvar == 'X') {
									if (abs(c1 - item1) <= tolX)
										extfile << std::format("{:.5e}", c2 * b) << "\t" << std::format("{:.5e}", c3 * h) << "\t"
										<< std::format("{:.5e}", D.Area[i][j].u[l]) << "\n";
								}
								else if (fixvar == 'Y') {
									if (abs(c2 - item1) <= tolY)
										extfile << std::format("{:.5e}", c1 * a) << "\t" << std::format("{:.5e}", c3 * h) << "\t"
										<< std::format("{:.5e}", D.Area[i][j].u[l]) << "\n";
								}
							}
							break;

						} /* kval case */

							break;
						}// kval case
						break; /* entnum case */

					case 9: /* VELOCITY (V-component) */

						switch (kval) {
						case 0: /* DL-values */ {

							break;
						}// kval case

						case 1: /* DM-values */ {

							h = (D.Area[i][j].z2 - D.Area[i][j].z1) +
								(D.Area[i][j].d2 + D.Area[i][j].d1) + D.h0;	/* m */

							for (l = 0; l < Nz_curr; ++l) {

								c3 = z1_curr + (1.0 * l / (Nz_curr - 1));

								if (fixvar == '-')
									if (D.dkey == 2) extfile << std::format("{:.5e}", c1 * a) << "\t" << std::format("{:.5e}", c2 * b) << "\t" << std::format("{:.5e}", c3 * h) << "\t"
										<< std::format("{:.5e}", D.Area[i][j].v[l]) << "\n";
									else extfile << std::format("{:.5e}", c1 * a) << "\t" << std::format("{:.5e}", c3 * h) << "\t"
										<< std::format("{:.5e}", D.Area[i][j].v[l]) << "\n";
								else if (fixvar == 'X') {
									if (abs(c1 - item1) <= tolX)
										extfile << std::format("{:.5e}", c2 * b) << "\t" << std::format("{:.5e}", c3 * h) << "\t"
										<< std::format("{:.5e}", D.Area[i][j].v[l]) << "\n";
								}
								else if (fixvar == 'Y') {
									if (abs(c2 - item1) <= tolY)
										extfile << std::format("{:.5e}", c1 * a) << "\t" << std::format("{:.5e}", c3 * h) << "\t"
										<< std::format("{:.5e}", D.Area[i][j].v[l]) << "\n";
								}
							}
							break;
						} // kval case
						} // kval case
						break; /* entnum case */

					case 10: /* VELOCITY (U,V) */

						switch (kval) {}
						break;

					case 11: /* TEMPERATURE */

						switch (kval) {}
						break;

					case 12: /* THERMAL CONDUCTIVITY */

						switch (kval) {}
						break;

					case 13: /* HEAT CAPACITY */
						switch (kval) {}
						break;
					} /* entnum switch */
						

				} /* enttype */

				if (enttype == 2) { // Friction results (distributed)
				
					switch (entnum) {
					
					case 0: break; /* reserved */
					
					case 1: { //  

													
						break; 
					}
								
					
					} // entnum switch
				
				
				}

			} /* j-loop */
	}

	
	extfile.close(); //std::cout << "we closed file: " << fname + pstfx + ".TXT\n";
			
} /* end proc */


