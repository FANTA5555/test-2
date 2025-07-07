#include <iostream>
#include "main.h"
#include "models.h"
#include "deform.h"
#include "output.h"
#include <format>
/* UMFPACK lib for solving the system of linear equations with sparse matrix */
#include "suitesparse/umfpack.h" 


/********************* COMMON PROCEDURES ******************************/

void dumpvar1D(Domain D, int kval, int dtype, int kout) {
	//****** write variable from domain to screen or external file ******
	// kval: 1 - pressure, 2 - thickness, 3 - temperature
	// dtype: 0 - as is; 1/2/3/4 - 2x/4x/8x/16x linear interpolate
	// kout: 0 - for initalization; 1 - for plotting profile

	int Nx, nNx, i, prv, nxt, ndiv, lx_end;
	double dx, stpX, cA, cB, delta, curr;
	double x_curr, y_curr;
	double deltaAB;
	string fname, dlm, pstfx, strrsx;

	Nx = D.Nx;
	dx = D.dx;
	
	std::fstream extfile;

	ndiv = static_cast <int>(1 / pow(0.5, dtype) - 1); // division number

	strrsx = "_r" + std::to_string(D.rsx);
	if (kout == 0) pstfx = strrsx + "_init"; else pstfx = strrsx + "_plot";

	switch (kval) {

	case 1: // pressure 
	{
		fname = "LC-PRS_DL" + pstfx + ".DAT";
		extfile.open(fname, std::ios::out);
		break;
	} // case 1

	case 2: // thickness 
	{
		fname = "LC-THK_DL" + pstfx + ".DAT";
		extfile.open(fname, std::ios::out);

	} // case 2
	break;
	} // switch

	stpX = dx * pow(0.5, dtype); // step X
	nNx = static_cast<int>(floor(2 * D.kx / stpX) + 1);
	
	if (kout == 0) {

		dlm = " ";
		extfile << D.kx << " " << nNx << " " << D.h0 << "\n";
	}

	else dlm = "\t";

	for (i = 0; i < (Nx - 1); ++i) {

		if (i == Nx - 2) lx_end = ndiv + 1; else lx_end = ndiv;
		for (int lx = 0; lx <= lx_end; ++lx) {

			x_curr = (-1.0 * D.kx) + i * D.dx + stpX * lx;

				switch (kval) {

				case 1: // pressure 
				{
					cA = D.Area[i    ][0].p;
					cB = D.Area[i + 1][0].p;
					break;
				} // case 1

				case 2: // thickness 
				{
					cA = D.Area[i    ][0].h;
					cB = D.Area[i + 1][0].h;
					break;
				} // case 2
				break;
				} // switch

				deltaAB = (cB - cA) / (ndiv + 1);
				
				extfile << x_curr << dlm << cA + (deltaAB * lx) << "\n";
				
		} // lx loop
	} // i loop

	extfile.close();

} // proc

void dumpvar2D(Domain D, int kval, int dtype, int kout) {
//****** write variable from domain to screen or external file ******
// kval: 1 - pressure, 2 - thickness, 3 - temperature
// dtype: 0 - as is; 1/2/3/4 - 2x/4x/8x/16x linear interpolate
// kout: 0 - for initalization; 1 - for plotting profile

	int Nx, Ny, nNx, nNy, i, j, prv, nxt, ndiv, lx_end, ly_end;
	double dx, dy, stpX, stpY, cA, cB, cC, cD, delta, curr;
	double x_curr, y_curr;
	double deltaAB, deltaAC, deltaCD, deltaAC_curr;
	string fname, dlm, pstfx,strrsx;

	Nx = D.Nx;
	Ny = D.Ny;
	dx = D.dx;
	dy = D.dy;
	
	std::fstream extfile;
	
	ndiv = static_cast <int>(1 / pow(0.5, dtype) - 1); // division number

	strrsx = "_r" + std::to_string(D.rsx);
	if (kout == 0) pstfx = strrsx + "_init"; else pstfx = strrsx + "_plot";
		
	switch (kval) {

	case 1: // pressure 
	{
		fname = "PC-PRS_DL" + pstfx + ".DAT";
		extfile.open(fname, std::ios::out);
		break;
	} // case 1

	case 2: // thickness 
	{
		fname = "PC-THK_DL" + pstfx + ".DAT";
		extfile.open(fname, std::ios::out);

	} // case 2
	break;
	} // switch

	stpX = dx * pow(0.5, dtype); // step X
	stpY = dy * pow(0.5, dtype); // step Y
	nNx = static_cast<int>(floor(2 * D.kx / stpX) + 1);
	nNy = static_cast<int>(floor(2 * D.ky / stpY) + 1);

	

	if (kout == 0) {

		dlm = " ";
		extfile << D.kx << " " << D.ky << " " << nNx * nNy << " " << D.h0 << "\n";
	}

	else dlm = "\t";
	
	for (i = 0; i < (Nx - 1); ++i) {
		
		if (i == Nx - 2) lx_end = ndiv + 1; else lx_end = ndiv;
		for (int lx = 0; lx <= lx_end; ++lx) {

			x_curr = (-1.0 * D.kx) + i * D.dx + stpX * lx;

			for (j = 0; j < (Ny - 1); ++j) {

				y_curr = (-1.0 * D.ky) + j * D.dy;

				switch (kval) {

				case 1: // pressure 
				{
					cA = D.Area[i][j].p;
					cB = D.Area[i + 1][j].p;
					cC = D.Area[i][j + 1].p;
					cD = D.Area[i + 1][j + 1].p;
					break;
				} // case 1

				case 2: // thickness 
				{
					cA = D.Area[i][j].h;
					cB = D.Area[i + 1][j].h;
					cC = D.Area[i][j + 1].h;
					cD = D.Area[i + 1][j + 1].h;
					break;
				} // case 2
				break;
				} // switch
				
				deltaAB = (cB - cA) / (ndiv + 1);
				deltaAC = (cC - cA) / (ndiv + 1);
				deltaCD = (cD - cC) / (ndiv + 1);
				
				deltaAC_curr = ((cC + (deltaCD * lx)) - (cA + (deltaAB * lx))) / (ndiv + 1);
				
				if (j == Ny - 2) ly_end = ndiv + 1; else ly_end = ndiv;
				for (int ly = 0; ly <= ly_end; ++ly) {
					extfile << x_curr << dlm << y_curr + stpY * ly << dlm << cA + (deltaAB * lx) + deltaAC * ly << "\n";
				}
				
			} //j loop

		} // lx loop
	} // i loop

	extfile.close();
	
} // proc

void calc_triplets(int n, int m, double* A, int*& Ti, int*& Tj, double*& Tx) {
/* optional procedure, since we construct the matrix in triplet form directly */
}

void calc_rhoeps(int i, int j, Domain D, Body B1, Body B2, double &eps_curr, double &rhox_curr, double &rhoy_curr)
//****** calculate epsilon, rhostar_x, rhostar_y in [i,j]-point ******
{
	dbvector inveta, rho, rhoeta, rhozeta;
	
	int		 l, ll, Nz;
	
	double	h,h1, RX, z1, z2, dz, z_curr, res;
	double	eta_e, rho_e, eta_e1, rho_1, rho_2;


	Nz = D.Area[i][j].Nz;
		
	//h = (D.Area[i][j].z2 - D.Area[i][j].z1) + (D.Area[i][j].d2 + D.Area[i][j].d1) + D.h0;	// m
	
	h = D.Area[i][j].h;

	// all parameters below is in dimensionless form
	z1 = (D.Area[i][j].z1 - D.Area[i][j].d1		  ) / h;	// start Z-coord
	z2 = (D.Area[i][j].z2 + D.Area[i][j].d2 + D.h0) / h;	// end Z-coord
	dz = abs(z2 - z1) / (Nz - 1); // each [i,j]-point has its own dz
	
	// calculate the vectors inveta, rhoeta, rhozeta
	for (l = 0; l < Nz; ++l) {
		rho.push_back	(	  D.Area[i][j].rho[l] / D.rhoR);
		inveta.push_back(1 / (D.Area[i][j].eta[l] / D.etaR));
		
		if (l > 0) { // integrate (inveta) and (z_inveta) 
			rhoeta.push_back ((D.Area[i][j].rho[l] / D.rhoR) * calc_integrals(inveta, z1, dz, l + 1, 1, 0));
			rhozeta.push_back((D.Area[i][j].rho[l] / D.rhoR) * calc_integrals(inveta, z1, dz, l + 1, 2, 0));
		}
		else { // first time the integration limits are equal
			rhoeta.push_back(0.0);
			rhozeta.push_back(0.0);
		}
	}

	eta_e	= 1 / calc_integrals(inveta, z1, dz, Nz, 1, 0);
	rho_e	=	  calc_integrals(rho,	 z1, dz, Nz, 1, 0);

	eta_e1  = 1 / calc_integrals(inveta, z1, dz, Nz, 2, 0);
	rho_1	=     calc_integrals(rhoeta, z1, dz, Nz, 1, 0);
	rho_2	=	  calc_integrals(rhozeta,z1, dz, Nz, 1, 0);

	// epsilon in dimensionless form
	eps_curr = (eta_e / eta_e1) * rho_1 - rho_2; 
	
	// rhostar_x and rhostar_y in dimensionless form
	// for steady state case we don't use a um or vm (see my theoretical background)
	rhox_curr = (rho_e * B1.vel[0] + eta_e * (B2.vel[0] - B1.vel[0]) * rho_1); 
	
	if (D.dkey == 2)
		 rhoy_curr = (rho_e * B1.vel[1] + eta_e * (B2.vel[1] - B1.vel[1]) * rho_1);
		
	else rhoy_curr = 0.0;
} // proc

double dFdP(int k, int l, int i, int j, double **rhox, double **rhoy, double dx, double dy, double theta) {
//****** calculate the derivative dF[k,l] by dP[i,j] ******
	
	double p1, p2;

	p1 = (2 / (pi * pi) * dx) * (rhox[k][l] * calc_IC_2D(k - i, l - j, dx) - rhox[k - 1][l] * calc_IC_2D(k - 1 - i, l - j, dx));

	p2 = (2 * theta / (pi * pi) * dy) * (rhoy[k][l] * calc_IC_2D(k - i, l - j, dx) - rhoy[k][l - 1] * calc_IC_2D(k - i, l - 1 - j, dx));

	return (p1 + p2);
} // proc


  
/********************* LINE CONTACT PROCEDURES ******************************/

void calc_coef1D(int i, Domain D, Body B1, Body B2, double* eps, double* rhox, double* H,
	double& cA, double& cB, double& cE, double& cF, double &cJ, int kval) {
//****** calculate the discretization coefficients in [i]-point (line contact) ******/

	double dx, alpha, um, RX;
	double eH3_i05, eH3_05i, A1;

	dx = D.dx;
	RX = pow(1 / B1.Rx + 1 / B2.Rx, -1);
	
	um = 1.0; // for steady state we don't use um though leave it in the formulas
	alpha = D.ph * pow(D.a, 3) / (um * D.etaR * pow(RX, 2));
			
	// (eH3)[i+1/2]
	eH3_i05 = 0.5 * (eps[i + 1] * pow(H[i + 1], 3) + eps[i] * pow(H[i], 3));
	// (eH3)[i-1/2]
	eH3_05i = 0.5 * (eps[i - 1] * pow(H[i - 1], 3) + eps[i] * pow(H[i], 3));
	
	// standart discretization scheme 
	cA = alpha / (dx * dx) * eH3_i05;	// P[i+1][j] coeff
	cB = alpha / (dx * dx) * eH3_05i;	// P[i-1][j] coeff

	cE = - alpha * (eH3_05i + eH3_i05) / (dx * dx); // P[i][j] coeff
	A1 =   alpha * (eH3_i05 + eH3_05i) / (dx * dx);
	
	// wedge term (or Couette term)
	if (kval==0)	
		// 1st order backward scheme
		cF = (1 / dx) * (rhox[i] * H[i] - rhox[i - 1] * H[i - 1]); 
	
	else {
		// 2nd order backward scheme 
		if (i == 1) // for first line we use 1st order scheme
			cF = (1 / (2 * dx)) * (rhox[i] * H[i] - rhox[i - 1] * H[i - 1]);
		else
			cF = (1 / (2 * dx)) * (3 * rhox[i] * H[i] - 4 * rhox[i - 1] * H[i - 1] + rhox[i - 2] * H[i - 2]);
	}
	
	// coefficient for Jacobi relaxation formula
	cJ = alpha * (-2 * eH3_05i - eH3_i05) / (dx * dx);

} //proc
 

void calc_pressure_field_1D_old(Domain& D, Body B1, Body B2)
{
	double
		* eps,		/* dimensionless epsilon		*/
		* rhox,		/* dimensionless rhostar_x		*/
		* rhoy,		/* dimensionless rhostar_y (for this case it equal to zero) */
		* H;		/* dimensionless thickness */
				
	int Nx, i;
	int lc_prs, lc_h0, maxlc_prs0, maxlc_prs1, maxlc_prs, maxlc_h0, wt_dst;

	double RX, RY, x_curr, eps_curr, rhox_curr, rhoy_curr, a;
	double dx, p1, um, prs, i_beg, i_end;
	double alpha, cA, cB, cE, cF, cJ;

	double h0_rlx, H_i, H_1i, H_i1, eH3_05i, eH3_i05;
	double T1, T2, gammai, derivi, ptol, swc, summ0, summ1, perr, perr_old, Hmin;
	double gs_rlx, jc_rlx, H0min, H_2i, rhox_2i, k, pint_curr, prs_errtol;
	bool f_prs, f_h0;
	
	double rs, summ2, rmsres, res_errtol;
	double H1, H1k, H2;
	
	Nx = D.Nx;
	dx = D.dx;
	a = D.a;
	RX = pow(1 / B1.Rx + 1 / B2.Rx, -1);
		
	um = (B1.vel[0] + B2.vel[0]) / 2.0; // m/s 
	alpha = D.ph * pow(D.a, 3) / (um * D.etaR * pow(RX, 2));
	H0min = -0.596574; // limiting value for H0 from Herzian dry contact theory (Venner, 1991, p.80)
	
	// relaxation coefficient for h0 (apply (Huang,2013)-method to Hamrock-formula for Hm)
	//h0_rlx = abs((-0.166) * RX * 2.922 * pow(D.W_HD, -1.166) * pow(D.U_HD, 0.692) * pow(D.G_HD, 0.47));
	
	/**** ADJUSTABLE PARAMETERS *******/
		
	h0_rlx = 0.004;		// H0 relaxation factor

	gs_rlx = 0.8;		// Gauss-Seidel relaxation factor
	jc_rlx = 0.8;		// Jacobi relaxation factor
		
	maxlc_h0  = 270;		// max H0 relaxation steps
	maxlc_prs0 = 28;	// max pressure sweeps before H0 updating
	prs_errtol = 1e-15;	// pressure error tolerance
	res_errtol  = 1e-1;	// pressure error tolerance
	
	H1 = 0.001;			// H-level for decrease the H0 relaxation faclor
	H1k = 0.1;			// decreasing factor (< 1)
	maxlc_prs1 = 4;		// max pressure sweeps before H0 updating (after decreasing)
	
	H2 = 0.05;			// H-level for freeze the H0 changes (h0_rlx = 0)
	
	/**********************************/

	/* boundary conditions (zero-pressure) */
	D.Area[0][0].p = 0.0; 
	D.Area[Nx - 1][0].p = 0.0;
	i_beg = 1;
	i_end = Nx - 1; //limits for i-index
	
	/* optional boundary conditions (zero pressure gradient) */
	//D.Area[1][0].p = 0.0;
	//D.Area[Nx - 2][0].p = 0.0;
	//i_beg = 1; i_end = Nx - 2; //limits for i-index
		
	eps  = new double[Nx];
	rhox = new double[Nx];
	H	 = new double[Nx];
		 
	lc_h0 = 1;
	f_h0 = false; // force balance flag
	maxlc_prs = maxlc_prs0;

	std::fstream extfile;
	extfile.open("lc_log.txt", std::ios::out);
	extfile << "PRS ITER\tPRS ERR\tPRS INT\tRMSRES\tHmin\n\n";

	// loop for force balancing
	while ((!f_h0) && (lc_h0 < maxlc_h0 + 1)) {
		 
		lc_prs = 1;
		f_prs = false; // pressure convergence flag
								
		// loop for pressure calculation
		while ((!f_prs) && (lc_prs < maxlc_prs + 1)) {
		
			// update film thickness and lubricant properties
			calc_deform(D, B1, B2);
			//calc_physical_fields(D, D.RHO, D.LSV, D.GNV, D.THC, D.HCP, B1.vel, B2.vel);
			calc_densviscNB(D, D.RHO, D.LSV);
			summ0 = 0.0; summ1 = 0.0; summ2 = 0.0; Hmin = 1e10;

			// fill arrays for epsilon, rhostarx and thickness (dimensionless form)
			for (i = 0; i < Nx; ++i) {

				x_curr = (-1.0 * D.kx) + i * D.dx;

				eps_curr = 0; rhox_curr = 0; rhoy_curr = 0;
				calc_rhoeps(i, 0, D, B1, B2, eps_curr, rhox_curr, rhoy_curr);

				eps[i]  = eps_curr;
				rhox[i] = rhox_curr;
				
				p1 = (RX / (a * a)) * (D.Area[i][0].d1 + D.Area[i][0].d2);
				H[i] = (RX / (a * a)) * D.h0 + 0.5 * pow(x_curr, 2) + p1;
				D.Area[i][0].h = ((a * a) / RX) * H[i];

				if (H[i] < Hmin) Hmin = H[i]; // update the minimum film thickness
				
				if (Hmin < 0)
					extfile << "\n<<<< THE NEGATIVE HMIN WAS DETECTED, YOUR RESULTS MAY BE NOT CORRECT >>>>\n";

				//std::cout << "eps[" << i << "] = " << eps_curr << "\n";
				//std::cout << "rhox[" << i << "] = " << rhox_curr << "\n";

			} /* i loop */

			// calculate new pressure field (hybrid relaxation)
			wt_dst = 0; // type of wedge term discretization (0 - 1st order; 1 - 2nd order) 
			for (i = i_beg; i < i_end; ++i) {
								
				calc_coef1D(i, D, B1, B2, eps, rhox, H, cA, cB, cE, cF, cJ, wt_dst);
				
				if (alpha * eps[i] * pow(H[i], 3)/ (dx * dx) >= 0.3) // classic condition
				//if (alpha * eps[i] * pow(H[i], 3) / dx >= 0.01) // conservative condition
				{ 
				// Gauss-Seidel relaxation
					
					gammai = cF - cA * D.Area[i + 1][0].pold - cB * D.Area[i - 1][0].p - cE * D.Area[i][0].pold;
										
					if ((i > i_beg) && (wt_dst == 1))
					// use 2nd order wedge term in derivative
					derivi = cE + (1 / (pi * 2 * dx)) * (3 * rhox[i] * calc_IC_1D(0, dx) - 4 * rhox[i - 1] * calc_IC_1D(-1, dx) + 
							rhox[i - 2] * calc_IC_1D(-2, dx));
					else 
					// use 1st order wedge term in derivative
					derivi = cE + (1 / (pi * dx)) * (rhox[i] * calc_IC_1D(0, dx) - rhox[i - 1] * calc_IC_1D(-1, dx));
					
					D.Area[i][0].p = D.Area[i][0].pold + gs_rlx * (1 / derivi) * gammai;
				}
				else 
				{ 
				// Jacobi relaxation
					
					gammai = cF - cA * D.Area[i + 1][0].pold - cB * D.Area[i - 1][0].pold - cE * D.Area[i][0].pold;
															
					if ((i > i_beg) && (wt_dst ==1))
						// use 2nd order wedge term in derivative
						derivi = (cE - cB) + 1 / (pi * 2 * dx) * ((3 * rhox[i] + 4 * rhox[i - 1]) * (calc_IC_1D(0, dx) - calc_IC_1D(-1, dx)) - 
							      rhox[i - 2] * (calc_IC_1D(-1, dx) - calc_IC_1D(-2, dx)));
					else
						// use 1st order wedge term in derivative
						/* My term for derivative (only for 1st order discr. of wedge term) */
						//derivi = (cE - cB) + (calc_IC_1D(0, dx) - calc_IC_1D(1, dx)) / (pi * dx) * (rhox[i] + rhox[i - 1]);
						
						/* Venner term for derivative (only for 1st order discr. of wedge term) */
						derivi = cJ + (2 / pi) * (rhox[i] * calc_IC_1D(0, dx) - rhox[i - 1] * calc_IC_1D(-1, dx));
					
					D.Area[i    ][0].p = D.Area[i    ][0].pold + jc_rlx * (1 / derivi) * gammai;
					D.Area[i - 1][0].p = D.Area[i - 1][0].pold - jc_rlx * (1 / derivi) * gammai;
				}
				
				// cavitation condition
				if (D.Area[i][0].p < 0.0) D.Area[i][0].p = 0.0; 
				
				// residual at node i
				rs = abs(cF - cA * D.Area[i + 1][0].p - cB * D.Area[i - 1][0].p - cE * D.Area[i][0].p);
				
				// cumulative sum for pressure changes, pressure integral and residual
				summ0 = summ0 + abs(D.Area[i][0].pold - D.Area[i][0].p);
				summ1 = summ1 + D.Area[i][0].p;
				summ2 = summ2 + rs * rs;

				// save current pressure value for next iteration
				D.Area[i][0].pold = D.Area[i][0].p; 

			} // i loop

			perr = summ0 / summ1;	// pressure error on current iteration 
			pint_curr = summ1 * dx; // current value for pressure intergral
			rmsres = pow(summ2 / (i_end - i_beg), 0.5); // current RMS residual

			std::cout << "PRS ITER = " << lc_prs << "\tPRS ERR = " << std::format("{:.5e}", perr) << "\tPRS INT = " <<
			std::format("{:.5e}", pint_curr) << "\tRMSRES = " << std::format("{:.5e}", rmsres) << "\tHmin = " << std::format("{:.5e}", Hmin) << "\n";

			extfile << lc_prs << "\t" << std::format("{:.5e}", perr) << "\t" <<
				std::format("{:.5e}", pint_curr) << "\t" << std::format("{:.5e}", rmsres) << "\t" << std::format("{:.5e}", Hmin) << "\n";


			/* check the pressure convergence criterion */
			if ((perr < prs_errtol)) {
				std::cout << "\n << PRESSURES ARE CONVERGED << \n";
				extfile << "\n << PRESSURES ARE CONVERGED << \n";
				f_prs = true;
			}
			else lc_prs = lc_prs + 1;
						
		} // while (pressure convergence)

		// adjust h0 and check error for pressure integral
		if ((abs(pint_curr - (pi / 2)) / (pi / 2) < 1e-4) && (rmsres <= res_errtol)) {

			std::cout << " << FORCE BALANCE IS REACHED << \n";
			extfile << " << FORCE BALANCE IS REACHED << \n";
			f_h0 = true;
		}
		else {

			/* update h0 and continue calculations */
			D.h0 = D.h0 - h0_rlx * ((pi / 2.0) - pint_curr) * a * a / RX;
			/*
			{
				if (Hmin < H1) {
					if (Hmin > H2)
						D.h0 = D.h0 - H1k * h0_rlx * ((pi / 2.0) - pint_curr) * a * a / RX;

					maxlc_prs = maxlc_prs1;

				}
				else {
					D.h0 = D.h0 - h0_rlx * ((pi / 2.0) - pint_curr) * a * a / RX;
					maxlc_prs = maxlc_prs0;
				}
			}
			*/
			lc_h0 = lc_h0 + 1;
			if (lc_h0 < maxlc_h0 + 1) {
				std::cout << "\n<<<< HO ITER = " << lc_h0 << "\tNEW H0 = " << std::format("{:.5e}", RX * D.h0 / (a * a)) << "\n\n";
				extfile << "\n<<<<<< H0 ITER = " << lc_h0 << "\tNEW H0 = " << std::format("{:.5e}", RX * D.h0 / (a * a)) << "\n\n";
			}
		}
	} // while (force balance)

	/* free arrays for eps, rho */
	{
		delete[] eps;
		delete[] rhox;
	}
	
	extfile.close();

} /* proc */
void calc_pressure_field_1D_backup(Domain& D, Body B1, Body B2)
{
	double
		* eps,		/* dimensionless epsilon		*/
		* rhox,		/* dimensionless rhostar_x		*/
		* rhoy,		/* dimensionless rhostar_y (for this case it equal to zero) */
		* H;		/* dimensionless thickness */

	int Nx, i;
	int lc_prs, lc_h0, maxlc_prs0, maxlc_prs1, maxlc_prs, maxlc_h0, wt_dst;
	int nd_cnt;

	double RX, RY, x_curr, eps_curr, rhox_curr, rhoy_curr, a;
	double dx, p1, um, prs, i_beg, i_end;
	double alpha, cA, cB, cE, cF, cJ;

	double h0_rlx, H_i, H_1i, H_i1, eH3_05i, eH3_i05;
	double T1, T2, gammai, derivi, ptol, swc, summ0, summ1, summ3, perr, perr_old, Hmin;
	double gs_rlx, jc_rlx, H0min, H_2i, rhox_2i, k, pint_curr, prs_errtol;
	bool f_prs, f_h0, f_neg;

	double rs, summ2, rmsres, res_errtol, frc_errtol, A1;
	double H1, H1k, H2, delta, A2, H0lim;
	string gapstr;

	std::fstream extfile;
	std::fstream extfile1;
	std::fstream extfile2;

	Nx = D.Nx;
	dx = D.dx;
	a = D.a;
	RX = pow(1 / B1.Rx + 1 / B2.Rx, -1);

	um = (B1.vel[0] + B2.vel[0]) / 2.0; // m/s 
	alpha = D.ph * pow(D.a, 3) / (um * D.etaR * pow(RX, 2));

	H0lim = -0.596574; // limiting value for H0 from Herzian dry contact theory (Venner, 1991, p.80)

	// the minimum film thickness for line contact (empirical formula)
	H0min = 1.714 * pow(D.W_HD, -0.128) * pow(D.U_HD, 0.694) * pow(D.G_HD, 0.568);

	// relaxation coefficient for h0 (apply (Huang,2013)-method to Hamrock-formula for Hc)
	//h0_rlx = abs((-0.166) * RX * 2.922 * pow(D.W_HD, -1.166) * pow(D.U_HD, 0.692) * pow(D.G_HD, 0.47));

	/**** ADJUSTABLE PARAMETERS *******/

	h0_rlx = 0.05;		// H0 relaxation factor

	gs_rlx = 0.7;		// Gauss-Seidel relaxation factor
	jc_rlx = 0.4;		// Jacobi relaxation factor

	maxlc_h0 = 40;		// max H0 relaxation steps
	maxlc_prs0 = 6;		// max pressure sweeps before H0 updating

	prs_errtol = 1e-15;	// pressure error tolerance
	res_errtol = 1e-5;	// residual tolerance
	frc_errtol = 1e-4;  // force balance tolerance

	/**********************************/

	/* boundary conditions (zero-pressure) */
	D.Area[0][0].p = 0.0;
	D.Area[Nx - 1][0].p = 0.0;
	i_beg = 1;
	i_end = Nx - 1; //limits for i-index

	/* optional boundary conditions (zero pressure gradient) */
	//D.Area[1][0].p = 0.0;
	//D.Area[Nx - 2][0].p = 0.0;
	//i_beg = 1; i_end = Nx - 2; //limits for i-index

	eps = new double[Nx];
	rhox = new double[Nx];
	H = new double[Nx];

	lc_h0 = 2;
	f_h0 = false; // force balance flag
	maxlc_prs = maxlc_prs0;

	extfile.open("lc_log.txt", std::ios::out);
	extfile << "\t*** INITIAL VALUE FOR H0 = " << std::format("{:.5e}", D.h0 * RX / (a * a))
		<< " (" << std::format("{:.5e}", D.h0) << " m)\n\n";
	extfile << "\tH0 ITER = 1\n\n";


	// loop for force balancing
	while ((!f_h0) && (lc_h0 < maxlc_h0 + 2)) {

		lc_prs = 1;
		f_prs = false; // pressure convergence flag

		// loop for pressure calculation
		while ((!f_prs) && (lc_prs < maxlc_prs + 1)) {

			// update film thickness and lubricant properties
			calc_deform(D, B1, B2);
			calc_densviscNB(D, D.RHO, D.LSV);
			summ0 = 0.0; summ1 = 0.0; summ2 = 0.0; summ3 = 0.0; Hmin = 1e10;

			// fill arrays for epsilon, rhostarx and thickness (dimensionless form)
			f_neg = false; // flag for negative Hmin
			for (i = 0; i < Nx; ++i) {

				x_curr = (-1.0 * D.kx) + i * D.dx;

				// calculate film thickness
				p1 = (RX / (a * a)) * (D.Area[i][0].d1 + D.Area[i][0].d2);
				H[i] = (RX / (a * a)) * D.h0 + 0.5 * pow(x_curr, 2) + p1;

				D.Area[i][0].h = ((a * a) / RX) * H[i];

				if (H[i] < Hmin) Hmin = H[i]; // update the minimum film thickness

				if ((Hmin < 0) && (!f_neg)) {
					std::cout << "\n\t*** HMIN = " << Hmin << " WAS DETECTED AT X = " << x_curr << "\n";
					extfile << "\n\t*** HMIN = " << Hmin << " WAS DETECTED AT X = " << x_curr << "\n";
					f_neg = true;
				}

				// calculate epsilon and rhostarx (rhoy_curr remains but is not used here)
				eps_curr = 0; rhox_curr = 0; rhoy_curr = 0;
				calc_rhoeps(i, 0, D, B1, B2, eps_curr, rhox_curr, rhoy_curr);

				eps[i] = eps_curr;
				rhox[i] = rhox_curr;


			} /* i loop */

			// correct H0 if negative Hmin is found
			if (Hmin < 0) {
				delta = abs(Hmin) + H0min;
				D.h0 = D.h0 + delta * (a * a) / RX;

				std::cout << "\t*** CORRECT H0 BY " << std::format("{:.5e}", delta) << ", NEW VALUE IS " <<
					std::format("{:.5e}", D.h0 * RX / (a * a)) << " (" << std::format("{:.5e}", D.h0) << " m)\n";
				extfile << "\t*** CORRECT H0 BY " << std::format("{:.5e}", delta) << ", NEW VALUE IS " <<
					std::format("{:.5e}", D.h0 * RX / (a * a)) << " (" << std::format("{:.5e}", D.h0) << " m)\n";

				Hmin = 1e10;
				for (i = 0; i < Nx; ++i) {

					H[i] = H[i] + delta;

					if (H[i] < Hmin) Hmin = H[i];
				}

				std::cout << "\t*** AFTER CORRECTION HMIN = " << std::format("{:.5e}", Hmin) << "\n\n";
				extfile << "\t*** AFTER CORRECTION HMIN = " << std::format("{:.5e}", Hmin) << "\n\n";
			} // if


			// calculate new pressure field (hybrid relaxation)
			wt_dst = 0; // type of wedge term discretization (0 - 1st order; 1 - 2nd order) 
			nd_cnt = 0; // number of nodes in "contact" zone (for jacobi relaxations)
			for (i = i_beg; i < i_end; ++i) {

				calc_coef1D(i, D, B1, B2, eps, rhox, H, cA, cB, cE, cF, cJ, wt_dst);
				// term for Huang's conditions (if it used)
				A2 = (1 / (pi * dx)) * (rhox[i] * calc_IC_1D(0, dx) - rhox[i - 1] * calc_IC_1D(-1, dx));

				if (abs(cE) > 0.05 * abs(A2)) // Huang's, 2013 condition
					//if (alpha * eps[i] * pow(H[i], 3) / (dx * dx) >= 0.5) // classic condition (Venner, 1991)
					//if (alpha * eps[i] * pow(H[i], 3) / dx >= 0.01) // conservative condition (Venner, 1991)
				{
					// Gauss-Seidel relaxation

					gammai = cF - cA * D.Area[i + 1][0].pold - cB * D.Area[i - 1][0].p - cE * D.Area[i][0].pold;

					if ((i > i_beg) && (wt_dst == 1))
						// use 2nd order wedge term in derivative
						derivi = cE + (1 / (pi * 2 * dx)) * (3 * rhox[i] * calc_IC_1D(0, dx) - 4 * rhox[i - 1] * calc_IC_1D(-1, dx) +
							rhox[i - 2] * calc_IC_1D(-2, dx));
					else
						// use 1st order wedge term in derivative
						derivi = cE + (1 / (pi * dx)) * (rhox[i] * calc_IC_1D(0, dx) - rhox[i - 1] * calc_IC_1D(-1, dx));

					D.Area[i][0].p = D.Area[i][0].pold + gs_rlx * (1 / derivi) * gammai;
				}
				else
				{
					// Jacobi relaxation
					nd_cnt = nd_cnt + 1;
					gammai = cF - cA * D.Area[i + 1][0].pold - cB * D.Area[i - 1][0].pold - cE * D.Area[i][0].pold;

					if ((i > i_beg) && (wt_dst == 1))
						// use 2nd order wedge term in derivative
						derivi = (cE - cB) + 1 / (pi * 2 * dx) * ((3 * rhox[i] + 4 * rhox[i - 1]) * (calc_IC_1D(0, dx) - calc_IC_1D(-1, dx)) -
							rhox[i - 2] * (calc_IC_1D(-1, dx) - calc_IC_1D(-2, dx)));
					else
						// use 1st order wedge term in derivative
						/* My term for derivative (only for 1st order discr. of wedge term) */
						//derivi = (cE - cB) + (calc_IC_1D(0, dx) - calc_IC_1D(1, dx)) / (pi * dx) * (rhox[i] + rhox[i - 1]);

						/* Venner term for derivative (only for 1st order discr. of wedge term) */
						derivi = cJ + (2 / pi) * (rhox[i] * calc_IC_1D(0, dx) - rhox[i - 1] * calc_IC_1D(-1, dx));

					//dipole (or bipolar) Jacobi relaxation method
					D.Area[i][0].p = D.Area[i][0].pold + jc_rlx * (1 / derivi) * gammai;
					D.Area[i - 1][0].p = D.Area[i - 1][0].pold - jc_rlx * (1 / derivi) * gammai;
				}

				// cavitation condition
				if (D.Area[i][0].p < 0.0) D.Area[i][0].p = 0.0;

				// residual at node i
				if (D.Area[i][0].p > 1e-5)
					rs = abs(cF - cA * D.Area[i + 1][0].p - cB * D.Area[i - 1][0].p - cE * D.Area[i][0].p);
				else rs = 0.0;

				// cumulative sums
				summ0 = summ0 + abs(D.Area[i][0].pold - D.Area[i][0].p); // for pressure error
				summ1 = summ1 + (D.Area[i][0].p + D.Area[i + 1][0].p);   // for pressure integral
				summ2 = summ2 + rs * rs;								 // for residual
				summ3 = summ3 + D.Area[i][0].p;							 // for pressure error

				// save current pressure value for next iteration
				D.Area[i][0].pold = D.Area[i][0].p;

			} // i loop

			perr = summ0 / summ3;	// pressure error on current iteration 
			pint_curr = 0.5 * summ1 * dx; // current value for pressure intergral
			rmsres = pow(summ2 / Nx, 0.5); // current RMS residual

			std::cout << "PRS ITER = " << lc_prs << "\tPRS INT = " << std::format("{:.5e}", pint_curr) << "\tHMIN = " <<
				std::format("{:.5e}", Hmin) << "\tPRS ERR = " << std::format("{:.5e}", perr) << "\tRMS RES = " << std::format("{:.5e}", rmsres) <<
				"\tJRA = " << std::format("{:.2f}", 1.0 * nd_cnt / Nx) << "\n";
			extfile << "PRS ITER = " << lc_prs << "\tPRS INT = " << std::format("{:.5e}", pint_curr) << "\tHMIN = " <<
				std::format("{:.5e}", Hmin) << "\tPRS ERR = " << std::format("{:.5e}", perr) << "\tRMS RES = " << std::format("{:.5e}", rmsres) <<
				"\tJRA = " << std::format("{:.2f}", 1.0 * nd_cnt / Nx) << "\n";

			/* check the pressure convergence criterion */
			if ((perr < prs_errtol)) {
				std::cout << "\n<< PRESSURES ARE CONVERGED\n";
				extfile << "\n<< PRESSURES ARE CONVERGED\n";
				f_prs = true;
			}
			else lc_prs = lc_prs + 1;

		} // while (pressure sweeps)

		// check error for pressure integral and adjust h0 
		if ((abs(pint_curr - (pi / 2)) / (pi / 2) < frc_errtol) && (rmsres <= res_errtol)) {

			std::cout << "<< FORCE BALANCE IS REACHED\n";
			extfile << "<< FORCE BALANCE IS REACHED\n";
			f_h0 = true;
		}
		else {
			// update h0 and continue calculations
			p1 = h0_rlx * ((pi / 2) - pint_curr) * a * a / RX;
			if (p1 > 0) gapstr = "HAS DECREASED"; else gapstr = "HAS INCREASED";
			if (abs(p1) < 1e-12) gapstr = "HASN'T CHANGED";
			D.h0 = D.h0 - p1;

			p1 = RX * D.h0 / (a * a);
			if (lc_h0 < maxlc_h0 + 1) {
				std::cout << "\n\tHO ITER = " << lc_h0 << "\tNEW H0 = " << std::format("{:.5e}", p1)
					<< "(" << std::format("{:.5e}", D.h0) << " m)" << "\t*** THE GAP BETWEEN THE SURFACES " << gapstr << "\n\n";
				extfile << "\n\tH0 ITER = " << lc_h0 << "\tNEW H0 = " << std::format("{:.5e}", p1)
					<< "(" << std::format("{:.5e}", D.h0) << " m)" << "\t*** THE GAP BETWEEN THE SURFACES " << gapstr << "\n\n";
			}

			// check the H0lim
			if (abs(p1) >= abs(H0lim)) {

				D.h0 = (a * a) * H0lim / RX;

				std::cout << "\n\t*** H0 IS EQUAL OR GREATER THAN THE SURFACE'S PROXIMITY IN DRY HERZIAN CONTACT\n";
				std::cout << "\t*** NEW VALUE FOR H0 " << std::format("{:.5e}", H0lim) << "\n";
				extfile << "\n\t*** H0 IS EQUAL OR GREATER THAN THE SURFACE'S PROXIMITY IN DRY HERZIAN CONTACT\n";
				extfile << "\t*** NEW VALUE FOR H0 " << std::format("{:.5e}", H0lim) << "\n\n";
			}

			lc_h0 = lc_h0 + 1;
		} // else
	} // while (force balance)
	extfile.close();

	// write the service information

	/* write the residuals in external file */
	extfile.open("lc_NDRES.txt", std::ios::out);

	for (i = i_beg; i < i_end; ++i) {

		calc_coef1D(i, D, B1, B2, eps, rhox, H, cA, cB, cE, cF, cJ, 0);

		x_curr = (-1.0 * D.kx) + i * D.dx;

		// residual at node i
		if (D.Area[i][0].p > 1e-5)
			rs = abs(cF - cA * D.Area[i + 1][0].p - cB * D.Area[i - 1][0].p - cE * D.Area[i][0].p);
		else rs = 0.0;

		extfile << std::format("{:.5f}", x_curr) << "\t" << std::format("{:.5e}", rs) << "\n";

	}
	extfile.close();

	/* write the values for JRA-criterions (Venner's and Huang's) in external files */
	extfile.open("lc_EH3.txt", std::ios::out);
	extfile1.open("lc_A1.txt", std::ios::out);
	extfile2.open("lc_A2.txt", std::ios::out);

	for (i = i_beg; i < i_end; ++i) {

		calc_coef1D(i, D, B1, B2, eps, rhox, H, cA, cB, cE, cF, cJ, 0);

		x_curr = (-1.0 * D.kx) + i * D.dx;

		//if ((x_curr >= -1) && (x_curr <= 1)) p1 = alpha * eps[i] * pow(H[i], 3) / (dx * dx); else p1 = 0.0;

		// (eH3)[i+1/2]
		eH3_i05 = 0.5 * (eps[i + 1] * pow(H[i + 1], 3) + eps[i] * pow(H[i], 3));
		// (eH3)[i-1/2]
		eH3_05i = 0.5 * (eps[i - 1] * pow(H[i - 1], 3) + eps[i] * pow(H[i], 3));

		p1 = alpha * eps[i] * pow(H[i], 3) / (dx * dx);
		A1 = (eH3_i05 + eH3_05i) / (dx * dx);
		A2 = (1 / (pi * dx)) * (rhox[i] * calc_IC_1D(0, dx) - rhox[i - 1] * calc_IC_1D(-1, dx));

		extfile << std::format("{:.5f}", x_curr) << "\t" << std::format("{:.5e}", p1) << "\n";
		extfile1 << std::format("{:.5f}", x_curr) << "\t" << std::format("{:.5e}", log10(A1)) << "\n";
		extfile2 << std::format("{:.5f}", x_curr) << "\t" << std::format("{:.5e}", log10(abs(A2))) << "\n";

	}
	extfile.close();
	extfile1.close();
	extfile2.close();

	/* free arrays for eps, rho */
	{
		delete[] eps;
		delete[] rhox;
	}



} // proc

void calc_pressure_field_1D(Domain& D, Body B1, Body B2)
{
	double
		* eps,		/* dimensionless epsilon		*/
		* rhox,		/* dimensionless rhostar_x		*/
		* rhoy,		/* dimensionless rhostar_y (for this case it equal to zero) */
		* H;		/* dimensionless thickness */

	int Nx, i;
	int lc_prs, lc_h0, maxlc_prs0, maxlc_prs1, maxlc_prs, maxlc_h0, wt_dst;
	int nd_cnt;

	double RX, RY, x_curr, eps_curr, rhox_curr, rhoy_curr, a;
	double dx, p1, um, prs, i_beg, i_end;
	double alpha, cA, cB, cE, cF, cJ;

	double h0_rlx, H_i, H_1i, H_i1, eH3_05i, eH3_i05;
	double T1, T2, gammai, derivi, ptol, swc, summ0, summ1, summ3, perr, perr_old, Hmin;
	double gs_rlx, jc_rlx, HM, H_2i, rhox_2i, k, pint_curr, prs_errtol;
	bool f_prs, f_h0, f_neg;

	double rs, summ2, rmsres, res_errtol, frc_errtol, A1;
	double H0, H1, H1k, H2, delta, A2, H0lim;
	string gapstr;

	std::fstream extfile;
	std::fstream extfile1;
	std::fstream extfile2;

	Nx = D.Nx;
	dx = D.dx;
	a = D.a;
	RX = pow(1 / B1.Rx + 1 / B2.Rx, -1);

	um = 1.0; // for steady state we don't use um though leave it in the formulas 
	alpha = D.ph * pow(D.a, 3) / (um * D.etaR * pow(RX, 2));

	H0lim = -0.596574; // limiting value for H0 from Herzian dry contact theory (Venner, 1991, p.80)
	
	// the minimum film thickness for line contact (empirical formula)
	// we use a (RX*RX/(a*a)), since our dimensionless expression for H is different than Hamrock's
	//HM = 1.714 * pow(D.W_HD, -0.128) * pow(D.U_HD, 0.694) * pow(D.G_HD, 0.568) * pow((RX / a), 2); // so big, I don't know why
	
	// relaxation coefficient for h0 (apply (Huang,2013)-method to Hamrock-formula for Hc)
	//h0_rlx = abs((-0.166) * RX * 2.922 * pow(D.W_HD, -1.166) * pow(D.U_HD, 0.692) * pow(D.G_HD, 0.47));

	/**** ADJUSTABLE PARAMETERS *******/
	
	HM = 0.01;			// minimum for dimensionless film thickness
	H0 = -0.5;			// initial value for dimensionless H0

	h0_rlx = 0.02;		// H0 relaxation factor (0.1...0.001)

	gs_rlx = 0.008;		// Gauss-Seidel relaxation factor (0.5...1.0)
	jc_rlx = 0.002;		// Jacobi relaxation factor (0.1...0.6)

	maxlc_h0  = 20;	// max H0 relaxation steps
	maxlc_prs = 3;		// max pressure sweeps before H0 updating
	
	prs_errtol = 1e-5;	// pressure error tolerance
	res_errtol = 1e-3;	// residual tolerance
	frc_errtol = 1e-3;  // force balance tolerance

	/**********************************/

	if (D.h0 == 0.0) D.h0 = H0 * (D.a * D.a) / pow(1 / B1.Rx + 1 / B2.Rx, -1);

	/* boundary conditions (zero-pressure) */
	D.Area[0][0].p = 0.0;
	D.Area[Nx - 1][0].p = 0.0;
	i_beg = 1;
	i_end = Nx - 1; //limits for i-index

	/* optional boundary conditions (zero pressure gradient) */
	//D.Area[1][0].p = 0.0;
	//D.Area[Nx - 2][0].p = 0.0;
	//i_beg = 1; i_end = Nx - 2; //limits for i-index

	eps  = new double[Nx];
	rhox = new double[Nx];
	H    = new double[Nx];

	lc_h0 = 2;
	f_h0 = false; // force balance flag
	
	extfile.open("LC_LOG.TXT", std::ios::out);
	extfile << "\t*** INITIAL VALUE FOR H0 = " << std::format("{:.5e}", D.h0 * RX / (a * a))
		<< " (" << std::format("{:.5e}", D.h0) << " m)\n\n";
	extfile << "\tH0 ITER = 1\n\n";
		

	// loop for force balancing
	while ((!f_h0) && (lc_h0 < maxlc_h0 + 2)) {

		lc_prs = 1;
		f_prs = false; // pressure convergence flag

		// loop for pressure calculation
		while ((!f_prs) && (lc_prs < maxlc_prs + 1)) {

			// update film thickness and lubricant properties
			calc_deform(D, B1, B2);
			calc_densviscNB(D, D.RHO, D.LSV);
			summ0 = 0.0; summ1 = 0.0; summ2 = 0.0; summ3 = 0.0; Hmin = 1e10;

			// calculate film thickness (dimensionless form)
			f_neg = false; // flag for negative Hmin
			for (i = 0; i < Nx; ++i) {

				x_curr = (-1.0 * D.kx) + i * D.dx;

				p1 = (RX / (a * a)) * (D.Area[i][0].d1 + D.Area[i][0].d2);
				H[i] = (RX / (a * a)) * D.h0 + 0.5 * pow(x_curr, 2) + p1;

				D.Area[i][0].h = ((a * a) / RX) * H[i];

				if (H[i] < Hmin) Hmin = H[i]; // update the minimum film thickness

				if ((Hmin < 0) && (!f_neg)) {
					std::cout << "\n\t*** HMIN = " << Hmin << " WAS DETECTED AT X = " << x_curr << "\n";
					extfile << "\n\t*** HMIN = " << Hmin << " WAS DETECTED AT X = " << x_curr << "\n";
					f_neg = true;
				}

			} // i loop

			// correct H0 if negative Hmin is found
			if (Hmin < HM) {
				
				delta = abs(Hmin) + HM;
				D.h0 = D.h0 + delta * (a * a) / RX;
			
				std::cout << "\t*** CORRECT H0 BY " << std::format("{:.5e}", delta) << ", NEW VALUE IS " <<
				std::format("{:.5e}", D.h0 * RX / (a * a)) << " (" << std::format("{:.5e}", D.h0) << " m)\n";
					
				extfile << "\t*** CORRECT H0 BY " << std::format("{:.5e}", delta) << ", NEW VALUE IS " <<
				std::format("{:.5e}", D.h0 * RX / (a * a)) << " (" << std::format("{:.5e}", D.h0) << " m)\n";

				// define new Hmin
				Hmin = 1e10;
				for (i = 0; i < Nx; ++i) {

					H[i] = H[i] + delta;

					if (H[i] < Hmin) Hmin = H[i];
				}

				std::cout << "\t*** AFTER CORRECTION HMIN = " << std::format("{:.5e}", Hmin) << "\n\n";
				extfile << "\t*** AFTER CORRECTION HMIN = " << std::format("{:.5e}", Hmin) << "\n\n";
			} // if

			// fill arrays for epsilon, rhostarx (dimensionless form)
			
			for (i = 0; i < Nx; ++i) {
									
				// calculate epsilon and rhostarx (rhoy_curr remains but is not used here)
				eps_curr = 0; rhox_curr = 0; rhoy_curr = 0;
				calc_rhoeps(i, 0, D, B1, B2, eps_curr, rhox_curr, rhoy_curr);

				eps[i] = eps_curr;
				rhox[i] = rhox_curr;

			} /* i loop */

			// calculate new pressure field (hybrid relaxation)
			wt_dst = 0; // type of wedge term discretization (0 - 1st order; 1 - 2nd order) 
			nd_cnt = 0; // number of nodes in "contact" zone (for jacobi relaxations)
			for (i = i_beg; i < i_end; ++i) {

				calc_coef1D(i, D, B1, B2, eps, rhox, H, cA, cB, cE, cF, cJ, wt_dst);
				// terms for Huang's conditions (if it used)
				eH3_i05 = 0.5 * (eps[i + 1] * pow(H[i + 1], 3) + eps[i] * pow(H[i], 3));
				eH3_05i = 0.5 * (eps[i - 1] * pow(H[i - 1], 3) + eps[i] * pow(H[i], 3));
				A1 = alpha * (eH3_i05 + eH3_05i) / (dx * dx);
				A2 = (1 / (pi * dx)) * (rhox[i] * calc_IC_1D(0, dx) - rhox[i - 1] * calc_IC_1D(-1, dx));
				//std::cout << i << "\t" << A2 / A1 << "\n";
				//if (A1 > 0.1 * abs(A2)) // Huang's, 2013 condition
				if (alpha * eps[i] * pow(H[i], 3) / (dx * dx) >= 0.3) // classic condition (Venner, 1991)
				//if (alpha * eps[i] * pow(H[i], 3) / dx >= 0.01) // conservative condition (Venner, 1991)
				{
					// Gauss-Seidel relaxation

					gammai = cF - cA * D.Area[i + 1][0].pold - cB * D.Area[i - 1][0].p - cE * D.Area[i][0].pold;

					if ((i > i_beg) && (wt_dst == 1))
						// use 2nd order wedge term in derivative
						derivi = cE + (1 / (pi * 2 * dx)) * (3 * rhox[i] * calc_IC_1D(0, dx) - 4 * rhox[i - 1] * calc_IC_1D(-1, dx) +
							rhox[i - 2] * calc_IC_1D(-2, dx));
					else
						// use 1st order wedge term in derivative
						derivi = cE + (1 / (pi * dx)) * (rhox[i] * calc_IC_1D(0, dx) - rhox[i - 1] * calc_IC_1D(-1, dx));

					D.Area[i][0].p = D.Area[i][0].pold + gs_rlx * (1 / derivi) * gammai;
				}
				else
				{
					// Jacobi relaxation
					nd_cnt = nd_cnt + 1;
					gammai = cF - cA * D.Area[i + 1][0].pold - cB * D.Area[i - 1][0].pold - cE * D.Area[i][0].pold;

					if ((i > i_beg) && (wt_dst == 1))
						// use 2nd order wedge term in derivative
						derivi = (cE - cB) + 1 / (pi * 2 * dx) * ((3 * rhox[i] + 4 * rhox[i - 1]) * (calc_IC_1D(0, dx) - calc_IC_1D(-1, dx)) -
							rhox[i - 2] * (calc_IC_1D(-1, dx) - calc_IC_1D(-2, dx)));
					else
						// use 1st order wedge term in derivative
						/* My term for derivative (only for 1st order discr. of wedge term) */
						//derivi = (cE - cB) + (calc_IC_1D(0, dx) - calc_IC_1D(1, dx)) / (pi * dx) * (rhox[i] + rhox[i - 1]);

						/* Venner term for derivative (only for 1st order discr. of wedge term) */
						derivi = cJ + (2 / pi) * (rhox[i] * calc_IC_1D(0, dx) - rhox[i - 1] * calc_IC_1D(-1, dx));
					
					//dipole (or bipolar) Jacobi relaxation method
					D.Area[i    ][0].p = D.Area[i    ][0].pold + jc_rlx * (1 / derivi) * gammai;
					D.Area[i - 1][0].p = D.Area[i - 1][0].pold - jc_rlx * (1 / derivi) * gammai;
				}

				// cavitation condition
				if (D.Area[i][0].p < 0.0) D.Area[i][0].p = 0.0;
				//if (D.Area[i][0].p < 0.0) D.Area[i][0].p = abs(D.Area[i][0].p);

				// residual at node i
				//rs = abs(cF - cA * D.Area[i + 1][0].p - cB * D.Area[i - 1][0].p - cE * D.Area[i][0].p);
				
				if (D.Area[i][0].p > 0.01)
					rs = abs(cF - cA * D.Area[i + 1][0].p - cB * D.Area[i - 1][0].p - cE * D.Area[i][0].p);
				else rs = 0.0;
				
				// cumulative sums
				summ0 = summ0 + abs(D.Area[i][0].pold - D.Area[i][0].p); // for pressure error
				summ1 = summ1 + (D.Area[i][0].p + D.Area[i + 1][0].p);   // for pressure integral
				summ2 = summ2 + rs * rs;								 // for residual
				summ3 = summ3 + D.Area[i][0].p;							 // for pressure error

				// save current pressure value for next iteration
				D.Area[i][0].pold = D.Area[i][0].p;

			} // i loop

			perr = summ0 / summ3;	// pressure error on current iteration 
			pint_curr = 0.5 * summ1 * dx; // current value for pressure intergral
			rmsres = pow(summ2 / Nx, 0.5); // current RMS residual

			std::cout << "PRS ITER = " << lc_prs << "\tPRS INT = " << std::format("{:.5e}", pint_curr) << "\tHMIN = " <<
				std::format("{:.5e}", Hmin) << "\tPRS ERR = " << std::format("{:.5e}", perr) << "\tRMS RES = " << std::format("{:.5e}", rmsres) << 
				"\tJRA = " << std::format("{:.2f}",1.0 * nd_cnt / Nx) << "\n";
		 	  extfile << "PRS ITER = " << lc_prs << "\tPRS INT = " << std::format("{:.5e}", pint_curr) << "\tHMIN = " <<
				std::format("{:.5e}", Hmin) << "\tPRS ERR = " << std::format("{:.5e}", perr) << "\tRMS RES = " << std::format("{:.5e}", rmsres) << 
			    "\tJRA = " << std::format("{:.2f}",1.0 * nd_cnt / Nx) << "\n";

			/* check the pressure convergence criterion */
			if ((perr < prs_errtol)) {
				std::cout << "\n<< PRESSURES ARE CONVERGED\n";
				  extfile << "\n<< PRESSURES ARE CONVERGED\n";
				f_prs = true;
			}
			else lc_prs = lc_prs + 1;

		} // while (pressure sweeps)

		// check error for pressure integral and adjust h0 
		if ((abs(pint_curr - (pi / 2)) / (pi / 2) < frc_errtol) && (rmsres <= res_errtol)) {

			std::cout << "<< FORCE BALANCE IS REACHED\n";
			  extfile << "<< FORCE BALANCE IS REACHED\n";
			f_h0 = true;
		}
		else {
			// update h0 and continue calculations
			p1 = h0_rlx * ((pi / 2) - pint_curr) * a * a / RX;
			if (p1 > 0) gapstr = "HAS DECREASED"; else gapstr = "HAS INCREASED";
			if (abs(p1) < 1e-12) gapstr = "HASN'T CHANGED";
			D.h0 = D.h0 - p1;
			
			p1 = RX * D.h0 / (a * a); 
			if (lc_h0 < maxlc_h0 + 1) {
				std::cout << "\n\tHO ITER = " << lc_h0 << "\tNEW H0 = " << std::format("{:.5e}", p1) 
					<< "(" << std::format("{:.5e}", D.h0) << " m)" << "\t*** THE GAP BETWEEN THE SURFACES " << gapstr << "\n\n";
				  extfile << "\n\tH0 ITER = " << lc_h0 << "\tNEW H0 = " << std::format("{:.5e}", p1) 
				    << "(" << std::format("{:.5e}", D.h0) << " m)" << "\t*** THE GAP BETWEEN THE SURFACES " << gapstr << "\n\n";
			}
			
			// check the H0lim
			/*
			if (p1 <= H0lim) {
				
				D.h0 = (a * a) * H0lim / RX;
				
				std::cout << "\n\t*** H0 REACHED A LIMIT VALUE (" << H0lim << ") CALCULATED FROM DRY HERZIAN CONTACT\n";
				std::cout << "\t*** NEW VALUE FOR H0 " << std::format("{:.5e}", H0lim) << "\n";
				  extfile << "\n\t*** H0 REACHED A LIMIT VALUE (" << H0lim << ") CALCULATED FROM DRY HERZIAN CONTACT\n";
				  extfile << "\t*** NEW VALUE FOR H0 " << std::format("{:.5e}", H0lim) << "\n\n";
			}
			*/

			lc_h0 = lc_h0 + 1;
		} // else
	} // while (force balance)
	extfile.close();
	
	// write the service information

	/* write the residuals in external file */
	extfile.open("LC_NDRES.txt", std::ios::out);
	
	for (i = i_beg; i < i_end; ++i) {
	
		calc_coef1D(i, D, B1, B2, eps, rhox, H, cA, cB, cE, cF, cJ, 0);

		x_curr = (-1.0 * D.kx) + i * D.dx;
		
		// residual at node i
		//rs = abs(cF - cA * D.Area[i + 1][0].p - cB * D.Area[i - 1][0].p - cE * D.Area[i][0].p);
				
		if (D.Area[i][0].p > 0.01)
			rs = abs(cF - cA * D.Area[i + 1][0].p - cB * D.Area[i - 1][0].p - cE * D.Area[i][0].p);
		else rs = 0.0;
		
		extfile << std::format("{:.5f}",x_curr) << "\t" << std::format("{:.5e}",rs) << "\n";
		
	}
	extfile.close();

	/*
	//write the values for JRA-criterions (Venner's and Huang's) in external files
	extfile.open("lc_EH3.txt", std::ios::out);
	extfile1.open("lc_A1.txt", std::ios::out);
	extfile2.open("lc_A2.txt", std::ios::out);

	for (i = i_beg; i < i_end; ++i) {

		calc_coef1D(i, D, B1, B2, eps, rhox, H, cA, cB, cE, cF, cJ, 0);

		x_curr = (-1.0 * D.kx) + i * D.dx;
		
		//if ((x_curr >= -1) && (x_curr <= 1)) p1 = alpha * eps[i] * pow(H[i], 3) / (dx * dx); else p1 = 0.0;
		
		// (eH3)[i+1/2]
		eH3_i05 = 0.5 * (eps[i + 1] * pow(H[i + 1], 3) + eps[i] * pow(H[i], 3));
		// (eH3)[i-1/2]
		eH3_05i = 0.5 * (eps[i - 1] * pow(H[i - 1], 3) + eps[i] * pow(H[i], 3));
		
		p1 = alpha * eps[i] * pow(H[i], 3) / (dx * dx);
		A1 = alpha * (eH3_i05 + eH3_05i) / (dx * dx);
		A2 = (1 / (pi * dx)) * (rhox[i] * calc_IC_1D(0, dx) - rhox[i - 1] * calc_IC_1D(-1, dx));

		extfile  << std::format("{:.5f}", x_curr) << "\t" << std::format("{:.5e}", log10(abs(p1))) << "\n";
		extfile1 << std::format("{:.5f}", x_curr) << "\t" << std::format("{:.5e}", log10(A1)) << "\n";
		extfile2 << std::format("{:.5f}", x_curr) << "\t" << std::format("{:.5e}", log10(abs(A2))) << "\n";

	}
	extfile.close();
	extfile1.close();
	extfile2.close();
	*/
	
	/* free arrays for eps, rho */
	{
		delete[] eps;
		delete[] rhox;
	}

	

} // proc



/********************* POINT CONTACT PROCEDURES ******************************/
void calc_coef2D(int i, int j, Domain D, Body B1, Body B2, double** eps, double** rhox, double** rhoy, double** H,
	double& cA, double& cB, double& cC, double& cD, double& cE, double& cF, int kval) {
	//****** calculate the discretization coefficients in [i,j]-point (point contact) ******/

	double dx, dy, alpha, theta, um, RX;
	double eH3_i05_j, eH3_05i_j, eH3_i_j05, eH3_i_05j;

	dx = D.dx;
	dy = D.dy;
	RX = pow(1 / B1.Rx + 1 / B2.Rx, -1);

	um = 1.0; // for steady state we don't use um though leave it in the formulas
	alpha = D.ph * pow(D.a, 3) / (um * D.etaR * pow(RX, 2));
	theta = D.a / D.b;

	// (eH3)[i+1/2][j]
	eH3_i05_j = 0.5 * (eps[i + 1][j] * pow(H[i + 1][j], 3) + eps[i][j] * pow(H[i][j], 3));
	// (eH3)[i-1/2][j]
	eH3_05i_j = 0.5 * (eps[i - 1][j] * pow(H[i - 1][j], 3) + eps[i][j] * pow(H[i][j], 3));
	// (eH3)[i][j+1/2]
	eH3_i_j05 = 0.5 * (eps[i][j + 1] * pow(H[i][j + 1], 3) + eps[i][j] * pow(H[i][j], 3));
	// (eH3)[i][j-1/2]
	eH3_i_05j = 0.5 * (eps[i][j - 1] * pow(H[i][j - 1], 3) + eps[i][j] * pow(H[i][j], 3));
	
	// standart discretization scheme 
	cA = alpha / (dx * dx) * eH3_i05_j;	// P[i+1][j] coeff
	cB = alpha / (dx * dx) * eH3_05i_j;	// P[i-1][j] coeff

	cC = (alpha * theta * theta / (D.dy * D.dy)) * eH3_i_j05; // P[i][j+1] coeff
	cD = (alpha * theta * theta / (D.dy * D.dy)) * eH3_i_05j; // P[i][j-1] coeff

	cE = -alpha / (dx * dx) * (eH3_05i_j + eH3_i05_j)
		-(alpha * theta * theta) / (dy * dy) * (eH3_i_05j + eH3_i_j05); // P[i][j] coeff

	// wedge term (or Couette term)
	if (kval == 0)
		// 1st order backward scheme
		cF = (1 / dx) * (rhox[i][j] * H[i][j] - rhox[i - 1][j] * H[i - 1][j]) +
		 (theta / dy) * (rhoy[i][j] * H[i][j] - rhoy[i][j - 1] * H[i][j - 1]);

	else {
		// 2nd order backward scheme 
		if ((i == 1) or (j == 1)) // anyway, for first lines we use 1st order scheme
			cF = (1 / dx) * (rhox[i][j] * H[i][j] - rhox[i - 1][j] * H[i - 1][j]) +
			 (theta / dy) * (rhoy[i][j] * H[i][j] - rhoy[i][j - 1] * H[i][j - 1]);
		else
			cF = (1 / (2 * dx)) * (3 * rhox[i][j] * H[i][j] - 4 * rhox[i - 1][j] * H[i - 1][j] + rhox[i - 2][j] * H[i - 2][j]) +
			 (theta / (2 * dy)) * (3 * rhoy[i][j] * H[i][j] - 4 * rhoy[i - 1][j] * H[i - 1][j] + rhoy[i - 2][j] * H[i - 2][j]);
	}
} //proc

void calc_pressure_field_2D_old(Domain& D, Body B1, Body B2)
{
	double
		** eps,			/* dimensionless epsilon		*/
		** rhox,		/* dimensionless rhostar_x		*/
		** rhoy,		/* dimensionless rhostar_y		*/
		** H,
		** M,
		* RHS,
		* Tx, * Ax, * X;
	int
		* Ti, * Tj, * Ai, * Ap;

	int	Nx, Ny, i, j, nzv, vdim, tc;
	double RX, RY, x_curr, y_curr, eps_curr, rhox_curr, rhoy_curr, a, b;
	double dx, dy, p1, um, prs;
	double alpha, theta, cA, cB, cC, cD, cE, cF;

	double Hmin0, A1;
	double eH3_05i_j, eH3_i05_j, eH3_i_05j, eH3_i_j05;
	double T1, T2, gammai, derivi, ptol, swc, summ0, summ1, perr, perr_old, Hmin;
	double h0_rlx, gs_rlx, jc_rlx, k, pint_curr, prs_errtol;
	bool f_prs, f_h0;
	int lc_prs, lc_h0, maxlc_prs, maxlc_h0, cnt;
	double K_00, K_10, K_20, prs_int;

	dbvector			Tx_vec;
	std::vector<int>	Ti_vec, Tj_vec;

	int status;
	void* Symbolic, * Numeric;

	Nx = D.Nx;
	Ny = D.Ny;
	dx = D.dx;
	dy = D.dy;
	a = D.a;
	b = D.b;
	RX = pow(1 / B1.Rx + 1 / B2.Rx, -1); /* it needs for dimesionless eps */
	RY = pow(1 / B1.Ry + 1 / B2.Ry, -1);
	k = 1.0339 * pow(RY / RX, 0.636); //ellipticity ratio

	//um = (B1.vel[0] + B2.vel[0]) / 2.0; /* m/s */
	um = 1.0; // we have to throw out it from entire equation 
	alpha = D.ph * pow(D.a, 3) / (um * D.etaR * pow(RX, 2));
	theta = a / b;
	Hmin0 = -0.569; // limiting value for H0 from Herzian dry contact theory (Venner, 1991, p.80)

	// relaxation coefficient for h0 (apply (Huang,2013)-method to Hamrock-formula for Hm)
	h0_rlx = (-0.067) * RX * 2.69 * pow(D.W_HD, -1.067) * pow(D.U_HD, 0.67) * pow(D.G_HD, 0.53) *
		(1 - 0.61 * exp(-0.73 * k));

	gs_rlx = 0.3;  // [0.3...0.9], Nurgat 1997 (lower for high load)
	jc_rlx = 0.000001;   //  0.1 or lower for high load, Nurgat, 1997

	maxlc_h0 = 10; // max h0 relaxation steps
	maxlc_prs = 150; // max pressure sweeps before H0-updating
	prs_errtol = 8e-5;


	/* intervals for inner grid nodes (exclude nodes for BC) */
	int i_beg = 1;			/* zero pressure on (Y=-1) side */
	//int i_end = Nx - 1;		/* zero pressure on (Y=1) side */
	int i_end = Nx - 1;		/* zero pressure gradient on (Y=1) side */

	int j_beg = 1;			/* zero pressure on (X=-1) side */
	int j_end = Ny - 1;		/* zero pressure gradient on (X=1) side */
	//int j_end = Ny - 1;		/* zero pressure on (X=1) side */

	//int IGNN = (i_end - i_beg) * (j_end - j_beg); /* number of inner grid nodes */
	int nrow = j_end - j_beg;					  /* number of inner grid nodes in row */
	int ncol = i_end - i_beg;					  /* number of inner grid nodes in row */

	/* set boundary conditions - Y axis */
	for (j = 0; j < Ny; ++j) {

		for (i = 0; i < i_beg; ++i) D.Area[i][j].p = 0.0;

		for (i = i_end; i < Nx; ++i) D.Area[i][j].p = 0.0;
	}

	/* set boundary conditions - X axis */
	for (i = 0; i < Nx; ++i) {

		for (j = 0; j < j_beg; ++j) D.Area[i][j].p = 0.0;

		for (j = j_end; j < Ny; ++j) D.Area[i][j].p = 0.0;
	}

	for (i = 0; j < i_beg; ++i)
		for (j = 0; j < Ny; ++j) D.Area[i][j].p = 0.0;

	D.Area[Nx - 1][0].p = 0.0;

	/* initialize dynamic arrays */
	eps = new double* [Nx];
	rhox = new double* [Nx];
	rhoy = new double* [Nx];
	H = new double* [Nx];

	for (i = 0; i < Ny; ++i) {
		eps[i] = new double[Ny];
		rhox[i] = new double[Ny];
		rhoy[i] = new double[Ny];
		H[i] = new double[Ny];
	}

	M = new double* [ncol];
	for (i = 0; i < ncol; ++i) {
		M[i] = new double[ncol];
		for (j = 0; j < ncol; ++j) M[i][j] = 0.0;
	}

	RHS = new double[ncol];

	/* tridiagonal Jacobi matrix for Y-line relaxation */
	//Ti = new int[3 * ncol - 2];
	//Tj = new int[3 * ncol - 2];
	//Tx = new double[3 * ncol - 2];

	/* matrices for solving the system */
	Ai = new int[3 * ncol - 2];
	Ap = new int[ncol + 1];
	Ax = new double[3 * ncol - 2];
	X = new double[ncol];

	lc_h0 = 1;
	f_h0 = false; // force balance flag

	K_00 = calc_IC_2D(0, 0, dx);
	K_10 = calc_IC_2D(1, 0, dx);
	K_20 = calc_IC_2D(2, 0, dx);

	output_data(D, B1, B2, 1, 3, '-', 1, 0, "E:\\_work\\EHL-prg\\TEHL\\TEHL-v2\\_output\\PRS_DL");

	// loop for force balancing
	while ((lc_h0 < maxlc_h0)) {

		lc_prs = 1;
		f_prs = false; // pressure convergence flag

		// loop for pressure calculation
		while ((lc_prs < maxlc_prs)) {

			/* update film thickness and lubricant properties */
			calc_deform(D, B1, B2);
			calc_physical_fields(D, D.RHO, D.LSV, D.GNV, D.THC, D.HCP, B1.vel, B2.vel);

			/* fill arrays for epsilon, rho and thickness (dimensionless form) */
			for (i = 0; i < Nx; ++i) {

				for (j = 0; j < Ny; ++j) {

					/* dimensionless (i,j)-position */
					x_curr = (-1.0 * D.kx) + i * dx;
					y_curr = (-1.0 * D.ky) + j * dy;

					eps_curr = 0; rhox_curr = 0; rhoy_curr = 0;
					calc_rhoeps(i, j, D, B1, B2, eps_curr, rhox_curr, rhoy_curr);

					eps[i][j] = eps_curr;
					rhox[i][j] = rhox_curr;
					rhoy[i][j] = rhoy_curr;

					/* H[i][j] */
					p1 = (RX / (a * a)) * (D.Area[i][j].d1 + D.Area[i][j].d2);
					H[i][j] = (RX / (a * a)) * D.h0 + 0.5 * pow(x_curr, 2) +
						0.5 * (RX / RY) * (1 / pow(theta, 2)) * pow(y_curr, 2) + p1;
					D.Area[i][j].h = H[i][j] * ((a * a) / RX);

					//std::cout << "eps[" << i << "][" << j << "] = " << eps_curr << "\n";
					//std::cout << "rhox[" << i << "][" << j << "] = " << rhox_curr << "\n";
					//std::cout << "rhoy[" << i << "][" << j << "] = " << rhoy_curr << "\n";
					//std::cout << "H[" << i << "][" << j << "] = " << H[i][j] << "\n";

				} // j loop
			} // i loop

			/* calculate new pressure field (Y-line relaxation) */
			//output_data(D, B1, B2, 1, 3, '-', 1, 0, "E:\\_work\\EHL-prg\\TEHL\\TEHL-v2\\_output\\PRS_DL");


			/* main loop */
			for (j = j_beg; j < j_end; ++j) {

				cnt = 0; nzv = 0; tc = 0;

				for (i = i_beg; i < i_end; ++i) {

					/* (eH3)[i+1/2][j] */
					eH3_i05_j = 0.5 * (eps[i + 1][j] * pow(H[i + 1][j], 3) + eps[i][j] * pow(H[i][j], 3));
					/* (eH3)[i-1/2][j] */
					eH3_05i_j = 0.5 * (eps[i - 1][j] * pow(H[i - 1][j], 3) + eps[i][j] * pow(H[i][j], 3));
					/* (eH3)[i][j+1/2] */
					eH3_i_j05 = 0.5 * (eps[i][j + 1] * pow(H[i][j + 1], 3) + eps[i][j] * pow(H[i][j], 3));
					/* (eH3)[i][j-1/2] */
					eH3_i_05j = 0.5 * (eps[i][j - 1] * pow(H[i][j - 1], 3) + eps[i][j] * pow(H[i][j], 3));

					/* standart discretization schema */
					cA = alpha / (dx * dx) * eH3_i05_j;	// P[i+1][j] coeff
					cB = alpha / (dx * dx) * eH3_05i_j;	// P[i-1][j] coeff

					cC = (alpha * theta * theta / (D.dy * D.dy)) * eH3_i_j05; // P[i][j+1] coeff
					cD = (alpha * theta * theta / (D.dy * D.dy)) * eH3_i_05j; // P[i][j-1] coeff

					cE = -alpha / (dx * dx) * (eH3_05i_j + eH3_i05_j) -
						(alpha * theta * theta) / (dy * dy) * (eH3_i_05j + eH3_i_j05);	// P[i] coeff

					/* backward (1st order) */
					cF = (1 / dx) * (rhox[i][j] * H[i][j] - rhox[i - 1][j] * H[i - 1][j]) +
						(theta / dy) * (rhoy[i][j] * H[i][j] - rhox[i][j - 1] * H[i][j - 1]); // RHO-term

					/* write the Jacobi matrix */
					p1 = cE - (2 / (pi * pi * dx) * (rhox[i][j] * K_00 - rhox[i - 1][j] * K_10) +
						2 * theta / (pi * pi * dy) * (rhoy[i][j] * K_00 - rhoy[i - 1][j] * K_10));

					M[cnt][cnt] = p1;
					Ti_vec.push_back(cnt);
					Tj_vec.push_back(cnt);
					Tx_vec.push_back(p1);

					//Ti[tc] = cnt;
					//Tj[tc] = cnt;
					//Tx[tc] = p1;

					//std::cout << "Ti[" << tc << "] = " << Ti[tc] << "\n";
					//std::cout << "Tj[" << tc << "] = " << Tj[tc] << "\n";
					//std::cout << "Tx[" << tc << "] = " << Tx[tc] << "\n";

					//tc = tc + 1;
					nzv = nzv + 1;
					//std::cout << "M[" << cnt << "][" << cnt << "] = " << p1 << "\n";


					if (i < (i_end - 1)) {
						p1 = cA - (2 / (pi * pi * dx) * (rhox[i][j] * K_10 - rhox[i - 1][j] * K_20) +
							2 * theta / (pi * pi * dy) * (rhoy[i][j] * K_10 - rhoy[i - 1][j] * K_20));

						M[cnt][cnt + 1] = p1;
						Ti_vec.push_back(cnt);
						Tj_vec.push_back(cnt + 1);
						Tx_vec.push_back(p1);

						//Ti[tc] = cnt;
						//Tj[tc] = cnt + 1;
						//Tx[tc] = p1;

						//std::cout << "Ti[" << tc << "] = " << Ti[tc] << "\n";
						//std::cout << "Tj[" << tc << "] = " << Tj[tc] << "\n";
						//std::cout << "Tx[" << tc << "] = " << Tx[tc] << "\n";

						//tc = tc + 1;
						nzv = nzv + 1;
						//std::cout << "M[" << cnt << "][" << cnt+1 << "] = " << p1 << "\n";

					}

					if (i > i_beg) {
						p1 = cB - (2 / (pi * pi * dx) * (rhox[i][j] * K_10 - rhox[i - 1][j] * K_00) +
							2 * theta / (pi * pi * dy) * (rhoy[i][j] * K_10 - rhoy[i - 1][j] * K_00));

						M[cnt][cnt - 1] = p1;
						Ti_vec.push_back(cnt);
						Tj_vec.push_back(cnt - 1);
						Tx_vec.push_back(p1);
						//Ti[tc] = cnt;
						//Tj[tc] = cnt - 1;
						//Tx[tc] = p1;
						//std::cout << "Ti[" << tc << "] = " << Ti[tc] << "\n";
						//std::cout << "Tj[" << tc << "] = " << Tj[tc] << "\n";
						//std::cout << "Tx[" << tc << "] = " << Tx[tc] << "\n";

						//tc = tc + 1;
						nzv = nzv + 1;
						//std::cout << "M[" << cnt << "][" << cnt-1 << "] = " << p1 << "\n";

					}

					/* write the RHS vector */
					//A1 = cA * D.Area[i + 1][j].pold + cB * D.Area[i - 1][j].pold + cC * D.Area[i][j + 1].pold + cD * D.Area[i][j - 1].pold + cE * D.Area[i][j].pold;

					if (alpha * eps[i][j] * pow(H[i][j], 3) / (dx * dx) < 0.3) { // contact zone, Jacobi residual

						RHS[cnt] = cA * D.Area[i + 1][j].pold + cB * D.Area[i - 1][j].pold + cC * D.Area[i][j + 1].pold +
							cD * D.Area[i][j - 1].pold + cE * D.Area[i][j].pold - cF;
					}
					else { // non-contact zone, Gauss-Seidel residual

						RHS[cnt] = cA * D.Area[i + 1][j].pold + cB * D.Area[i - 1][j].pold + cC * D.Area[i][j + 1].pold +
							cD * D.Area[i][j - 1].p + cE * D.Area[i][j].pold - cF;

					}
					//std::cout << "RHS[" << cnt << "] = " << RHS[cnt] << "\n";
					//std::cout << "eps[" << i << "][" << j << "] / dxdx = " << alpha * eps[i][j] * pow(H[i][j],3) / (dx * dx) << "\t WEDGE TERM = " << cF << "\n";

					cnt = cnt + 1;
				} // i loop

				int vdim = Ti_vec.size();

				int* Ti, * Tj;
				double* Tx;

				Ti = new int[vdim];
				Tj = new int[vdim];
				Tx = new double[vdim];

				for (i = 0; i < vdim; ++i) {

					Ti[i] = Ti_vec[i];
					Tj[i] = Tj_vec[i];
					Tx[i] = Tx_vec[i];
				}

				Ti_vec.clear();
				Tj_vec.clear();
				Tx_vec.clear();

				/* convert triplet form of matrix to compressed sparse column form */
				status = umfpack_di_triplet_to_col(ncol, ncol, nzv, Ti, Tj, Tx, Ap, Ai, Ax, NULL);

				/* free Ti, Tj, Tx */
				delete[] Ti; 	delete[] Tj; 	delete[] Tx;

				/* symbolic and numeric factorization */
				status = umfpack_di_symbolic(ncol, ncol, Ap, Ai, Ax, &Symbolic, NULL, NULL);

				status = umfpack_di_numeric(Ap, Ai, Ax, Symbolic, &Numeric, NULL, NULL);

				umfpack_di_free_symbolic(&Symbolic);

				/* solve the linear system */
				status = umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, X, RHS, Numeric, NULL, NULL);

				/* free the memory associated with the numeric factorization */
				umfpack_di_free_numeric(&Numeric);

				/* update the pressure field in non-contact zone */
				for (i = i_beg; i < i_end; ++i) {

					if (alpha * eps[i][j] * pow(H[i][j], 3) / (dx * dx) > 0.3) { // for node in non-contact zone

						D.Area[i][j].p = D.Area[i][j].pold + gs_rlx * X[i - i_beg];

						if (D.Area[i][j].p < 0.0) D.Area[i][j].p = 0.0; // cavitation condition	
					}

					else // for node in contact zone we just save the addition
						D.Area[i][j].dp = jc_rlx * X[i - i_beg];

					//std::cout << "X[" << i-i_beg << "] = " << X[i-i_beg] << "\n";
				}

			} // j loop

			/* update pressure for all nodes in contact zone */
			prs_int = 0.0;
			for (i = i_beg; i < i_end; ++i)
				for (j = j_beg; j < j_end; ++j) {

					if (alpha * eps[i][j] * pow(H[i][j], 3) / (dx * dx) <= 0.3) {

						D.Area[i][j].p = D.Area[i][j].pold + D.Area[i][j].dp;
						D.Area[i][j].dp = 0.0;
						if (D.Area[i][j].p < 0.0) D.Area[i][j].p = 0.0; // cavitation condition
					}


					D.Area[i][j].pold = D.Area[i][j].p;
					prs_int = prs_int + D.Area[i][j].p;
				}

			lc_prs = lc_prs + 1;
			output_data(D, B1, B2, 1, 3, '-', 1, 0, "E:\\_work\\EHL-prg\\TEHL\\TEHL-v2\\_output\\PRS_DL");
		} //while prs


		D.h0 = D.h0 + 0.00000001 * (2.0944 - prs_int * dx * dy);
		//D.h0 = 0.0;
		std::cout << "ITER = " << lc_h0 << "\tUPDATE H0, NEW VALUE " << D.h0 <<
			"\tPRS_INT = " << prs_int * dx * dy << "\n";
		lc_h0 = lc_h0 + 1;
	} // while h0

	  /* free arrays for eps, rho, H */

	for (i = 0; i < Nx; ++i) {
		delete[] eps[i];
		delete[] rhox[i];
		delete[] rhoy[i];
		delete[] H[i];

	}

	delete[] eps;
	delete[] rhox;
	delete[] rhoy;
	delete[] H;

	for (i = 0; i < ncol; ++i) delete[] M[i];
	delete[] M;

	delete[] RHS;

	output_data(D, B1, B2, 1, 3, '-', 1, 0, "E:\\_work\\EHL-prg\\TEHL\\TEHL-v2\\_output\\PRS_DL");
} /* proc */

void calc_pressure_field_2D_LNRLX(Domain& D, Body B1, Body B2)
{
	double
		** eps,			/* dimensionless epsilon		*/
		** rhox,		/* dimensionless rhostar_x		*/
		** rhoy,		/* dimensionless rhostar_y		*/
		**H,
		**M,
		*RHS,
		* Tx, *Ax, *X;
	int
		* Ti, * Tj, * Ai, * Ap;

	int	Nx, Ny, i, j, nzv, vdim, tc;
	double RX, RY, x_curr, y_curr, eps_curr, rhox_curr, rhoy_curr, a, b;
	double dx, dy, p1, um, prs;
	double alpha, theta, cA, cB, cC, cD, cE, cF;

	double Hmin0, A1;
	double eH3_05i_j, eH3_i05_j, eH3_i_05j, eH3_i_j05;
	double T1, T2, gammai, derivi, ptol, swc, summ0, summ1, perr, perr_old, Hmin;
	double h0_rlx, gs_rlx, jc_rlx, k, pint_curr, prs_errtol;
	bool f_prs, f_h0;
	int lc_prs, lc_h0, maxlc_prs0, maxlc_prs1, maxlc_prs, maxlc_h0, cnt;
	double K_00, K_10, K_20, prs_int, rmsres, prs_err;
	double H0, H1, H1k, H2, res_errtol;

	dbvector			Tx_vec;
	std::vector<int>	Ti_vec, Tj_vec;

	int status;
	void* Symbolic, * Numeric;

	Nx = D.Nx;
	Ny = D.Ny;
	dx = D.dx;
	dy = D.dy;
	a = D.a;
	b = D.b;
	RX = pow(1 / B1.Rx + 1 / B2.Rx, -1); /* it needs for dimesionless eps */
	RY = pow(1 / B1.Ry + 1 / B2.Ry, -1);
	k = 1.0339 * pow(RY / RX, 0.636); //ellipticity ratio

	//um = (B1.vel[0] + B2.vel[0]) / 2.0; /* m/s */
	um = 1.0; // we have to throw out it from entire equation 
	alpha = D.ph * pow(D.a, 3) / (um * D.etaR * pow(RX, 2));
	theta = a / b;
	Hmin0 = -0.569; // limiting value for H0 from Herzian dry contact theory (Venner, 1991, p.80)

	// relaxation coefficient for h0 (apply (Huang,2013)-method to Hamrock-formula for Hm)
	//h0_rlx = (-0.067) * RX * 2.69 * pow(D.W_HD, -1.067) * pow(D.U_HD, 0.67) * pow(D.G_HD, 0.53) * (1 - 0.61 * exp(-0.73 * k));

			 std::fstream extfile1;
			 extfile1.open("pc_log.txt", std::ios::out);
			 extfile1 << "PRS ITER\tPRS ERR\tPRS INT\tRMSRES\tHmin\n\n";


	/**** ADJUSTABLE PARAMETERS *******/

	h0_rlx = 0.5;		// H0 relaxation factor

	gs_rlx = 0.07;		// Gauss-Seidel relaxation factor
	jc_rlx = 0.03;		// Jacobi relaxation factor

	maxlc_h0   = 35;		// max H0 relaxation steps
	maxlc_prs0 = 3;	// max pressure sweeps before H0 updating
	prs_errtol = 1e-5;	// pressure error tolerance
	res_errtol = 1e-1;	// pressure error tolerance

	H1 = 0.4;			// H-level for decrease the H0 relaxation faclor
	H1k = 0.1;			// decreasing factor (< 1)
	maxlc_prs1 = 3;		// max pressure sweeps before H0 updating (after decreasing)

	H2 = 0.05;			// H-level for freeze the H0 changes (h0_rlx = 0)

	/**********************************/
			
	/* intervals for inner grid nodes (exclude nodes for BC) */
	int i_beg = 1;			/* zero pressure on (X=-1) side */
	int i_end = Nx - 1;		/* zero pressure on (X=1) side */
	//int i_end = Nx - 2;	/* zero pressure gradient on (X=1) side */

	int j_beg = 1;			/* zero pressure on (Y=-1) side */
	//int j_end = Ny - 2;	/* zero pressure gradient on (Y=1) side */
	int j_end = Ny - 1;		/* zero pressure on (Y=1) side */

	int IGNN = (i_end - i_beg) * (j_end - j_beg); /* number of inner grid nodes */
	int nrow = i_end - i_beg;					  /* number of inner grid nodes along X-axis */
	int ncol = j_end - j_beg;					  /* number of inner grid nodes along Y-axis */

	/* set boundary conditions - Y axis */
	for (j = 0; j < Ny; ++j) {
		
		for (i = 0; i < i_beg; ++i) D.Area[i][j].p = 0.0;

		for (i = i_end; i < Nx; ++i) D.Area[i][j].p = 0.0;
	}
	
	/* set boundary conditions - X axis */
	for (i = 0; i < Nx; ++i) {

		for (j = 0; j < j_beg; ++j) D.Area[i][j].p = 0.0;

		for (j = j_end; j < Ny; ++j) D.Area[i][j].p = 0.0;
	}
	
	for (i = 0; j < i_beg; ++i)
		for (j = 0; j < Ny; ++j) D.Area[i][j].p = 0.0;

	/* initialize dynamic arrays */
	{
		eps  = new double* [Nx];
		rhox = new double* [Nx];
		rhoy = new double* [Nx];
		H	 = new double* [Nx];
		for (i = 0; i < Ny; ++i) {
			eps[i]  = new double[Ny];
			rhox[i] = new double[Ny];
			rhoy[i] = new double[Ny];
			H[i]	= new double[Ny];
		}

		// tridiagonal Jacobi matrix for nodes on one X-line
		M = new double* [nrow];
		for (i = 0; i < nrow; ++i) {
			M[i] = new double[nrow];
			for (j = 0; j < nrow; ++j) M[i][j] = 0.0;
		}

		RHS = new double[nrow];

		// triplets of Jacobi matrix
		Ti = new int[3 * nrow - 2];
		Tj = new int[3 * nrow - 2];
		Tx = new double[3 * nrow - 2];

		// matrices and vector for UMFPACK solving procedure
		Ai = new int[3 * nrow - 2];
		Ap = new int[nrow + 1];
		Ax = new double[3 * nrow - 2];
		X  = new double[nrow];
	}
	
	lc_h0 = 1;
	f_h0 = false; // force balance flag

	K_00 = calc_IC_2D(0, 0, dx);
	K_10 = calc_IC_2D(1, 0, dx);
	K_20 = calc_IC_2D(2, 0, dx);

	//output_data(D, B1, B2, 1, 3, '-', 1, 0, "E:\\_work\\EHL-prg\\TEHL\\TEHL-v2\\_output\\PRS_DL");

	maxlc_prs = maxlc_prs0;
	// loop for force balancing
	while ((lc_h0 < maxlc_h0)) {

		lc_prs = 1;
		f_prs = false; // pressure convergence flag
				
		// loop for pressure calculation
		while ((lc_prs < maxlc_prs)) {
			
			/* update film thickness and lubricant properties */
			calc_deform(D, B1, B2);
			//calc_physical_fields(D, D.RHO, D.LSV, D.GNV, D.THC, D.HCP, B1.vel, B2.vel);
			calc_densviscNB(D, D.RHO, D.LSV);

			// fill arrays for epsilon, rho and thickness (dimensionless form)
			for (i = 0; i < Nx; ++i) {

				for (j = 0; j < Ny; ++j) {

					/* dimensionless (i,j)-position */
					x_curr = (-1.0 * D.kx) + i * dx;
					y_curr = (-1.0 * D.ky) + j * dy;

					eps_curr = 0; rhox_curr = 0; rhoy_curr = 0;
					calc_rhoeps(i, j, D, B1, B2, eps_curr, rhox_curr, rhoy_curr);

					eps[i][j]  = eps_curr;
					rhox[i][j] = rhox_curr;
					rhoy[i][j] = rhoy_curr;

					/* H[i][j] */
					p1 = (RX / (a * a)) * (D.Area[i][j].d1 + D.Area[i][j].d2);
					H[i][j] = (RX / (a * a)) * D.h0 + 0.5 * pow(x_curr, 2) +
						0.5 * (RX / RY) * (1 / pow(theta, 2)) * pow(y_curr, 2) + p1;
					
					D.Area[i][j].h = H[i][j] * ((a * a) / RX);
					
					//std::cout << "eps[" << i << "][" << j << "] = " << eps_curr << "\n";
					//std::cout << "rhox[" << i << "][" << j << "] = " << rhox_curr << "\n";
					//std::cout << "rhoy[" << i << "][" << j << "] = " << rhoy_curr << "\n";
					//std::cout << "H[" << i << "][" << j << "] = " << H[i][j] << "\n";

				} // j loop
			} // i loop

			/* calculate new pressure field (line relaxation) */
			
			//output_data(D, B1, B2, 1, 3, '-', 1, 0, "E:\\_work\\EHL-prg\\TEHL\\TEHL-v2\\_output\\PRS_DL");
			
			/* main loop (relaxation along X-axis) */
			for (j = j_beg; j < j_end; ++j) {
				
				cnt = 0; nzv = 0; tc = 0;
				
				for (i = i_beg; i < i_end; ++i) {

					calc_coef2D(i, j, D, B1, B2, eps, rhox, rhoy, H, cA, cB, cC, cD, cE, cF, 0);
					
					// adding wedge term contributons (semi-system approach, Liu, 2006)
					 //cE, cA, cB is needed for line relaxation along X */
					//cE = cE - (1 / dx) * (rhox[i][j] * K_00 - rhox[i - 1][j] * K_10);
					//cA = cA - (1 / dx) * (rhox[i][j] * K_10 - rhox[i - 1][j] * K_20);
					//cB = cB - (1 / dx) * (rhox[i][j] * K_10 - rhox[i - 1][j] * K_00);
					//cF = -(cC + cD)  
					//	+ (1 / dx) * (rhox[i][j] * (H[i][j] - (K_10 * D.Area[i - 1][j].pold + K_00 * D.Area[i][j].pold + K_10 * D.Area[i + 1][j].pold)) - 
					//	rhox[i - 1][j] * (H[i - 1][j] - (K_00 * D.Area[i - 1][j].pold + K_10 * D.Area[i][j].pold + K_20 * D.Area[i + 1][j].pold)));

					// write the Jacobi matrix 
					p1 = cE - (2  / (pi * pi * dx) * (rhox[i][j] * K_00 - rhox[i - 1][j] * K_10) +
					    2 * theta / (pi * pi * dy) * (rhoy[i][j] * K_00 - rhoy[i - 1][j] * K_10));
					
					M[cnt][cnt] = p1;
					Ti[tc] = cnt;
					Tj[tc] = cnt;
					Tx[tc] = p1;
					
					tc = tc + 1;
					nzv = nzv + 1;
					
					//std::cout << "Ti[" << tc << "] = " << Ti[tc] << "\n";
					//std::cout << "Tj[" << tc << "] = " << Tj[tc] << "\n";
					//std::cout << "Tx[" << tc << "] = " << Tx[tc] << "\n";
					//std::cout << "M[" << cnt << "][" << cnt << "] = " << p1 << "\n";
					

					if (i < (i_end - 1)) {
						p1 = cA - (2 / (pi * pi * dx) * (rhox[i][j] * K_10 - rhox[i - 1][j] * K_20) +
						   2 * theta / (pi * pi * dy) * (rhoy[i][j] * K_10 - rhoy[i - 1][j] * K_20));
						
						M[cnt][cnt + 1] = p1;
						Ti[tc] = cnt;
						Tj[tc] = cnt + 1;
						Tx[tc] = p1;
						
						tc = tc + 1;
						nzv = nzv + 1;
						
						//std::cout << "Ti[" << tc << "] = " << Ti[tc] << "\n";
						//std::cout << "Tj[" << tc << "] = " << Tj[tc] << "\n";
						//std::cout << "Tx[" << tc << "] = " << Tx[tc] << "\n";
						//std::cout << "M[" << cnt << "][" << cnt+1 << "] = " << p1 << "\n";
						
					}
					
					if (i > i_beg) {
						p1 = cB - (2 / (pi * pi * dx) * (rhox[i][j] * K_10 - rhox[i - 1][j] * K_00) +
						   2 * theta / (pi * pi * dy) * (rhoy[i][j] * K_10 - rhoy[i - 1][j] * K_00));
						
						M[cnt][cnt - 1] = p1;
						Ti[tc] = cnt;
						Tj[tc] = cnt - 1;
						Tx[tc] = p1;
						
						tc = tc + 1;
						nzv = nzv + 1;
						
						//std::cout << "Ti[" << tc << "] = " << Ti[tc] << "\n";
						//std::cout << "Tj[" << tc << "] = " << Tj[tc] << "\n";
						//std::cout << "Tx[" << tc << "] = " << Tx[tc] << "\n";
						//std::cout << "M[" << cnt << "][" << cnt-1 << "] = " << p1 << "\n";
						
					}
					
					// write the RHS vector
					
					if (alpha * eps[i][j] * pow(H[i][j],3) / (dx * dx) < 0.3) { // contact zone, Jacobi residual
					//if (alpha * eps[i][j] * pow(H[i][j], 3) / dx < 0.01) { // contact zone, Jacobi residual
						RHS[cnt] = -cF + (cA * D.Area[i + 1][j].pold + cB * D.Area[i - 1][j].pold + cC * D.Area[i][j + 1].pold +
							cD * D.Area[i][j - 1].pold + cE * D.Area[i][j].pold);
					}
					else { // non-contact zone, Gauss-Seidel residual

						RHS[cnt] = -cF + (cA * D.Area[i + 1][j].pold + cB * D.Area[i - 1][j].pold + cC * D.Area[i][j + 1].pold +
							cD * D.Area[i][j - 1].p + cE * D.Area[i][j].pold);

					}
					//std::cout << "RHS[" << cnt << "] = " << RHS[cnt] << "\n";
					//std::cout << "eps[" << i << "][" << j << "]/dxdx = " << alpha * eps[i][j] * pow(H[i][j],3) / (dx * dx) << "\t WEDGE TERM = " << cF << "\n";
				
					cnt = cnt + 1;
				} // i loop
				
				// write Jacobi matrix to the text file (if needed)
				std::fstream extfile;
				extfile.open("_M.TXT", std::ios::out);
				for (int k = 0; k < nrow; ++k)
					for (int l = 0; l < nrow; ++l) {
						if (l != (nrow - 1)) extfile << M[k][l] << "\t"; else extfile << M[k][l] << "\n";
					}

				extfile.close();

				/* convert triplet form of matrix to compressed sparse column form */
				status = umfpack_di_triplet_to_col(nrow, nrow, nzv, Ti, Tj, Tx, Ap, Ai, Ax, NULL);

				/* symbolic and numeric factorization */
				status = umfpack_di_symbolic(nrow, nrow, Ap, Ai, Ax, &Symbolic, NULL, NULL);
				
				status = umfpack_di_numeric(Ap, Ai, Ax, Symbolic, &Numeric, NULL, NULL);
				
				umfpack_di_free_symbolic(&Symbolic);

				/* solve the linear system */
				status = umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, X, RHS, Numeric, NULL, NULL);

				/* free the memory associated with the numeric factorization */
				umfpack_di_free_numeric(&Numeric);
	
				/* update the pressure field in non-contact zone */
				for (int ii = i_beg; ii < i_end; ++ii) {
										 
					if (alpha * eps[ii][j] * pow(H[ii][j],3) / (dx * dx) > 0.3) { // for node in non-contact zone
					//if (alpha * eps[ii][j] * pow(H[ii][j], 3) / dx > 0.01) { // for node in non-contact zone

						D.Area[ii][j].p = D.Area[ii][j].pold + gs_rlx * X[ii - i_beg];
						
						if (D.Area[ii][j].p < 0.0) D.Area[ii][j].p = 0.0; // cavitation condition	
					}

					else // for node in contact zone we just save the addition
						D.Area[ii][j].dp = jc_rlx * X[ii - i_beg];
										
					//std::cout << "X[" << ii-i_beg << "] = " << X[ii-i_beg] << "\n";
				}
			
			} // j loop

			// update pressure for all nodes in contact zone
			prs_int = 0.0;
			for (i = i_beg; i < i_end; ++i) 
				for (j = j_beg; j < j_end; ++j) {
					
					if (alpha * eps[i][j] * pow(H[i][j], 3) / (dx * dx) <= 0.3) {
					//if (alpha * eps[i][j] * pow(H[i][j], 3) / dx <= 0.01) {
					
						D.Area[i][j].p = D.Area[i][j].pold + D.Area[i][j].dp;
						D.Area[i][j].dp = 0.0;
						
						if (D.Area[i][j].p < 0.0) D.Area[i][j].p = 0.0; // cavitation condition	
					}
					
					prs_int = prs_int + D.Area[i][j].p;
					
				}
				
			prs_int = prs_int * dx * dy; // dimensionless pressure integral
			
			// calculate the pressure error and RMS-residual 
			rmsres = 0.0; Hmin = 1e6; summ0 = 0.0;
			for (i = i_beg; i < i_end; ++i)
				for (j = j_beg; j < j_end; ++j) {
					
					calc_coef2D(i, j, D, B1, B2, eps, rhox, rhoy, H, cA, cB, cC, cD, cE, cF, 0);

					p1 = cA * D.Area[i + 1][j].p + cB * D.Area[i - 1][j].p + cC * D.Area[i][j + 1].p +
						cD * D.Area[i][j - 1].p + cE * D.Area[i][j].p - cF;
					rmsres = rmsres + p1 * p1;
					
									
					if (H[i][j] < Hmin) {
					
						Hmin = H[i][j]; // update the minimum film thickness
						//std::cout << "new Hmin detected " << Hmin << " at point i = " << i << ", j = " << j << "\n";
						
					}
					summ0 = summ0 + abs(D.Area[i][j].pold - D.Area[i][j].p);
					
					// save current pressure as old for next sweep
					D.Area[i][j].pold = D.Area[i][j].p;
				}
			
			// RMSRES at current iteration
			rmsres = pow(rmsres / IGNN, 0.5);
			
			// pressure error at current iteration
			prs_err = summ0 * dx * dy / prs_int; 
			
			std::cout << "PRS_ITER = " << lc_prs << "\tPRS_INT = " << prs_int <<
				"\tPRS_ERR = " << std::format("{:.5e}",prs_err) << "\tRMSRES = " << std::format("{:.5e}",rmsres) << "\tHMIN = " << Hmin << "\n";
			extfile1  << "PRS_ITER = " << lc_prs << "\tPRS_INT = " << prs_int <<
				"\tPRS_ERR = " << std::format("{:.5e}", prs_err) << "\tRMSRES = " << std::format("{:.5e}", rmsres) << "\tHMIN = " << Hmin << "\n";

			lc_prs = lc_prs + 1;
 			//output_data(D, B1, B2, 1, 3, '-', 1, 0, "E:\\_work\\EHL-prg\\TEHL\\TEHL-v2\\_output\\PRS_DL");
			
		} //while prs
		
		// adjust h0 and check error for pressure integral
		if ((abs(prs_int - (2 * pi / 3)) / (2 * pi / 3) < 1e-4) && (rmsres <= res_errtol)) {

			std::cout << " << FORCE BALANCE IS REACHED << \n";
			extfile1  << " << FORCE BALANCE IS REACHED << \n";
			f_h0 = true;
		}
		else {
			// update h0 and continue calculations
			//newH0(D, B1, B2, prs_int, h0_rlx, Hmin, H1, H2, maxlc_prs0, maxlc_prs1, H0, maxlc_prs);
//			D.h0 = H0 * (a * a) / RX;
			
			lc_h0 = lc_h0 + 1;
			if (lc_h0 < maxlc_h0 + 1) {
				std::cout << "\nTHK ITER = " << lc_h0 << "\tNEW H0 = " << std::format("{:.5e}", RX * D.h0 / (a * a)) << "\n\n";
				extfile1  << "\nTHK ITER = " << lc_h0 << "\tNEW H0 = " << std::format("{:.5e}", RX * D.h0 / (a * a)) << "\n\n";
			}
		}
	} // while (force balance)

	extfile1.close();
	
	  /* free arrays for eps, rho, H */
		
		for (i = 0; i < Nx; ++i) {
			delete[] eps[i];
			delete[] rhox[i];
			delete[] rhoy[i];
			delete[] H[i];

		}
		
		delete[] eps;
		delete[] rhox;
		delete[] rhoy;
		delete[] H;

		for (i = 0; i < ncol; ++i) delete[] M[i];
		delete[] M;

		delete[] RHS;

		/* free Ti, Tj, Tx */
		delete[] Ti; 	delete[] Tj; 	delete[] Tx;

		/* free Ai, Ap, Ax */
		delete[] Ai; 	delete[] Ap; 	delete[] Ax;

		output_data(D, B1, B2, 1, 3, '-', 1, 0, "E:\\_work\\EHL-prg\\TEHL\\TEHL-v2\\_output\\PRS_DL");
} /* proc */

void calc_pressure_field_2D(Domain& D, Body B1, Body B2)
{
	double
		** eps,			/* dimensionless epsilon		*/
		** rhox,		/* dimensionless rhostar_x		*/
		** rhoy,		/* dimensionless rhostar_y		*/
		** H,
		** M,
		* RHS,
		* Tx, * Ax, * X;
	int
		* Ti, * Tj, * Ai, * Ap;

	int	Nx, Ny, i, j, kk, nzv, vdim, tc;
	double RX, RY, x_curr, y_curr, eps_curr, rhox_curr, rhoy_curr, a, b;
	double dx, dy, p1, um, prs;
	double alpha, theta, cA, cB, cC, cD, cE, cF;

	double H0, Hmin0, A1;
	double eH3_05i_j, eH3_i05_j, eH3_i_05j, eH3_i_j05;
	double T1, T2, gammai, derivi, ptol, swc, Xmax, summ1, perr, perr_old, Hmin;
	double h0_rlx, gs_rlx, jc_rlx, k, pint_curr, prs_errtol, rms_errtol;
	bool f_prs, f_h0, f_neg, h0rlx_flg;
	int lc_prs, lc_h0, maxlc_prs, maxlc_h0, cnt, strtIRj;
	double K_00, K_10, K_20, prs_int, Xstep, rms, g_int;
	double delta, H0min, dplim, Hmin_prv;
	string gapstr;

	dbvector			Tx_vec;
	std::vector<int>	Ti_vec, Tj_vec;

	int status;
	void* Symbolic, * Numeric;

	Nx = D.Nx;
	Ny = D.Ny;
	dx = D.dx;
	dy = D.dy;
	a = D.a;
	b = D.b;
	RX = pow(1 / B1.Rx + 1 / B2.Rx, -1); /* it needs for dimesionless eps */
	RY = pow(1 / B1.Ry + 1 / B2.Ry, -1);
	k = 1.0339 * pow(RY / RX, 0.636); //ellipticity ratio

	//um = (B1.vel[0] + B2.vel[0]) / 2.0; /* m/s */
	um = 1.0; // we have to throw out it from entire equation 
	alpha = D.ph * pow(D.a, 3) / (um * D.etaR * pow(RX, 2));
	theta = a / b;
	
	// relaxation coefficient for h0 (apply (Huang,2013)-method to Hamrock-formula for Hm)
	//h0_rlx = (-0.067) * RX * 2.69 * pow(D.W_HD, -1.067) * pow(D.U_HD, 0.67) * pow(D.G_HD, 0.53) *
		(1 - 0.61 * exp(-0.73 * k));

	//D.Hc = 2.69 * pow(D.W_HD, -0.067) * pow(D.U_HD, 0.67) * pow(D.G_HD, 0.53) *
	//	(1 - 0.61 * exp(-0.73 * k)) * pow((RX / a),2);
	
	//D.Hm = 3.63 * pow(D.W_HD, -0.073) * pow(D.U_HD, 0.68) * pow(D.G_HD, 0.49) *
	//	(1 - exp(-0.68 * k)) * pow((RX / a), 2); 

	//H0min = 0.001;
	H0min = 0.99 * D.Hm;

	//H0 = 0.1;
	H0 = -0.99 * D.delta * RX / (a * a);
	
	h0_rlx = 0.05;
	h0rlx_flg = true;

	maxlc_h0 = 150; // max h0 relaxation steps
	maxlc_prs = 3; // max pressure sweeps before H0-updating
	//Xstep = 0.003; // the "stepsize" in pressures update 
	dplim = 0.003;

	prs_errtol = 1e-3;
	rms_errtol = 1e-3;

	if (D.h0 == 0.0) D.h0 = H0 * (D.a * D.a) / pow(1 / B1.Rx + 1 / B2.Rx, -1);

	/* intervals for inner grid nodes (exclude nodes for BC) */
	int i_beg = 1;			/* zero pressure on (Y=-1) side */
	//int i_end = Nx - 1;		/* zero pressure on (Y=1) side */
	int i_end = Nx - 1;		/* zero pressure gradient on (Y=1) side */

	int j_beg = 1;			/* zero pressure on (X=-1) side */
	int j_end = Ny - 1;		/* zero pressure gradient on (X=1) side */
	//int j_end = Ny - 1;		/* zero pressure on (X=1) side */

	int IGNN = (i_end - i_beg) * (j_end - j_beg); /* number of inner grid nodes */
	int nrow = j_end - j_beg;					  /* number of inner grid nodes in row */
	int ncol = i_end - i_beg;					  /* number of inner grid nodes in row */

	/* set boundary conditions - Y axis */
	for (j = 0; j < Ny; ++j) {

		for (i = 0; i < i_beg; ++i) D.Area[i][j].p = 0.0;

		for (i = i_end; i < Nx; ++i) D.Area[i][j].p = 0.0;
	}

	/* set boundary conditions - X axis */
	for (i = 0; i < Nx; ++i) {

		for (j = 0; j < j_beg; ++j) D.Area[i][j].p = 0.0;

		for (j = j_end; j < Ny; ++j) D.Area[i][j].p = 0.0;
	}

	/* initialize dynamic arrays */
	{
		eps = new double* [Nx];
		rhox = new double* [Nx];
		rhoy = new double* [Nx];
		H = new double* [Nx];

		for (i = 0; i < Ny; ++i) {
			eps[i] = new double[Ny];
			rhox[i] = new double[Ny];
			rhoy[i] = new double[Ny];
			H[i] = new double[Ny];
		}

		M = new double* [IGNN];
		for (i = 0; i < IGNN; ++i) {
			M[i] = new double[IGNN];
			for (j = 0; j < IGNN; ++j) M[i][j] = 0.0;
		}

		RHS = new double[IGNN];
	}

	/* open external file for logging */
	std::fstream extfile1;
	extfile1.open("pc_log.txt", std::ios::out);
	extfile1 << "\t*** INITIAL VALUE FOR H0 = " << std::format("{:.5e}", D.h0 * RX / (a * a))
		<< " (" << std::format("{:.5e}", D.h0) << " m)\n\n";
	extfile1 << "\tH0 ITER = 1\n\n";
		
	lc_h0 = 1;
	f_h0 = false; // force balance flag

	//output_data(D, B1, B2, 1, 3, '-', 1, 0, "E:\\_work\\EHL-prg\\TEHL\\TEHL-v2\\_output\\PRS_DL");

	/* set the influence region */
			
	//int maxIRi = static_cast<int>(floor(ncol / 8)); // for full Jacobian (maxIRi = ncol)
	//int maxIRj = static_cast<int>(floor(nrow / 8)); // for full Jacobian (maxIRj = nrow)
			
	int maxIRi = 2; // for full Jacobian (maxIRi = ncol)
	int maxIRj = 2; // for full Jacobian (maxIRj = nrow)

	std::fstream extfile2;
	extfile2.open("pc_DPLAST.txt", std::ios::out);

	// loop for force balancing
	while ((!f_h0) && (lc_h0 < maxlc_h0 + 2)) {

		lc_prs = 1;
		f_prs = false; // pressure convergence flag

		// loop for pressure calculation
		while ((!f_prs) && (lc_prs < maxlc_prs + 1)) {
			
			calc_deform(D, B1, B2);
			calc_densviscNB(D, D.RHO, D.LSV);
			
			// calculate film thickness (dimensionless form)
			f_neg = false; // flag for negative Hmin
			Hmin = 1e3;
			for (i = 0; i < Nx; ++i) {
				for (j = 0; j < Ny; ++j) {

					// dimensionless (i,j)-position
					x_curr = (-1.0 * D.kx) + i * dx;
					y_curr = (-1.0 * D.ky) + j * dy;

					// H[i][j]
					p1 = (RX / (a * a)) * (D.Area[i][j].d1 + D.Area[i][j].d2);
					H[i][j] = (RX / (a * a)) * D.h0 + 0.5 * pow(x_curr, 2) +
						0.5 * (RX / RY) * (1 / pow(theta, 2)) * pow(y_curr, 2) + p1;

					if (H[i][j] < Hmin) Hmin = H[i][j]; // update the minimum film thickness

					if ((Hmin < 0) && (!f_neg)) {
						std::cout << "\n\t*** HMIN = " << Hmin << " WAS DETECTED AT POINT (" << x_curr << ", " << y_curr << ")\n";
						 extfile1 << "\n\t*** HMIN = " << Hmin << " WAS DETECTED AT POINT (" << x_curr << ", " << y_curr << ")\n";
						f_neg = true;
					}
				} // j loop
			} // i loop	
			
			if (Hmin > 0.0) Hmin_prv = Hmin;

			//output_data(D, B1, B2, 0, 5, 'Y', 0.0, 0, "E:/_work/EHL-prg/TEHL/TEHL-v2/_OUTPUT/PC_THK-XZ_DL_pos=0.000.TXT");
			//output_data(D, B1, B2, 0, 5, 'Y', 0.0, 0, "E:/_work/EHL-prg/TEHL/TEHL-v2/_OUTPUT/PC_PRS-XZ_DL_pos=0.000.TXT");
			
			//output_data(D, B1, B2, 0, 5, '-', 1, 0, "E:/_work/EHL-prg/TEHL/TEHL-v2/_OUTPUT/PC_THK_DL.TXT");
			//output_data(D, B1, B2, 0, 5, '-', 1, 0, "E:/_work/EHL-prg/TEHL/TEHL-v2/_OUTPUT/PC_PRS_DL.TXT");


			// correct H0 if negative Hmin is found
			if (Hmin < 0) {

				delta = abs(Hmin) + H0min;
				//delta = abs(Hmin) + Hmin_prv;
				D.h0 = D.h0 + delta * (a * a) / RX;

				//h0_rlx = h0_rlx / 10;
				//extfile1 << "\n\n h0_rlx = " << h0_rlx << "\n\n";

				std::cout << "\t*** CORRECT H0 BY " << std::format("{:.5e}", delta) << ", NEW VALUE IS " <<
					std::format("{:.5e}", D.h0 * RX / (a * a)) << " (" << std::format("{:.5e}", D.h0) << " m)\n";

				extfile1 << "\t*** CORRECT H0 BY " << std::format("{:.5e}", delta) << ", NEW VALUE IS " <<
					std::format("{:.5e}", D.h0 * RX / (a * a)) << " (" << std::format("{:.5e}", D.h0) << " m)\n";
				
				// write new nodal thickness and define new Hmin
				Hmin = 1e10;
				for (i = 0; i < Nx; ++i) {
					for (j = 0; j < Ny; ++j) {
						H[i][j] = H[i][j] + delta;
						if (H[i][j] < Hmin) Hmin = H[i][j];
						D.Area[i][j].h = H[i][j] * ((a * a) / RX);
					}
				}

				std::cout << "\t*** AFTER CORRECTION HMIN = " << std::format("{:.5e}", Hmin) << "\n\n";
				extfile1 << "\t*** AFTER CORRECTION HMIN = " << std::format("{:.5e}", Hmin) << "\n\n";
			} // if
			
			/* fill arrays for epsilon, rho and thickness (dimensionless form) */
			Hmin = 1e6;
			for (i = 0; i < Nx; ++i) {

				for (j = 0; j < Ny; ++j) {

					/* dimensionless (i,j)-position */
					x_curr = (-1.0 * D.kx) + i * dx;
					y_curr = (-1.0 * D.ky) + j * dy;
					
					/* H[i][j] */
					p1 = (RX / (a * a)) * (D.Area[i][j].d1 + D.Area[i][j].d2);
					H[i][j] = (RX / (a * a)) * D.h0 + 0.5 * pow(x_curr, 2) +
						0.5 * (RX / RY) * (1 / pow(theta, 2)) * pow(y_curr, 2) + p1;
					if (H[i][j] < 0.0) H[i][j] = 5e-5;
					if (H[i][j] < Hmin) Hmin = H[i][j];
					D.Area[i][j].h = H[i][j] * ((a * a) / RX);

					eps_curr = 0; rhox_curr = 0; rhoy_curr = 0;
					calc_rhoeps(i, j, D, B1, B2, eps_curr, rhox_curr, rhoy_curr);

					eps[i][j] = eps_curr;
					rhox[i][j] = rhox_curr;
					rhoy[i][j] = rhoy_curr;	
					//std::cout << "eps[" << i << "][" << j << "] = " << eps_curr << "\n";
					//std::cout << "rhox[" << i << "][" << j << "] = " << rhox_curr << "\n";
					//std::cout << "rhoy[" << i << "][" << j << "] = " << rhoy_curr << "\n";
					//std::cout << "H[" << i << "][" << j << "] = " << H[i][j] << "\n";

				} // j loop
			} // i loop
			
			/* main loop */
			cnt = 0; nzv = 0;
			for (i = i_beg; i < i_end; ++i) {
				for (j = j_beg; j < j_end; ++j) {

					/* calculate the nodal coefficients */
					calc_coef2D(i,j, D, B1, B2, eps, rhox, rhoy, H, cA, cB, cC, cD, cE, cF, 0);

					/* building the Jacobian */
					
					// (1) Write derivatives by nodal pressure: [i,j], [i,j+1], [i,j-1], [i+1,j], [i-1,j]
					p1 = cE - dFdP(i, j, i, j, rhox, rhoy, dx, dy, theta);
					M[cnt][cnt] = p1;
					if ((p1 < 0.0f) || (p1 > 0.0f)) {
						Ti_vec.push_back(cnt);
						Tj_vec.push_back(cnt);
						Tx_vec.push_back(p1);
						nzv = nzv + 1;
					}

					if (j < nrow) {

						p1 = cC - dFdP(i, j, i, j + 1, rhox, rhoy, dx, dy, theta);
						M[cnt][cnt + 1] = p1;
						if ((p1 < 0.0f) || (p1 > 0.0f)) {
							Ti_vec.push_back(cnt);
							Tj_vec.push_back(cnt + 1);
							Tx_vec.push_back(p1);
							nzv = nzv + 1;
						}
					}

					if (j != j_beg) {

						p1 = cD - dFdP(i, j, i, j - 1, rhox, rhoy, dx, dy, theta);
						M[cnt][cnt - 1] = p1;
						if ((p1 < 0.0f) || (p1 > 0.0f)) {
							Ti_vec.push_back(cnt);
							Tj_vec.push_back(cnt - 1);
							Tx_vec.push_back(p1);
							nzv = nzv + 1;
						}
					}

					if (i != i_end - 1) {
						p1 = cA - dFdP(i, j, i + 1, j, rhox, rhoy, dx, dy, theta);
						M[cnt][cnt + nrow] = p1;
						if ((p1 < 0.0f) || (p1 > 0.0f)) {
							Ti_vec.push_back(cnt);
							Tj_vec.push_back(cnt + nrow);
							Tx_vec.push_back(p1);
							nzv = nzv + 1;
						}
					}

					if (cnt >= nrow) {
						p1 = cB - dFdP(i, j, i - 1, j, rhox, rhoy, dx, dy, theta);
						M[cnt][cnt - nrow] = p1;
						if ((p1 < 0.0f) || (p1 > 0.0f)) {
							Ti_vec.push_back(cnt);
							Tj_vec.push_back(cnt - nrow);
							Tx_vec.push_back(p1);
							nzv = nzv + 1;
						}
					}

					// (2) Write derivatives by other nodal pressures from the influence region */
					for (int IRi = 0; IRi <= maxIRi; ++IRi) {

						if (IRi == 0) strtIRj = 2; else strtIRj = 1;

						if ((i + IRi) <= (i_end - 1)) {

							if (IRi >= 2) {
								p1 = dFdP(i, j, (i + IRi), j, rhox, rhoy, dx, dy, theta);
								M[cnt][(i - i_beg + IRi) * nrow + j - j_beg] = p1;
								Ti_vec.push_back(cnt);
								Tj_vec.push_back((i - i_beg + IRi) * nrow + j - j_beg);
								Tx_vec.push_back(p1);
								nzv = nzv + 1;
							}

							for (int IRj = strtIRj; IRj <= maxIRj; ++IRj) {

								if ((j - IRj) >= j_beg) {
									p1 = dFdP(i, j, (i + IRi), (j - IRj), rhox, rhoy, dx, dy, theta);
									M[cnt][(i - i_beg + IRi) * nrow + j - j_beg - IRj] = p1;
									Ti_vec.push_back(cnt);
									Tj_vec.push_back((i - i_beg + IRi) * nrow + j - j_beg - IRj);
									Tx_vec.push_back(p1);
									nzv = nzv + 1;
								}

								if ((j + IRj) <= (j_end - 1)) {
									p1 = dFdP(i, j, (i + IRi), (j + IRj), rhox, rhoy, dx, dy, theta);
									M[cnt][(i - i_beg + IRi) * nrow + j - j_beg + IRj] = p1;
									Ti_vec.push_back(cnt);
									Tj_vec.push_back((i - i_beg + IRi) * nrow + j - j_beg + IRj);
									Tx_vec.push_back(p1);
									nzv = nzv + 1;
								}

							} // IRj loop

						} // if

						if ((i - IRi) >= i_beg) {

							if (IRi >= 2) {
								p1 = dFdP(i, j, (i - IRi), j, rhox, rhoy, dx, dy, theta);
								M[cnt][(i - i_beg - IRi) * nrow + j - j_beg] = p1;
								Ti_vec.push_back(cnt);
								Tj_vec.push_back((i - i_beg - IRi) * nrow + j - j_beg);
								Tx_vec.push_back(p1);
								nzv = nzv + 1;
							}

							for (int IRj = strtIRj; IRj <= maxIRj; ++IRj) {

								if ((j - IRj) >= j_beg) {
									p1 = dFdP(i, j, (i - IRi), (j - IRj), rhox, rhoy, dx, dy, theta);
									M[cnt][(i - i_beg - IRi) * nrow + j - j_beg - IRj] = p1;
									Ti_vec.push_back(cnt);
									Tj_vec.push_back((i - i_beg - IRi) * nrow + j - j_beg - IRj);
									Tx_vec.push_back(p1);
									nzv = nzv + 1;
								}

								if ((j + IRj) <= (j_end - 1)) {
									p1 = dFdP(i, j, (i - IRi), (j + IRj), rhox, rhoy, dx, dy, theta);
									M[cnt][(i - i_beg - IRi) * nrow + j - j_beg + IRj] = p1;
									Ti_vec.push_back(cnt);
									Tj_vec.push_back((i - i_beg - IRi) * nrow + j - j_beg + IRj);
									Tx_vec.push_back(p1);
									nzv = nzv + 1;
								}

							} // IRJ loop
						} // if
					} // IRi loop

					/* calculate residuals for RHS */
					RHS[cnt] = -(cA * D.Area[i + 1][j].p + cB * D.Area[i - 1][j].p + cC * D.Area[i][j + 1].p +
						cD * D.Area[i][j - 1].p + cE * D.Area[i][j].p - cF);
					
					cnt = cnt + 1;
				} // j loop
			} // i loop

			/* write Jacobian in textfile */
			/*
			std::fstream extfile;
			extfile.open("_M.TXT", std::ios::out);
			for (i = 0; i < IGNN; ++i)
				for (j = 0; j < IGNN; ++j) {
					if (j != (IGNN - 1)) extfile << M[i][j] << "\t"; else extfile << M[i][j] << "\n";
				}

			extfile.close();
			*/

			/* solving the linear system */

				/* write data from vectors to arrays */
				int vdim = Ti_vec.size();
				int* Ti, * Tj;
				double* Tx;

				Ti = new int[vdim];
				Tj = new int[vdim];
				Tx = new double[vdim];

				for (i = 0; i < vdim; ++i) {

					Ti[i] = Ti_vec[i];
					Tj[i] = Tj_vec[i];
					Tx[i] = Tx_vec[i];
				}

				Ti_vec.clear();
				Tj_vec.clear();
				Tx_vec.clear();

				int* Ai, * Ap;
				double* Ax, * X;

				Ai = new int[nzv];
				Ap = new int[IGNN + 1];
				Ax = new double[nzv];
				X = new double[IGNN];

				int status;
				void* Symbolic, * Numeric;

				/* convert triplet form of matrix to compressed sparse column form */
				status = umfpack_di_triplet_to_col(IGNN, IGNN, nzv, Ti, Tj, Tx, Ap, Ai, Ax, NULL);

				/* free Ti, Tj, Tx */
				delete[] Ti; 	delete[] Tj; 	delete[] Tx;

				/* symbolic and numeric factorization */
				status = umfpack_di_symbolic(IGNN, IGNN, Ap, Ai, Ax, &Symbolic, NULL, NULL);

				status = umfpack_di_numeric(Ap, Ai, Ax, Symbolic, &Numeric, NULL, NULL);

				umfpack_di_free_symbolic(&Symbolic);

				/* solve the linear system */
				status = umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, X, RHS, Numeric, NULL, NULL);

				/* free the memory associated with the numeric factorization */
				umfpack_di_free_numeric(&Numeric);

			// calculate the norm of solution vector
			//find max X[i]
			Xmax = 0.0;
			for (i = 0; i < IGNN; ++i) {
				if (Xmax < abs(X[i])) Xmax = abs(X[i]);
				//summ0 = summ0 + X[i] * X[i];
				//std::cout << X[i] << "\n";
			}
					
			//summ0 = pow(summ0, 0.5);

			/* write the pressure field in domain */
			
			kk = 0; // index in the result vector X 
			prs_int = 0.0; // summ for pressure integral 
			summ1 = 0.0; // summ for pressure error
			for (i = i_beg; i < i_end; ++i)
				for (j = j_beg; j < j_end; ++j) {
					
					//D.Area[i][j].p = D.Area[i][j].pold + Xstep * X[kk] / summ0; 
					
					if (Xmax < 0.05) Xstep = 1.0;
						else Xstep = abs(Xmax / dplim);
					if (Xstep != 0.0) D.Area[i][j].p = D.Area[i][j].pold + X[kk] / Xstep;
					else D.Area[i][j].p = D.Area[i][j].pold;
					
					//D.Area[i][j].p = D.Area[i][j].pold + Xstep * X[kk];
					
					if (D.Area[i][j].p < 0.0) D.Area[i][j].p = 0.0; //caviation condition
					
					summ1 = summ1 + abs(D.Area[i][j].pold - D.Area[i][j].p);
					D.Area[i][j].pold = D.Area[i][j].p;

					prs_int = prs_int + D.Area[i][j].p;
					
					if ((lc_h0 >= maxlc_h0) && (lc_prs >= maxlc_prs)) {

						/* dimensionless (i,j)-position */
						x_curr = (-1.0 * D.kx) + i * dx;
						y_curr = (-1.0 * D.ky) + j * dy;
						//extfile2 << x_curr << "\t" << y_curr << "\t" << Xstep * X[kk] / summ0 << "\n";
						extfile2 << x_curr << "\t" << y_curr << "\t" << Xstep * X[kk] << "\n";
					}
					
					kk = kk + 1;
				}
		
			/* calculate the pressure error and the RMS-residual */
			rms = 0.0; //Hmin = 1e6;
			for (i = i_beg; i < i_end; ++i)
				for (j = j_beg; j < j_end; ++j) {
					
					if (D.Area[i][j].p > 1e-5) {
						calc_coef2D(i, j, D, B1, B2, eps, rhox, rhoy, H, cA, cB, cC, cD, cE, cF, 0);
						p1 = cA * D.Area[i + 1][j].p + cB * D.Area[i - 1][j].p + cC * D.Area[i][j + 1].p +
							cD * D.Area[i][j - 1].p + cE * D.Area[i][j].p - cF;
						rms = rms + p1 * p1;
					}
					
				}
			rms = pow(rms / IGNN, 0.5); // RMSRES in current iteration
			prs_int = prs_int * dx * dy; // pressure integral

			std::cout << "PRS ITER = " << lc_prs << "\tPRS INT = " << std::format("{:.5e}", prs_int) <<
				"\tHMIN = " << std::format("{:.5e}", Hmin) << "\tPRS ERR = " << std::format("{:.5e}", summ1 / prs_int) <<
				"\tRMS RES = " << std::format("{:.5e}", rms) << "\n";

			extfile1 << "PRS ITER = " << lc_prs << "\tPRS INT = " << std::format("{:.5e}", prs_int) <<
				"\tHMIN = " << std::format("{:.5e}", Hmin) << "\tPRS ERR = " << std::format("{:.5e}", summ1 / prs_int) <<
				"\tRMS RES = " << std::format("{:.5e}", rms) << "\n";
			
			lc_prs = lc_prs + 1;
			
		} // prs loop

		// adjust H0 and check error for pressure integral
		g_int = 2 * pi / 3;
		if ((abs(prs_int - g_int) / g_int < 1e-4) && (rms <= rms_errtol)) {

			std::cout << " << FORCE BALANCE IS REACHED << \n";
			 extfile1 << " << FORCE BALANCE IS REACHED << \n";
			f_h0 = true;
		}
		else {
			if ((Hmin < 0.01) && (!h0rlx_flg)) {
				h0_rlx = 0.1 * h0_rlx; 
				h0rlx_flg = true;
			}
			p1 = h0_rlx * (g_int - prs_int) * a * a / RX;
			if (p1 > 0) gapstr = "HAS DECREASED"; else gapstr = "HAS INCREASED";
			if (abs(p1) < 1e-12) gapstr = "HASN'T CHANGED";
			D.h0 = D.h0 - p1;

			p1 = RX * D.h0 / (a * a);
			if (lc_h0 < maxlc_h0 + 1) {
				std::cout << "\n\tHO ITER = " << lc_h0 << "\tNEW H0 = " << std::format("{:.5e}", p1)
					<< "(" << std::format("{:.5e}", D.h0) << " m)" << "\t*** THE GAP BETWEEN THE SURFACES " << gapstr << "\n\n";
				extfile1 << "\n\tH0 ITER = " << lc_h0 << "\tNEW H0 = " << std::format("{:.5e}", p1)
					<< "(" << std::format("{:.5e}", D.h0) << " m)" << "\t*** THE GAP BETWEEN THE SURFACES " << gapstr << "\n\n";
			}

			//D.h0 = D.h0 - h0_rlx * (g_int - prs_int) * (a * a) / RX;

			lc_h0 = lc_h0 + 1;

			
		} // else
					
	} // while h0

	extfile1.close();
	 
	// write nodal residuals in file
	extfile1.open("pc_NDRES.txt", std::ios::out);
	rms = 0.0; //Hmin = 1e6;
	for (i = i_beg; i < i_end; ++i)
		for (j = j_beg; j < j_end; ++j) {
			
			/* dimensionless (i,j)-position */
			x_curr = (-1.0 * D.kx) + i * dx;
			y_curr = (-1.0 * D.ky) + j * dy;
			
			if (D.Area[i][j].p > 1e-5) {
				calc_coef2D(i, j, D, B1, B2, eps, rhox, rhoy, H, cA, cB, cC, cD, cE, cF, 0);
				p1 = cA * D.Area[i + 1][j].p + cB * D.Area[i - 1][j].p + cC * D.Area[i][j + 1].p +
					cD * D.Area[i][j - 1].p + cE * D.Area[i][j].p - cF;
			}
			else p1 = 0.0;

			extfile1 << std::format("{:.5e}", x_curr) << "\t" << std::format("{:.5e}",y_curr) 
				<< "\t" << std::format("{:.5e}", p1) << "\n";
						
		} // j and i loop
	extfile1.close();
	
	extfile2.close();
	
	/* free arrays for eps, rho, H */

	for (i = 0; i < Nx; ++i) {
		delete[] eps[i];
		delete[] rhox[i];
		delete[] rhoy[i];
		delete[] H[i];

	}

	delete[] eps;
	delete[] rhox;
	delete[] rhoy;
	delete[] H;

	for (i = 0; i < ncol; ++i) delete[] M[i];
	delete[] M;

	delete[] RHS;
	

	//output_data(D, B1, B2, 1, 3, '-', 1, 0, "E:\\_work\\EHL-prg\\TEHL\\TEHL-v2\\_output\\PRS_DL");
} /* proc */


/******************** MULTIGRID PROCEDURES ***************************/
void coarsenLHS()
{

}

void coarsenRHS()
{

}

void relax()
{

}

void refine()
{

}


void cycle()
{

}

void fmg_interpolate()
{

}

void fmg()
{

}


void calc_pressure_field_2D_FMG(Domain& D, Body B1, Body B2, int kmax) {
	
	MGData		MG{};
	int lc_lev, arrbnd;

	MG.level = new LData[kmax];
	/* initialize MGData structures */
	lc_lev = 1;
	while (lc_lev <= kmax) {
		arrbnd = pow(2, lc_lev + 1) + 1;
		MG.level[lc_lev - 1].P = new double* [arrbnd];
		MG.level[lc_lev - 1].H = new double* [arrbnd];

		for (int i = 0; i < arrbnd; ++i) {
			MG.level[lc_lev - 1].P[i] = new double[arrbnd];
			MG.level[lc_lev - 1].H[i] = new double[arrbnd];
		}
		lc_lev = lc_lev + 1;
	}

	MG.level[0].P[0][0] = 1.0;
	MG.level[0].P[4][4] = 1.0;

	MG.level[1].P[0][0] = 1.0;
	MG.level[1].P[8][8] = 1.0;

	MG.level[2].P[0][0] = 1.0;
	MG.level[2].P[16][16] = 1.0;

	MG.level[3].P[0][0] = 1.0;
	MG.level[3].P[32][32] = 1.0;

	

	//static_cast<int>((D.Nx - 1) / 2);
	//MG.level

	/* free all MG-structures */
	while (lc_lev <= kmax) {
		arrbnd = pow(2, lc_lev + 1) + 1;
		
		for (int i = 0; i < arrbnd; ++i) {
			delete [] MG.level[lc_lev - 1].P[i];
			delete [] MG.level[lc_lev - 1].H[i];
		}
		delete[] MG.level[lc_lev - 1].P;
		delete[] MG.level[lc_lev - 1].H;
		
		lc_lev = lc_lev + 1;
	}


} // proc

